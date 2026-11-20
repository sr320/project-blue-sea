"""Readers for the tool outputs the cascade consumes.

Everything here parses a file produced by an external tool on a compute node.
The parsers are deliberately tolerant about columns they do not need and
strict about the ones they do, because a silently mis-parsed SignalP table
would propagate into the ranking without ever raising.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import yaml

from .records import Epitope, HostHit, Ortholog, Topology


def _open(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path) as handle:
        return yaml.safe_load(handle)


def read_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield (identifier, description, sequence) from a protein FASTA."""
    header: str | None = None
    chunks: list[str] = []
    with _open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield (*_split_header(header), "".join(chunks))
                header = line[1:]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield (*_split_header(header), "".join(chunks))


def _split_header(header: str) -> tuple[str, str]:
    parts = header.split(None, 1)
    return parts[0], (parts[1] if len(parts) > 1 else "")


def read_tsv(path: str | Path) -> Iterator[dict[str, str]]:
    """Tab-separated with a header row. Lines starting with '#' are comments,
    including before the header -- the control tables carry a note at the top
    explaining what the file is for, and that must not become the header."""
    with _open(path) as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    yield from csv.DictReader(lines, delimiter="\t")


def read_signalp(path: str | Path) -> dict[str, dict[str, Any]]:
    """SignalP 6 eukaryotic summary: ID, Prediction, OTHER, SP, CS Position."""
    out = {}
    for fields in _iter_rows(path):
        if len(fields) < 4:
            raise ValueError(f"Malformed SignalP row: {fields}")
        pid, prediction = fields[:2]
        if pid in out or prediction not in {"OTHER", "SP"}:
            raise ValueError(f"Duplicate ID or unsupported SignalP class: {pid}")
        probability = float(fields[3])
        if not 0 <= probability <= 1:
            raise ValueError(f"Invalid SignalP probability for {pid}")
        cleavage = re.search(r"CS pos:\s*(\d+)-", " ".join(fields[4:]))
        out[pid] = {
            "prediction": prediction,
            "probability": probability,
            "signal_peptide_end": int(cleavage[1]) if cleavage else None,
        }
    return out


def read_deeptmhmm(path: str | Path) -> dict[str, dict[str, Any]]:
    """Strict DeepTMHMM 1 three-line parser; preserves residue coordinates."""
    with _open(path) as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if len(lines) % 3:
        raise ValueError("DeepTMHMM output must contain complete three-line records")
    out = {}
    for index in range(0, len(lines), 3):
        header, sequence, topology = lines[index : index + 3]
        if not header.startswith(">") or len(sequence) != len(topology):
            raise ValueError(f"Malformed DeepTMHMM record at line {index + 1}")
        pid = header[1:].split()[0]
        if pid in out or set(topology) - set("SMIOB"):
            raise ValueError(f"Duplicate ID or unknown topology labels: {pid}")
        out[pid] = {
            "sequence": sequence,
            "topology_string": topology,
            "n_transmembrane": _count_runs(topology, "M") + _count_runs(topology, "B"),
            "extracellular_aa": topology.count("O"),
            "has_signal": "S" in topology,
        }
    return out


def _count_runs(string: str, char: str) -> int:
    runs, previous = 0, ""
    for current in string:
        if current == char and previous != char:
            runs += 1
        previous = current
    return runs


def read_blast_tabular(path: str | Path) -> dict[str, list[dict[str, Any]]]:
    """BLAST/DIAMOND outfmt 6 with the standard 12 columns.

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send
    evalue bitscore
    """
    hits: dict[str, list[dict[str, Any]]] = {}
    for row in _iter_rows(path):
        if len(row) < 12:
            raise ValueError(f"Malformed BLAST row: {row}")
        hits.setdefault(row[0], []).append(
            {
                "subject": row[1],
                "identity_pct": float(row[2]),
                "length": int(row[3]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "evalue": float(row[10]),
                "bitscore": float(row[11]),
                "sstart": int(row[8]),
                "send": int(row[9]),
                "qseq": row[12] if len(row) > 12 else None,
                "sseq": row[13] if len(row) > 13 else None,
            }
        )
    return hits


def _iter_rows(path: str | Path) -> Iterator[list[str]]:
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def host_alignments(hits: list[dict[str, Any]], query_length: int) -> list[HostHit]:
    return [
        HostHit(
            accession=h["subject"],
            identity_pct=h["identity_pct"],
            evalue=h["evalue"],
            query_coverage=(abs(h["qend"] - h["qstart"]) + 1) / max(query_length, 1),
        )
        for h in hits
    ]


def best_host_hit(
    hits: list[dict[str, Any]],
    query_length: int,
    cfg: dict[str, Any] | None = None,
) -> HostHit | None:
    """Select among significant, sufficiently covered HSPs before identity.

    All HSPs are retained separately for review of short conserved patches.
    """
    candidates = host_alignments(hits, query_length)
    if cfg:
        candidates = [
            h
            for h in candidates
            if h.evalue <= cfg["max_evalue"]
            and h.query_coverage >= cfg["min_query_coverage"]
        ]
    return max(candidates, key=lambda h: (h.identity_pct, h.query_coverage), default=None)


def read_epitopes(path: str | Path) -> dict[str, list[Epitope]]:
    """BepiPred/DiscoTope-style table: protein_id, start, end, score, kind."""
    out: dict[str, list[Epitope]] = {}
    for row in read_tsv(path):
        out.setdefault(row["protein_id"], []).append(
            Epitope(
                start=int(row["start"]),
                end=int(row["end"]),
                score=float(row["score"]),
                kind=row.get("kind", "linear"),
                sequence=row.get("sequence", ""),
            )
        )
    return out


def read_protein_gene_map(path: str | Path) -> dict[str, tuple[str, str]]:
    """protein accession -> (gene_id, gene_symbol), from a genome GFF.

    This is the only authoritative link between the two. NCBI `protein.faa`
    headers carry no gene tag at all, so without the GFF every protein looks
    like its own gene and nothing keyed on genes can join to it.

    The two annotations here do not agree on how to say it, and a parser
    written against one silently returns nothing for the other:

        RefSeq   protein_id=XP_040583074.1 gene=LOC121131766
                 Dbxref=GeneID:121131766
        GenBank  protein_id=QQP55611.1     locus_tag=FKW44_000001

    So gene_id prefers the stable numeric GeneID, falls back to the locus tag,
    then to the gene name; the symbol prefers the readable name. CDS features
    repeat per exon, and the first mention of a protein wins.
    """
    out: dict[str, tuple[str, str]] = {}
    for row in _iter_rows(path):
        if len(row) < 9 or row[2] != "CDS":
            continue
        attributes = _gff_attributes(row[8])
        protein_id = attributes.get("protein_id")
        if not protein_id or protein_id in out:
            continue
        gene_id = (
            _gene_dbxref(attributes.get("Dbxref", ""))
            or attributes.get("locus_tag")
            or attributes.get("gene")
            or ""
        )
        symbol = attributes.get("gene") or attributes.get("locus_tag") or ""
        if gene_id:
            out[protein_id] = (gene_id, symbol)
    return out


def _gff_attributes(field: str) -> dict[str, str]:
    attributes = {}
    for item in field.split(";"):
        key, _, value = item.partition("=")
        if value:
            attributes[key.strip()] = value.strip()
    return attributes


def _gene_dbxref(dbxref: str) -> str:
    for item in dbxref.split(","):
        if item.startswith("GeneID:"):
            return item
    return ""


def read_expression(path: str | Path, stages: list[str]) -> dict[str, dict[str, float]]:
    """Gene-by-stage TPM matrix. Columns absent from the file are left out of
    the per-gene dict rather than zero-filled, so "not assayed" stays visible."""
    out: dict[str, dict[str, float]] = {}
    expression_stages(path, stages)
    for row in read_tsv(path):
        gene_id = row.get("gene_id") or row.get("gene")
        if not gene_id:
            continue
        out[gene_id] = {
            stage: float(row[stage])
            for stage in stages
            if stage in row and row[stage] not in ("", "NA")
        }
    if any(
        not math.isfinite(v) or v < 0 for values in out.values() for v in values.values()
    ):
        raise ValueError(f"Invalid TPM in {path}")
    return out


def species_stages(config: dict[str, Any], species: str) -> list[str]:
    return config["species"][species].get("life_stages", config.get("life_stages", []))


def expression_stages(path: str | Path, stages: list[str]) -> list[str]:
    with _open(path) as handle:
        header = next((line for line in handle if not line.startswith("#")), "")
    columns = header.rstrip().split("\t")
    unknown = set(columns) - {"gene_id", "gene", *stages}
    if unknown:
        raise ValueError(f"Unknown expression stages in {path}: {sorted(unknown)}")
    return [stage for stage in stages if stage in columns]


def read_rnai(path: str | Path) -> dict[str, dict[str, str]]:
    """LiceBase knockdown phenotypes, keyed by gene."""
    out: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        gene_id = row.get("gene_id") or row.get("gene")
        if not gene_id:
            continue
        out[gene_id] = {
            "phenotype": row.get("phenotype", "untested"),
            "source": row.get("source", ""),
            "note": row.get("note", ""),
        }
    return out


def read_rbh(path: str | Path, species: str) -> dict[str, Ortholog]:
    """Reciprocal-best-hit table produced by `bluesea.orthology`."""
    out: dict[str, Ortholog] = {}
    for row in read_tsv(path):
        out[row["query"]] = Ortholog(
            species=species,
            protein_id=row["subject"],
            identity_pct=float(row["identity_pct"]),
            evalue=float(row["evalue"]),
            reciprocal_best=row.get("reciprocal", "true").lower() == "true",
            epitope_region_identity_pct=_maybe_float(
                row.get("epitope_region_identity_pct")
            ),
        )
    return out


def _maybe_float(value: str | None) -> float | None:
    if value in (None, "", "NA", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def topology_from_tools(
    protein_id: str,
    signalp: dict[str, dict[str, Any]],
    tmhmm: dict[str, dict[str, Any]],
    secretome: set[str] | None = None,
    gpi: set[str] | None = None,
) -> Topology:
    sp = signalp.get(protein_id, {})
    tm = tmhmm.get(protein_id, {})
    prediction = sp.get("prediction", "OTHER")
    return Topology(
        topology_string=tm.get("topology_string", ""),
        signal_peptide_end=sp.get("signal_peptide_end"),
        has_signal_peptide=prediction not in ("OTHER", ""),
        signalp_prob=float(sp.get("probability", 0.0)),
        n_transmembrane=int(tm.get("n_transmembrane", 0)),
        gpi_anchored=protein_id in (gpi or set()),
        extracellular_aa=int(tm.get("extracellular_aa", 0)),
        secretome_confirmed=protein_id in (secretome or set()),
    )


def write_tsv(path: str | Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
