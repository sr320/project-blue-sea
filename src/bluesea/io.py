"""Readers for the tool outputs the cascade consumes.

Everything here parses a file produced by an external tool on a compute node.
The parsers are deliberately tolerant about columns they do not need and
strict about the ones they do, because a silently mis-parsed SignalP table
would propagate into the ranking without ever raising.
"""

from __future__ import annotations

import csv
import gzip
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
    """SignalP 6 summary table: ID, Prediction, probabilities, CS position."""
    out: dict[str, dict[str, Any]] = {}
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            protein_id, prediction = fields[0], fields[1]
            # Column layout varies with the model used; the signal probability
            # is the largest non-"OTHER" probability in the numeric columns.
            probs = [_maybe_float(f) for f in fields[2:]]
            probs = [p for p in probs if p is not None]
            out[protein_id] = {
                "prediction": prediction,
                "probability": max(probs) if prediction != "OTHER" and probs else 0.0,
            }
    return out


def read_deeptmhmm(path: str | Path) -> dict[str, dict[str, Any]]:
    """DeepTMHMM 3-line output: per-residue topology string per protein.

    Labels: S signal, M transmembrane, I inside (cytoplasmic), O outside,
    B beta-barrel. The outside count is what the exposure filter cares about.
    """
    out: dict[str, dict[str, Any]] = {}
    protein_id: str | None = None
    with _open(path) as handle:
        lines = [ln.rstrip("\n") for ln in handle if ln.strip()]
    for _index, line in enumerate(lines):
        if line.startswith(">"):
            protein_id = line[1:].split("|")[0].split()[0]
            continue
        if protein_id is None or set(line) - set("SMIOB"):
            continue
        out[protein_id] = {
            "topology_string": line,
            "n_transmembrane": _count_runs(line, "M") + _count_runs(line, "B"),
            "extracellular_aa": line.count("O"),
            "has_signal": "S" in line,
        }
        protein_id = None
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
            continue
        hits.setdefault(row[0], []).append(
            {
                "subject": row[1],
                "identity_pct": float(row[2]),
                "length": int(row[3]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "evalue": float(row[10]),
                "bitscore": float(row[11]),
            }
        )
    return hits


def _iter_rows(path: str | Path) -> Iterator[list[str]]:
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            yield line.rstrip("\n").split("\t")


def best_host_hit(
    hits: list[dict[str, Any]], query_length: int
) -> HostHit | None:
    """The salmon hit that matters is the most similar one, not the top-scoring
    one: cross-reactivity risk is set by the closest match at any coverage."""
    if not hits:
        return None
    best = max(hits, key=lambda h: h["identity_pct"])
    coverage = (best["qend"] - best["qstart"] + 1) / max(query_length, 1)
    return HostHit(
        accession=best["subject"],
        identity_pct=best["identity_pct"],
        evalue=best["evalue"],
        query_coverage=coverage,
    )


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


def read_expression(path: str | Path, stages: list[str]) -> dict[str, dict[str, float]]:
    """Gene-by-stage TPM matrix. Columns absent from the file are left out of
    the per-gene dict rather than zero-filled, so "not assayed" stays visible."""
    out: dict[str, dict[str, float]] = {}
    for row in read_tsv(path):
        gene_id = row.get("gene_id") or row.get("gene")
        if not gene_id:
            continue
        out[gene_id] = {
            stage: float(row[stage])
            for stage in stages
            if stage in row and row[stage] not in ("", "NA")
        }
    return out


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
