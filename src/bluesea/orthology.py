"""Reciprocal best hits between the two louse proteomes.

The cascade runs once per species; this is what lets the two shortlists be
intersected afterwards. A protein surviving in both genera is the smaller,
more interesting set -- conservation across *Lepeophtheirus* and *Caligus* is
evidence of shared sequence that merits experimental follow-up.
It does not establish functional essentiality or cross-protection.

Conservation is checked where it has to hold: across the epitope-bearing
region. Two orthologs at 60% identity overall are no use if every predicted
epitope sits in the divergent 40%.
"""

from __future__ import annotations

from typing import Any

from .records import Epitope, Ortholog, ProteinRecord


def reciprocal_best_hits(
    forward: dict[str, list[dict[str, Any]]],
    reverse: dict[str, list[dict[str, Any]]],
    cfg: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Pair A->B and B->A best hits that agree.

    `forward` and `reverse` are parsed BLAST tables (see io.read_blast_tabular)
    for A-against-B and B-against-A respectively.
    """
    best_reverse = {
        query: max(hits, key=lambda h: h["bitscore"])["subject"]
        for query, hits in reverse.items()
        if hits
    }

    pairs: dict[str, dict[str, Any]] = {}
    for query, hits in forward.items():
        eligible = [
            hit
            for hit in hits
            if hit["evalue"] <= cfg["max_evalue"]
            and hit["identity_pct"] >= cfg["min_identity_pct"]
        ]
        if not eligible:
            continue
        best = max(hits, key=lambda h: h["bitscore"])
        if best not in eligible:
            continue
        pairs[query] = {
            "subject": best["subject"],
            "identity_pct": best["identity_pct"],
            "evalue": best["evalue"],
            "reciprocal": best_reverse.get(best["subject"]) == query,
            "alignment_hit": best,
        }
    return pairs


def epitope_region_identity(
    query_sequence: str,
    subject_sequence: str,
    epitopes: list[Epitope],
    alignment: list[tuple[int, int]] | None = None,
) -> float | None:
    """Percent identity restricted to the epitope-bearing residues.

    `alignment` is a list of (query_index, subject_index) 0-based residue
    pairs from the pairwise alignment. Without one -- for instance when only a
    tabular BLAST summary is on disk -- this returns None and the caller
    flags the record rather than guessing, since an ungapped index-for-index
    comparison would silently mis-measure any ortholog pair with an indel.
    """
    if not epitopes or alignment is None:
        return None

    epitope_positions = {
        position
        for epitope in epitopes
        for position in range(epitope.start - 1, epitope.end)
    }
    compared = matched = 0
    for query_index, subject_index in alignment:
        if query_index not in epitope_positions:
            continue
        if query_index >= len(query_sequence) or subject_index >= len(subject_sequence):
            continue
        compared += 1
        matched += query_sequence[query_index] == subject_sequence[subject_index]
    if compared == 0:
        return None
    # Unaligned/deleted epitope residues count as nonmatches.
    return 100.0 * matched / len(epitope_positions)


def attach_orthologs(
    records: dict[str, ProteinRecord],
    pairs: dict[str, dict[str, Any]],
    other_species: str,
    cfg: dict[str, Any],
    other_sequences: dict[str, str] | None = None,
    alignments: dict[str, list[tuple[int, int]]] | None = None,
) -> None:
    for protein_id, pair in pairs.items():
        record = records.get(protein_id)
        if record is None:
            continue
        region_identity = None
        coverage = None
        alignment = (alignments or {}).get(protein_id)
        hit = pair.get("alignment_hit", {})
        if alignment is None and hit.get("qseq") and hit.get("sseq"):
            alignment = alignment_positions(hit)
        if other_sequences and alignment is not None:
            region_identity = epitope_region_identity(
                record.sequence,
                other_sequences.get(pair["subject"], ""),
                record.accessible_epitopes(),
                alignment,
            )
            positions = {
                p for e in record.accessible_epitopes() for p in range(e.start - 1, e.end)
            }
            covered = {q for q, _ in alignment} & positions
            coverage = len(covered) / len(positions) if positions else None
        record.orthologs.append(
            Ortholog(
                species=other_species,
                protein_id=pair["subject"],
                identity_pct=pair["identity_pct"],
                evalue=pair["evalue"],
                reciprocal_best=pair["reciprocal"],
                epitope_region_identity_pct=region_identity,
                epitope_region_coverage=coverage,
            )
        )
        if (
            region_identity is not None
            and region_identity < cfg["epitope_region_min_identity_pct"]
        ):
            record.flag("epitopes_not_conserved_across_species")


def intersect_shortlists(
    shortlist_a: list[ProteinRecord],
    shortlist_b: list[ProteinRecord],
    cfg: dict[str, Any] | None = None,
) -> list[tuple[ProteinRecord, ProteinRecord]]:
    """Candidates that survived the cascade in both species.

    Matched through the recorded orthology, not through name similarity --
    annotation names across these two genomes are not reliable enough to join
    on, and a wrong join here would be invisible downstream.
    """
    cfg = cfg or {}
    minimum = cfg.get("epitope_region_min_identity_pct", 50.0)
    coverage = cfg.get("epitope_region_min_coverage", 0.8)

    def supported(o: Ortholog) -> bool:
        return (
            o.reciprocal_best
            and o.epitope_region_identity_pct is not None
            and o.epitope_region_identity_pct >= minimum
            and o.epitope_region_coverage is not None
            and o.epitope_region_coverage >= coverage
        )

    by_id_b = {(record.species, record.protein_id): record for record in shortlist_b}
    for record in [*shortlist_a, *shortlist_b]:
        record.flags = [f for f in record.flags if f != "cross_species_candidate"]
    out: list[tuple[ProteinRecord, ProteinRecord]] = []
    for record in shortlist_a:
        for ortholog in record.orthologs:
            partner = by_id_b.get((ortholog.species, ortholog.protein_id))
            reverse = partner.best_ortholog(record.species) if partner else None
            if (
                partner is not None
                and supported(ortholog)
                and reverse is not None
                and reverse.protein_id == record.protein_id
                and supported(reverse)
            ):
                record.flag("cross_species_candidate")
                partner.flag("cross_species_candidate")
                out.append((record, partner))
                break
    return out


def alignment_positions(hit: dict[str, Any]) -> list[tuple[int, int]]:
    """Map the gapped DIAMOND qseq/sseq columns to zero-based coordinates."""
    qseq, sseq = hit["qseq"], hit["sseq"]
    if len(qseq) != len(sseq):
        raise ValueError("Unequal gapped alignment lengths")
    q, s = hit["qstart"] - 1, hit["sstart"] - 1
    pairs = []
    for qa, sa in zip(qseq, sseq, strict=True):
        if qa != "-" and sa != "-":
            pairs.append((q, s))
        q += qa != "-"
        s += sa != "-"
    if q != hit["qend"] or s != hit["send"]:
        raise ValueError("Alignment sequence does not match reported coordinates")
    return pairs
