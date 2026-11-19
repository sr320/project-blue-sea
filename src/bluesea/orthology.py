"""Reciprocal best hits between the two louse proteomes.

The cascade runs once per species; this is what lets the two shortlists be
intersected afterwards. A protein surviving in both genera is the smaller,
more interesting set -- conservation across *Lepeophtheirus* and *Caligus* is
independent evidence of functional essentiality that no single-species screen
can produce.

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
        best = max(eligible, key=lambda h: h["bitscore"])
        pairs[query] = {
            "subject": best["subject"],
            "identity_pct": best["identity_pct"],
            "evalue": best["evalue"],
            "reciprocal": best_reverse.get(best["subject"]) == query,
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
    return 100.0 * matched / compared


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
        if other_sequences and alignments:
            region_identity = epitope_region_identity(
                record.sequence,
                other_sequences.get(pair["subject"], ""),
                record.epitopes,
                alignments.get(protein_id),
            )
        record.orthologs.append(
            Ortholog(
                species=other_species,
                protein_id=pair["subject"],
                identity_pct=pair["identity_pct"],
                evalue=pair["evalue"],
                reciprocal_best=pair["reciprocal"],
                epitope_region_identity_pct=region_identity,
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
) -> list[tuple[ProteinRecord, ProteinRecord]]:
    """Candidates that survived the cascade in both species.

    Matched through the recorded orthology, not through name similarity --
    annotation names across these two genomes are not reliable enough to join
    on, and a wrong join here would be invisible downstream.
    """
    by_id_b = {record.protein_id: record for record in shortlist_b}
    out: list[tuple[ProteinRecord, ProteinRecord]] = []
    for record in shortlist_a:
        for ortholog in record.orthologs:
            partner = by_id_b.get(ortholog.protein_id)
            if partner is not None:
                record.flag("cross_species_candidate")
                partner.flag("cross_species_candidate")
                out.append((record, partner))
                break
    return out
