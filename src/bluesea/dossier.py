"""One-page candidate dossiers.

Written for someone who was not involved in building the pipeline: what the
protein is, why it survived, every score with its justification, and the full
provenance trail. The uncertainty is on the page rather than in a footnote,
because a dossier that reads as more confident than the evidence is worse
than no dossier.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from .records import ProteinRecord
from .scoring import CRITERIA

FLAG_NOTES = {
    "antibody_access_unconfirmed": (
        "Secretion/topology is predicted. Tissue localization and contact with "
        "host antibody have not been established; O can include organelle lumens."
    ),
    "short_host_similar_patch": (
        "A short high-identity host alignment needs epitope-level review. The "
        "whole-query coverage heuristic did not exclude this patch."
    ),
    "weak_predicted_epitope_support": (
        "Few predicted accessible linear patches. Prediction is advisory; "
        "confirm binding experimentally with salmon antibodies."
    ),
    "epitope_accessibility_unknown": (
        "Residue topology is missing; epitope access is unknown."
    ),
    "expression_partial_coverage": "Some species-valid life stages were not assayed.",
    "expression_insufficient_stages": "Too few stages were assayed to test breadth.",
    "orthology_not_detected": "No qualifying cross-species match was detected.",
    "literature_sources_missing": (
        "No curated source packet; no model assessment was made."
    ),
    "literature_response_unsupported": (
        "The model returned unsupported citations or scores."
    ),
    "essentiality_unknown": (
        "No RNAi knockdown has been reported. RNAi coverage of this genome is "
        "partial, so this is absence of evidence, not evidence that the gene "
        "is dispensable. Knockdown is the cheapest way to settle it."
    ),
    "expression_unknown": "No stage-resolved expression data covers this gene.",
    "nonclassical_secretion_only": (
        "Secretion is predicted by a non-classical route, which is the weakest "
        "exposure evidence the cascade accepts."
    ),
    "epitope_prediction_mammal_calibrated": (
        "B-cell epitope predictors are trained largely on human and mouse data. "
        "Salmon IgM/IgT recognition differs, so treat the epitope map as weak "
        "evidence and test binding using salmon antibodies."
    ),
    "epitopes_not_conserved_across_species": (
        "The ortholog is conserved overall but the predicted epitopes are not, "
        "so a single construct is unlikely to cover both louse species."
    ),
    "epitope_region_identity_unknown": (
        "Epitope-region identity was not computed; only whole-protein identity "
        "supports the conservation score."
    ),
    "ortholog_not_reciprocal_best": (
        "The cross-species match is not a reciprocal best hit, so it may be a "
        "paralog rather than a true ortholog."
    ),
    "difficult_recombinant_expression": (
        "Multiple transmembrane domains; expect difficulty producing the full "
        "protein recombinantly. Consider expressing an extracellular domain."
    ),
    "literature_not_assessed": "The literature axis did not run for this protein.",
    "literature_score_uncited": (
        "The literature score came back without citations and should be "
        "treated as unsupported."
    ),
    "literature_skipped_budget": "Literature triage halted at its spend ceiling.",
    "cross_species_candidate": (
        "Reciprocal matches survive in both species with measured epitope "
        "conservation. Shared efficacy and essentiality remain untested."
    ),
    "single_source_dependent": (
        "Top-tier rank depends on one evidence source. Fragile: a change in "
        "that source would move this candidate substantially."
    ),
}


def construct_boundary(record: ProteinRecord) -> tuple[int, int] | None:
    """Suggest only an epitope-bearing, contiguous predicted outside domain."""
    topology = record.topology.topology_string
    if not topology or len(topology) != record.length:
        return None
    cleavage = record.topology.signal_peptide_end
    if record.topology.has_signal_peptide and cleavage is None and "S" not in topology:
        return None
    regions = []
    start = None
    for index, label in enumerate(topology + "I", start=1):
        if label == "O" and start is None:
            start = index
        if label != "O" and start is not None:
            end = index - 1
            start = max(start, (cleavage or 0) + 1)
            if any(
                start <= e.start <= e.end <= end for e in record.accessible_epitopes()
            ):
                regions.append((start, end))
            start = None
    return max(regions, key=lambda r: r[1] - r[0], default=None)


def render(
    record: ProteinRecord,
    weights: dict[str, int],
    provenance: dict[str, Any],
) -> str:
    lines: list[str] = []
    add = lines.append

    add(f"# {record.protein_id} — {record.description or 'hypothetical protein'}")
    add("")
    add(f"**Species** {record.species}  ")
    add(f"**Gene** {record.gene_id}  ")
    add(f"**Length** {record.length} aa  ")
    add(f"**Total score** {record.total_score:.1f} / {sum(weights.values())}")
    add("")

    add("## Scores")
    add("")
    add("| Criterion | Weight | Score | Weighted |")
    add("| --- | --- | --- | --- |")
    for criterion in CRITERIA:
        value = record.scores.get(criterion)
        # NaN marks a criterion dropped by a leave-one-source-out run.
        if value is None or math.isnan(value):
            add(f"| {criterion} | {weights[criterion]} | dropped | — |")
            continue
        add(
            f"| {criterion} | {weights[criterion]} | {value:.2f} | "
            f"{value * weights[criterion]:.1f} |"
        )
    add("")

    add("## Why it survived")
    add("")
    for verdict in record.verdicts:
        mark = "pass" if verdict["passed"] else "FAIL"
        add(f"- **{verdict['step']}** ({mark}): {verdict['reason']}")
    add("")

    add("## Localisation and topology")
    add("")
    topology = record.topology
    add(f"- Exposure class: `{topology.exposure_class}`")
    add(
        f"- Signal peptide: {'yes' if topology.has_signal_peptide else 'no'} "
        f"(p = {topology.signalp_prob:.2f})"
    )
    add(f"- Transmembrane domains: {topology.n_transmembrane}")
    add(f"- Extracellular residues: {topology.extracellular_aa}")
    add(
        "- Secretome proteomics: "
        + ("confirmed" if topology.secretome_confirmed else "not confirmed")
    )
    add("")

    add("## Epitope map")
    add("")
    if record.epitopes:
        add("| Start | End | Length | Score | Kind |")
        add("| --- | --- | --- | --- | --- |")
        for epitope in sorted(record.epitopes, key=lambda e: e.start):
            add(
                f"| {epitope.start} | {epitope.end} | {epitope.length} | "
                f"{epitope.score:.2f} | {epitope.kind} |"
            )
    else:
        add("No epitopes predicted above threshold.")
    add("")
    boundary = construct_boundary(record)
    if boundary:
        start, end = boundary
        add(
            f"**Suggested construct**: residues {start}–{end} "
            "(predicted outside domain; validate localization and folding)."
        )
    else:
        add("**Suggested construct**: withheld; no supported epitope-bearing domain.")
    add("")

    add("## Expression across life stages")
    add("")
    if record.expression:
        add("| Stage | TPM |")
        add("| --- | --- |")
        for stage, tpm in sorted(record.expression.items()):
            add(f"| {stage} | {tpm:.1f} |")
    else:
        add("No stage-resolved expression data for this gene.")
    add("")

    add("## Essentiality")
    add("")
    add(f"- RNAi phenotype: **{record.rnai.phenotype}**")
    if record.rnai.source:
        add(f"- Source: {record.rnai.source}")
    if record.rnai.note:
        add(f"- Note: {record.rnai.note}")
    add("")

    add("## Host cross-reactivity")
    add("")
    if record.host_hit:
        hit = record.host_hit
        add(
            f"- Closest *Salmo salar* protein: `{hit.accession}` at "
            f"{hit.identity_pct:.1f}% identity "
            f"(E = {hit.evalue:.1e}, coverage {hit.query_coverage:.0%})"
        )
    else:
        add("- No salmon homolog detected.")
    add("")

    add("## Cross-species orthology")
    add("")
    if record.orthologs:
        for ortholog in record.orthologs:
            region = (
                f", epitope region {ortholog.epitope_region_identity_pct:.1f}%"
                if ortholog.epitope_region_identity_pct is not None
                else ""
            )
            add(
                f"- {ortholog.species} `{ortholog.protein_id}`: "
                f"{ortholog.identity_pct:.1f}% identity"
                + (", reciprocal best hit" if ortholog.reciprocal_best else "")
                + region
            )
    else:
        add("- No ortholog found in the other louse species.")
    add("")

    add("## Literature")
    add("")
    if record.literature:
        evidence = record.literature
        blinded = " (blinded run)" if evidence.blinded else ""
        add(f"**Score** {evidence.score:.2f}{blinded}")
        add("")
        add(evidence.justification or "_No justification returned._")
        add("")
        if evidence.citations:
            add("Citations: " + ", ".join(f"`{c}`" for c in evidence.citations))
        else:
            add("_No citations returned; treat this score as unsupported._")
    else:
        add("_Literature triage did not run for this protein._")
    add("")

    if record.flags:
        add("## Caveats")
        add("")
        for flag in record.flags:
            add(f"- **{flag}** — {FLAG_NOTES.get(flag, 'see pipeline log')}")
        add("")

    add("## Provenance")
    add("")
    for key, value in sorted(provenance.items()):
        add(f"- {key}: `{value}`")
    add("")
    add(
        "_This dossier is a hypothesis for the bench. The pipeline does not "
        "predict efficacy, and nothing here substitutes for a challenge trial._"
    )
    return "\n".join(lines) + "\n"


def write_dossiers(
    records: list[ProteinRecord],
    out_dir: str | Path,
    weights: dict[str, int],
    provenance: dict[str, Any],
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for position, record in enumerate(records, start=1):
        path = out_dir / f"{position:03d}_{record.protein_id}.md"
        path.write_text(render(record, weights, provenance))
        written.append(path)
    return written
