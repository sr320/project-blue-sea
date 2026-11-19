"""The filter cascade: ~14,000 proteins per species down to a few dozen.

Each step is a pure function of one record plus its config block, returning
(passed, reason, detail). The driver records every verdict whether or not the
protein survived, so a dropped candidate can be traced to the exact rule that
dropped it and that rule can be revisited without rerunning the predictors.

Order is cost-ascending. The exposure step does most of the cutting and is
also the least certain -- predicted secretion is not confirmed secretion --
so it is the first place to look when the survivor count comes out wrong.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from typing import Any

from .records import ProteinRecord

Verdict = tuple[bool, str, dict[str, Any]]
Step = Callable[[ProteinRecord, dict[str, Any]], Verdict]

_LOW_COMPLEXITY_WINDOW = 20


def filter_length(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Fragments carry too little sequence to judge; very long proteins are
    rarely tractable as recombinant antigens."""
    low, high = cfg["min_aa"], cfg["max_aa"]
    length = record.length
    if length < low:
        return False, f"length {length} aa below minimum {low}", {"length": length}
    if length > high:
        return False, f"length {length} aa above maximum {high}", {"length": length}
    return True, f"length {length} aa within [{low}, {high}]", {"length": length}


def filter_low_complexity(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Repeat-rich sequences produce spurious epitope and homology signal."""
    fraction = record.low_complexity_fraction or estimate_low_complexity(
        record.sequence
    )
    record.low_complexity_fraction = fraction
    limit = cfg["max_fraction"]
    passed = fraction <= limit
    return (
        passed,
        (
            f"low-complexity fraction {fraction:.2f} "
            f"{'within' if passed else 'above'} limit {limit}"
        ),
        {"low_complexity_fraction": round(fraction, 4)},
    )


def estimate_low_complexity(sequence: str) -> float:
    """Fraction of residues in windows dominated by one or two residue types.

    A cheap stand-in for SEG that needs no external call; the threshold is
    loose enough that the approximation does not decide borderline cases.
    """
    if len(sequence) < _LOW_COMPLEXITY_WINDOW:
        return 0.0
    flagged = 0
    for start in range(len(sequence) - _LOW_COMPLEXITY_WINDOW + 1):
        window = sequence[start : start + _LOW_COMPLEXITY_WINDOW]
        top_two = sum(count for _, count in Counter(window).most_common(2))
        if top_two / _LOW_COMPLEXITY_WINDOW >= 0.7:
            flagged += 1
    return flagged / (len(sequence) - _LOW_COMPLEXITY_WINDOW + 1)


def filter_exposure(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Keep proteins antibody could actually reach.

    Antibody arrives with the host material the louse ingests, so a target
    must be secreted, surface-anchored, or present enough extracellular or
    gut-lumen-facing sequence to be bound. This cuts roughly eight proteins
    in nine and is the single largest source of false negatives in the run.
    """
    topology = record.topology
    exposure_class = topology.exposure_class
    detail = {
        "exposure_class": exposure_class,
        "signalp_prob": round(topology.signalp_prob, 3),
        "n_transmembrane": topology.n_transmembrane,
        "extracellular_aa": topology.extracellular_aa,
    }

    if topology.secretome_confirmed:
        return True, "confirmed in secretome proteomics", detail

    if topology.has_signal_peptide and topology.signalp_prob >= cfg["signalp_min_prob"]:
        return True, f"signal peptide (p={topology.signalp_prob:.2f})", detail

    if topology.gpi_anchored and cfg.get("allow_gpi_anchored", True):
        return True, "GPI-anchored surface protein", detail

    if (
        topology.n_transmembrane >= 1
        and topology.extracellular_aa >= cfg["min_extracellular_aa"]
    ):
        return (
            True,
            f"membrane protein with {topology.extracellular_aa} aa extracellular",
            detail,
        )

    if topology.nonclassical_secretion and cfg.get("accept_nonclassical", True):
        record.flag("nonclassical_secretion_only")
        return True, "predicted non-classical secretion (weak evidence)", detail

    return False, "no predicted route to host-antibody contact", detail


def filter_host_similarity(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Hard exclusion on salmon cross-reactivity.

    A vaccine antigen too close to a host protein risks raising autoantibodies
    in the fish being protected. No score compensates for this, so it is a
    filter rather than a penalty.
    """
    hit = record.host_hit
    if hit is None:
        return True, "no salmon homolog detected", {"host_identity_pct": None}

    detail = {
        "host_identity_pct": hit.identity_pct,
        "host_evalue": hit.evalue,
        "host_coverage": round(hit.query_coverage, 3),
        "host_accession": hit.accession,
    }
    significant = (
        hit.evalue <= cfg["max_evalue"]
        and hit.query_coverage >= cfg["min_query_coverage"]
    )
    if significant and hit.identity_pct > cfg["max_identity_pct"]:
        return (
            False,
            (
                f"{hit.identity_pct:.1f}% identity to salmon {hit.accession} "
                f"exceeds {cfg['max_identity_pct']}%"
            ),
            detail,
        )
    return True, f"closest salmon homolog at {hit.identity_pct:.1f}% identity", detail


def filter_epitopes(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Require predicted B-cell epitopes on the accessible surface.

    Weak evidence by construction: the predictors are calibrated on human and
    mouse immunology, and salmon IgM/IgT recognition differs. The threshold is
    set low on purpose so this step removes only the clearly epitope-free.
    """
    good = [
        epitope
        for epitope in record.epitopes
        if epitope.score >= cfg["min_epitope_score"]
        and epitope.length >= cfg.get("window", 15) // 2
    ]
    detail = {
        "n_epitopes": len(good),
        "mean_epitope_score": (
            round(sum(e.score for e in good) / len(good), 3) if good else 0.0
        ),
    }
    minimum = cfg["min_linear_epitopes"]
    if len(good) < minimum:
        return False, f"{len(good)} qualifying epitopes, need {minimum}", detail
    record.flag("epitope_prediction_mammal_calibrated")
    return True, f"{len(good)} qualifying epitopes", detail


def filter_expression(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Drop targets not expressed in the parasitic stages.

    Stages missing from the matrix were not assayed and are not counted
    against the protein; the filter asks only that the assayed stages show
    expression somewhere.
    """
    if not record.expression:
        record.flag("expression_unknown")
        return True, "no expression data available; not penalised", {}

    expressed = [
        stage for stage, tpm in record.expression.items() if tpm >= cfg["min_tpm"]
    ]
    detail = {
        "stages_expressed": sorted(expressed),
        "n_stages_expressed": len(expressed),
        "max_tpm": round(max(record.expression.values()), 2),
    }
    minimum = cfg["min_stages_expressed"]
    if len(expressed) < minimum:
        return (
            False,
            f"expressed in {len(expressed)} stages, need {minimum}",
            detail,
        )
    return True, f"expressed in {len(expressed)} stages", detail


def annotate_expressibility(record: ProteinRecord, cfg: dict[str, Any]) -> Verdict:
    """Never excludes. Flags proteins likely to fight back in expression."""
    if record.topology.n_transmembrane > cfg["warn_above"]:
        record.flag("difficult_recombinant_expression")
    return (
        True,
        f"{record.topology.n_transmembrane} transmembrane domains",
        {"n_transmembrane": record.topology.n_transmembrane},
    )


CASCADE: list[tuple[str, str, Step]] = [
    ("length", "length", filter_length),
    ("low_complexity", "low_complexity", filter_low_complexity),
    ("exposure", "exposure", filter_exposure),
    ("host_similarity", "host_similarity", filter_host_similarity),
    ("epitope", "epitope", filter_epitopes),
    ("expression", "expression", filter_expression),
    ("transmembrane", "transmembrane", annotate_expressibility),
]


def run_cascade(
    record: ProteinRecord, cascade_cfg: dict[str, Any], stop_on_fail: bool = True
) -> ProteinRecord:
    """Apply every step, recording a verdict for each.

    With `stop_on_fail` the remaining steps are skipped once one fails, which
    is what the production run wants: the expensive predictors never see a
    protein that is already excluded. Set it False to profile the cascade and
    see how much each step would have cut independently.
    """
    for step_name, cfg_key, step in CASCADE:
        if stop_on_fail and not record.survived:
            break
        passed, reason, detail = step(record, cascade_cfg[cfg_key])
        record.record_verdict(step_name, passed, reason, **detail)
    return record


def cascade_summary(records: list[ProteinRecord]) -> list[dict[str, Any]]:
    """Per-step attrition table: how many entered, how many each step dropped.

    Written beside every run. A step that drops nothing is either misconfigured
    or redundant, and both are worth seeing before reading the ranked list.
    """
    rows: list[dict[str, Any]] = []
    remaining = len(records)
    for step_name, _, _ in CASCADE:
        dropped = sum(1 for r in records if r.dropped_at == step_name)
        rows.append(
            {
                "step": step_name,
                "entered": remaining,
                "dropped": dropped,
                "survived": remaining - dropped,
                "drop_rate": round(dropped / remaining, 4) if remaining else 0.0,
            }
        )
        remaining -= dropped
    return rows
