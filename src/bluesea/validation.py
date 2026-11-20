"""Retrospective benchmark. Build it before the first full run, not after.

A ranked list is worthless until someone has tested whether the ranking means
anything. Three checks, and a gate:

* positive controls -- antigens with measurable effects in published trials
  must surface in the top tier on deterministic evidence alone, which is why
  the literature stage is blinded for this run;
* negative controls -- salmon proteins must be rejected by cross-reactivity,
  intracellular housekeeping proteins by exposure, and randomly drawn louse
  proteins should distribute around the middle;
* robustness -- re-rank under perturbed weights and under leave-one-source-out,
  and label candidates that depend on a single evidence source as fragile
  rather than ranking them quietly alongside the robust ones.

The gate is the six-week decision: if the known-working antigens land mid-pack,
the rubric is wrong and gets rebuilt before anyone reads a novel candidate.
"""

from __future__ import annotations

import copy
import random
import statistics
from collections.abc import Mapping
from typing import Any

from .records import ProteinRecord
from .scoring import CRITERIA, gene_representatives, rank

# Which criteria each evidence source feeds, for leave-one-source-out.
SOURCE_CRITERIA = {
    "rnai": {"essentiality"},
    "expression": {"expression_breadth"},
    "literature": {"literature"},
    "epitope": {"epitope"},
    "orthology": {"conservation"},
}


def top_tier(ranked: list[ProteinRecord], fraction: float) -> list[ProteinRecord]:
    cutoff = max(1, round(len(ranked) * fraction))
    return ranked[:cutoff]


def as_groups(controls: Any) -> dict[str, set[str]]:
    """Normalise controls to {label: accessions}.

    A control is a gene; the ranking holds proteins, and RefSeq carries one
    entry per isoform. Scoring a control against a single accession means that
    when a sibling isoform reaches the top tier instead, a correct result is
    recorded as a miss. So a control may name every accession that represents
    it, and a plain set of ids still works -- each id standing for itself.
    """
    if isinstance(controls, Mapping):
        return {label: set(ids) for label, ids in controls.items()}
    return {ident: {ident} for ident in controls}


def positive_control_report(
    ranked: list[ProteinRecord],
    controls: Any,
    fraction: float,
    minimum_recovered: float,
) -> dict[str, Any]:
    """Where did the known-working antigens land?

    A control counts as recovered when any accession representing it reaches
    the top tier, and its reported position is the best its isoforms achieved.
    """
    groups = as_groups(controls)
    tier = {record.protein_id for record in top_tier(ranked, fraction)}
    position_of = {record.protein_id: i for i, record in enumerate(ranked, start=1)}

    positions: dict[str, int] = {}
    matched: dict[str, str] = {}
    in_tier: set[str] = set()
    for label, accessions in groups.items():
        ranked_here = {a: position_of[a] for a in accessions if a in position_of}
        if ranked_here:
            best = min(ranked_here, key=lambda a: ranked_here[a])
            positions[label] = ranked_here[best]
            matched[label] = best
        if accessions & tier:
            in_tier.add(label)

    missing = set(groups) - set(positions)
    recovered = len(in_tier) / len(groups) if groups else 0.0

    return {
        "n_controls": len(groups),
        "n_in_ranking": len(positions),
        "n_in_top_tier": len(in_tier),
        "recovered_fraction": round(recovered, 3),
        "threshold": minimum_recovered,
        "passed": recovered >= minimum_recovered,
        "positions": dict(sorted(positions.items(), key=lambda kv: kv[1])),
        # Which isoform actually carried the control, so a result that depends
        # on a sibling accession is visible rather than implicit.
        "matched_accession": {k: matched[k] for k in sorted(matched)},
        # Controls absent from the ranking were dropped by the cascade; that
        # is a louder failure than ranking them low, so they are listed apart.
        "dropped_by_cascade": sorted(missing),
    }


def as_negative_groups(expected: Mapping[str, Any]) -> dict[str, tuple[str, set[str]]]:
    """Normalise to {label: (expected_step, accessions)}.

    A value that is a bare step keeps the old one-accession-per-control
    meaning; a (step, accessions) pair names every accession representing the
    control, as `as_groups` does for positives.
    """
    groups: dict[str, tuple[str, set[str]]] = {}
    for label, value in expected.items():
        if isinstance(value, str):
            groups[label] = (value, {label})
        else:
            step, accessions = value
            groups[label] = (step, set(accessions))
    return groups


def negative_control_report(
    records: list[ProteinRecord],
    expected: Mapping[str, Any],
) -> dict[str, Any]:
    """Each negative control must be dropped, and by the right step.

    Being rejected for the wrong reason counts as a failure: a salmon protein
    that fails the exposure filter tells us nothing about whether the
    cross-reactivity filter works.

    Where a control names several accessions, every one that appears in the
    records must agree. Isoforms of one gene diverging on the step that
    rejected them is itself worth seeing, so it fails rather than being
    resolved by majority.
    """
    by_id = {record.protein_id: record for record in records}
    rows = []
    for label, (want, accessions) in as_negative_groups(expected).items():
        present = [by_id[a] for a in sorted(accessions) if a in by_id]
        if not present:
            rows.append(
                {
                    "protein_id": label,
                    "n_accessions": 0,
                    "expected_drop_step": want,
                    "actual_drop_step": "missing",
                    "passed": False,
                }
            )
            continue
        got = sorted({record.dropped_at or "survived" for record in present})
        rows.append(
            {
                "protein_id": label,
                "n_accessions": len(present),
                "expected_drop_step": want,
                "actual_drop_step": got[0] if len(got) == 1 else got,
                "passed": got == [want],
            }
        )
    return {
        "n": len(rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "passed": all(row["passed"] for row in rows) if rows else False,
        "details": rows,
    }


def random_control_report(
    ranked: list[ProteinRecord], random_ids: set[str]
) -> dict[str, Any]:
    """Randomly drawn proteins should sit around the middle of the ranking.

    A pipeline that ranks random proteins highly is measuring nothing, so this
    checks the median percentile rather than any individual position.
    """
    total = len(ranked)
    percentiles = [
        index / total
        for index, record in enumerate(ranked, start=1)
        if record.protein_id in random_ids
    ]
    if not percentiles:
        return {"n": 0, "passed": False, "reason": "no random controls in ranking"}
    median = statistics.median(percentiles)
    return {
        "n": len(percentiles),
        "median_percentile": round(median, 3),
        "fraction_in_top_decile": round(
            sum(1 for p in percentiles if p <= 0.1) / len(percentiles), 3
        ),
        # Generous band: the point is to catch a ranking that loves everything,
        # not to demand the random set be exactly central.
        "passed": 0.3 <= median <= 0.7,
    }


def weight_perturbation(
    records: list[ProteinRecord],
    weights_cfg: dict[str, Any],
    pct: float,
    n: int,
    fraction: float,
    seed: int,
) -> dict[str, Any]:
    """How often does each candidate hold its top-tier place under jitter?

    Weights are a judgment, not a measurement. If the top ten reshuffles when
    essentiality moves a few points, that instability is itself a finding.
    """
    if n < 1:
        raise ValueError("n_perturbations must be positive")
    rng = random.Random(seed)
    baseline = rank(copy.deepcopy(records), weights_cfg)
    baseline_tier = [r.protein_id for r in top_tier(baseline, fraction)]
    appearances = dict.fromkeys(baseline_tier, 0)

    for _ in range(n):
        perturbed = {
            "weights": {
                criterion: max(
                    0.0,
                    weight * (1 + rng.uniform(-pct, pct) / 100),
                )
                for criterion, weight in weights_cfg["weights"].items()
            },
            **{k: v for k, v in weights_cfg.items() if k != "weights"},
        }
        ranked = rank(copy.deepcopy(records), perturbed)
        for record in top_tier(ranked, fraction):
            if record.protein_id in appearances:
                appearances[record.protein_id] += 1

    stability = {pid: round(count / n, 3) for pid, count in appearances.items()}
    return {
        "n_perturbations": n,
        "perturbation_pct": pct,
        "baseline_top_tier": baseline_tier,
        "stability": stability,
        "mean_stability": (
            round(statistics.mean(stability.values()), 3) if stability else 0.0
        ),
    }


def leave_one_source_out(
    records: list[ProteinRecord],
    weights_cfg: dict[str, Any],
    sources: list[str],
    fraction: float,
) -> dict[str, Any]:
    """Which candidates survive dropping each evidence source in turn?

    A candidate that holds its place across every variant is one to take to
    the bench. One that depends entirely on a single source is interesting but
    fragile, and gets labelled that way in its dossier.
    """
    baseline = rank(copy.deepcopy(records), weights_cfg)
    baseline_tier = {r.protein_id for r in top_tier(baseline, fraction)}

    per_source: dict[str, list[str]] = {}
    for source in sources:
        dropped = SOURCE_CRITERIA.get(source, set())
        ranked = rank(copy.deepcopy(records), weights_cfg, drop_sources=dropped)
        per_source[source] = [r.protein_id for r in top_tier(ranked, fraction)]

    robust = set(baseline_tier)
    for tier in per_source.values():
        robust &= set(tier)

    fragile = sorted(baseline_tier - robust)
    by_id = {record.protein_id: record for record in records}
    for protein_id in fragile:
        by_id[protein_id].flag("single_source_dependent")

    return {
        "sources": sources,
        "baseline_top_tier": sorted(baseline_tier),
        "per_source_top_tier": per_source,
        "robust": sorted(robust),
        "fragile": fragile,
    }


def run_benchmark(
    records: list[ProteinRecord],
    weights_cfg: dict[str, Any],
    validation_cfg: dict[str, Any],
    positive_ids: Any,
    negative_expected: Mapping[str, Any],
    random_ids: set[str],
) -> dict[str, Any]:
    """The full benchmark, with the gate verdict at the top level.

    `passed` here is the six-week gate: a False means rebuild the rubric
    before running the full screen, not tune it until this turns True.
    """
    fraction = validation_cfg["top_tier_fraction"]
    # Deterministic benchmark: remove literature evidence rather than asking
    # a model to forget known trial results. Production has separate sensitivity.
    records = copy.deepcopy(records)
    for record in records:
        record.literature = None
    survivors = gene_representatives(
        rank([r for r in records if r.survived and not r.benchmark_only], weights_cfg)
    )
    ranked = rank(list(survivors), weights_cfg)

    positives = positive_control_report(
        ranked,
        positive_ids,
        fraction,
        validation_cfg["min_positive_controls_in_top_tier"],
    )
    negatives = negative_control_report(records, negative_expected)
    randoms = random_control_report(ranked, random_ids)
    robustness_cfg = validation_cfg["robustness"]
    perturbation = weight_perturbation(
        survivors,
        weights_cfg,
        robustness_cfg["weight_perturbation_pct"],
        robustness_cfg["n_perturbations"],
        fraction,
        validation_cfg["random_seed"],
    )
    loso = leave_one_source_out(
        survivors, weights_cfg, robustness_cfg["leave_one_source_out"], fraction
    )

    sensitivity = sensitivity_verdict(perturbation, loso, robustness_cfg)
    return {
        "gate_passed": positives["passed"]
        and negatives["passed"]
        and sensitivity["passed"],
        "sensitivity_validation": sensitivity,
        "positive_controls": positives,
        "negative_controls": negatives,
        "random_controls": randoms,
        "weight_perturbation": perturbation,
        "leave_one_source_out": loso,
        "criteria": list(CRITERIA),
        "literature_policy": "excluded; constant missing-evidence score for all records",
    }


def sensitivity_report(records, weights_cfg, validation_cfg):
    """Sensitivity of the actual ranking, preserving baseline scores in place."""
    cfg = validation_cfg["robustness"]
    for record in records:
        record.flags = [f for f in record.flags if f != "single_source_dependent"]
    report = {
        "population": "one representative per species/gene",
        "weight_perturbation": weight_perturbation(
            records,
            weights_cfg,
            cfg["weight_perturbation_pct"],
            cfg["n_perturbations"],
            validation_cfg["top_tier_fraction"],
            validation_cfg["random_seed"],
        ),
        "leave_one_source_out": leave_one_source_out(
            records,
            weights_cfg,
            cfg["leave_one_source_out"],
            validation_cfg["top_tier_fraction"],
        ),
    }
    report["validation"] = sensitivity_verdict(
        report["weight_perturbation"], report["leave_one_source_out"], cfg
    )
    return report


def sensitivity_verdict(perturbation, loso, cfg):
    size = len(loso["baseline_top_tier"])
    robust_fraction = len(loso["robust"]) / size if size else 0.0
    minimum_stability = cfg.get("min_mean_stability", 0.8)
    minimum_robust = cfg.get("min_robust_fraction", 0.5)
    return {
        "passed": bool(size)
        and perturbation["mean_stability"] >= minimum_stability
        and robust_fraction >= minimum_robust,
        "mean_stability": perturbation["mean_stability"],
        "min_mean_stability": minimum_stability,
        "robust_top_tier_fraction": robust_fraction,
        "min_robust_fraction": minimum_robust,
        "scope": "Weight and score-source sensitivity among survivors; "
        "does not test filter false negatives.",
    }


def construct_report(records, weights_cfg, validation_cfg, registry):
    """Trace experimental constructs without assigning efficacy to components.

    Recovery is diagnostic until a complete independently sourced benchmark
    is curated. A multi-component construct requires every louse component.
    """
    from pathlib import Path

    from .io import read_fasta
    from .provenance import file_digest

    population = copy.deepcopy(records)
    for record in population:
        record.literature = None
    ranked = gene_representatives(
        rank([r for r in population if r.survived and not r.benchmark_only], weights_cfg)
    )
    tier = {r.protein_id for r in top_tier(ranked, validation_cfg["top_tier_fraction"])}
    positions = {r.protein_id: i for i, r in enumerate(ranked, 1)}
    by_id = {r.protein_id: r for r in population}
    rows = []
    studies = set()
    for construct in registry["constructs"]:
        components = []
        for component in construct["louse_components"]:
            ids = set(component["mapped_accessions"])
            present = [by_id[pid] for pid in ids if pid in by_id]
            components.append(
                {
                    "name": component["name"],
                    "source_accession": component["source_accession"],
                    "source_region": component["source_region"],
                    "mapping_status": component["mapping_status"],
                    "present": bool(present),
                    "in_top_tier": bool(ids & tier),
                    "best_rank": min(
                        (positions[p] for p in ids if p in positions), default=None
                    ),
                    "drop_steps": sorted({r.dropped_at or "survived" for r in present}),
                }
            )
        sequence_file = construct.get("sequence_file")
        sequence_verified = bool(
            sequence_file
            and Path(sequence_file).is_file()
            and file_digest(sequence_file) == construct.get("sequence_sha256")
            and len(list(read_fasta(sequence_file))) == 1
        )
        eligible = (
            construct["curation_status"] == "verified"
            and all(c["mapping_status"] == "verified" for c in components)
            and sequence_verified
        )
        if eligible:
            studies.add(construct["study_doi"])
        rows.append(
            {
                **construct,
                "component_recovery": components,
                "all_components_in_top_tier": bool(components)
                and all(c["in_top_tier"] for c in components),
                "eligible_for_gate": eligible,
            }
        )
    recovered_by_study = {}
    for row in rows:
        recovered_by_study.setdefault(row["study_doi"], []).append(
            row["all_components_in_top_tier"]
        )
    # A study counts once regardless of how many constructs share its controls.
    recovered = (
        sum(all(values) for values in recovered_by_study.values())
        / len(recovered_by_study)
        if recovered_by_study
        else 0.0
    )
    return {
        "recovered_fraction": recovered,
        "recovery_passed": recovered
        >= validation_cfg["min_positive_controls_in_top_tier"],
        "ready": bool(rows)
        and all(r["eligible_for_gate"] for r in rows)
        and len(studies) >= validation_cfg.get("min_independent_studies", 2),
        "independent_studies": len(studies),
        "constructs": rows,
        "interpretation": "Component recovery does not establish full-protein efficacy. "
        "Pending sequence/mapping curation blocks readiness.",
    }
