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

import random
import statistics
from typing import Any

from .records import ProteinRecord
from .scoring import CRITERIA, rank

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


def positive_control_report(
    ranked: list[ProteinRecord],
    control_ids: set[str],
    fraction: float,
    minimum_recovered: float,
) -> dict[str, Any]:
    """Where did the known-working antigens land?"""
    tier = {record.protein_id for record in top_tier(ranked, fraction)}
    positions = {
        record.protein_id: index
        for index, record in enumerate(ranked, start=1)
        if record.protein_id in control_ids
    }
    found = set(positions)
    missing = control_ids - found
    recovered = len(found & tier) / len(control_ids) if control_ids else 0.0

    return {
        "n_controls": len(control_ids),
        "n_in_ranking": len(found),
        "n_in_top_tier": len(found & tier),
        "recovered_fraction": round(recovered, 3),
        "threshold": minimum_recovered,
        "passed": recovered >= minimum_recovered,
        "positions": dict(sorted(positions.items(), key=lambda kv: kv[1])),
        # Controls absent from the ranking were dropped by the cascade; that
        # is a louder failure than ranking them low, so they are listed apart.
        "dropped_by_cascade": sorted(missing),
    }


def negative_control_report(
    records: list[ProteinRecord],
    expected: dict[str, str],
) -> dict[str, Any]:
    """Each negative control must be dropped, and by the right step.

    Being rejected for the wrong reason counts as a failure: a salmon protein
    that fails the exposure filter tells us nothing about whether the
    cross-reactivity filter works.
    """
    rows = []
    for record in records:
        want = expected.get(record.protein_id)
        if want is None:
            continue
        got = record.dropped_at
        rows.append(
            {
                "protein_id": record.protein_id,
                "expected_drop_step": want,
                "actual_drop_step": got or "survived",
                "passed": got == want,
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
    rng = random.Random(seed)
    baseline = rank(list(records), weights_cfg)
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
        ranked = rank(list(records), perturbed)
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
    baseline = rank(list(records), weights_cfg)
    baseline_tier = {r.protein_id for r in top_tier(baseline, fraction)}

    per_source: dict[str, list[str]] = {}
    for source in sources:
        dropped = SOURCE_CRITERIA.get(source, set())
        ranked = rank(list(records), weights_cfg, drop_sources=dropped)
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
    positive_ids: set[str],
    negative_expected: dict[str, str],
    random_ids: set[str],
) -> dict[str, Any]:
    """The full benchmark, with the gate verdict at the top level.

    `passed` here is the six-week gate: a False means rebuild the rubric
    before running the full screen, not tune it until this turns True.
    """
    fraction = validation_cfg["top_tier_fraction"]
    survivors = [record for record in records if record.survived]
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

    return {
        "gate_passed": positives["passed"] and negatives["passed"],
        "positive_controls": positives,
        "negative_controls": negatives,
        "random_controls": randoms,
        "weight_perturbation": perturbation,
        "leave_one_source_out": loso,
        "criteria": list(CRITERIA),
    }
