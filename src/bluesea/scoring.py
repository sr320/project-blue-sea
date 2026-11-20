"""Weighted scoring over the cascade survivors.

Seven of the eight criteria are deterministic functions of tool output. Only
the literature axis comes from an API call, and even there the score arrives
with a written justification so a reviewer can check it.

Every criterion returns a value in [0, 1]; the weighted sum is the rank key.
Criteria whose evidence is missing return a configured neutral value and set
a flag, rather than scoring zero -- absence of evidence is not evidence of
absence, and silently scoring zero would bias the screen toward the
already-studied proteins the project exists to look past.
"""

from __future__ import annotations

from typing import Any

from .records import ProteinRecord

CRITERIA = (
    "essentiality",
    "exposure",
    "epitope",
    "conservation",
    "expression_breadth",
    "host_dissimilarity",
    "literature",
    "expressibility",
)


def _scale(value: float, low: float, high: float) -> float:
    """Linear ramp from low->0 to high->1, clamped."""
    if high <= low:
        return 0.0
    return max(0.0, min(1.0, (value - low) / (high - low)))


def score_essentiality(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """RNAi phenotype, backed by orthology to essential copepod genes.

    RNAi coverage is partial, so "untested" maps to a mid-range unknown and
    raises a flag; it must not read as "dispensable".
    """
    table = cfg["phenotype_scores"]
    phenotype = record.rnai.phenotype or "untested"
    score = table.get(phenotype, table["untested"])
    if phenotype == "untested" or phenotype not in table:
        record.flag(cfg.get("untested_flag", "essentiality_unknown"))

    ortholog_evidence = record.scores.get("_ortholog_essentiality")
    if ortholog_evidence is not None:
        weight = cfg.get("ortholog_essentiality_weight", 0.5)
        score = (1 - weight) * score + weight * ortholog_evidence
    return max(0.0, min(1.0, score))


def score_exposure(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    return cfg["class_scores"].get(record.topology.exposure_class, 0.0)


def score_epitope(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """Count and quality of predicted epitopes, with count saturating.

    Saturation keeps a 1,500-residue protein from outscoring a good small one
    purely on surface area.
    """
    epitopes = record.accessible_epitopes()
    if not epitopes:
        return 0.0
    mean_score = sum(e.score for e in epitopes) / len(epitopes)
    count_term = min(1.0, len(epitopes) / cfg["count_saturation"])
    return (
        cfg["mean_score_weight"] * min(1.0, mean_score) + cfg["count_weight"] * count_term
    )


def score_conservation(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """Cross-species conservation between the two louse genera.

    Conservation supports a shared-sequence hypothesis. It does not establish
    functional essentiality, antibody binding, or protection in either species.
    """
    ortholog = record.best_ortholog()
    if ortholog is None:
        record.flag("orthology_not_detected")
        return cfg.get("default_when_missing", 0.4)

    identity = (
        ortholog.epitope_region_identity_pct
        if ortholog.epitope_region_identity_pct is not None
        else ortholog.identity_pct
    )
    if ortholog.epitope_region_identity_pct is None:
        record.flag("epitope_region_identity_unknown")

    score = _scale(identity, cfg["identity_floor_pct"], cfg["identity_ceiling_pct"])
    if ortholog.reciprocal_best:
        score = min(1.0, score + cfg["rbh_bonus"] * (1 - score))
    else:
        record.flag("ortholog_not_reciprocal_best")
    return score


def score_expression_breadth(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """Constitutive beats peak.

    A stage-restricted target gives ingested antibody a narrow window to act
    in, which is the opposite of what the durability problem needs. Early
    stages get a bonus because interrupting attachment prevents the
    infestation rather than reducing it.
    """
    if not record.expression:
        record.flag("expression_unknown")
        return 0.5

    expressed = {
        stage for stage, tpm in record.expression.items() if tpm >= cfg["min_tpm"]
    }
    denominator = len(record.expected_stages) or len(record.expression)
    if len(record.expression) < denominator:
        record.flag("expression_partial_coverage")
    # Unassayed stages contribute neutral evidence instead of perfect breadth.
    breadth = (
        len(expressed) + 0.5 * (denominator - len(record.expression))
    ) / denominator
    if expressed & set(record.early_stages or cfg.get("early_stages", [])):
        breadth += cfg.get("early_stage_bonus", 0.0)
    return min(1.0, breadth)


def score_host_dissimilarity(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """Reward being well below the hard cross-reactivity threshold.

    No salmon homolog at all is the ideal, so it scores 1.0.
    """
    if record.host_hit is None:
        return 1.0
    identity = record.host_hit.identity_pct
    return 1.0 - _scale(identity, cfg["identity_floor_pct"], cfg["identity_ceiling_pct"])


def score_literature(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    if record.literature is None:
        record.flag("literature_not_assessed")
        return cfg["default_when_missing"]
    return max(0.0, min(1.0, record.literature.score))


def score_expressibility(record: ProteinRecord, cfg: dict[str, Any]) -> float:
    """Proxies for how much trouble the protein will be at the bench."""
    low, high = cfg["ideal_length"]
    if low <= record.length <= high:
        score = 1.0
    elif record.length < low:
        score = _scale(record.length, low / 3, low)
    else:
        score = 1.0 - _scale(record.length, high, high * 3)

    score -= min(
        cfg["max_disulfide_penalty"],
        record.n_disulfide_pairs * cfg["disulfide_penalty_per_pair"],
    )
    score -= min(
        cfg["max_tm_penalty"],
        record.topology.n_transmembrane * cfg["tm_penalty_per_domain"],
    )
    score -= min(
        cfg["max_glycosylation_penalty"],
        record.n_glycosylation_sites * cfg["glycosylation_penalty_per_site"],
    )
    return max(0.0, min(1.0, score))


SCORERS = {
    "essentiality": score_essentiality,
    "exposure": score_exposure,
    "epitope": score_epitope,
    "conservation": score_conservation,
    "expression_breadth": score_expression_breadth,
    "host_dissimilarity": score_host_dissimilarity,
    "literature": score_literature,
    "expressibility": score_expressibility,
}


def score_record(
    record: ProteinRecord,
    weights_cfg: dict[str, Any],
    drop_sources: set[str] | None = None,
) -> ProteinRecord:
    """Score one record in place.

    `drop_sources` supports the leave-one-source-out robustness check: a
    dropped criterion contributes nothing and its weight leaves the
    denominator, so totals stay comparable across runs.
    """
    drop_sources = drop_sources or set()
    weights = weights_cfg["weights"]
    total = 0.0
    weight_sum = 0.0

    for criterion in CRITERIA:
        if criterion in drop_sources:
            record.scores[criterion] = float("nan")
            continue
        value = SCORERS[criterion](record, weights_cfg.get(criterion, {}))
        record.scores[criterion] = round(value, 4)
        total += value * weights[criterion]
        weight_sum += weights[criterion]

    # Normalise to the configured scale so a 100-point run and a
    # leave-one-out run can be read on the same axis.
    full_scale = sum(weights[c] for c in CRITERIA)
    record.total_score = round(total / weight_sum * full_scale, 3) if weight_sum else 0.0
    return record


def rank(
    records: list[ProteinRecord],
    weights_cfg: dict[str, Any],
    drop_sources: set[str] | None = None,
) -> list[ProteinRecord]:
    for record in records:
        score_record(record, weights_cfg, drop_sources=drop_sources)
    return sorted(
        records, key=lambda r: (-r.total_score, r.species, r.gene_id, r.protein_id)
    )


def ranking_rows(records: list[ProteinRecord]) -> list[dict[str, Any]]:
    rows = []
    for position, record in enumerate(records, start=1):
        row: dict[str, Any] = {
            "rank": position,
            "protein_id": record.protein_id,
            "gene_id": record.gene_id,
            "species": record.species,
            "total_score": record.total_score,
            "description": record.description,
            "exposure_class": record.topology.exposure_class,
            "rnai_phenotype": record.rnai.phenotype,
            "host_identity_pct": (
                record.host_hit.identity_pct if record.host_hit else ""
            ),
            "cross_species": "cross_species_candidate" in record.flags,
            "flags": ";".join(record.flags),
        }
        row.update({f"score_{c}": record.scores.get(c, "") for c in CRITERIA})
        rows.append(row)
    return rows


RANKING_COLUMNS = [
    "rank",
    "protein_id",
    "gene_id",
    "species",
    "total_score",
    *[f"score_{c}" for c in CRITERIA],
    "exposure_class",
    "rnai_phenotype",
    "host_identity_pct",
    "cross_species",
    "flags",
    "description",
]


def gene_representatives(records: list[ProteinRecord]) -> list[ProteinRecord]:
    """One best-scoring isoform per species/gene; keep all isoforms separately."""
    seen = set()
    result = []
    for record in sorted(
        records, key=lambda r: (-r.total_score, r.species, r.protein_id)
    ):
        key = (record.species, record.gene_id or record.protein_id)
        if key not in seen:
            seen.add(key)
            result.append(record)
    return result
