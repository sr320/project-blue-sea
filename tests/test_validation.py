from bluesea.cascade import run_cascade
from bluesea.scoring import rank
from bluesea.symbiont import MAG, gate
from bluesea.validation import (
    leave_one_source_out,
    negative_control_report,
    positive_control_report,
    random_control_report,
    run_benchmark,
    weight_perturbation,
)

from .conftest import make_record


def build_population():
    """Eight strong candidates, twenty mediocre ones. The strong set stands in
    for the published antigens the rubric must recover."""
    strong = [
        make_record(
            f"POS_{i}",
            seed=100 + i,
            rnai="lethal",
            host_identity=18.0,
            ortholog_identity=85.0,
            n_epitopes=7,
        )
        for i in range(8)
    ]
    weak = [
        make_record(
            f"BG_{i}",
            seed=200 + i,
            rnai="no_phenotype",
            host_identity=38.0,
            n_epitopes=2,
            epitope_score=0.55,
            stages={"adult_female": 4.0, "adult_male": 0.2},
        )
        for i in range(20)
    ]
    return strong + weak


def test_rubric_recovers_the_positive_controls(weights, config):
    records = build_population()
    ranked = rank(records, weights)
    report = positive_control_report(
        ranked,
        {f"POS_{i}" for i in range(8)},
        fraction=0.3,
        minimum_recovered=0.6,
    )
    assert report["passed"]
    assert report["dropped_by_cascade"] == []


def test_controls_dropped_by_the_cascade_are_reported_separately(weights):
    """A control missing from the ranking was excluded outright, which is a
    louder failure than ranking it low and must not be hidden."""
    ranked = rank(build_population(), weights)
    report = positive_control_report(
        ranked, {"POS_0", "NEVER_SEEN"}, fraction=0.3, minimum_recovered=0.6
    )
    assert report["dropped_by_cascade"] == ["NEVER_SEEN"]


def test_negative_controls_must_fail_at_the_right_step(config):
    host_like = make_record("SSAL_ACTB", host_identity=95.0)
    housekeeping = make_record("LSAL_RPL7", signal_peptide=False, seed=5)
    for record in (host_like, housekeeping):
        run_cascade(record, config["cascade"])

    report = negative_control_report(
        [host_like, housekeeping],
        {"SSAL_ACTB": "host_similarity", "LSAL_RPL7": "exposure"},
    )
    assert report["passed"]

    # Right verdict, wrong reason: still a failure.
    wrong = negative_control_report(
        [host_like], {"SSAL_ACTB": "exposure"}
    )
    assert not wrong["passed"]


def test_a_sibling_isoform_in_the_top_tier_counts_as_recovery(weights):
    """The control is a gene; the ranking holds proteins. RefSeq carries one
    entry per isoform, so scoring against a single accession records a correct
    result as a miss whenever a sibling is what ranked."""
    records = build_population()
    # The canonical accession is nowhere in the ranking; its isoform is first.
    ranked = rank(records, weights)
    report = positive_control_report(
        ranked,
        {"P0_GENE": {"P0_CANONICAL", "POS_0"}},
        fraction=0.3,
        minimum_recovered=0.6,
    )
    assert report["n_controls"] == 1
    assert report["passed"]
    assert report["dropped_by_cascade"] == []
    # And it says which isoform carried it, so the dependency is visible.
    assert report["matched_accession"]["P0_GENE"] == "POS_0"


def test_a_control_is_missing_only_when_every_isoform_is(weights):
    ranked = rank(build_population(), weights)
    report = positive_control_report(
        ranked,
        {"GHOST": {"NOT_HERE", "ALSO_NOT_HERE"}},
        fraction=0.3,
        minimum_recovered=0.6,
    )
    assert report["dropped_by_cascade"] == ["GHOST"]
    assert report["n_in_ranking"] == 0


def test_isoforms_disagreeing_on_the_drop_step_is_a_failure(config):
    """Identical paralogs rejected by different rules is worth seeing, not
    worth resolving by majority."""
    agree_a = make_record("H3_A", signal_peptide=False, seed=11)
    agree_b = make_record("H3_B", signal_peptide=False, seed=12)
    for record in (agree_a, agree_b):
        run_cascade(record, config["cascade"])
    report = negative_control_report(
        [agree_a, agree_b], {"H3": ("exposure", {"H3_A", "H3_B"})}
    )
    assert report["passed"]
    assert report["details"][0]["n_accessions"] == 2

    # One sibling dropped by length instead: the control no longer tests what
    # it was chosen to test.
    odd = make_record("H3_C", length=30, seed=13)
    run_cascade(odd, config["cascade"])
    split = negative_control_report(
        [agree_a, odd], {"H3": ("exposure", {"H3_A", "H3_C"})}
    )
    assert not split["passed"]
    assert split["details"][0]["actual_drop_step"] == ["exposure", "length"]


def test_random_controls_land_mid_pack(weights):
    records = build_population()
    ranked = rank(records, weights)
    report = random_control_report(ranked, {f"BG_{i}" for i in range(4, 16)})
    assert report["n"] == 12
    assert report["fraction_in_top_decile"] == 0.0


def test_weight_perturbation_reports_stability(weights):
    records = build_population()
    report = weight_perturbation(
        records, weights, pct=20, n=25, fraction=0.3, seed=1
    )
    assert report["n_perturbations"] == 25
    # Candidates this separated from the background should be stable under
    # jitter; if they are not, the rubric is riding on one weight.
    assert report["mean_stability"] > 0.8


def test_leave_one_source_out_flags_fragile_candidates(weights):
    records = build_population()
    # A candidate carried entirely by its RNAi phenotype.
    fragile = make_record(
        "FRAGILE",
        seed=999,
        rnai="lethal",
        host_identity=38.0,
        n_epitopes=2,
        epitope_score=0.5,
        stages={"adult_female": 3.0, "adult_male": 0.1},
    )
    records.append(fragile)
    report = leave_one_source_out(
        records, weights, ["rnai", "expression", "literature"], fraction=0.5
    )
    assert set(report["robust"]) <= set(report["baseline_top_tier"])
    for protein_id in report["fragile"]:
        record = next(r for r in records if r.protein_id == protein_id)
        assert "single_source_dependent" in record.flags


def test_full_benchmark_produces_a_gate_verdict(weights, config):
    records = build_population()
    for record in records:
        run_cascade(record, config["cascade"])
    # Fewer perturbations than the real run; the shape of the report is what
    # this test is about, not the stability estimate.
    validation_cfg = {
        **config["validation"],
        "robustness": {**config["validation"]["robustness"], "n_perturbations": 10},
    }
    report = run_benchmark(
        records,
        weights,
        validation_cfg,
        {f"POS_{i}" for i in range(8)},
        {},
        {f"BG_{i}" for i in range(20)},
    )
    assert "gate_passed" in report
    assert report["positive_controls"]["n_controls"] == 8


# -- symbiont track --------------------------------------------------------


def symbiont_cfg(config):
    return config["symbiont"]


def test_environmental_taxa_are_excluded_before_scoring(config):
    mags = [
        MAG(
            mag_id="MAG1",
            taxon="Vibrio splendidus",
            genome_size_mb=4.8,
            completeness=92,
            contamination=2,
            prevalence_by_population={"NO": 1.0, "SCO": 1.0, "CA": 1.0},
        )
    ]
    report = gate(mags, symbiont_cfg(config))
    assert report["n_excluded_environmental"] == 1
    assert not report["proceed"]


def test_no_signature_is_a_clean_negative_not_an_error(config):
    mags = [
        MAG(
            mag_id="MAG2",
            taxon="Moritella viscosa",
            genome_size_mb=4.2,
            completeness=88,
            contamination=3,
            prevalence_by_population={"NO": 0.4},
        )
    ]
    report = gate(mags, symbiont_cfg(config))
    assert not report["proceed"]
    assert "stops here" in report["conclusion"]


def test_full_signature_opens_the_gate(config):
    mags = [
        MAG(
            mag_id="MAG3",
            taxon="Candidatus Lepeophtheirella",
            genome_size_mb=0.7,
            completeness=95,
            contamination=1,
            prevalence_by_population={"NO": 1.0, "SCO": 0.98, "CL": 0.99},
            lost_pathways=["lysine", "biotin", "riboflavin"],
            host_species=["lsal", "crog"],
            phylogenetic_congruence=0.86,
        )
    ]
    report = gate(mags, symbiont_cfg(config))
    assert report["proceed"]
    assert report["qualifying"][0]["both_louse_species"]
