from bluesea.cascade import cascade_summary, estimate_low_complexity, run_cascade

from .conftest import make_record


def test_good_candidate_survives(config, good_candidate):
    run_cascade(good_candidate, config["cascade"])
    assert good_candidate.survived
    assert good_candidate.dropped_at is None


def test_intracellular_dropped_at_exposure(config, intracellular):
    run_cascade(intracellular, config["cascade"])
    assert intracellular.dropped_at == "exposure"


def test_host_like_dropped_at_host_similarity(config, host_like):
    run_cascade(host_like, config["cascade"])
    assert host_like.dropped_at == "host_similarity"


def test_host_hit_below_coverage_threshold_is_not_an_exclusion(config):
    """A short, high-identity patch is not grounds to exclude: the filter asks
    for a significant alignment, not any alignment."""
    record = make_record("LSAL_PATCH", host_identity=95.0)
    record.host_hit.query_coverage = 0.05
    run_cascade(record, config["cascade"])
    assert record.survived


def test_every_step_records_a_verdict_when_profiling(config, intracellular):
    run_cascade(intracellular, config["cascade"], stop_on_fail=False)
    steps = [verdict["step"] for verdict in intracellular.verdicts]
    assert steps == [
        "length",
        "low_complexity",
        "exposure",
        "host_similarity",
        "epitope",
        "expression",
        "transmembrane",
    ]


def test_cascade_stops_after_first_failure(config, intracellular):
    run_cascade(intracellular, config["cascade"])
    steps = [verdict["step"] for verdict in intracellular.verdicts]
    assert steps[-1] == "exposure"


def test_membrane_protein_with_large_ectodomain_passes_exposure(config):
    record = make_record("LSAL_TM", signal_peptide=False, n_transmembrane=1)
    run_cascade(record, config["cascade"])
    assert record.survived
    assert record.topology.exposure_class == "type_I_membrane"


def test_missing_expression_data_does_not_exclude(config):
    """Absence of evidence must not read as evidence of absence."""
    record = make_record("LSAL_NOEXPR", stages={})
    run_cascade(record, config["cascade"])
    assert record.survived
    assert "expression_unknown" in record.flags


def test_low_complexity_detects_repeats():
    assert estimate_low_complexity("A" * 100) > 0.9
    assert estimate_low_complexity("ACDEFGHIKLMNPQRSTVWY" * 5) < 0.1


def test_summary_accounts_for_every_protein(config):
    records = [
        make_record("A", rnai="lethal", host_identity=20.0),
        make_record("B", signal_peptide=False, seed=2),
        make_record("C", host_identity=99.0, seed=3),
    ]
    for record in records:
        run_cascade(record, config["cascade"])
    summary = cascade_summary(records)
    assert sum(row["dropped"] for row in summary) == 2
    assert summary[-1]["survived"] == 1
