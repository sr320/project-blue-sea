import math

from bluesea.records import LiteratureEvidence
from bluesea.scoring import CRITERIA, rank, score_record

from .conftest import make_record


def test_untested_rnai_scores_above_no_phenotype(weights):
    """Partial RNAi coverage means untested must outscore tested-and-null, or
    the screen quietly punishes every gene nobody has knocked down yet."""
    untested = score_record(make_record("A", rnai="untested"), weights)
    null = score_record(make_record("B", rnai="no_phenotype"), weights)
    assert untested.scores["essentiality"] > null.scores["essentiality"]
    assert "essentiality_unknown" in untested.flags


def test_lethal_phenotype_scores_highest(weights):
    lethal = score_record(make_record("A", rnai="lethal"), weights)
    assert lethal.scores["essentiality"] == 1.0


def test_host_dissimilarity_rewards_distance(weights):
    near = score_record(make_record("A", host_identity=39.0), weights)
    far = score_record(make_record("B", host_identity=16.0), weights)
    none = score_record(make_record("C", host_identity=None), weights)
    assert none.scores["host_dissimilarity"] == 1.0
    assert far.scores["host_dissimilarity"] > near.scores["host_dissimilarity"]


def test_epitope_count_saturates(weights):
    """A long protein must not win on surface area alone."""
    many = score_record(make_record("A", length=1800, n_epitopes=40), weights)
    enough = score_record(make_record("B", n_epitopes=8), weights)
    assert math.isclose(many.scores["epitope"], enough.scores["epitope"], abs_tol=1e-6)


def test_conservation_needs_an_ortholog(weights):
    with_ortholog = score_record(make_record("A", ortholog_identity=80.0), weights)
    without = score_record(make_record("B"), weights)
    assert with_ortholog.scores["conservation"] > without.scores["conservation"]


def test_expressibility_penalises_multipass_membrane(weights):
    simple = score_record(make_record("A", length=300), weights)
    multipass = score_record(make_record("B", length=300, n_transmembrane=7), weights)
    assert multipass.scores["expressibility"] < simple.scores["expressibility"]


def test_missing_literature_uses_neutral_default_not_zero(weights):
    record = score_record(make_record("A"), weights)
    assert record.scores["literature"] == weights["literature"]["default_when_missing"]
    assert "literature_not_assessed" in record.flags


def test_literature_score_is_used_when_present(weights):
    record = make_record("A")
    record.literature = LiteratureEvidence(score=0.95, justification="x", citations=["1"])
    score_record(record, weights)
    assert record.scores["literature"] == 0.95


def test_dropping_a_source_keeps_totals_comparable(weights):
    """Leave-one-source-out renormalises, so a dropped axis does not simply
    look like a universally worse run."""
    record = make_record("A", rnai="lethal", host_identity=20.0)
    full = score_record(make_record("A", rnai="lethal", host_identity=20.0), weights)
    dropped = score_record(record, weights, drop_sources={"literature"})
    assert math.isnan(dropped.scores["literature"])
    assert dropped.total_score > full.total_score  # literature was the weak axis


def test_ranking_is_descending_and_complete(weights):
    records = [
        make_record("GOOD", rnai="lethal", host_identity=18.0, ortholog_identity=85.0),
        make_record("MID", rnai="fecundity_reduced", host_identity=30.0),
        make_record("WEAK", rnai="no_phenotype", host_identity=39.0, n_epitopes=2),
    ]
    ranked = rank(records, weights)
    assert [r.protein_id for r in ranked] == ["GOOD", "MID", "WEAK"]
    assert all(set(CRITERIA) <= set(r.scores) for r in ranked)
    scores = [r.total_score for r in ranked]
    assert scores == sorted(scores, reverse=True)
