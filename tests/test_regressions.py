"""Scientific and ranking regressions found in the repository review."""

import copy
import json

import pytest

from bluesea import io, orthology, validation
from bluesea.cascade import filter_epitopes, filter_expression, filter_host_similarity
from bluesea.dossier import construct_boundary
from bluesea.records import Epitope, Ortholog, Topology
from bluesea.scoring import gene_representatives, rank
from bluesea.screening import require_pilot, signature

from .conftest import make_record


def test_short_identity_hit_cannot_mask_long_host_homolog(config):
    hits = [
        dict(subject="patch", identity_pct=95, evalue=1e-10, qstart=1, qend=100),
        dict(subject="homolog", identity_pct=80, evalue=1e-100, qstart=1, qend=900),
    ]
    cfg = config["cascade"]["host_similarity"]
    record = make_record(length=1000)
    record.host_hits = io.host_alignments(hits, record.length)
    record.host_hit = io.best_host_hit(hits, record.length, cfg)
    assert record.host_hit.accession == "homolog"
    assert not filter_host_similarity(record, cfg)[0]
    assert "short_host_similar_patch" in record.flags


def test_host_identity_at_threshold_is_excluded(config):
    record = make_record(host_identity=40)
    assert not filter_host_similarity(record, config["cascade"]["host_similarity"])[0]


def test_missing_negative_control_fails_even_when_another_passes():
    record = make_record("present")
    record.record_verdict("exposure", False, "test")
    report = validation.negative_control_report(
        [record], {"present": "exposure", "missing": "host_similarity"}
    )
    assert not report["passed"]
    assert report["n"] == 2
    assert report["details"][1]["actual_drop_step"] == "missing"


def test_caligus_stages_are_preserved_and_wrong_stages_rejected(config, tmp_path):
    path = tmp_path / "expression.tsv"
    path.write_text("gene_id\tchalimus_III\tchalimus_IV\ng1\t4\t8\n")
    stages = io.species_stages(config, "crog")
    assert io.read_expression(path, stages)["g1"] == {"chalimus_III": 4, "chalimus_IV": 8}
    with pytest.raises(ValueError, match="Unknown expression stages"):
        io.read_expression(path, io.species_stages(config, "lsal"))


def test_single_assayed_stage_is_insufficient_to_exclude(config):
    record = make_record(stages={"chalimus_I": 1.1})
    assert filter_expression(record, config["cascade"]["expression"])[0]
    assert "expression_insufficient_stages" in record.flags


def test_epitope_accessibility_and_advisory_policy(config):
    record = make_record(length=100)
    record.topology.topology_string = "I" * 50 + "O" * 50
    record.epitopes = [Epitope(2, 20, 0.9), Epitope(60, 80, 0.9)]
    assert [e.start for e in record.accessible_epitopes()] == [60]
    cfg = config["cascade"]["epitope"]
    assert filter_epitopes(record, cfg)[0]
    assert not filter_epitopes(record, {**cfg, "hard_filter": True})[0]


def test_construct_boundary_uses_coordinates_and_withholds_unknown():
    record = make_record(length=100)
    record.topology = Topology(
        has_signal_peptide=True,
        n_transmembrane=1,
        topology_string="S" * 20 + "O" * 50 + "M" * 20 + "I" * 10,
    )
    record.epitopes = [Epitope(25, 40, 0.9)]
    assert construct_boundary(record) == (21, 70)
    record.topology.topology_string = ""
    assert construct_boundary(record) is None


def test_partial_alignment_does_not_claim_perfect_epitope_conservation():
    epitope = [Epitope(1, 10, 0.9)]
    assert orthology.epitope_region_identity("A" * 10, "A" * 10, epitope, [(0, 0)]) == 10


def test_intersection_requires_reciprocity_and_bilateral_epitope_support(config):
    a, b = make_record("a"), make_record("b")
    b.species = "crog"
    a.orthologs = [Ortholog("crog", "b", 80, 1e-40, False, 10, 1)]
    b.orthologs = [Ortholog("lsal", "a", 80, 1e-40, True, 80, 1)]
    assert not orthology.intersect_shortlists([a], [b], config["orthology"])
    a.orthologs[0].reciprocal_best = True
    a.orthologs[0].epitope_region_identity_pct = 80
    assert len(orthology.intersect_shortlists([a], [b], config["orthology"])) == 1
    b.orthologs[0].epitope_region_coverage = None
    assert not orthology.intersect_shortlists([a], [b], config["orthology"])
    assert "cross_species_candidate" not in a.flags


def test_gapped_alignment_maps_residues():
    hit = dict(qseq="AC-DE", sseq="A-CDE", qstart=2, qend=5, sstart=3, send=6)
    assert orthology.alignment_positions(hit) == [(1, 2), (3, 4), (4, 5)]


def test_sensitivity_preserves_scores_and_uses_gene_representatives(config, weights):
    records = [make_record("a", rnai="lethal"), make_record("a_iso"), make_record("b")]
    records[1].gene_id = records[0].gene_id
    ranked = rank(records, weights)
    before = [copy.deepcopy(r.scores) for r in ranked]
    representatives = gene_representatives(ranked)
    assert len(representatives) == 2
    validation.sensitivity_report(representatives, weights, config["validation"])
    assert [r.scores for r in ranked] == before


def test_pending_constructs_cannot_pass_gate(config, weights):
    with open(config["validation"]["construct_controls"]) as handle:
        registry = json.load(handle)
    report = validation.construct_report([], weights, config["validation"], registry)
    assert not report["ready"]
    assert report["independent_studies"] == 0
    assert all(not c["eligible_for_gate"] for c in report["constructs"])


def test_full_run_rejects_stale_or_failed_pilot(config, weights, tmp_path):
    cfg = copy.deepcopy(config)
    cfg["screening"]["mode"] = "full"
    path = tmp_path / "benchmark.json"
    cfg["screening"]["pilot_report"] = str(path)
    path.write_text(
        json.dumps(
            {
                "gate_passed": True,
                "screening_mode": "pilot",
                "screening_signature": signature(cfg, weights),
            }
        )
    )
    require_pilot(cfg, weights)
    cfg["cascade"]["host_similarity"]["max_identity_pct"] = 90
    with pytest.raises(ValueError, match="matching"):
        require_pilot(cfg, weights)


def test_signalp_cleavage_is_respected_when_topology_disagrees():
    record = make_record(length=100)
    record.topology = Topology(
        has_signal_peptide=True,
        signal_peptide_end=25,
        topology_string="S" * 20 + "O" * 80,
    )
    record.epitopes = [Epitope(21, 30, 0.9), Epitope(40, 55, 0.9)]
    assert [e.start for e in record.accessible_epitopes()] == [40]
    assert construct_boundary(record) == (26, 100)
