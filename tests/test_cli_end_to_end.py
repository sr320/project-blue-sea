"""End-to-end run on a synthetic nine-protein proteome.

Covers what the unit tests cannot: the tool-output parsers, the JSON handoff
between stages, and whether the whole chain produces a ranking and dossiers
without a cluster or an API key.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from bluesea.cli import main

from .conftest import make_record, random_sequence

# protein_id, signal peptide, n_TM, salmon identity, n epitopes, rnai, expressed
PROTEOME = [
    ("LSAL_P0", True, 0, 28.0, 6, "lethal", True),
    ("LSAL_BM86", True, 1, 21.0, 7, "fecundity_reduced", True),
    ("LSAL_HPX1", True, 0, 33.0, 5, "untested", True),
    ("LSAL_NOVEL1", True, 0, 19.0, 4, "untested", True),
    ("LSAL_NOVEL2", True, 2, 25.0, 3, "moulting_defect", True),
    ("LSAL_RPL7", False, 0, 35.0, 4, "lethal", True),        # intracellular
    ("SSAL_ACTB", True, 0, 97.0, 5, "untested", True),       # host-identical
    ("LSAL_SHORT", True, 0, 20.0, 3, "untested", True),      # too short
    ("LSAL_SILENT", True, 0, 22.0, 4, "untested", False),    # not expressed
]

STAGES = ["copepodid", "chalimus_I", "preadult_I", "adult_female"]


@pytest.fixture
def workspace(tmp_path: Path, config: dict) -> Path:
    reference = tmp_path / "resources" / "reference" / "lsal"
    scratch = tmp_path / "scratch" / "lsal"
    reference.mkdir(parents=True)
    scratch.mkdir(parents=True)

    signalp, tmhmm, epitopes, host, expression, rnai = [], [], [], [], [], []
    fasta = []

    for index, (pid, sp, n_tm, identity, n_ep, phenotype, expressed) in enumerate(
        PROTEOME
    ):
        length = 40 if pid == "LSAL_SHORT" else 400
        sequence = random_sequence(length, seed=index)
        fasta.append(f">{pid} [gene={pid}_g] test protein\n{sequence}")

        signalp.append(
            f"{pid}\t{'SP' if sp else 'OTHER'}\t{0.96 if sp else 0.02}\t"
            f"{0.04 if sp else 0.98}"
        )
        # 3-line DeepTMHMM: header, sequence, topology string.
        inside = "I" * (length - 40 * n_tm - (10 if sp else 0))
        topology = (
            ("S" * 10 if sp else "")
            + ("M" * 20 + "O" * 20) * n_tm
            + inside[: max(0, length - (10 if sp else 0) - 40 * n_tm)]
        )
        topology = (topology + "O" * length)[:length]
        tmhmm.extend([f">{pid}", sequence, topology])

        for epitope_index in range(n_ep):
            start = 10 + epitope_index * 30
            epitopes.append(f"{pid}\t{start}\t{start + 14}\t0.82\tlinear")

        host.append(
            f"{pid}\tXP_salmon_{index}\t{identity}\t300\t0\t0\t1\t300\t1\t300\t1e-40\t400"
        )
        tpm = [50.0, 42.0, 31.0, 28.0] if expressed else [0.0, 0.0, 0.0, 0.1]
        expression.append(f"{pid}_g\t" + "\t".join(str(v) for v in tpm))
        rnai.append(f"{pid}_g\t{phenotype}\tLiceBase\t")

    (reference / "proteins.faa").write_text("\n".join(fasta) + "\n")
    (reference / "stage_tpm.tsv").write_text(
        "gene_id\t" + "\t".join(STAGES) + "\n" + "\n".join(expression) + "\n"
    )
    (reference / "rnai_phenotypes.tsv").write_text(
        "gene_id\tphenotype\tsource\tnote\n" + "\n".join(rnai) + "\n"
    )
    (scratch / "signalp.tsv").write_text("\n".join(signalp) + "\n")
    (scratch / "deeptmhmm.gff3").write_text("\n".join(tmhmm) + "\n")
    (scratch / "epitopes.tsv").write_text(
        "protein_id\tstart\tend\tscore\tkind\n" + "\n".join(epitopes) + "\n"
    )
    (scratch / "vs_host.tsv").write_text("\n".join(host) + "\n")

    # Single-species config; orthology needs two proteomes and is covered by
    # its own unit tests.
    local = dict(config)
    local["workdir"] = str(tmp_path)
    local["species"] = {"lsal": config["species"]["lsal"]}
    local["life_stages"] = STAGES
    local["validation"] = {
        **config["validation"],
        "positive_controls": str(tmp_path / "positive.tsv"),
        "negative_controls": str(tmp_path / "negative.tsv"),
        "random_control_n": 2,
        "min_positive_controls_in_top_tier": 0.5,
        "top_tier_fraction": 0.5,
        "robustness": {**config["validation"]["robustness"], "n_perturbations": 5},
    }
    (tmp_path / "positive.tsv").write_text(
        "protein_id\tgene_symbol\n LSAL_P0\tP0\n".replace(" ", "")
    )
    (tmp_path / "negative.tsv").write_text(
        "protein_id\texpected_drop_step\n"
        "SSAL_ACTB\thost_similarity\n"
        "LSAL_RPL7\texposure\n"
    )

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(local))
    return tmp_path


def run(workspace: Path, *args: str) -> int:
    weights = Path(__file__).resolve().parents[1] / "config" / "weights.yaml"
    return main(
        [
            "--config",
            str(workspace / "config.yaml"),
            "--weights",
            str(weights),
            *args,
        ]
    )


def test_pipeline_runs_end_to_end(workspace: Path, capsys):
    assert run(workspace, "prefilter", "--species", "lsal") == 0
    assert run(workspace, "assemble", "--species", "lsal") == 0
    assert run(workspace, "cascade", "--species", "lsal") == 0
    # Dry run: no API key, no network, no spend.
    assert run(workspace, "literature", "--species", "lsal") == 0
    assert run(workspace, "rank") == 0
    assert run(workspace, "dossiers", "--top", "5") == 0

    results = workspace / "results"
    ranked = (results / "ranked_all.tsv").read_text().splitlines()
    header, rows = ranked[0].split("\t"), ranked[1:]
    assert "total_score" in header

    survivors = {row.split("\t")[1] for row in rows}
    assert "LSAL_P0" in survivors
    assert "SSAL_ACTB" not in survivors    # host-identical
    assert "LSAL_RPL7" not in survivors    # intracellular
    assert "LSAL_SHORT" not in survivors   # below length floor
    assert "LSAL_SILENT" not in survivors  # not expressed

    dossiers = sorted((results / "dossiers").glob("*.md"))
    assert dossiers
    text = dossiers[0].read_text()
    assert "## Scores" in text
    assert "## Why it survived" in text
    assert "challenge trial" in text   # the standing caveat


def test_cascade_audit_traces_every_drop(workspace: Path):
    run(workspace, "prefilter", "--species", "lsal")
    run(workspace, "assemble", "--species", "lsal")
    run(workspace, "cascade", "--species", "lsal")

    audit = workspace / "logs" / "lsal" / "cascade_audit.jsonl"
    entries = {
        json.loads(line)["protein_id"]: json.loads(line)
        for line in audit.read_text().splitlines()
    }
    assert len(entries) == len(PROTEOME)
    assert entries["SSAL_ACTB"]["dropped_at"] == "host_similarity"
    assert entries["LSAL_RPL7"]["dropped_at"] == "exposure"
    assert entries["LSAL_SHORT"]["dropped_at"] == "length"
    assert entries["LSAL_P0"]["dropped_at"] is None
    # Every drop carries the rule that caused it, not just the step name.
    failing = [v for v in entries["SSAL_ACTB"]["verdicts"] if not v["passed"]]
    assert "identity to salmon" in failing[0]["reason"]


def test_prefilter_forwards_only_reachable_proteins(workspace: Path):
    """The point of the split: the costly predictors never see a protein that
    the cheap filters already excluded."""
    assert run(workspace, "prefilter", "--species", "lsal") == 0

    candidates = (workspace / "scratch" / "lsal" / "candidates.faa").read_text()
    forwarded = {
        line[1:].split()[0] for line in candidates.splitlines() if line.startswith(">")
    }
    assert "LSAL_RPL7" not in forwarded     # intracellular: no route to antibody
    assert "LSAL_SHORT" not in forwarded    # below the length floor
    assert "LSAL_P0" in forwarded
    # Still judged on cross-reactivity, so it has to reach the salmon search.
    assert "SSAL_ACTB" in forwarded
    assert len(forwarded) < len(PROTEOME)

    # Dropped proteins stay in the record set, or the audit would lose them.
    records = json.loads(
        (workspace / "scratch" / "lsal" / "prefiltered.json").read_text()
    )
    assert len(records) == len(PROTEOME)


def test_disabling_the_prefilter_forwards_everything(workspace: Path):
    """The escape hatch for measuring what the exposure filter costs: every
    protein goes to the predictors, and the verdicts are unchanged."""
    config_path = workspace / "config.yaml"
    local = yaml.safe_load(config_path.read_text())
    local["cascade"]["prefilter"]["enabled"] = False
    config_path.write_text(yaml.safe_dump(local))

    assert run(workspace, "prefilter", "--species", "lsal") == 0
    candidates = (workspace / "scratch" / "lsal" / "candidates.faa").read_text()
    forwarded = {
        line[1:].split()[0] for line in candidates.splitlines() if line.startswith(">")
    }
    assert len(forwarded) == len(PROTEOME)

    run(workspace, "assemble", "--species", "lsal")
    run(workspace, "cascade", "--species", "lsal", "--profile")

    audit = workspace / "logs" / "lsal" / "cascade_audit.jsonl"
    entries = {
        json.loads(line)["protein_id"]: json.loads(line)
        for line in audit.read_text().splitlines()
    }
    # Drop attribution is unchanged: the first failing step still owns the drop.
    assert entries["LSAL_RPL7"]["dropped_at"] == "exposure"
    assert entries["SSAL_ACTB"]["dropped_at"] == "host_similarity"
    # But now an exposure-dropped protein also carries the later verdicts,
    # which is the whole point of running this way.
    steps = [v["step"] for v in entries["LSAL_RPL7"]["verdicts"]]
    assert "epitope" in steps and "expression" in steps


def test_random_controls_never_include_a_control_isoform(config):
    """Drawing a control's own sibling as a 'random' protein would put a
    known antigen in the set that is supposed to represent the background."""
    from bluesea.cli import _sample_random_controls

    records = [make_record(f"P{i:03d}", seed=i) for i in range(40)]
    positives = {"CTRL_A": {"P000", "P001"}}
    negatives = {"CTRL_B": ("exposure", {"P002"})}
    drawn = _sample_random_controls(
        records,
        {**config["validation"], "random_control_n": 30},
        positives,
        negatives,
    )
    assert drawn
    assert not drawn & {"P000", "P001", "P002"}


def test_isoforms_column_is_optional(config):
    from bluesea.cli import _control_accessions

    assert _control_accessions({"protein_id": "XP_1"}) == {"XP_1"}
    assert _control_accessions({"protein_id": "XP_1", "isoforms": ""}) == {"XP_1"}
    assert _control_accessions(
        {"protein_id": "XP_1", "isoforms": "XP_2; XP_3"}
    ) == {"XP_1", "XP_2", "XP_3"}


def test_a_broken_gene_join_is_loud(capsys):
    """Expression and RNAi join on gene_id. NCBI protein.faa carries no
    [gene=] tag, so a table keyed on anything but the protein accession joins
    nothing -- and essentiality, the heaviest criterion, silently goes flat."""
    from bluesea.cli import _report_joins

    _report_joins(10, [("rnai", {"LOC121131766": {}}, 0)])
    captured = capsys.readouterr()
    assert "rnai: 0/10 joined (0.0%)" in captured.out
    assert "no record joined to the rnai table" in captured.err


def test_a_working_gene_join_is_quiet(capsys):
    from bluesea.cli import _report_joins

    _report_joins(10, [("expression", {"g1": {}}, 7)])
    captured = capsys.readouterr()
    assert "expression: 7/10 joined (70.0%)" in captured.out
    assert captured.err == ""


def test_an_absent_table_is_not_a_broken_join(capsys):
    """crog has `rnai: null` -- no resource published. That is a documented
    absence, not a join that failed, and must not raise the alarm."""
    from bluesea.cli import _report_joins

    _report_joins(10, [("rnai", {}, 0)])
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_dry_run_reports_projected_spend_without_calling(workspace: Path, capsys):
    run(workspace, "prefilter", "--species", "lsal")
    run(workspace, "assemble", "--species", "lsal")
    run(workspace, "cascade", "--species", "lsal")
    run(workspace, "literature", "--species", "lsal")

    report = json.loads((workspace / "logs" / "llm_usage_lsal.json").read_text())
    assert report["dry_run"] is True
    assert report["stages"]["literature"]["calls"] > 0
    assert report["total_cost_usd"] > 0
    # Nothing was cached, because nothing was called.
    assert not list((workspace / "cache" / "llm").rglob("*.json"))


def test_benchmark_gate_reports_a_verdict(workspace: Path):
    run(workspace, "prefilter", "--species", "lsal")
    run(workspace, "assemble", "--species", "lsal")
    run(workspace, "cascade", "--species", "lsal")
    exit_code = run(workspace, "benchmark")

    report = json.loads((workspace / "results" / "benchmark.json").read_text())
    assert report["negative_controls"]["passed"]
    assert report["positive_controls"]["n_controls"] == 1
    # Exit code is the gate: non-zero stops the workflow before the full run.
    assert exit_code == (0 if report["gate_passed"] else 1)
