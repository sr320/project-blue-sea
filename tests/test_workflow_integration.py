"""Execute the real Snakemake DAG with native-format stand-in executables.

This tests scheduling, handoffs and reruns. It does not validate model accuracy.
"""

import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from .conftest import REPO, random_sequence

DRIVER = r"""
import csv
import pathlib
import shutil
import sys

args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]
def fasta(path):
    result = []
    for line in pathlib.Path(path).read_text().splitlines():
        if line.startswith(">"):
            result.append([line[1:].split()[0], ""])
        else:
            result[-1][1] += line
    return result
if args[0] == "makedb":
    shutil.copyfile(value("--in"), value("-d") + ".dmnd")
elif args[0] == "blastp":
    with open(value("-o"), "w") as handle:
        for qid, qseq in fasta(value("-q")):
            for sid, sseq in fasta(value("-d")):
                length = min(len(qseq), len(sseq))
                identity = 100 if qseq == sseq else 20
                row = [qid, sid, identity, length, 0, 0, 1, length, 1, length,
                       "1e-30", identity * length]
                if "qseq" in args:
                    row += [qseq[:length], sseq[:length]]
                handle.write("\t".join(map(str, row)) + "\n")
else:
    tool, source, destination = args
    destination = pathlib.Path(destination)
    destination.mkdir()
    records = fasta(source)
    if tool == "signalp":
        with open(destination / "prediction_results.txt", "w") as handle:
            for pid, seq in records:
                sp = "NEG" not in pid
                row = [pid, "SP" if sp else "OTHER", 0.1 if sp else 0.9,
                       0.9 if sp else 0.1]
                handle.write("\t".join(map(str, row)) + "\n")
    elif tool == "deeptmhmm":
        with open(destination / "predicted_topologies.3line", "w") as handle:
            for pid, seq in records:
                topology = ("I" * len(seq) if "NEG" in pid
                            else "S" * 20 + "O" * (len(seq) - 20))
                handle.write(f">{pid} | SP\n{seq}\n{topology}\n")
    else:
        with open(destination / "raw_output.csv", "w") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Accession", "Residue", "BepiPred-3.0 score",
                             "BepiPred-3.0 linear epitope score"])
            for pid, seq in records:
                for index, residue in enumerate(seq):
                    writer.writerow([pid, residue, 0.8, 0.8 if index >= 20 else 0.1])
"""


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="requires workflow extra")
def test_real_dag_runs_with_isolated_native_predictors(tmp_path, config):
    (tmp_path / "config").mkdir()
    shutil.copyfile(REPO / "config/weights.yaml", tmp_path / "config/weights.yaml")
    (tmp_path / "src").symlink_to(REPO / "src", target_is_directory=True)
    driver = tmp_path / "driver.py"
    driver.write_text(DRIVER)
    binaries = tmp_path / "bin"
    binaries.mkdir()
    diamond = binaries / "diamond"
    diamond.write_text(f"#!{sys.executable}\n" + DRIVER)
    diamond.chmod(0o755)
    cfg = copy.deepcopy(config)
    cfg["paths"].update(scratch="scratch", results="results", logs="logs")
    cfg["n_shards"] = cfg["n_shards_candidates"] = 2
    cfg["validation"].update(
        positive_controls="positive.tsv",
        negative_controls="negative.tsv",
        construct_controls=None,
        host_spike_n=2,
        require_host_spikes=True,
        top_tier_fraction=0.5,
    )
    cfg["validation"]["robustness"]["n_perturbations"] = 5
    (tmp_path / "positive.tsv").write_text("protein_id\tlsal\nlsal_POS\tP0\n")
    (tmp_path / "negative.tsv").write_text(
        "protein_id\texpected_drop_step\nlsal_NEG\texposure\ncrog_NEG\texposure\n"
    )
    for name in cfg["species"]:
        directory = tmp_path / "resources/reference" / name
        directory.mkdir(parents=True)
        (directory / "proteins.faa").write_text(
            "".join(
                f">{name}_{pid}\n{random_sequence(100, seed=i)}\n"
                for i, pid in enumerate(["POS", "NOVEL", "NEG"])
            )
        )
        (directory / "annotation.gff3").write_text("##gff-version 3\n")
    host = tmp_path / "resources/reference/ssal"
    host.mkdir(parents=True)
    (host / "proteins.faa").write_text(
        "".join(f">HOST_{i}\n{random_sequence(110, seed=i + 100)}\n" for i in range(2))
    )
    for name, predictor in cfg["predictors"].items():
        predictor.update(
            version="native-format-fixture",
            required_files=[str(driver)],
            command=[sys.executable, str(driver), name, "{input}", "{output_dir}"],
        )
    (tmp_path / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    env = {
        **os.environ,
        "PATH": f"{binaries}:{Path(sys.executable).parent}:{os.environ['PATH']}",
        "XDG_CACHE_HOME": str(tmp_path / "cache-home"),
    }
    command = [
        shutil.which("snakemake"),
        "-s",
        str(REPO / "workflow/Snakefile"),
        "all",
        "--cores",
        "4",
        "--resources",
        "network=1",
    ]
    result = subprocess.run(
        command, cwd=tmp_path, env=env, text=True, capture_output=True
    )
    (tmp_path / "workflow.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, (result.stdout + result.stderr)[-10000:]
    report = json.loads((tmp_path / "results/benchmark.json").read_text())
    assert report["negative_controls"]["passed"]
    assert report["host_spikes_present"]
    for species in cfg["species"]:
        original = json.loads((tmp_path / f"scratch/{species}/cascaded.json").read_text())
        enriched = json.loads(
            (tmp_path / f"scratch/{species}/orthologized.json").read_text()
        )
        assert not any(r["orthologs"] for r in original)
        assert any(r["orthologs"] for r in enriched)
        assert not any(r["literature"] for r in enriched)
    assert "BLUESEA_HOST_SPIKE_" not in (tmp_path / "results/ranked_all.tsv").read_text()
    assert list((tmp_path / "scratch/lsal/raw/signalp").glob("*/manifest.json"))
    sensitivity = json.loads((tmp_path / "results/ranking_sensitivity.json").read_text())
    assert sensitivity["weight_perturbation"]["n_perturbations"] == 5
    again = subprocess.run(command, cwd=tmp_path, env=env, text=True, capture_output=True)
    assert again.returncode == 0
    assert "Nothing to be done" in again.stdout + again.stderr
    # A changed weight is a declared dependency, so rank/sensitivity rerun.
    with open(tmp_path / "config/weights.yaml", "a") as handle:
        handle.write("\n# dependency change\n")
    rerun = subprocess.run(
        [*command, "-n"], cwd=tmp_path, env=env, text=True, capture_output=True
    )
    assert rerun.returncode == 0
    assert "rule rank:" in rerun.stdout + rerun.stderr

    override = tmp_path / "override.yaml"
    override.write_text(
        yaml.safe_dump(
            {
                "paths": {
                    **cfg["paths"],
                    "scratch": "other-scratch",
                    "results": "other-results",
                    "logs": "other-logs",
                }
            }
        )
    )
    overridden = subprocess.run(
        [*command, "-n", "--printshellcmds", "--configfile", str(override)],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
    )
    assert overridden.returncode == 0, overridden.stdout + overridden.stderr
    assert "other-scratch/resolved_config.yaml" in overridden.stdout + overridden.stderr
