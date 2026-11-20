"""Native-format adapters, isolated execution, and coverage validation."""

import copy
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor

import pytest

from bluesea import io
from bluesea.bepipred import tiles
from bluesea.predictors import bepipred_intervals, run_predictor


def native_csv(path, pid="p", sequence="ACDEFGHIKLMNPQRSTVWY", scores=None):
    with open(path, "w") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Accession",
                "Residue",
                "BepiPred-3.0 score",
                "BepiPred-3.0 linear epitope score",
            ]
        )
        for residue, score in zip(sequence, scores or [0.8] * len(sequence), strict=True):
            writer.writerow([pid, residue, score, score])


def test_bepipred_native_scores_become_nonoverlapping_intervals(tmp_path):
    path = tmp_path / "raw.csv"
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    native_csv(path, scores=[0.8] * 7 + [0.1] * 3 + [0.9] * 10)
    rows = bepipred_intervals(path, {"p": sequence}, 0.5, 7)
    assert [(r["start"], r["end"]) for r in rows] == [(1, 7), (11, 20)]
    assert rows[1]["sequence"] == sequence[10:]


def test_bepipred_no_epitopes_is_valid_but_truncation_is_not(tmp_path):
    path = tmp_path / "raw.csv"
    sequence = "ACDEFGHIKLMNPQRSTVWY"
    native_csv(path, scores=[0.1] * len(sequence))
    assert bepipred_intervals(path, {"p": sequence}, 0.5, 7) == []
    with pytest.raises(ValueError, match="Truncated"):
        bepipred_intervals(path, {"p": sequence + "A"}, 0.5, 7)


def test_signalp_uses_sp_column_and_preserves_cleavage(tmp_path):
    path = tmp_path / "signalp.tsv"
    path.write_text(
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "p\tSP\t0.1\t0.9\tCS pos: 23-24. Pr: 0.95\n"
    )
    result = io.read_signalp(path)["p"]
    assert result["probability"] == 0.9
    assert result["signal_peptide_end"] == 23


def test_three_line_parser_preserves_pipe_identifiers_and_rejects_truncation(tmp_path):
    path = tmp_path / "tm.3line"
    path.write_text(">sp|p|name | SP\nACDE\nSSOO\n")
    assert "sp|p|name" in io.read_deeptmhmm(path)
    path.write_text(">p\nACDE\nSSO\n")
    with pytest.raises(ValueError, match="Malformed"):
        io.read_deeptmhmm(path)


def test_long_sequences_tile_with_complete_coverage():
    result = tiles({"long": "A" * 2000})
    covered = {
        p for _, start, seq in result.values() for p in range(start, start + len(seq))
    }
    assert covered == set(range(2000))
    assert max(len(seq) for _, _, seq in result.values()) <= 1000


def test_concurrent_predictors_have_separate_native_files(config, tmp_path):
    driver = tmp_path / "driver.py"
    driver.write_text("""import pathlib, sys
source, destination = map(pathlib.Path, sys.argv[1:])
header, sequence = source.read_text().splitlines()
destination.mkdir()
(destination / "prediction_results.txt").write_text(header[1:] + "\\tSP\\t0.1\\t0.9\\n")
""")
    cfg = copy.deepcopy(config)
    cfg["predictors"]["signalp"].update(
        version="fixture-1",
        required_files=[str(driver)],
        command=[sys.executable, str(driver), "{input}", "{output_dir}"],
    )
    for pid in ("one", "two"):
        (tmp_path / f"{pid}.faa").write_text(f">{pid}\nACDEFGHIK\n")

    def run(pid):
        run_predictor(
            cfg,
            "signalp",
            tmp_path / f"{pid}.faa",
            tmp_path / f"{pid}.tsv",
            tmp_path / f"raw-{pid}",
        )

    with ThreadPoolExecutor(2) as executor:
        list(executor.map(run, ("one", "two")))
    for pid in ("one", "two"):
        assert set(io.read_signalp(tmp_path / f"{pid}.tsv")) == {pid}
        assert (
            json.loads((tmp_path / f"raw-{pid}" / "manifest.json").read_text())[
                "n_sequences"
            ]
            == 1
        )
