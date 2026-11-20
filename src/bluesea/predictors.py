"""Isolated local predictor execution and strict native-output normalization."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from . import io
from .provenance import file_digest

EPITOPE_COLUMNS = ["protein_id", "start", "end", "score", "kind", "sequence"]


def fasta_sequences(path: str | Path) -> dict[str, str]:
    records = list(io.read_fasta(path))
    sequences = {pid: seq for pid, _, seq in records}
    if not sequences or len(sequences) != len(records):
        raise ValueError("Predictor input must be nonempty with unique FASTA identifiers")
    return sequences


def require_ids(actual: set[str], expected: set[str], name: str) -> None:
    if actual != expected:
        raise ValueError(
            f"{name} ID mismatch: missing={sorted(expected - actual)[:10]}, "
            f"unexpected={sorted(actual - expected)[:10]}"
        )


def bepipred_intervals(
    path: str | Path,
    sequences: dict[str, str],
    threshold: float,
    min_length: int,
) -> list[dict[str, Any]]:
    """Native BP3 raw_output.csv -> contiguous above-threshold linear patches.

    Positions are implicit per accession in the native CSV. Full sequence
    coverage and exact residues are required, including proteins with no hits.
    """
    if not 0 <= threshold <= 1 or min_length < 1:
        raise ValueError("Invalid epitope threshold or minimum length")
    scores: dict[str, list[float]] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle)
        required = {"Accession", "Residue", "BepiPred-3.0 linear epitope score"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Unsupported BepiPred raw CSV schema")
        for row in reader:
            pid = row["Accession"]
            if pid not in sequences:
                raise ValueError(f"Unknown BepiPred accession {pid}")
            values = scores.setdefault(pid, [])
            position = len(values)
            if (
                position >= len(sequences[pid])
                or row["Residue"].upper() != sequences[pid][position].upper()
            ):
                raise ValueError(f"BepiPred residue mismatch for {pid} at {position + 1}")
            score = float(row["BepiPred-3.0 linear epitope score"])
            if not math.isfinite(score) or not 0 <= score <= 1:
                raise ValueError(f"Invalid BepiPred score for {pid}")
            values.append(score)
    require_ids(set(scores), set(sequences), "BepiPred")
    rows = []
    for pid, sequence in sequences.items():
        values = scores[pid]
        if len(values) != len(sequence):
            raise ValueError(
                f"Truncated BepiPred output for {pid}: {len(values)}/{len(sequence)}"
            )
        start = None
        for index, score in enumerate([*values, -1.0]):
            if score >= threshold and start is None:
                start = index
            if score < threshold and start is not None:
                if index - start >= min_length:
                    rows.append(
                        dict(
                            protein_id=pid,
                            start=start + 1,
                            end=index,
                            score=sum(values[start:index]) / (index - start),
                            kind="linear",
                            sequence=sequence[start:index],
                        )
                    )
                start = None
    return rows


def check_tools(config: dict[str, Any], hash_files: bool = False) -> dict[str, Any]:
    """Check the configured local installation without loading model weights."""
    report = {}
    for name in ("signalp", "deeptmhmm", "bepipred"):
        cfg = config["predictors"][name]
        command = cfg.get("command", [])
        errors = []
        if not command or not shutil.which(command[0]):
            errors.append("configure an installed local command")
        if not cfg.get("version"):
            errors.append("record the installed tool/source version")
        files = cfg.get("required_files", [])
        if not files:
            errors.append("declare local model/script files in required_files")
        errors.extend(f"missing file: {p}" for p in files if not Path(p).is_file())
        report[name] = {
            "configuration": cfg,
            "ready": not errors,
            "errors": errors,
            "version": cfg.get("version"),
        }
        if hash_files and not errors:
            report[name]["file_sha256"] = {p: file_digest(p) for p in files}
    return report


def run_predictor(
    config: dict[str, Any],
    name: str,
    fasta: str | Path,
    output: str | Path,
    raw_dir: str | Path,
    threads: int = 1,
) -> None:
    """Each job gets its own cwd, output tree, logs, and checksum manifest.

    Commands are argument arrays, never evaluated by a shell. They must run
    locally with staged models; no cloud runner or model download is invoked.
    """
    cfg = config["predictors"][name]
    status = check_tools(config)[name]
    if not status["ready"]:
        raise ValueError(f"{name} installation incomplete: {status['errors']}")
    fasta, output, raw_dir = (
        Path(fasta).resolve(),
        Path(output).resolve(),
        Path(raw_dir).resolve(),
    )
    sequences = fasta_sequences(fasta)
    raw_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f"{raw_dir.name}-", dir=raw_dir.parent))
    substitutions = {
        "input": str(fasta),
        "output_dir": str(work / "output"),
        "work_dir": str(work),
        "threads": threads,
    }
    command = [str(arg).format(**substitutions) for arg in cfg["command"]]
    # Drivers create output directories themselves; some reject an existing one.
    with (
        open(work / "stdout.log", "w") as stdout,
        open(work / "stderr.log", "w") as stderr,
    ):
        subprocess.run(command, cwd=work, stdout=stdout, stderr=stderr, check=True)
    native = work / cfg["native_output"].format(**substitutions)
    if name == "signalp":
        parsed = io.read_signalp(native)
        require_ids(set(parsed), set(sequences), name)
    elif name == "deeptmhmm":
        parsed = io.read_deeptmhmm(native)
        require_ids(set(parsed), set(sequences), name)
        if any(
            parsed[pid]["sequence"].upper() != seq.upper()
            for pid, seq in sequences.items()
        ):
            raise ValueError("DeepTMHMM sequence mismatch")
    else:
        epitope_cfg = config["cascade"]["epitope"]
        rows = bepipred_intervals(
            native,
            sequences,
            epitope_cfg["min_epitope_score"],
            epitope_cfg.get("min_length", 7),
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = work / "normalized"
    if name == "bepipred":
        io.write_tsv(temp, rows, EPITOPE_COLUMNS)
    else:
        shutil.copyfile(native, temp)
    installation = (
        Path(config.get("workdir", ".")) / config["paths"]["scratch"] / "tools.ready"
    )
    if installation.is_file():
        entry = json.loads(installation.read_text())[name]
        if entry["configuration"] != cfg or not entry["ready"]:
            raise ValueError(
                "Predictor installation manifest is stale; rerun tools_check"
            )
        model_hashes = entry["file_sha256"]
    else:
        model_hashes = {p: file_digest(p) for p in cfg["required_files"]}
    manifest = {
        "tool": name,
        "version": cfg["version"],
        "command": command,
        "input_sha256": file_digest(fasta),
        "native_sha256": file_digest(native),
        "normalized_sha256": file_digest(temp),
        "model_files": model_hashes,
        "n_sequences": len(sequences),
    }
    (work / "manifest.json").write_text(json.dumps(manifest, indent=2))
    # Keep raw outputs on success and on failure; never share another job's cwd.
    shutil.copyfile(temp, output.with_suffix(output.suffix + ".tmp"))
    output.with_suffix(output.suffix + ".tmp").replace(output)
    if raw_dir.exists():
        shutil.rmtree(raw_dir)
    work.rename(raw_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--tool", choices=["signalp", "deeptmhmm", "bepipred"])
    parser.add_argument("--fasta")
    parser.add_argument("--output")
    parser.add_argument("--raw-dir")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    config = io.load_config(args.config)
    if args.check:
        report = check_tools(config, hash_files=bool(args.output))
        if args.output:
            destination = Path(args.output)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        raise SystemExit(0 if all(r["ready"] for r in report.values()) else 1)
    if not all((args.tool, args.fasta, args.output, args.raw_dir)):
        parser.error("--tool, --fasta, --output and --raw-dir are required")
    run_predictor(config, args.tool, args.fasta, args.output, args.raw_dir, args.threads)


if __name__ == "__main__":
    main()
