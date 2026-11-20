"""Extract small control sets and diagnose exposure-filter recovery.

This is a prefilter diagnostic. Pending construct sequence/mapping curation
cannot be cleared by a successful predictor run; use the main pilot/gate for
validation. See docs/validation.md.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bluesea import io
from bluesea.cascade import run_cascade, steps_for
from bluesea.records import ProteinRecord


def load_controls(config: dict) -> list[tuple[str, str, str]]:
    """(protein_id, label, expectation) for every control that resolved."""
    validation = config["validation"]
    controls = []
    for row in io.read_tsv(validation["positive_controls"]):
        controls.append((row["protein_id"], row.get("gene_symbol", ""), "survive"))
    construct_path = validation.get("construct_controls")
    if construct_path:
        registry = json.loads(Path(construct_path).read_text())
        seen = {pid for pid, _, _ in controls}
        for construct in registry["constructs"]:
            for component in construct["louse_components"]:
                for pid in component["mapped_accessions"]:
                    if pid not in seen:
                        controls.append((pid, component["name"], "survive"))
                        seen.add(pid)
    for row in io.read_tsv(validation["negative_controls"]):
        controls.append(
            (row["protein_id"], row.get("class", ""), row["expected_drop_step"])
        )
    return controls


def extract(proteome: Path, wanted: set[str], destination: Path) -> dict[str, str]:
    found = {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w") as handle:
        for protein_id, description, sequence in io.read_fasta(proteome):
            if protein_id not in wanted:
                continue
            found[protein_id] = sequence
            handle.write(f">{protein_id} {description}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    return found


def evaluate(config: dict, controls, sequences, signalp_path, tmhmm_path) -> int:
    signalp = io.read_signalp(signalp_path)
    tmhmm = io.read_deeptmhmm(tmhmm_path)
    cascade_cfg = config["cascade"]
    prefilter_steps = [name for name, _, _ in steps_for("prefilter")]

    print(f"\nPrefilter steps applied: {', '.join(prefilter_steps)}")
    print(
        "'ok' means the control did what it was chosen to do: positives "
        "survive, negatives drop at their expected step.\n"
    )
    header = f"{'protein':18} {'label':14} {'exposure class':22} {'verdict':9} why"
    print(header)
    print("-" * len(header))

    survived_positive = positives = 0
    wrong_step = []
    for protein_id, label, expectation in controls:
        sequence = sequences.get(protein_id)
        if sequence is None:
            print(f"{protein_id:18} {label:14} {'NOT IN PROTEOME':22} {'ERROR':9}")
            wrong_step.append((protein_id, expectation, "missing"))
            continue
        record = ProteinRecord(
            protein_id=protein_id,
            gene_id=protein_id,
            species="",
            sequence=sequence,
            topology=io.topology_from_tools(protein_id, signalp, tmhmm),
        )
        run_cascade(record, cascade_cfg, phase="prefilter")
        dropped = record.dropped_at
        reason = next(
            (v["reason"] for v in record.verdicts if not v["passed"]), "passed all"
        )

        if expectation == "survive":
            positives += 1
            ok = dropped is None
            survived_positive += ok
        else:
            # A negative dropped by the wrong rule proves nothing about the
            # rule it was chosen to test.
            ok = dropped == expectation
            if dropped is not None and dropped != expectation:
                wrong_step.append((protein_id, expectation, dropped))
            if dropped is None and expectation in prefilter_steps:
                wrong_step.append((protein_id, expectation, "survived"))

        print(
            f"{protein_id:18} {label:14} {record.topology.exposure_class:22} "
            f"{'ok' if ok else 'NO':9} {reason[:44]}"
        )

    print()
    needed = config["validation"]["min_positive_controls_in_top_tier"]
    # positive_control_report passes when recovered/n >= needed, so the bar is
    # the smallest integer count that clears it -- a ceiling, not a trunc.
    required = math.ceil(needed * positives) if positives else 0
    print(f"positives surviving the prefilter: {survived_positive}/{positives}")
    print(
        f"gate needs {needed:.0%} of {positives} in the top tier, so a dropped "
        f"positive is unrecoverable"
    )
    for protein_id, expected, actual in wrong_step:
        print(f"  note: {protein_id} expected to drop at {expected}, dropped at {actual}")

    if survived_positive < required:
        # The table above is the evidence for what follows; keep them in order
        # when stdout is a pipe and stderr is not.
        sys.stdout.flush()
        print(
            f"\nThe gate cannot pass: {positives - survived_positive} positive "
            f"control(s) never reach the ranking, leaving at most "
            f"{survived_positive}/{positives} recoverable against a "
            f"{needed:.0%} threshold.\n"
            f"This is a finding about the rubric, not a crash. The exposure "
            f"filter and the published antigens disagree; see README 'Known "
            f"limits' on the concealed-antigen mechanism.",
            file=sys.stderr,
        )
        return 1
    path = config["validation"].get("construct_controls")
    pending = not positives or (
        path
        and any(
            c["curation_status"] != "verified"
            for c in json.loads(Path(path).read_text())["constructs"]
        )
    )
    if pending or wrong_step:
        print(
            "\nPrefilter diagnostic incomplete: inspect missing/failed controls "
            "and construct curation."
        )
        return 1
    print("\nPrefilter diagnostic complete; the full pilot benchmark is still required.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--species", default="lsal")
    parser.add_argument("--outdir", default="scratch/smoke")
    args = parser.parse_args(argv)

    config = io.load_config(args.config)
    root = Path(config.get("workdir", "."))
    reference = root / config["paths"]["reference"]
    proteome = reference / config["species"][args.species]["proteome"]
    outdir = Path(args.outdir)
    fasta = outdir / "controls.faa"

    if not proteome.exists():
        print(
            f"proteome not found: {proteome}\nRun `snakemake provision` first.",
            file=sys.stderr,
        )
        return 2

    controls = load_controls(config)
    wanted = {protein_id for protein_id, _, _ in controls}
    sequences = extract(proteome, wanted, fasta)
    missing = wanted - set(sequences)

    print(f"{len(sequences)}/{len(wanted)} controls extracted -> {fasta}")
    for protein_id in sorted(missing):
        print(f"  MISSING from the proteome: {protein_id}", file=sys.stderr)

    signalp_path = outdir / "signalp.tsv"
    tmhmm_path = outdir / "deeptmhmm.gff3"
    if not (signalp_path.exists() and tmhmm_path.exists()):
        print(
            f"\nNow run the two predictors on {fasta}, writing:\n"
            f"  {signalp_path}\n  {tmhmm_path}\n\n"
            "Use the configured local predictor adapter with separate raw directories;\n"
            f"see docs/tools.md. Cloud runners are not configured for compute nodes.\n\n"
            f"Then run this script again to get the verdicts."
        )
        return 0

    return evaluate(config, controls, sequences, signalp_path, tmhmm_path)


if __name__ == "__main__":
    raise SystemExit(main())
