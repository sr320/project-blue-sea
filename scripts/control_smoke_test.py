"""Ask whether the gate can pass, for the cost of nine proteins.

    python scripts/control_smoke_test.py            # step 1: extract
    # ...run the two predictors on the FASTA it wrote...
    python scripts/control_smoke_test.py            # step 2: evaluate

The positive controls are antigens with measurable effects in published
trials, so the rubric has to rank them highly or it is not measuring what it
claims to. But three of them -- P0, akirin and peroxiredoxin -- are
intracellular proteins, and the exposure filter exists to drop proteins with
no route to host antibody. If it drops them they can never reach the top tier,
`recovered_fraction` is capped below its threshold, and the gate fails no
matter how good the ranking is.

That outcome is decided by the first three cascade steps, which need only
SignalP and DeepTMHMM. Running them on the controls alone answers the question
in minutes instead of an allocation. A failure here is not a bug: it is the
rubric telling you its central filter is incompatible with the only antigens
known to work in this system, which is worth knowing before the full run
rather than after.

The verdicts come from `bluesea.cascade` itself, in the same phase the
workflow uses, so this cannot drift from what the pipeline would decide.
"""

from __future__ import annotations

import argparse
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
    print("\nEvery positive control survives the prefilter; the gate is reachable.")
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
            f"  signalp6 --fastafile {fasta} --organism eukarya --format none \\\n"
            f"      --output_dir {outdir}\n"
            f"  mv {outdir}/prediction_results.txt {signalp_path}\n\n"
            f"  biolib run DTU/DeepTMHMM --fasta {fasta} --out-dir {outdir}\n"
            f"  mv {outdir}/predicted_topologies.3line {tmhmm_path}\n\n"
            f"Then run this script again to get the verdicts."
        )
        return 0

    return evaluate(config, controls, sequences, signalp_path, tmhmm_path)


if __name__ == "__main__":
    raise SystemExit(main())
