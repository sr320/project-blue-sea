"""Run the official BepiPred Python API with local ESM weights, without downloads.

The BP3 checkout and its model dependencies are provisioned separately. Long
sequences are tiled with overlap; overlapping residue scores are averaged.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .predictors import fasta_sequences


def tiles(sequences: dict[str, str], size: int = 1000, overlap: int = 100):
    if not 0 <= overlap < size:
        raise ValueError("Invalid tile overlap")
    result = {}
    for index, (pid, sequence) in enumerate(sequences.items()):
        for start in range(0, len(sequence), size - overlap):
            tile_id = f"tile_{index}_{start}"
            result[tile_id] = (pid, start, sequence[start : start + size])
            if start + size >= len(sequence):
                break
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--esm-model", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    model = args.esm_model.resolve()
    regression = model.with_name(model.stem + "-contact-regression.pt")
    if not model.is_file() or not regression.is_file():
        raise ValueError("Stage ESM model and contact-regression weights before running")
    sys.path.insert(0, str(args.repo.resolve()))
    from bp3 import bepipred3

    sequences = fasta_sequences(args.fasta)
    chunks = tiles(sequences)
    args.output_dir.mkdir(parents=True)
    tiled = args.output_dir / "tiles.faa"
    tiled.write_text("".join(f">{key}\n{seq}\n" for key, (_, _, seq) in chunks.items()))
    antigens = bepipred3.Antigens(
        tiled, args.output_dir / "embeddings", run_esm_model_local=str(model)
    )
    predictor = bepipred3.BP3EnsemblePredict(
        antigens, rolling_window_size=9, top_pred_pct=0.2
    )
    predictor.run_bp3_ensemble()
    native_dir = args.output_dir / "tiles"
    predictor.create_csvfile(native_dir)
    sums = {pid: [0.0] * len(seq) for pid, seq in sequences.items()}
    counts = {pid: [0] * len(seq) for pid, seq in sequences.items()}
    positions = dict.fromkeys(chunks, 0)
    with open(native_dir / "raw_output.csv") as handle:
        for row in csv.DictReader(handle):
            tile = row["Accession"]
            pid, start, seq = chunks[tile]
            offset = positions[tile]
            if offset >= len(seq) or row["Residue"].upper() != seq[offset].upper():
                raise ValueError("Tiled BepiPred residue mismatch")
            position = start + offset
            sums[pid][position] += float(row["BepiPred-3.0 linear epitope score"])
            counts[pid][position] += 1
            positions[tile] += 1
    if any(positions[t] != len(seq) for t, (_, _, seq) in chunks.items()):
        raise ValueError("Truncated BepiPred tile output")
    with open(args.output_dir / "raw_output.csv", "w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Accession", "Residue", "BepiPred-3.0 linear epitope score"])
        for pid, seq in sequences.items():
            for index, residue in enumerate(seq):
                writer.writerow([pid, residue, sums[pid][index] / counts[pid][index]])


if __name__ == "__main__":
    main()
