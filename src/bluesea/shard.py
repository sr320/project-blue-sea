"""Deterministic FASTA sharding that handles fewer proteins than requested shards."""

import argparse
from pathlib import Path

from . import io


def shard(source, destination, n):
    records = list(io.read_fasta(source))
    if n < 1 or not records:
        raise ValueError("Sharding needs a positive count and nonempty FASTA")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for old in destination.glob("part_*.faa"):
        old.unlink()
    n = min(n, len(records))
    for index in range(n):
        (destination / f"part_{index + 1:03d}.faa").write_text(
            "".join(f">{pid} {desc}\n{seq}\n" for pid, desc, seq in records[index::n])
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n", type=int, required=True)
    args = parser.parse_args()
    shard(args.input, args.output_dir, args.n)


if __name__ == "__main__":
    main()
