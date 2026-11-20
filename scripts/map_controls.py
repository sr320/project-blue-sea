"""Establish control accessions by alignment, not by annotation text.

    python scripts/map_controls.py resources/reference/lsal/proteins.faa

Annotation text does not work for this. RefSeq names P0 "large ribosomal
subunit protein uL10" and L7 "uL30", so neither is findable by the name the
literature uses; meanwhile two unrelated proteins are both called
"peroxiredoxin 2" and only alignment says which one the published antigen is.

For each control this fetches a UniProt reference sequence, finds the best
local alignment against the proteome, and prints identity and coverage so the
assignment can be checked or rejected. It decides nothing on its own: a hit
with low coverage is a non-result and is reported as such.

Needs biopython:  pip install -e '.[curation]'
"""

from __future__ import annotations

import argparse
import sys
import urllib.request

from Bio.Align import PairwiseAligner, substitution_matrices

# (label, UniProt accession). Chosen from UniProt entries for Lepeophtheirus
# salmonis (taxon 72036), except Bm86, which has no louse entry and is
# referenced against the canonical Rhipicephalus microplus antigen.
CONTROLS = [
    ("P0", "C1BTD7"),
    ("my32/akirin", "D3PIL3"),
    ("prx2", "C1BU78"),
    ("bm86-like", "P20736"),
    ("RPL7", "C1BTB6"),
    ("EF1A", "A7TZ68"),
    ("histone H3", "C1BSB6"),
    ("proteasome b1", "C1BVW5"),
]

K = 5
CANDIDATES = 25
# Below this, a best hit is reported but should not be treated as an
# assignment: a short perfect match to a shared domain is not orthology.
MIN_COVERAGE = 0.50
MIN_IDENTITY = 50.0

MATRIX = substitution_matrices.load("BLOSUM62")
ALPHABET = set(MATRIX.alphabet)


def read_fasta(path: str):
    header, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def clean(sequence: str) -> str:
    """BLOSUM62 has no column for '*' or unusual letters; fold them to X."""
    return "".join(c if c in ALPHABET else "X" for c in sequence.upper())


def kmers(sequence: str) -> set[str]:
    return {sequence[i : i + K] for i in range(len(sequence) - K + 1)}


def fetch_reference(accession: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=60) as response:
        body = response.read().decode()
    return clean("".join(body.split("\n")[1:]))


def align(aligner: PairwiseAligner, reference: str, subject: str):
    """Identity across the aligned region, and the fraction of the reference
    that region covers. Coverage is what separates an ortholog from a protein
    that merely shares a domain."""
    alignment = aligner.align(reference, subject)[0]
    top, bottom = alignment[0], alignment[1]
    # Both rows of an alignment are the same length by construction; strict
    # makes a biopython change that broke that assumption fail loudly.
    pairs = [
        (a, b) for a, b in zip(top, bottom, strict=True) if a != "-" and b != "-"
    ]
    if not pairs:
        return 0.0, 0.0
    identity = 100.0 * sum(1 for a, b in pairs if a == b) / len(pairs)
    return identity, len(pairs) / len(reference)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proteome", help="target proteome FASTA")
    parser.add_argument("--top", type=int, default=3, help="hits to show each")
    args = parser.parse_args(argv)

    entries = [(h.split()[0], h, clean(s)) for h, s in read_fasta(args.proteome)]
    sequences = {ident: seq for ident, _, seq in entries}
    descriptions = {
        ident: " ".join(header.split()[1:]).split(" [")[0]
        for ident, header, _ in entries
    }

    references = {}
    for label, accession in CONTROLS:
        print(f"fetching {accession} ({label})", file=sys.stderr)
        references[label] = (accession, fetch_reference(accession))

    # One pass over the proteome, scoring every protein against the union of
    # reference k-mers, so the expensive alignment only sees plausible hits.
    pool: set[str] = set()
    for _, sequence in references.values():
        pool |= kmers(sequence)
    shared = {
        ident: kmers(seq) & pool for ident, _, seq in entries if kmers(seq) & pool
    }

    aligner = PairwiseAligner()
    aligner.substitution_matrix = MATRIX
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"

    for label, (accession, reference) in references.items():
        reference_kmers = kmers(reference)
        shortlist = sorted(shared.items(), key=lambda kv: -len(kv[1] & reference_kmers))
        scored = []
        for ident, _ in shortlist[:CANDIDATES]:
            identity, coverage = align(aligner, reference, sequences[ident])
            scored.append((identity, coverage, ident))
        scored.sort(key=lambda row: -(row[0] / 100.0 * row[1]))

        print(f"\n### {label}  (reference {accession}, {len(reference)} aa)")
        if not scored:
            print("    no candidate shared a k-mer with the reference")
            continue
        for identity, coverage, ident in scored[: args.top]:
            confident = identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE
            print(
                f"    {'OK ' if confident else '   '} {ident:18} "
                f"{identity:5.1f}% id  cov {coverage:4.0%}  "
                f"{descriptions.get(ident, '')[:50]}"
            )
        best = scored[0]
        if not (best[0] >= MIN_IDENTITY and best[1] >= MIN_COVERAGE):
            print(
                f"    -> UNRESOLVED: best hit is {best[0]:.1f}% over {best[1]:.0%} "
                f"of the reference; do not record this as an assignment"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
