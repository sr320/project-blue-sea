"""Synthetic fixtures. No real proteome is needed to test the logic."""

from __future__ import annotations

import random
from pathlib import Path

import pytest
import yaml

from bluesea.records import Epitope, HostHit, Ortholog, ProteinRecord, Topology

REPO = Path(__file__).resolve().parents[1]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture
def config() -> dict:
    with open(REPO / "config" / "config.yaml") as handle:
        return yaml.safe_load(handle)


@pytest.fixture
def weights() -> dict:
    with open(REPO / "config" / "weights.yaml") as handle:
        return yaml.safe_load(handle)


def random_sequence(length: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(length))


def make_record(
    protein_id: str = "LSAL_0001",
    length: int = 400,
    seed: int = 1,
    *,
    signal_peptide: bool = True,
    n_transmembrane: int = 0,
    host_identity: float | None = None,
    n_epitopes: int = 4,
    epitope_score: float = 0.8,
    rnai: str = "untested",
    stages: dict[str, float] | None = None,
    ortholog_identity: float | None = None,
    sequence: str | None = None,
) -> ProteinRecord:
    record = ProteinRecord(
        protein_id=protein_id,
        gene_id=protein_id.replace("_P", "_G"),
        species="lsal",
        sequence=sequence or random_sequence(length, seed),
        description="test protein",
        topology=Topology(
            has_signal_peptide=signal_peptide,
            signalp_prob=0.95 if signal_peptide else 0.02,
            n_transmembrane=n_transmembrane,
            extracellular_aa=120 if n_transmembrane else 0,
        ),
        epitopes=[
            Epitope(start=10 + i * 30, end=24 + i * 30, score=epitope_score)
            for i in range(n_epitopes)
        ],
        expression=stages
        if stages is not None
        else {"copepodid": 25.0, "chalimus_I": 30.0, "adult_female": 18.0},
    )
    record.rnai.phenotype = rnai
    if host_identity is not None:
        record.host_hit = HostHit(
            accession="XP_salmon_1",
            identity_pct=host_identity,
            evalue=1e-40,
            query_coverage=0.8,
        )
    if ortholog_identity is not None:
        record.orthologs.append(
            Ortholog(
                species="crog",
                protein_id="CROG_0001",
                identity_pct=ortholog_identity,
                evalue=1e-50,
                reciprocal_best=True,
                epitope_region_identity_pct=ortholog_identity,
            )
        )
    return record


@pytest.fixture
def good_candidate() -> ProteinRecord:
    """Secreted, essential, conserved, salmon-distant: should rank high."""
    return make_record(
        "LSAL_GOOD",
        rnai="lethal",
        host_identity=22.0,
        ortholog_identity=78.0,
        n_epitopes=6,
    )


@pytest.fixture
def intracellular() -> ProteinRecord:
    """No route to host antibody: must fail the exposure filter."""
    return make_record("LSAL_RPL7", signal_peptide=False, seed=7)


@pytest.fixture
def host_like() -> ProteinRecord:
    """Too close to salmon: must fail cross-reactivity."""
    return make_record("SSAL_ACTB", host_identity=94.0, seed=9)
