"""Record types carried through the cascade.

One `ProteinRecord` per protein, accumulating evidence as it passes each
stage. Evidence that was never gathered is `None`, which is distinct from
evidence that was gathered and came back negative -- the scoring code relies
on that distinction (an untested RNAi target is not a dispensable one).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class Epitope:
    start: int
    end: int
    score: float
    kind: str = "linear"          # linear | conformational
    sequence: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class HostHit:
    """Best alignment of this protein against the salmon proteome."""

    accession: str
    identity_pct: float
    evalue: float
    query_coverage: float
    description: str = ""


@dataclass
class Ortholog:
    species: str
    protein_id: str
    identity_pct: float
    evalue: float
    reciprocal_best: bool = False
    epitope_region_identity_pct: float | None = None


@dataclass
class Topology:
    has_signal_peptide: bool = False
    signalp_prob: float = 0.0
    n_transmembrane: int = 0
    gpi_anchored: bool = False
    extracellular_aa: int = 0
    nonclassical_secretion: bool = False
    secretome_confirmed: bool = False
    gut_lumen_facing: bool = False

    @property
    def exposure_class(self) -> str:
        """The strongest exposure claim the evidence supports.

        Ordered by how well it predicts antibody contact, not by how confident
        the prediction is -- experimental secretome data outranks everything
        because it is the only non-predicted entry here.
        """
        if self.secretome_confirmed:
            return "secretome_confirmed"
        if self.gut_lumen_facing:
            return "gut_lumen_facing"
        if self.has_signal_peptide and self.n_transmembrane == 0:
            return "signal_peptide"
        if self.gpi_anchored:
            return "gpi_anchored"
        if self.n_transmembrane == 1:
            return "type_I_membrane"
        if self.n_transmembrane > 1:
            return "multipass_membrane"
        if self.nonclassical_secretion:
            return "nonclassical_secretion"
        return "intracellular"


@dataclass
class RNAiEvidence:
    phenotype: str = "untested"
    source: str = ""
    note: str = ""


@dataclass
class LiteratureEvidence:
    score: float
    justification: str
    citations: list[str] = field(default_factory=list)
    model: str = ""
    blinded: bool = False


@dataclass
class ProteinRecord:
    protein_id: str
    gene_id: str
    species: str
    sequence: str
    description: str = ""
    # The readable gene name where the annotation has one, kept beside the
    # stable gene_id because evidence tables may be keyed on either.
    gene_symbol: str = ""

    topology: Topology = field(default_factory=Topology)
    epitopes: list[Epitope] = field(default_factory=list)
    host_hit: HostHit | None = None
    orthologs: list[Ortholog] = field(default_factory=list)
    rnai: RNAiEvidence = field(default_factory=RNAiEvidence)
    literature: LiteratureEvidence | None = None

    # gene_id -> TPM per life stage; absent stages mean "not assayed".
    expression: dict[str, float] = field(default_factory=dict)

    n_disulfide_pairs: int = 0
    n_glycosylation_sites: int = 0
    low_complexity_fraction: float = 0.0

    # Populated by the cascade and the scorer.
    verdicts: list[dict[str, Any]] = field(default_factory=list)
    scores: dict[str, float] = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)
    total_score: float = 0.0

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def dropped_at(self) -> str | None:
        for verdict in self.verdicts:
            if not verdict["passed"]:
                return verdict["step"]
        return None

    @property
    def survived(self) -> bool:
        return self.dropped_at is None

    def record_verdict(self, step: str, passed: bool, reason: str, **detail: Any) -> None:
        self.verdicts.append(
            {"step": step, "passed": passed, "reason": reason, "detail": detail}
        )

    def flag(self, name: str) -> None:
        if name not in self.flags:
            self.flags.append(name)

    def best_ortholog(self, species: str | None = None) -> Ortholog | None:
        candidates = [
            o for o in self.orthologs if species is None or o.species == species
        ]
        if not candidates:
            return None
        return max(candidates, key=lambda o: (o.reciprocal_best, o.identity_pct))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, sort_keys=True)


def from_dict(payload: dict[str, Any]) -> ProteinRecord:
    """Rebuild a record from its JSON form, so stages can hand off via disk."""
    payload = dict(payload)
    payload["topology"] = Topology(**payload.get("topology", {}) or {})
    payload["epitopes"] = [Epitope(**e) for e in payload.get("epitopes", [])]
    payload["orthologs"] = [Ortholog(**o) for o in payload.get("orthologs", [])]
    payload["rnai"] = RNAiEvidence(**(payload.get("rnai") or {}))
    host_hit = payload.get("host_hit")
    payload["host_hit"] = HostHit(**host_hit) if host_hit else None
    literature = payload.get("literature")
    payload["literature"] = LiteratureEvidence(**literature) if literature else None
    return ProteinRecord(**payload)
