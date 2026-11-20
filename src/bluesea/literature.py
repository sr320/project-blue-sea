"""Optional model-assisted triage over curated, checksummed source excerpts.

No sources means no model call. Benchmark scoring excludes this evidence
entirely; prompt-based redaction is retained only as a presentation helper.
Returned identifiers must refer to supplied sources, whose claims still need
human review. Source-ID validation cannot establish entailment on its own.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

from .llm import CachedClient
from .records import LiteratureEvidence, ProteinRecord

STAGE = "literature"

SYSTEM = """\
You are assisting a reverse-vaccinology screen of the salmon louse proteome.
Use ONLY the supplied source excerpts. Do not supply facts or identifiers
from model memory. Cite only source IDs present in the supplied evidence.
For one protein, judge how much the published literature supports it as a
vaccine antigen candidate against a parasitic arthropod.

Score 0.0 to 1.0:
  0.9-1.0  an ortholog gave durable, replicated protection in a published trial
  0.7-0.9  an ortholog gave a measurable but partial or short-lived effect
  0.4-0.6  the protein family is functionally characterised but untested as an
           antigen, or published evidence is mixed
  0.2-0.4  little published work; judged on function alone
  0.0-0.2  published evidence argues against it (host-similar, poorly exposed,
           redundant, or a failed trial)

An uncharacterised protein is not a low-scoring protein. Score it in the
0.2-0.4 band and say the evidence is thin, rather than penalising it for being
unstudied -- this screen exists partly to find candidates nobody has looked at.

Do not infer results that the sources do not state. If you do not find
relevant work, say so and score accordingly. Cite every claim with a PMID,
DOI or accession. Reply with a single JSON object and nothing else:

{"score": <float>, "justification": "<2-4 sentences>",
 "citations": ["<PMID/DOI>", ...], "evidence_strength": "strong|moderate|weak|none"}
"""

PROMPT = """\
Protein: {protein_id} ({gene_id}), {species}
Annotation: {description}
Length: {length} aa
Predicted localisation: {exposure_class}
Closest salmon homolog: {host_summary}
Ortholog in the other louse species: {ortholog_summary}
RNAi knockdown phenotype: {rnai_summary}
Expressed in: {stages}

Sequence:
{sequence}

{blinding_note}
Assess the published evidence for this protein and its orthologs in parasitic
arthropods, then return the JSON object.
"""

BLINDING_NOTE = """\
IMPORTANT -- this is a blinded benchmark run. Disregard and do not search for
any vaccine, immunization, antigen, adjuvant or challenge-trial literature
about this protein. Judge it only on its molecular function, localisation and
essentiality. If you already know this protein was used as a vaccine antigen,
that knowledge must not influence the score.
"""


def redact(text: str, terms: list[str]) -> str:
    """Strip blinded terms from text handed to the model.

    An annotation reading "vaccine antigen Bm86 homolog" would defeat the
    blinding on its own, so the terms go from the prompt as well as from the
    retrieval query.
    """
    for term in terms:
        pattern = rf"\b{re.escape(term)}\w*\b"
        text = re.sub(pattern, "[redacted]", text, flags=re.IGNORECASE)
    return text


def build_prompt(record: ProteinRecord, blinded: bool, blinded_terms: list[str]) -> str:
    host = record.host_hit
    ortholog = record.best_ortholog()
    description = record.description
    if blinded:
        description = redact(description, blinded_terms)

    return PROMPT.format(
        protein_id=record.protein_id,
        gene_id=record.gene_id,
        species=record.species,
        description=description or "unannotated / hypothetical protein",
        length=record.length,
        exposure_class=record.topology.exposure_class,
        host_summary=(
            f"{host.accession} at {host.identity_pct:.1f}% identity"
            if host
            else "none detected"
        ),
        ortholog_summary=(
            f"{ortholog.protein_id} ({ortholog.species}) at "
            f"{ortholog.identity_pct:.1f}% identity"
            + (", reciprocal best hit" if ortholog.reciprocal_best else "")
            if ortholog
            else "none detected"
        ),
        rnai_summary=record.rnai.phenotype or "untested",
        stages=", ".join(sorted(record.expression)) or "no expression data",
        sequence=_wrap(record.sequence),
        blinding_note=BLINDING_NOTE if blinded else "",
    )


def _wrap(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def assess(
    record: ProteinRecord,
    client: CachedClient,
    blinded: bool = False,
    blinded_terms: list[str] | None = None,
    sources: list[dict[str, Any]] | None = None,
) -> ProteinRecord:
    record.literature = None
    if blinded:
        record.flag("literature_excluded_benchmark")
        return record
    if not sources:
        record.flag("literature_sources_missing")
        return record
    prompt = build_prompt(record, False, blinded_terms or [])
    prompt += "\n\nSource excerpts (cite source_id only):\n" + json.dumps(sources)
    payload = client.complete_json(STAGE, SYSTEM, prompt)

    if not payload:
        # Dry run, or a response that would not parse. Either way the record
        # keeps no literature evidence and the scorer falls back to its
        # configured neutral value with a flag, rather than a fabricated score.
        record.flag("literature_not_assessed")
        return record

    try:
        score = float(payload["score"])
    except (KeyError, TypeError, ValueError):
        record.flag("literature_response_malformed")
        return record

    citations = payload.get("citations")
    allowed = {source["source_id"] for source in sources}
    if (
        not math.isfinite(score)
        or not 0 <= score <= 1
        or not isinstance(citations, list)
        or not citations
        or any(not isinstance(c, str) or c not in allowed for c in citations)
    ):
        record.flag("literature_response_unsupported")
        return record
    record.literature = LiteratureEvidence(
        score=max(0.0, min(1.0, score)),
        justification=str(payload.get("justification", "")).strip(),
        citations=[str(c) for c in payload.get("citations", [])],
        model=client.model,
        blinded=blinded,
    )
    if not record.literature.citations:
        record.flag("literature_score_uncited")
    if str(payload.get("evidence_strength", "")).lower() in ("weak", "none"):
        record.flag("literature_evidence_thin")
    return record


def assess_all(
    records: list[ProteinRecord],
    client: CachedClient,
    blinded: bool = False,
    blinded_terms: list[str] | None = None,
    sources_by_protein: dict[str, list[dict[str, Any]]] | None = None,
) -> list[ProteinRecord]:
    """Triage the whole shortlist, stopping cleanly if the ceiling is hit.

    A halted stage leaves the records it already assessed intact and flags the
    rest as unassessed, so the run still produces a ranked list rather than
    nothing.
    """
    from .llm import BudgetExceeded

    for index, record in enumerate(records):
        try:
            assess(
                record,
                client,
                blinded=blinded,
                blinded_terms=blinded_terms,
                sources=(sources_by_protein or {}).get(record.protein_id),
            )
        except BudgetExceeded:
            for remaining in records[index:]:
                remaining.flag("literature_skipped_budget")
            break
    return records


def load_sources(paths: list[str]) -> dict[str, list[dict[str, Any]]]:
    """Read curated excerpts backed by checksummed local source snapshots.

    Retrieval happens before this stage; the packet and source snapshot are
    independently inspectable and are workflow dependencies.
    """
    from .provenance import file_digest

    result: dict[str, list[dict[str, Any]]] = {}
    for path in paths:
        for source in json.loads(Path(path).read_text()):
            snapshot = Path(path).parent / source["snapshot_path"]
            if (
                file_digest(snapshot) != source["snapshot_sha256"]
                or not source["excerpt"]
                or source["excerpt"] not in snapshot.read_text()
            ):
                raise ValueError(
                    f"Source excerpt/snapshot mismatch: {source['source_id']}"
                )
            if not source["url"].startswith("https://") or not source["source_id"]:
                raise ValueError("Source needs a stable identifier and HTTPS source URL")
            for pid in source["protein_ids"]:
                result.setdefault(pid, []).append(
                    {
                        k: source[k]
                        for k in ("source_id", "url", "excerpt", "snapshot_sha256")
                    }
                )
    return result
