"""Literature triage: the one Claude-judged axis.

What the model is asked for is a judgment no tool emits as a number -- has
this protein, or an ortholog of it, been tried as an antigen in any arthropod,
and what happened. The score comes back with a written justification and
citations so a human can check the reasoning rather than trust the number.

Blinding matters for the benchmark. If the model can read that P0 was used in
a vaccine trial, a high score for P0 proves nothing about the rubric. During
positive-control runs the retrieval excludes vaccine and immunization terms,
so the ranking has to rest on the deterministic axes alone.
"""

from __future__ import annotations

import re

from .llm import CachedClient
from .records import LiteratureEvidence, ProteinRecord

STAGE = "literature"

SYSTEM = """\
You are assisting a reverse-vaccinology screen of the salmon louse proteome.
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


def build_prompt(
    record: ProteinRecord, blinded: bool, blinded_terms: list[str]
) -> str:
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
) -> ProteinRecord:
    prompt = build_prompt(record, blinded, blinded_terms or [])
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
) -> list[ProteinRecord]:
    """Triage the whole shortlist, stopping cleanly if the ceiling is hit.

    A halted stage leaves the records it already assessed intact and flags the
    rest as unassessed, so the run still produces a ranked list rather than
    nothing.
    """
    from .llm import BudgetExceeded

    for index, record in enumerate(records):
        try:
            assess(record, client, blinded=blinded, blinded_terms=blinded_terms)
        except BudgetExceeded:
            for remaining in records[index:]:
                remaining.flag("literature_skipped_budget")
            break
    return records
