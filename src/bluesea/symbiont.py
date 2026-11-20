"""Symbiont track: a question before it is a screen.

The logic is sound -- if the louse depends on a bacterial symbiont, antibody
against that symbiont's surface proteins could damage the louse indirectly,
and a bacterial proteome is far more divergent from salmon than the louse
proteome is. But the precedent is in insects (*Wigglesworthia* in tsetse,
*Buchnera* in aphids), and no obligate endosymbiont has been established for
either louse species.

So this module does not screen. It evaluates one testable claim from
metagenome-assembled genomes: **does a consistent, vertically transmitted
symbiont exist at all?** Only if the signature appears do the MAG's surface
proteins enter the main cascade, with host dissimilarity relaxed.

The likely outcome is a negative, which is worth writing up on its own. The
track is capped at 15% of API spend and gated here so it cannot quietly
consume the budget.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class MAG:
    """One metagenome-assembled genome from a louse microbiome dataset."""

    mag_id: str
    taxon: str
    genome_size_mb: float
    completeness: float
    contamination: float
    # Population/site -> fraction of individuals the MAG was recovered from.
    prevalence_by_population: dict[str, float] = field(default_factory=dict)
    # Biosynthetic pathways the genome has lost (amino acids, cofactors).
    lost_pathways: list[str] = field(default_factory=list)
    host_species: list[str] = field(default_factory=list)
    # Correlation between symbiont and host phylogenies, if computed.
    phylogenetic_congruence: float | None = None


def is_environmental(mag: MAG, blocklist: list[str]) -> bool:
    """Exclude taxa routinely acquired from seawater or salmon mucus.

    The main failure mode of this whole track is mistaking environmental
    contamination for symbiosis, so presence alone never qualifies.
    """
    return any(genus.lower() in mag.taxon.lower() for genus in blocklist)


def symbiosis_signature(mag: MAG, cfg: dict[str, Any]) -> dict[str, Any]:
    """Score a MAG against the signature of obligate symbiosis.

    Four marks, all of which the insect endosymbionts show and transient gut
    flora do not: present in every population regardless of host and site,
    strongly reduced genome, lost biosynthetic pathways the host would supply,
    and phylogeny tracking the louse's own.
    """
    populations = [
        population
        for population, prevalence in mag.prevalence_by_population.items()
        if prevalence >= cfg["min_prevalence"]
    ]
    checks = {
        "mag_quality": mag.completeness >= cfg.get("min_completeness", 90)
        and mag.contamination <= cfg.get("max_contamination", 5),
        "present_across_populations": len(populations) >= cfg["min_populations_present"],
        "reduced_genome": mag.genome_size_mb <= cfg["max_genome_size_mb"],
        "pathway_loss": (
            bool(mag.lost_pathways) if cfg.get("require_pathway_loss", True) else True
        ),
        "phylogenetic_congruence": (
            mag.phylogenetic_congruence is not None and mag.phylogenetic_congruence >= 0.7
        ),
    }
    # Two-species presence is the strongest form of the evidence: a bacterium
    # in L. salmonis in the North Atlantic and in C. rogercresseyi off Chile is
    # hard to explain as shared seawater or shared host mucus.
    both_species = len(set(mag.host_species)) >= 2

    return {
        "mag_id": mag.mag_id,
        "taxon": mag.taxon,
        "checks": checks,
        "n_populations": len(populations),
        "both_louse_species": both_species,
        "score": sum(checks.values()) + (1 if both_species else 0),
        "max_score": len(checks) + 1,
    }


def gate(mags: list[MAG], cfg: dict[str, Any]) -> dict[str, Any]:
    """Go/no-go after MAG assembly.

    `proceed` False stops the track and the negative result gets written up.
    That is the expected outcome and is not a pipeline failure.
    """
    blocklist = cfg.get("environmental_blocklist", [])
    signature_cfg = cfg["signature"]

    assessments = []
    excluded = []
    for mag in mags:
        if is_environmental(mag, blocklist):
            excluded.append({"mag_id": mag.mag_id, "taxon": mag.taxon})
            continue
        assessments.append(symbiosis_signature(mag, signature_cfg))

    qualifying = [
        assessment for assessment in assessments if all(assessment["checks"].values())
    ]
    return {
        "proceed": bool(qualifying),
        "n_mags": len(mags),
        "n_excluded_environmental": len(excluded),
        "excluded": excluded,
        "assessments": sorted(assessments, key=lambda a: a["score"], reverse=True),
        "qualifying": qualifying,
        "conclusion": (
            f"{len(qualifying)} MAG(s) show the obligate-symbiosis signature; "
            "their surface proteins proceed to the cascade."
            if qualifying
            else "No MAG shows the obligate-symbiosis signature. The track "
            "stops here; the negative result is the output."
        ),
    }
