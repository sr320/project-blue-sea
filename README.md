# project-blue-sea

An automated reverse-vaccinology screen of the salmon louse proteome, looking
for **protein vaccine antigens** against *Lepeophtheirus salmonis* and
*Caligus rogercresseyi*.

It produces a ranked, cited candidate dossier: predicted localisation and
topology, B-cell epitope map, expression across life stages, salmon
cross-reactivity risk, orthology to antigens tried in other parasitic
arthropods, and a literature summary. **It produces hypotheses for the bench.
It does not predict efficacy, and nothing it outputs substitutes for a
challenge trial.**

## Why

[Bm86](https://pmc.ncbi.nlm.nih.gov/articles/PMC10295593) remains the only
clearly effective commercial ectoparasite vaccine, working on the
concealed-antigen principle: ingested antibody lyses tick gut epithelium.
Sea lice work has followed that logic with mixed results — the best current
data, from [chimeric P0 constructs](https://pmc.ncbi.nlm.nih.gov/articles/PMC12133460),
report 88–91% lice reduction at two to three months decaying to 20% and ~4%
by eleven months.

Three things follow, and the rubric encodes them:

1. **Durability is the binding constraint, not initial titer.** Essentiality
   and constitutive expression outweigh predicted immunogenicity.
2. **Antibody titer does not predict protection in salmon.** Functional
   essentiality has to be an independent axis.
3. **The search has been narrow.** Most tested antigens were picked because
   they were already characterised in ticks. The 14,014-gene proteome has
   never been screened end to end against explicit criteria.

## Layout

```
config/config.yaml      cascade thresholds, paths, budgets
config/weights.yaml     scoring weights (version-controlled on purpose)
src/bluesea/
  records.py            the record carried through every stage
  io.py                 parsers for SignalP, DeepTMHMM, BLAST, epitopes
  cascade.py            the filter cascade + per-protein audit
  orthology.py          reciprocal best hits, cross-species intersection
  scoring.py            the eight weighted criteria
  llm.py                cached, budgeted API client with dry-run
  literature.py         the one Claude-judged axis, blindable
  dossier.py            one-page candidate write-ups
  validation.py         positive/negative controls, robustness, the gate
  symbiont.py           the exploratory symbiont arm and its go/no-go
workflow/Snakefile      the three-stage workflow
resources/controls/     positive and negative control sets
```

## Running it

```bash
pip install -e '.[dev]'
pytest
```

On the cluster, run the pre-flight checks first — they confirm the account,
the bind path, and whether compute nodes really lack outbound access, which
is the assumption the whole architecture rests on:

```bash
bash scripts/hyak_preflight.sh --account srlab --partition compute
```

See [docs/hyak.md](docs/hyak.md) for what each check means.

The workflow is three stages, split because Hyak compute nodes cannot be
assumed to have outbound network access:

```bash
# Stage 1 and 3: deterministic, no network, SLURM array jobs
snakemake --profile workflow/profiles/slurm --use-apptainer all

# Stage 2: API calls, from a node with outbound access. Dry run by default.
BLUESEA_LIVE=1 bluesea literature --species lsal --live
```

Run the gate before trusting anything:

```bash
snakemake gate   # non-zero exit means rebuild the rubric, not tune it
```

## The cascade

Order is cost-ascending, so the cheap filters run first and the expensive
predictors only ever see survivors.

| Step | What it asks | Typical attrition |
| --- | --- | --- |
| length | 50–2000 aa | small |
| low_complexity | not repeat-dominated | small |
| **exposure** | **is there any route to host antibody?** | **~8 in 9** |
| host_similarity | <40% identity to any salmon protein (hard) | moderate |
| epitope | ≥2 predicted B-cell epitopes | moderate |
| expression | expressed in ≥2 parasitic stages | small |
| transmembrane | expressibility warning only, never excludes | none |

The exposure filter does most of the work and is also the least certain:
predicted secretion is not confirmed secretion. Every verdict is logged per
protein to `logs/<species>/cascade_audit.jsonl`, so a dropped candidate traces
to the exact rule that dropped it.

The cascade runs once per species; the two shortlists are then intersected.
A protein conserved across two genera is a candidate for a vaccine serving
both industries, and that conservation is independent evidence of functional
essentiality.

## Scoring

Seven of eight criteria are deterministic code over tool output. Only the
literature axis calls the API, and it returns a score with a written
justification and citations so a reviewer can check the reasoning. Weights
live in `config/weights.yaml`; re-ranking under different weights is cheap and
should be routine. If the top ten changes completely when essentiality drops
from 20 to 15, that instability is a finding and belongs in the output.

Missing evidence scores as a configured neutral value and raises a flag —
never zero. RNAi coverage is partial, so absence of a knockdown phenotype is
absence of evidence, and scoring it zero would bias the screen toward the
already-studied proteins the project exists to look past.

## Cost controls

Across both species the starting set is ~37,700 proteins, so the most
important budget decision is **not calling the API on all of them**. Three
controls:

- **per-stage spend ceilings** that halt the stage rather than overrun it;
- a **disk cache** keyed on a hash of the exact prompt, so reruns re-spend
  nothing — this is where debugging budgets actually go;
- **dry-run by default**, reporting projected calls and tokens without
  calling. Any config change should get a dry run first.

## Validation

Build the benchmark before the first full run.

- **Positive controls** (`resources/controls/positive_controls.tsv`): P0,
  my32/akirin, peroxiredoxin, two heme peroxidases, the Bm86 ortholog. These
  must surface in the top tier. **The literature stage is blinded for this
  test** — if the model can read that P0 was used in a trial, a high score
  proves nothing.
- **Negative controls**: salmon proteins must be rejected by cross-reactivity
  and intracellular housekeeping proteins by exposure. Being rejected for the
  *wrong* reason counts as a failure.
- **Robustness**: weight perturbation and leave-one-source-out. Candidates
  that hold across all variants go to the bench; ones depending entirely on a
  single source are labelled `single_source_dependent` in their dossier
  rather than ranked quietly alongside the robust ones.

## Known limits

- **Epitope prediction is calibrated on mammals.** Salmon IgM/IgT biology
  differs enough that these scores are weak evidence. This is the softest
  axis in the rubric and the one most likely to be wrong.
- **The concealed-antigen mechanism may not transfer.** Ticks ingest a large
  blood meal; lice feed on mucus, skin and blood in smaller continuous
  amounts, and how much intact host immunoglobulin reaches the louse gut is
  not well established.
- **Durability may be a delivery problem, not an antigen problem.** If decay
  is driven by adjuvant and memory rather than antigen choice, a
  better-chosen antigen decays the same way and the screen optimised the
  wrong variable.
- **Annotation quality limits everything downstream.** A large fraction of
  gene models are hypothetical proteins. Uncertainty is carried forward as
  flags rather than silently scored as zero.

## Open decisions

See `docs/decisions.md`. Two of them change the build.
