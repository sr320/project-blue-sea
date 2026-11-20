# project-blue-sea

A reverse-vaccinology research pipeline for prioritizing protein antigen
hypotheses against *Lepeophtheirus salmonis* and *Caligus rogercresseyi*.
It produces auditable filtering decisions, gene-level rankings, sensitivity
reports and candidate dossiers. Ranking scores are heuristic priorities;
challenge trials are required to establish efficacy.

The default configuration runs a **bounded pilot**: mapped experimental
components and negative controls, 200 reproducibly sampled background proteins
per species, and three real salmon-proteome sequence spikes. A full predictor
screen requires a passing pilot report under the same biology, tools, controls
and scoring settings.

## Current readiness

The Python logic and workflow handoffs can be tested locally. The checked-in
configuration is **not biologically validated**. Licensed local predictors,
model weights, exact experimental construct sequences and additional
independent benchmark studies still need provisioning/curation. These are
explicit readiness failures. A diagnostic ranking can still be generated and
inspected when a biological benchmark fails.

## Development

```bash
pip install -e '.[dev]'
pytest -q
ruff check .
# Python >=3.11, for workflow integration tests and Snakemake:
pip install -e '.[dev,workflow]'
pytest -q tests/test_workflow_integration.py
```

The workflow integration test executes the actual DAG with stand-in tools
that emit native predictor schemas. It tests concurrency, parsers, checkpoints,
immutable handoffs, host spikes, ranking, sensitivity and a clean second run.
It does not evaluate predictor accuracy.

## Running a pilot

1. Install local predictors and weights; set commands, versions and required
   files in `config/config.yaml`. See [docs/tools.md](docs/tools.md).
2. Build and check the container; provision NCBI references. See
   [docs/reference-data.md](docs/reference-data.md) and [docs/hyak.md](docs/hyak.md).
3. Run the pilot, inspect the diagnostic outputs, then evaluate its gate.

```bash
bash scripts/build_container.sh
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm provision
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm all
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm gate
```

`all` preserves `benchmark.json` even when biological validation fails; `gate`
returns nonzero on a failed report. Before submitting jobs, inspect the DAG:

```bash
snakemake -s workflow/Snakefile --profile workflow/profiles/local -n all
```

Full-mode configuration uses a complete copy of the pilot configuration with
`screening.mode: full` and separate scratch/results/log directories. Supply it
with `--configfile`; the workflow checks the pilot inside the configured container before preparing full predictor inputs.
See [docs/validation.md](docs/validation.md) for readiness and curation rules.

## Evidence and biological assumptions

- **Exposure:** SignalP/DeepTMHMM provide predicted secretion and topology,
  not tissue localization or demonstrated host-antibody contact. Predicted
  routes carry an uncertainty flag. Intracellular proteins are excluded by
  the current exposure policy; loss of experimental components is reported
  rather than repaired by tuning controls into the shortlist.
- **Host similarity:** every returned significant HSP with at least 30% query
  coverage is considered before the 40% identity heuristic is applied.
  Short high-identity patches remain flagged for review. This screen neither
  establishes cross-reactivity nor certifies safety.
- **Epitopes:** native BepiPred residue scores become contiguous intervals.
  Only linear intervals wholly within predicted outside regions contribute
  when residue topology is available. This is advisory by default because
  mammalian antibody-complex training does not calibrate salmon protection.
- **Expression:** stages are species-specific. *Caligus* has four chalimus
  stages and no preadult stages; *L. salmonis* has two chalimus and two
  preadult stages. Unassayed stages are distinguished from zero expression.
  [Primary *Caligus* genome study](https://www.nature.com/articles/s41597-021-00842-w).
- **Conservation:** cross-species status requires bilateral reciprocal matches
  and measured epitope-region identity and coverage. Unaligned epitope residues
  count against identity. Shared sequence does not establish essentiality or
  cross-protection.
- **Essentiality:** RNAi evidence is scored separately; untested is distinct
  from a reported null phenotype. RNAi knockdown and antibody perturbation
  need not produce the same effect.
- **Literature:** optional model judgments use curated excerpts backed by
  checksummed local snapshots. Missing sources yield no model call and a
  neutral missing-evidence score. Returned citations must refer to supplied
  sources. The deterministic benchmark excludes literature evidence entirely.

The [2025 P0 study](https://doi.org/10.3389/fimmu.2025.1570948) evaluated chimeric
constructs with adjuvant and prime-boost delivery. At the first main count,
reductions were 88–91% against non-injected controls and 77–84% against adjuvant
controls. Effects differed at later challenges. These results motivate
longitudinal testing but do not establish that constitutive expression or
essentiality predicts protection duration. The study's limited antibody
sampling also does not settle whether titer predicts protection.

## Outputs and layout

The default paths are `scratch/pilot`, `logs/pilot`, and `results/pilot`.

| Output | Purpose |
| --- | --- |
| `ranked_all.tsv`, `<species>/ranked.tsv` | One best-scoring isoform per species/gene |
| `ranking_sensitivity.json` | Weight perturbations and removal of each score source |
| `benchmark.json` | Deterministic recovery, control completeness and readiness |
| `dossiers/` | Scores, evidence, caveats and supported domain suggestions |
| `provenance.json` | Code/configuration hashes and run context |
| `scratch/<species>/raw/<tool>/<shard>/` | Native predictions, logs and checksum manifests |
| `logs/<species>/cascade_audit.jsonl` | Every protein's filter decisions |

Records progress through separate `prefiltered`, `assembled`, `cascaded`,
`orthologized`, `assessed`, and `ranked` JSON files. Enrichment stages do not
overwrite their inputs. All isoforms remain in the record files for inspection.
Sensitivity flags are persisted into ranked records and dossiers.

Weights and sensitivity thresholds are explicit judgments in `config/`.
Sensitivity analyzes score dependence among survivors; it does not test
false negatives introduced by hard filters. Disable the prefilter's predictor
restriction and profile the cascade only when undertaking that separate,
more expensive analysis.

The symbiont arm remains exploratory. Its MAG gate checks quality alongside
prevalence and genomic signatures; those signatures do not prove vertical
transmission or obligate dependence. MAG assembly is not implemented.
