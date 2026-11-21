# Benchmark and validation policy

The default pilot evaluates a bounded sample before a full predictor run.
It includes all mapped construct components and negative controls, a seeded
background sample per species, and actual sequences from the host proteome.
Host spikes carry reserved identifiers and cannot enter candidate rankings.
Their exposure step explicitly records a benchmark bypass so the experiment
isolates host-similarity filtering instead of being preempted by localization.

## Traceable experimental constructs

`resources/controls/experimental_constructs.json` records a primary DOI and
source location, construct description, louse components, source accessions
and regions, reference mappings, formulation, comparator-specific endpoints,
timing, sequence file/checksum and curation status.

The current entries describe TT-P0 and P0-my32 from one primary study.
They are diagnostic and pending the exact administered sequence. Both share
the study and the identical P0 fragment, from one template, so they are not
independent positive protein controls.

Fragment mapping is partly settled. A component carries a `derived_fragment`
or `derived_insert` block when its sequence was reconstructed here rather
than read from a paper: what it was derived from, how it matches the
configured annotation, and what independently confirms it. A component
reaches `verified` only on an exact annotation match whose boundaries hold
up against the cloning primers. Deriving a component never advances a
construct's own `curation_status`, which tracks the administered sequence.

That study prints no construct sequence: only PCR primers, and MS/MS identity
confirmations rather than sequences. So both records stay `sequence_pending`
on evidence, not on unfinished reading. Each carries its design references,
expression details, cloning primers with restriction sites, and a blockers
and next-steps list naming what would actually close it.

Study-level facts shared by both constructs live under `studies`, keyed by
DOI: the timeline in study weeks and the unresolved timing conflicts. Two are
recorded. The reported "2, 5, 11 months" labels match no challenge date under
either the prime or the boost anchor, and Figure 5's 32 dpc contradicts the
week 44-to-51 interval. Arithmetic that would explain each is recorded as a
derived hypothesis, marked as not stated by the paper, and is not used to
populate any field. Endpoints carry their comparator, p-value and the section
they came from; one abstract-only figure with no Results counterpart is
labelled as such rather than promoted to a result.

Before marking a construct verified:

- Obtain the exact administered sequence, including carriers, linkers and tags,
  from the primary paper, supplements or deposited sequence. Save one FASTA
  record and its SHA-256 file checksum.
- Align each louse fragment to the configured annotation; record source
  coordinates and verify its mapping. Whole-protein family similarity is
  insufficient to establish the exact vaccine fragment.
- Check comparator, endpoint, timing, design and limitations against the
  primary source. Preserve inconsistent or missing details explicitly.
- Add independently sourced studies. The default policy requires at least two;
  that is a declared readiness policy, not a power calculation.

All louse components of a construct must reach the top tier to count as
recovered. Studies count once, regardless of how many constructs they tested.
Component recovery diagnoses a whole-proteome ranking; it does not predict
how a fusion or carrier will perform.

The previous proposed protein controls are retained in
`legacy_candidate_controls.tsv` for curation. The active standalone-positive
file has no entries because their exact construct-level provenance is not yet
verified. Do not reactivate entries on annotation names or incomplete review
citations alone.

## Deterministic gate

`benchmark` reads the immutable `orthologized.json` records, including dropped
proteins, and explicitly removes all literature evidence. It uses the same
orthology and deterministic evidence available to production ranking. Missing
negative controls fail completeness. Every observed negative isoform must
fail at the expected step.

The report includes construct readiness/recovery, negative controls, host
spikes, weight perturbations and score-source removal. The default sensitivity
policy requires mean top-tier stability >=0.8 and at least half the baseline
top tier retained under every configured source removal. Thresholds live in
`config.yaml`; set them before inspecting the outcome. Robustness is calculated
with copied records, so the final scores remain the baseline scores.

Leave-one-source-out only tests a source that varies. A source whose table is
unconfigured scores every record identically, so removing it reorders nothing
and hands back the baseline tier intact, which the policy above would read as
flawless stability -- most convincingly where there is no evidence at all. The
report therefore counts the distinct scores each configured source produces
and lists any that are constant under
`sensitivity_validation.uninformative_sources`. Those fail the verdict instead
of passing it vacuously, so an absent evidence table surfaces as absence. An
unrecognised source name drops no criteria and is caught the same way.

Random-control percentiles describe the surviving background sample. A random
sample's middle rank is not evidence of biological predictive accuracy and
is not used as an efficacy gate.

```bash
bluesea benchmark --report-only   # preserve a diagnostic report even on failure
bluesea check-gate --report results/pilot/benchmark.json
```

Snakemake `all` writes diagnostics; `gate` requires the report to pass. A failed
gate leaves no `gate.passed` marker and preserves its report.

## Production ranking and sensitivity

Ranked TSVs and dossiers use one best-scoring representative per species/gene.
Every isoform remains in `ranked.json`. The final, optionally literature-enriched
ranking receives its own `ranking_sensitivity.json`; source-dependent flags
are saved in the ranked records before dossiers are generated.

Source removal changes score contributions among survivors. It does not rerun
filters. Expression or exposure false negatives require a separate ablation
run, with predictor coverage expanded accordingly.

## Full mode

Use separate output directories and a complete configuration file. Set
`screening.mode: full` and point `screening.pilot_report` to a passed pilot.
The pipeline checks a signature of biological configuration, tool settings,
controls and weights and input/code checksums before preparing full predictor inputs. Changing those invalidates approval.
A failed/pending benchmark must be addressed with evidence or an explicitly
revised scientific policy and a new pilot.

A passing software test suite does not supply the missing benchmark sequences,
independent experiments, localization evidence or challenge-trial results.

Snakemake writes its effective, merged configuration to
`scratch/.../resolved_config.yaml`. Every CLI stage consumes that snapshot,
including configfile and command-line overrides. Changing a configured evidence
file, weights, code or a present container image is a declared workflow input.
