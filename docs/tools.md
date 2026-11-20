# Local predictor integration

The workflow invokes `python -m bluesea.predictors`. Each shard gets a unique
working directory, its own native output directory, stdout/stderr logs and a
checksum manifest. Successful outputs are normalized only after their IDs,
coverage and residue sequences have been checked. Failed work directories
remain beside the raw output directories for diagnosis.

Configure each `predictors` entry with an argument-array `command`, a recorded
`version`, a `native_output` path relative to the working directory, and
`required_files` listing the installed scripts and model files. Commands support
`{input}`, `{output_dir}`, `{work_dir}` and `{threads}` placeholders. Use
absolute paths for scripts and interpreters; commands execute in an isolated
working directory. Native output paths can also be absolute.

```bash
python -m bluesea.predictors --config config/config.yaml --check
```

This installation check does not execute a model. A real control pilot is
still required. The checked-in configuration deliberately fails this check
until local tools, weights and versions have been supplied.

## SignalP 6

The configured example invokes `signalp6` with eukaryotic predictions. Its
native output must be `output/prediction_results.txt`, with ID, Prediction,
OTHER probability, SP probability and optional cleavage-site field. The parser
uses the SP column explicitly and retains the cleavage coordinate.

Install a licensed SignalP 6 distribution, its dependencies and weights; add
its script/model paths to `required_files` and set the actual version.

## DeepTMHMM 1

Supply a local driver command that consumes `{input}` and writes
`{output_dir}/predicted_topologies.3line`. The file must contain the original
FASTA identifier, exact protein sequence, and an equal-length topology string
using `S`, `M`, `I`, `O`, `B`. The parser consumes exactly three lines per
record, retaining pipe characters in identifiers. A driver using another
schema/version needs an explicit adapter rather than guessing its labels.

No default cloud command is supplied: a `biolib run` command is insufficient
for an offline compute-node deployment. Obtain and provision the permitted
local distribution/driver for the site and pin its version and weights.

`O` means the non-cytoplasmic side of a membrane. It can include organelle
lumens and does not establish exposure to host antibody. See the
[authors' topology explanation](https://github.com/fteufel/DeepTMHMM2).

## BepiPred 3

The bundled `bluesea.bepipred` driver calls the
[official Python API](https://github.com/UberClifford/BepiPred-3.0) with a local
ESM-2 model path. Provision a compatible pinned checkout, its BP3 ensemble
weights, Python dependencies, `esm2_t33_650M_UR50D.pt`, and the matching
`esm2_t33_650M_UR50D-contact-regression.pt`. It does not request model downloads.
Dependencies belong in the predictor's installation environment; the baseline
container does not install PyTorch/ESM or BP3.

The driver uses 1,000-residue tiles with 100-residue overlap for long proteins,
then averages overlapping smoothed scores. This is an explicit approximation;
boundary predictions need validation. Native tile outputs remain archived.

The adapter reads native `raw_output.csv`: `Accession`, `Residue`, and
`BepiPred-3.0 linear epitope score`. It turns contiguous above-threshold runs
into one-based inclusive intervals, retaining mean score and sequence. It
rejects missing proteins, truncated predictions, invalid scores and residue
mismatches. Proteins with no above-threshold run are valid zero-hit results.
The threshold and minimum run length are configurable; they are not calibrated
salmon efficacy thresholds.

## Container deployment

Stage installations, local weights and an explicit `install.sh` in
`workflow/predictors/`. The definition copies that directory to `/opt/predictors`
and executes the installer after installing the package. Configure commands
to use the resulting local interpreter/executable paths. Do not commit the
licensed software or weights.

```bash
bash scripts/build_container.sh
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm tools_check
```

The build script checks the configured tools inside the built image and
returns a failure if they are incomplete. CI exercises adapters and the real
DAG with native-format stand-ins; it does not certify the licensed models.
