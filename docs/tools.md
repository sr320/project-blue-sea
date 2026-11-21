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

## Availability, checked 2026-11-21

Distribution terms move; the notes below record what the primary pages said on
that date, not a standing guarantee. Re-check before provisioning, and treat
anything not quoted here as unverified.

None of the three vendors states whether their software may be baked into a
container image on shared institutional infrastructure. Staging into an image
for one site is ordinary practice, but if the image store is readable beyond
the group, confirm with the vendor rather than assuming.

## SignalP 6

The configured example invokes `signalp6` with eukaryotic predictions. Its
native output must be `output/prediction_results.txt`, with ID, Prediction,
OTHER probability, SP probability and optional cleavage-site field. The parser
uses the SP column explicitly and retains the cleavage coordinate.

Install a licensed SignalP 6 distribution, its dependencies and weights; add
its script/model paths to `required_files` and set the actual version.

The current DTU build is 6.0i. There is a download page for academic users;
others are directed to the DTU Health Tech Software Package Manager. The
licence text itself is not published on the service page, so read it at
download time and record which version and mode the site is licensed for.

`fast` and `slow_sequential` are separate downloads, and an installation only
carries the mode that was fetched; adding another means installing its weights
too. `fast` is a distilled model and is the appropriate choice for a
proteome-scale screen. `slow_sequential` takes six times as long and buys only
more accurate region borders. The bare `slow` mode has been discontinued.

Version constraints are tight and easy to miss: Python up to 3.10, PyTorch
strictly below 2.0, and `numpy<2`. SignalP runs on CPU by default; GPU use
needs `signalp6_convert_models gpu`, and GPU-converted weights then fail on
CPU. Weights live in the installed package directory or wherever
`signalp6 --model_dir` points.

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

Academic local use needs no licence negotiation. The DTU service page states
the software "is made freely available here for both academic and commercial
purposes" and sends users to BioLib for download instructions; only commercial
users self-hosting it need a licence, from `licensing@biolib.com`. What the
BioLib download actually delivers, and whether its client needs network access
at run time, is unconfirmed here and has to be settled against the obtained
package. An offline compute node is the requirement to test it against.

DeepTMHMM2 is a different tool, not a newer drop-in. It writes a file of the
same name, but its topology alphabet is wider than this parser accepts: `R`
(reentrant), `F` (interfacial) and `>` (transit peptide) survive even under
`--simplify-io`, and its default mode emits membrane-type-specific letters
instead of `I`/`O` -- `E`/`e` for eukaryotic plasma membrane, `H`/`h` for ER.
That would not merely trip the label check in `io.read_deeptmhmm`; it would
silently zero `extracellular_aa`, which the exposure filter thresholds on.
Adopting it means deciding what reentrant and interfacial residues mean for
antibody access, which is a scientific judgement and belongs in
[decisions.md](decisions.md), not a character substitution in the parser. Note
also that it ships no LICENSE file and fetches unpinned, unchecksummed weights
from a branch head.

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

The BP3 ensemble weights are committed to that repository under
`bp3/BP3Models/`, so the checkout supplies them; only the ESM-2 weights are
fetched separately. The repository carries no LICENSE file. Its README grants
free use "by academic groups for non-commercial purposes" and directs for-profit
use to a separate licence from Morten Nielsen at DTU. Record that the
installation is academic-only alongside the version.

Its pins are old enough to need their own environment rather than the shared
one: `numpy==1.22`, `pandas==1.4.3`, `fair-esm==1.0.3`, `bp3==0.0.12.7`, with
the authors reporting Python 3.8.8 and leaving PyTorch unpinned.

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

### ESM-2 weights

Both files are needed, in one directory, and the driver is given that path so
nothing is fetched at run time:

| file | size |
| --- | --- |
| `esm2_t33_650M_UR50D.pt` | 2,604,537,549 bytes |
| `esm2_t33_650M_UR50D-contact-regression.pt` | 3,687 bytes |

They are served from `dl.fbaipublicfiles.com/fair-esm/` under `models/` and
`regression/` respectively. The ESM code is MIT licensed, so mirroring them is
permitted, and mirroring is now the prudent course: the upstream repository was
archived read-only on 2026-11-20. Both URLs still resolved on 2026-11-21, but an
archived repository is no promise that its asset hosting persists. Fetch them
to site storage, record their SHA-256, and provision from the mirror. Those
digests enter every prediction manifest through `required_files`, so a later
substitution invalidates a full-mode approval signature rather than passing
unnoticed.

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
