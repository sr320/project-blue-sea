# First run on Hyak

Run the checks rather than reading them:

```bash
bash scripts/hyak_preflight.sh --account srlab --partition compute
```

It exits with the number of failed checks, so it can gate a first run. Every
check is read-only apart from two probe jobs (2 minutes, 1 core, 1 GB).

What follows is what each check is actually asking, and what to do when it
comes back the other way.

## Account and partition

`workflow/profiles/slurm/config.yaml` ships with `srlab` and `compute`. These
are placeholders written from outside the cluster and should be confirmed
before the first submission, not after.

The reason they are checked separately is that `sbatch` reports a bad account
and a bad account/partition *combination* with the same message, which does
not say which half is wrong. The script reads your associations with
`sacctmgr` and the partition list with `sinfo` so the answer is unambiguous.

To override without editing the file:

```bash
snakemake --profile workflow/profiles/slurm \
  --default-resources slurm_account=otheracct slurm_partition=ckpt
```

Note that Snakemake does **not** expand shell or environment syntax inside
`default-resources`. A `"${SLURM_ACCOUNT:-srlab}"` there reaches `sbatch`
verbatim and the submission fails. That is why the profile holds literals.

If you move to `ckpt` for throughput, it is preemptible: the profile's
`retries: 2` covers an ordinary requeue, but a long DIAMOND shard that keeps
getting preempted will stall. Watch the first sweep before trusting it
overnight.

## The apptainer bind path

The profile passes `--bind /gscratch:/gscratch`. A bind only helps if the
workflow's own paths sit under it — a container that cannot see `scratch/`
fails at the first tool invocation, and the error points at the tool rather
than at the bind, which wastes an afternoon.

So the check is two-part: does `/gscratch` exist, and is the working
directory under it. If you run from somewhere else, either move the checkout
or add that path to `apptainer-args`.

Build the image before the first run:

```bash
apptainer build containers/bluesea.sif workflow/envs/bluesea.def
```

Record its digest so provenance can pick it up:

```bash
export BLUESEA_CONTAINER_DIGEST=$(apptainer inspect --json containers/bluesea.sif \
  | python -c 'import json,sys; print(json.load(sys.stdin)["data"]["attributes"]["labels"].get("org.label-schema.build-id","unknown"))')
```

SignalP 6 and DeepTMHMM are licence-restricted and are not in
`environment.yml`. They have to be added to the image at build time from the
DTU download. Until they are, stage 1 cannot run even if everything else
checks out.

## Outbound access — the one that matters

The three-stage architecture exists because Hyak compute nodes cannot be
assumed to have outbound network access. That assumption has never been
tested here; the script tests it.

It probes `api.anthropic.com` from the login node, then submits a short job
that probes it again from a compute node. An HTTP 401 or 405 means the route
works (reachable, unauthenticated) — no API key is needed to get the answer.
A `000` means no route at all.

Three outcomes:

- **Login yes, compute no.** The expected result, and the one the design
  assumes. It confirms the `literature` rule must stay in `localrules`, where
  it runs on the submitting login node instead of being handed to SLURM.
- **Both yes.** The split is not strictly required. Keep it anyway — it makes
  the expensive stage separately restartable and cacheable, which is worth
  having on its own — but you have more freedom than the design assumed.
- **Login no.** The API stage has nowhere to run at all. Check whether klone
  wants an HTTPS proxy set (`https_proxy`); if so, the literature stage needs
  that in its environment and the finding belongs in this file.

The compute probe is skipped when an earlier check already failed, because a
bad account would make a `000` ambiguous — it would mean "job never ran",
not "no network". Force it with `FORCE_PROBE=1` once you know why.

## Node-local scratch

The profile sets `tmpdir: /scr` so per-shard temporary files stay off the
shared filesystem. If compute nodes have no writable `/scr`, every job's
temp writes go to a path that does not exist. The script checks this from
inside a probe job rather than from the login node, because the login node's
`/scr` tells you nothing about a compute node's.

## Before spending an allocation

```bash
# Resolves the DAG without submitting anything.
snakemake --profile workflow/profiles/local --dry-run all

# The gate: positive controls must reach the top tier on deterministic
# evidence alone. Non-zero exit means rebuild the rubric, not tune it.
snakemake --profile workflow/profiles/slurm --use-apptainer gate
```

Then, and only then:

```bash
snakemake --profile workflow/profiles/slurm --use-apptainer all
```

The literature stage stays a dry run until `BLUESEA_LIVE=1` is set. Leave it
that way through the first full sweep — the deterministic stages are what
you are debugging, and a dry run reports projected call count and token
volume without spending anything.
