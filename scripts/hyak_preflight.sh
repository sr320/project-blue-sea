#!/usr/bin/env bash
# Pre-flight checks before the first real run on Hyak (klone).
#
#   bash scripts/hyak_preflight.sh
#   bash scripts/hyak_preflight.sh --account srlab --partition ckpt
#
# Every check is read-only. Nothing is submitted except two short probe jobs
# (2 minutes, 1 core, 1 GB) that exist to answer the one question the whole
# architecture rests on: do compute nodes have outbound network access?
#
# Exit status is the number of failed checks, so this can gate a first run.

set -uo pipefail

ACCOUNT="srlab"
PARTITION="compute"
BIND_PATH="/gscratch"
CONTAINER="containers/bluesea.sif"
PROBE_URL="https://api.anthropic.com/v1/messages"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --account)   ACCOUNT="$2";   shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --bind)      BIND_PATH="$2"; shift 2 ;;
        --container) CONTAINER="$2"; shift 2 ;;
        -h|--help)   sed -n '2,12p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

FAILURES=0
WARNINGS=0

pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
warn() { printf '  \033[33mWARN\033[0m  %s\n' "$1"; WARNINGS=$((WARNINGS + 1)); }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAILURES=$((FAILURES + 1)); }
note() { printf '        %s\n' "$1"; }
head2() { printf '\n\033[1m%s\033[0m\n' "$1"; }

head2 "Where am I"
note "host:      $(hostname)"
note "user:      ${USER:-unknown}"
note "account:   $ACCOUNT"
note "partition: $PARTITION"

# ---------------------------------------------------------------------------
head2 "1. SLURM account and partition"
# ---------------------------------------------------------------------------
# A wrong account is the most common first-run failure and the error sbatch
# gives for it ("Invalid account or account/partition combination") does not
# say which of the two is wrong, so they are checked separately.

if ! command -v sinfo >/dev/null 2>&1; then
    fail "sinfo not found -- this does not look like a SLURM cluster"
    note "Run this on a klone login node, not on your laptop."
else
    ASSOCS=$(sacctmgr --noheader --parsable2 show assoc \
                 user="${USER}" format=account 2>/dev/null | sort -u)
    if [[ -z "$ASSOCS" ]]; then
        warn "could not read your SLURM associations"
        note "sacctmgr may be restricted here; check the account by hand."
    elif grep -qx "$ACCOUNT" <<<"$ASSOCS"; then
        pass "account '$ACCOUNT' is one of yours"
    else
        fail "account '$ACCOUNT' is not in your associations"
        note "Yours: $(tr '\n' ' ' <<<"$ASSOCS")"
        note "Set it in workflow/profiles/slurm/config.yaml:default-resources"
    fi

    if sinfo -h -p "$PARTITION" -o '%P' 2>/dev/null | grep -q .; then
        STATES=$(sinfo -h -p "$PARTITION" -o '%a %D %t' | tr '\n' ';')
        pass "partition '$PARTITION' exists ($STATES)"
    else
        fail "partition '$PARTITION' not found"
        note "Available: $(sinfo -h -o '%P' | sort -u | tr '\n' ' ')"
    fi

    # ckpt is preemptible. The profile retries twice, which is usually enough,
    # but a long DIAMOND shard that gets preempted repeatedly will stall.
    if [[ "$PARTITION" == ckpt* ]]; then
        warn "'$PARTITION' is preemptible; expect requeues on long jobs"
        note "retries: 2 in the profile covers the usual case, not a bad day."
    fi
fi

# ---------------------------------------------------------------------------
head2 "2. Apptainer and the bind path"
# ---------------------------------------------------------------------------

if command -v apptainer >/dev/null 2>&1; then
    pass "apptainer present ($(apptainer --version 2>&1))"
elif command -v singularity >/dev/null 2>&1; then
    warn "only 'singularity' found, not 'apptainer'"
    note "Snakemake's --use-apptainer wants the apptainer binary; load a module."
else
    fail "neither apptainer nor singularity found"
    note "module load apptainer   (or whatever klone calls it now)"
fi

if [[ -d "$BIND_PATH" ]]; then
    pass "bind path '$BIND_PATH' exists"
    # The bind only helps if the workflow's own paths live under it. A
    # container that cannot see scratch/ fails at the first tool invocation.
    WORKDIR=$(pwd -P)
    if [[ "$WORKDIR" == "$BIND_PATH"/* ]]; then
        pass "working directory is under the bind path"
    else
        warn "working directory is NOT under '$BIND_PATH'"
        note "here: $WORKDIR"
        note "Either run from under $BIND_PATH, or add this path to"
        note "apptainer-args in workflow/profiles/slurm/config.yaml."
    fi
else
    fail "bind path '$BIND_PATH' does not exist"
    note "Fix apptainer-args in workflow/profiles/slurm/config.yaml."
fi

if [[ -f "$CONTAINER" ]]; then
    pass "container image present ($(du -h "$CONTAINER" | cut -f1))"
else
    warn "container image '$CONTAINER' not built yet"
    note "apptainer build $CONTAINER workflow/envs/bluesea.def"
fi

# ---------------------------------------------------------------------------
head2 "3. Outbound network -- the assumption the architecture rests on"
# ---------------------------------------------------------------------------
# If compute nodes turn out to have outbound access, nothing breaks: the
# three-stage split is still good practice. If they do not -- the expected
# case -- then this is what makes the split mandatory, and it is worth having
# the evidence on record rather than the assumption.

probe_url() {
    # A 401 means reachable and unauthenticated, which is the answer we want:
    # it proves the route without needing a key. 000 means no route at all.
    curl -s -o /dev/null -w '%{http_code}' --max-time 15 "$PROBE_URL" 2>/dev/null
}

if ! command -v curl >/dev/null 2>&1; then
    fail "curl not found; cannot test outbound access"
else
    LOGIN_CODE=$(probe_url)
    if [[ "$LOGIN_CODE" =~ ^(200|401|403|405)$ ]]; then
        pass "login node reaches the API (HTTP $LOGIN_CODE)"
    else
        fail "login node cannot reach the API (code '$LOGIN_CODE')"
        note "The API stage has nowhere to run. Check proxy settings:"
        note "  https_proxy=${https_proxy:-unset} HTTPS_PROXY=${HTTPS_PROXY:-unset}"
    fi

    if command -v srun >/dev/null 2>&1 && [[ $FAILURES -eq 0 || -n "${FORCE_PROBE:-}" ]]; then
        note "submitting a 2-minute probe job to '$PARTITION'..."
        COMPUTE_CODE=$(srun --account="$ACCOUNT" --partition="$PARTITION" \
            --time=2 --cpus-per-task=1 --mem=1G --quiet \
            curl -s -o /dev/null -w '%{http_code}' --max-time 15 "$PROBE_URL" \
            2>/dev/null)
        COMPUTE_CODE=${COMPUTE_CODE:-000}
        if [[ "$COMPUTE_CODE" =~ ^(200|401|403|405)$ ]]; then
            warn "compute nodes DO reach the API (HTTP $COMPUTE_CODE)"
            note "The three-stage split is not strictly required here. Keep it"
            note "anyway -- it makes the expensive stage separately"
            note "restartable and cacheable -- but you have more freedom than"
            note "the design assumed."
        else
            pass "compute nodes have no outbound access (code '$COMPUTE_CODE')"
            note "This is the expected result and confirms the design: the"
            note "literature rule must stay a localrule."
        fi
    else
        warn "skipped the compute-node probe"
        note "Earlier failures make the result ambiguous. Re-run with"
        note "FORCE_PROBE=1 to submit it anyway."
    fi
fi

# ---------------------------------------------------------------------------
head2 "4. Node-local scratch"
# ---------------------------------------------------------------------------
# The profile sets tmpdir=/scr so per-shard temporary files stay off the
# shared filesystem. If compute nodes have no /scr, that setting sends every
# job's temp writes to a path that does not exist.

if command -v srun >/dev/null 2>&1 && [[ -n "${FORCE_PROBE:-}" || $FAILURES -eq 0 ]]; then
    SCR=$(srun --account="$ACCOUNT" --partition="$PARTITION" --time=2 \
        --cpus-per-task=1 --mem=1G --quiet \
        bash -c 'if [[ -d /scr && -w /scr ]]; then echo ok; else echo missing; fi' \
        2>/dev/null)
    case "${SCR:-}" in
        ok)      pass "/scr exists and is writable on compute nodes" ;;
        missing) fail "/scr is missing or not writable on compute nodes"
                 note "Change tmpdir in workflow/profiles/slurm/config.yaml." ;;
        *)       warn "could not determine whether /scr exists" ;;
    esac
else
    warn "skipped the /scr probe (see above)"
fi

# ---------------------------------------------------------------------------
head2 "5. Workflow tooling"
# ---------------------------------------------------------------------------

if command -v snakemake >/dev/null 2>&1; then
    pass "snakemake present ($(snakemake --version 2>&1))"
    # The SLURM executor is a separate package in Snakemake 8+. Without it,
    # `executor: slurm` in the profile is an invalid choice and the run dies
    # before submitting anything.
    if snakemake --executor slurm --version >/dev/null 2>&1; then
        pass "slurm executor plugin installed"
    else
        fail "snakemake-executor-plugin-slurm is not installed"
        note "pip install snakemake-executor-plugin-slurm"
    fi
else
    fail "snakemake not found"
fi

if command -v bluesea >/dev/null 2>&1; then
    pass "bluesea CLI on PATH"
else
    warn "bluesea CLI not on PATH (fine if it only lives in the container)"
fi

if [[ -n "${ANTHROPIC_API_KEY:-}" ]]; then
    pass "ANTHROPIC_API_KEY is set"
else
    warn "ANTHROPIC_API_KEY is not set"
    note "Only needed for a live literature run; dry runs do not use it."
fi

# ---------------------------------------------------------------------------
head2 "6. Does the DAG resolve here"
# ---------------------------------------------------------------------------
# Catches missing reference downloads before an allocation is spent on them.

if command -v snakemake >/dev/null 2>&1; then
    if DAG=$(snakemake --snakefile workflow/Snakefile \
                 --profile workflow/profiles/local --dry-run all 2>&1); then
        pass "dry run resolves ($(grep -c '^rule ' <<<"$DAG") rule instances)"
    else
        fail "dry run failed"
        grep -E 'Error|Exception|Missing input' <<<"$DAG" | head -5 \
            | while read -r line; do note "$line"; done
    fi
fi

# ---------------------------------------------------------------------------
head2 "Summary"
# ---------------------------------------------------------------------------
printf '  %d failed, %d warnings\n\n' "$FAILURES" "$WARNINGS"
if [[ $FAILURES -gt 0 ]]; then
    echo "Fix the failures before the first real run."
else
    echo "Clear to run:"
    echo "  snakemake --profile workflow/profiles/slurm --use-apptainer gate"
fi
exit "$FAILURES"
