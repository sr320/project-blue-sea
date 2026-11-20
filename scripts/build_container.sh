#!/usr/bin/env bash
# Build the Apptainer image every workflow rule runs inside.
#
#   bash scripts/build_container.sh
#   bash scripts/build_container.sh --def workflow/envs/bluesea.def --out containers/bluesea.sif
#
# This is a script rather than a Snakemake rule because the rule would have to
# run inside the image it is building. It is a prerequisite of the workflow,
# like scripts/hyak_preflight.sh, not a step in it.
#
# Stage licensed local tools, model weights and install.sh in
# workflow/predictors/ before building. See docs/tools.md. The base image
# alone is insufficient; a post-build predictor check must pass.

set -uo pipefail

DEF="workflow/envs/bluesea.def"
OUT="containers/bluesea.sif"
FAKEROOT=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --def)       DEF="$2"; shift 2 ;;
        --out)       OUT="$2"; shift 2 ;;
        --no-fakeroot) FAKEROOT=0; shift ;;
        -h|--help)   sed -n '2,15p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; }
note() { printf '        %s\n' "$1"; }
head2() { printf '\n\033[1m%s\033[0m\n' "$1"; }

head2 "Checks"

if command -v apptainer >/dev/null 2>&1; then
    pass "apptainer $(apptainer --version 2>/dev/null | awk '{print $NF}')"
    BUILDER=apptainer
elif command -v singularity >/dev/null 2>&1; then
    pass "singularity $(singularity --version 2>/dev/null | awk '{print $NF}')"
    BUILDER=singularity
else
    fail "no apptainer or singularity on PATH"
    note "on Hyak: module load apptainer"
    exit 1
fi

if [[ ! -f "$DEF" ]]; then
    fail "definition file not found: $DEF"
    exit 1
fi
pass "definition file: $DEF"

if [[ -f "$OUT" ]]; then
    note "$OUT already exists; it will be overwritten"
fi

mkdir -p "$(dirname "$OUT")"

head2 "Building"
note "this takes several minutes and writes a few GB to \$APPTAINER_CACHEDIR"

BUILD_ARGS=()
# Unprivileged build. Root works too, but nobody has it on a cluster login
# node, and --fakeroot is the supported path there.
if [[ "$FAKEROOT" -eq 1 && "$(id -u)" -ne 0 ]]; then
    BUILD_ARGS+=(--fakeroot)
    note "building with --fakeroot (pass --no-fakeroot to disable)"
fi

if "$BUILDER" build "${BUILD_ARGS[@]}" "$OUT" "$DEF"; then
    head2 "Done"
    pass "$OUT ($(du -h "$OUT" | cut -f1))"
    if ! "$BUILDER" exec "$OUT" python -m bluesea.predictors --config config/config.yaml --check; then
        fail "image built, but predictor installation/configuration is incomplete"
        note "configure the local tools and weights as described in docs/tools.md"
        exit 1
    fi
    note "next: snakemake --profile workflow/profiles/slurm --use-apptainer provision"
    exit 0
fi

head2 "Failed"
fail "build did not complete"
note "if --fakeroot was refused, your site may require a remote builder:"
note "  $BUILDER build --remote $OUT $DEF"
exit 1
