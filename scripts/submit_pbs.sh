#!/usr/bin/env bash
#
# Submits one PROfit subcommand to PBS as a single-node job.
#
# Usage: scripts/submit_pbs.sh <jobs/NAME.ini> <subcommand> [walltime] [queue] [-- extra PROfit args...]
#
#   scripts/submit_pbs.sh jobs/varset1_mc.ini process 06:00:00 capacity
#   scripts/submit_pbs.sh jobs/varset1_mc.ini plot 00:20:00
#
# Walltime defaults to 01:00:00 and queue to whatever input/site.conf says.
# Mind each queue's limits: `debug` is 1-2 nodes and caps at 01:00:00 (and
# allows only one job per user at a time); `capacity` is 1-16 nodes and allows
# up to 7 days. Anything longer than an hour needs capacity or the scheduler
# rejects the job outright.
#
# The generated job script is written into the run directory and submitted from
# there, so the submission directory is on the project filesystem rather than
# home (per ALCF guidance) and there is a durable record of what was run.
#
# This is deliberately one node and one task: multi-task fan-out is a separate
# concern and does not belong here yet.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_CONF="$REPO_ROOT/input/site.conf"

if [ "$#" -lt 2 ]; then
    echo "usage: scripts/submit_pbs.sh <jobs/NAME.ini> <subcommand> [walltime] [queue] [-- extra args...]" >&2
    exit 2
fi

INI="$1"; shift
SUBCOMMAND="$1"; shift

WALLTIME="01:00:00"
QUEUE_OVERRIDE=""
if [ "$#" -gt 0 ] && [ "$1" != "--" ]; then
    WALLTIME="$1"; shift
fi
if [ "$#" -gt 0 ] && [ "$1" != "--" ]; then
    QUEUE_OVERRIDE="$1"; shift
fi
[ "${1:-}" = "--" ] && shift
EXTRA=("$@")

if [ ! -f "$SITE_CONF" ]; then
    echo "error: $SITE_CONF not found." >&2
    echo "       cp input/site.conf.example input/site.conf and edit it for this machine." >&2
    exit 1
fi

get_setting() {
    local key="$1" line value
    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%%#*}"
        line="$(echo "$line" | xargs)"
        [ -z "$line" ] && continue
        [ "${line%%=*}" = "$key" ] || continue
        value="${line#*=}"
    done < "$SITE_CONF"
    echo "${value:-}"
}

ACCOUNT="$(get_setting pbs_account)"
QUEUE="$(get_setting pbs_queue)"
FILESYSTEMS="$(get_setting pbs_filesystems)"
OUTPUT_ROOT="$(get_setting output_root)"
NCPUS="$(get_setting pbs_ncpus)"
NCPUS="${NCPUS:-104}"
QUEUE="${QUEUE_OVERRIDE:-$QUEUE}"

for name in ACCOUNT QUEUE FILESYSTEMS OUTPUT_ROOT; do
    if [ -z "${!name}" ]; then
        echo "error: ${name,,} is not set in $SITE_CONF" >&2
        exit 1
    fi
done

JOB_NAME="$(basename "$INI" .ini)_${SUBCOMMAND}"
RUN_DIR="$OUTPUT_ROOT/$(basename "$INI" .ini)"
mkdir -p "$RUN_DIR"
JOB_SCRIPT="$RUN_DIR/${JOB_NAME}.pbs"

cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash -l
#PBS -N $JOB_NAME
#PBS -A $ACCOUNT
#PBS -l select=1
#PBS -l place=scatter
#PBS -l walltime=$WALLTIME
#PBS -l filesystems=$FILESYSTEMS
#PBS -q $QUEUE
#PBS -k doe

export OMP_NUM_THREADS=$NCPUS
export TMPDIR=/tmp

cd "$RUN_DIR"
exec "$REPO_ROOT/scripts/profit_run.sh" "$INI" "$SUBCOMMAND" ${EXTRA[@]+"${EXTRA[@]}"}
EOF

echo "job script: $JOB_SCRIPT"
echo "run dir:    $RUN_DIR"
echo "submitting: $JOB_NAME  (walltime $WALLTIME, queue $QUEUE, 1 node)"
cd "$RUN_DIR"
qsub "$JOB_SCRIPT"
