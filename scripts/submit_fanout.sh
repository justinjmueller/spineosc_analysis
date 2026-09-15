#!/usr/bin/env bash
#
# Runs several job configs concurrently as one multi-node PBS job: one config
# per node, one rank per node.
#
# Usage:
#   scripts/submit_fanout.sh <a.ini> <b.ini> ... -- <subcommand> [args...]
#   scripts/submit_fanout.sh --walltime 04:00:00 --queue capacity <inis> -- surface -g 30
#
# Example:
#   scripts/submit_fanout.sh jobs/contours/*.ini -- surface -g 30 --xlo 1e-2 --ylo 1e-1
#
# WHY PACK RATHER THAN SUBMIT SEPARATELY: capacity allows 16 nodes per job but
# only 2 running jobs per user. Six single-node jobs would run two at a time and
# serialise three deep; one six-node job runs all six at once. The node ceiling
# is not the binding constraint, the per-user job limit is.
#
# Ranks are entirely independent - PROfit is not MPI-parallel, mpiexec is only
# placing one process per node. Each rank writes its own cache, so a node
# failure costs only the incomplete ranks; re-run those configs.
#
# Each config must have a DISTINCT tag, since profit_run.sh keys the run
# directory on the tag and two ranks sharing one would collide.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_CONF="$REPO_ROOT/input/site.conf"

WALLTIME="04:00:00"
QUEUE_OVERRIDE=""
INIS=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        --walltime) WALLTIME="$2"; shift 2 ;;
        --queue)    QUEUE_OVERRIDE="$2"; shift 2 ;;
        --)         shift; break ;;
        *)          INIS+=("$1"); shift ;;
    esac
done
SUBCOMMAND="${1:-}"; [ "$#" -gt 0 ] && shift || true
EXTRA=("$@")

if [ "${#INIS[@]}" -eq 0 ] || [ -z "$SUBCOMMAND" ]; then
    echo "usage: scripts/submit_fanout.sh [--walltime HH:MM:SS] [--queue Q] <a.ini> ... -- <subcommand> [args...]" >&2
    exit 2
fi

if [ ! -f "$SITE_CONF" ]; then
    echo "error: $SITE_CONF not found (cp input/site.conf.example input/site.conf)" >&2
    exit 1
fi

get_setting() {
    local key="$1" line value
    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%%#*}"; line="$(echo "$line" | xargs)"
        [ -z "$line" ] && continue
        [ "${line%%=*}" = "$key" ] || continue
        value="${line#*=}"
    done < "$SITE_CONF"
    echo "${value:-}"
}

ACCOUNT="$(get_setting pbs_account)"
QUEUE="${QUEUE_OVERRIDE:-$(get_setting pbs_queue)}"
FILESYSTEMS="$(get_setting pbs_filesystems)"
OUTPUT_ROOT="$(get_setting output_root)"
NCPUS="$(get_setting pbs_ncpus)"; NCPUS="${NCPUS:-104}"

N="${#INIS[@]}"

# capacity tops out at 16 nodes; debug-scaling allows more but caps walltime at
# one hour, which a process run does not fit inside.
if [ "$N" -gt 16 ]; then
    echo "error: $N configs exceeds the 16-node ceiling on capacity. Split the batch." >&2
    exit 1
fi

# The site.conf default queue suits single-node interactive-scale work; a
# fan-out is neither. Catch the mismatch here rather than letting the scheduler
# reject the job after the batch has been staged.
wall_secs() { awk -F: '{print $1*3600 + $2*60 + $3}' <<< "$1"; }
case "$QUEUE" in
    debug)
        [ "$N" -gt 2 ] && { echo "error: queue 'debug' allows at most 2 nodes, this batch needs $N. Use --queue capacity." >&2; exit 1; }
        ;&
    debug-scaling)
        if [ "$(wall_secs "$WALLTIME")" -gt 3600 ]; then
            echo "error: queue '$QUEUE' caps walltime at 01:00:00, requested $WALLTIME. Use --queue capacity (up to 7 days)." >&2
            exit 1
        fi
        ;;
esac

declare -a TAGS=()
for ini in "${INIS[@]}"; do
    [ -f "$ini" ] || [ -f "$REPO_ROOT/$ini" ] || { echo "error: no such config: $ini" >&2; exit 1; }
    path="$ini"; [ -f "$path" ] || path="$REPO_ROOT/$ini"
    tag="$(grep -oP '^\s*tag\s*=\s*\K\S+' "$path" | tail -1)"
    [ -n "$tag" ] || { echo "error: $ini declares no tag" >&2; exit 1; }
    for t in ${TAGS[@]+"${TAGS[@]}"}; do
        [ "$t" = "$tag" ] && { echo "error: duplicate tag '$tag' - ranks would share a run directory" >&2; exit 1; }
    done
    TAGS+=("$tag")
done

BATCH="$OUTPUT_ROOT/_batches/${SUBCOMMAND}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BATCH"
TASKS="$BATCH/tasks.txt"
: > "$TASKS"
for ini in "${INIS[@]}"; do
    path="$ini"; [ -f "$path" ] || path="$REPO_ROOT/$ini"
    printf '%s %s %s %s\n' "$REPO_ROOT/scripts/profit_run.sh" "$path" "$SUBCOMMAND" "${EXTRA[*]}" >> "$TASKS"
done

JOB_SCRIPT="$BATCH/fanout.pbs"
cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash -l
#PBS -N fanout_${SUBCOMMAND}
#PBS -A $ACCOUNT
#PBS -l select=$N
#PBS -l place=scatter
#PBS -l walltime=$WALLTIME
#PBS -l filesystems=$FILESYSTEMS
#PBS -q $QUEUE
#PBS -k doe

cd "$BATCH"

# One rank per node. Cores 0 and 52 are left to the OS, so 102 are bound and
# OMP_NUM_THREADS matches that rather than the nominal $NCPUS.
mpiexec -n $N --ppn 1 --cpu-bind list:1-51,53-103 \\
  --env OMP_NUM_THREADS=102 \\
  --env TMPDIR=/tmp \\
  bash -c '
    RANK="\${PALS_RANKID:-0}"
    LINE=\$(( RANK + 1 ))
    CMD=\$(sed -n "\${LINE}p" "$TASKS")
    if [ -z "\$CMD" ]; then
      echo "rank \$RANK: no task on line \$LINE of $TASKS" >&2
      exit 1
    fi
    echo "rank \$RANK START: \$CMD"
    eval "\$CMD"
    rc=\$?
    echo "rank \$RANK DONE rc=\$rc: \$CMD"
    # Per-rank marker so a partial failure is easy to identify and resume.
    echo "\$rc" > "$BATCH/rank_\${RANK}.rc"
    # Always exit 0. The ranks are independent, but mpiexec tears down the whole
    # job the moment any rank exits non-zero - one bad config would otherwise
    # SIGTERM every other rank mid-fit and throw away their work. Failures are
    # recorded in the .rc markers instead; check those, not the job exit status.
    exit 0
  '
EOF

echo "batch dir:  $BATCH"
echo "tasks:      $N config(s), subcommand '$SUBCOMMAND' ${EXTRA[*]:-}"
for i in "${!INIS[@]}"; do printf '  rank %-2s %-46s tag=%s\n' "$i" "${INIS[$i]}" "${TAGS[$i]}"; done
echo "submitting: $N node(s), walltime $WALLTIME, queue $QUEUE"
cd "$BATCH"
qsub "$JOB_SCRIPT"
