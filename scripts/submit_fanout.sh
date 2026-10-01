#!/usr/bin/env bash
#
# Runs several job configs concurrently as one multi-node PBS job: one config
# per node, one rank per node.
#
# Usage:
#   scripts/submit_fanout.sh [OPTIONS] <a.ini> <b.ini> ... -- <subcommand> [args...]
#   scripts/submit_fanout.sh [OPTIONS] --variants <file> <base.ini> -- <subcommand> [args...]
#
# Options:
#   --walltime HH:MM:SS   default 04:00:00
#   --queue NAME          default from input/site.conf
#   --shared-cache        permit configs that share a tag (see below)
#   --variants FILE       expand ONE base config into many runs (implies --shared-cache).
#                         Each line is `label` plus exclude / inject-systs /
#                         inject-osc / syst-only
# Examples:
#   scripts/submit_fanout.sh jobs/contours/*.ini -- surface -g 30 --xlo 1e-2 --ylo 1e-1
#   scripts/submit_fanout.sh --shared-cache --variants jobs/contours/scan_detsyst.variants \
#       jobs/contours/contours_full.ini -- surface -g 30 --xlo 1e-2 --ylo 1e-1
#
# WHY PACK RATHER THAN SUBMIT SEPARATELY: capacity allows 16 nodes per job but
# only 2 running jobs per user. Six single-node jobs would run two at a time and
# serialise three deep; one six-node job runs all six at once. The node ceiling
# is not the binding constraint, the per-user job limit is.
#
# Ranks are entirely independent - PROfit is not MPI-parallel, mpiexec is only
# placing one process per node. Failures are recorded in per-rank .rc markers
# and every rank exits 0, so one bad config cannot tear down the others.
#
# SHARED CACHES. Normally each config needs a DISTINCT tag, because
# profit_run.sh keys the run directory on the tag. --shared-cache lifts that for
# the case where several runs are different *views* of one cache - excluding a
# systematic, fitting a different variable - which cost only the fit.
#
# It is gated on the cache already existing, and that gate is the whole point:
# `surface` and `plot` auto-process when no cache is present, so N ranks sharing
# a tag with no cache would all build <tag>_prop.bin into the same directory at
# once and corrupt each other. The adjacent case is already safe without help -
# if the cache exists but the XML changed, PROfit refuses to load on hash
# mismatch rather than silently rebuilding, so you get N clean refusals.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_CONF="$REPO_ROOT/input/site.conf"

WALLTIME="04:00:00"
QUEUE_OVERRIDE=""
SHARED_CACHE=0
VARIANTS=""
INIS=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        --walltime)     WALLTIME="$2"; shift 2 ;;
        --queue)        QUEUE_OVERRIDE="$2"; shift 2 ;;
        --shared-cache) SHARED_CACHE=1; shift ;;
        --variants)     VARIANTS="$2"; SHARED_CACHE=1; shift 2 ;;
        --)             shift; break ;;
        *)              INIS+=("$1"); shift ;;
    esac
done
SUBCOMMAND="${1:-}"; [ "$#" -gt 0 ] && shift || true
EXTRA=("$@")

if [ "${#INIS[@]}" -eq 0 ] || [ -z "$SUBCOMMAND" ]; then
    echo "usage: scripts/submit_fanout.sh [--walltime T] [--queue Q] [--shared-cache]" >&2
    echo "                                [--variants FILE] <ini>... -- <subcommand> [args...]" >&2
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

resolve() { local p="$1"; [ -f "$p" ] || p="$REPO_ROOT/$1"; [ -f "$p" ] || { echo "error: no such file: $1" >&2; exit 1; }; echo "$p"; }
tag_of()  { grep -oP '^\s*tag\s*=\s*\K\S+' "$1" | tail -1; }

BATCH="$OUTPUT_ROOT/_batches/${SUBCOMMAND}_$(date +%Y%m%d_%H%M%S)"

# ---- variant expansion: one base config becomes N per-variant configs --------
# Globals like -o and --exclude-systs belong in the config, not after the
# subcommand on the command line, so each variant gets a generated .ini in the
# batch directory. Those files are also the provenance record for the batch.
#
# A variant line is a label followed by one or more keyword sections, in any
# order. Every section is explicit - there is no default action:
#
#   exclude       <name> [<name>...]        -> exclude-systs = ...
#   inject-systs  <name> <sigma> [...]      -> inject-systs  = ...
#   inject-osc    <param> <value> [...]     -> inject        = ...
#   syst-only                               -> syst-only = true
#
#   nowmxtxw      exclude      WireModxThetaXW_sbnd_Run1
#   inj_ffqe_p1   inject-systs VecFFCCQEshape 1
#   qe_p1_noqe    inject-systs VecFFCCQEshape 1   exclude Cross-Section-QE
#   osc_1_0p5     inject-osc   dmsq 1 sinsq2thmm 0.5
#
# There is deliberately no bare `inject`: PROfit's `inject` key means physics
# parameters while `inject-systs` means systematics, and a variant keyword
# named `inject` would read naturally as either. A bare `inject` is rejected.
#
# A line whose second token is not a keyword is an error rather than being
# guessed at, so an old-style `label name...` line fails loudly at submit
# time instead of silently doing something else.
#
# OVERRIDE, NOT MERGE. A key the variant sets replaces the base config's value
# for that key; a key it does not set is inherited from the base unchanged.
#
# SECTIONS IN THE BASE ARE DROPPED. A `[section]` activates that subcommand,
# so a base carrying `[plot]` would make every rank run plot on top of the
# subcommand given after `--`. Each generated .ini keeps only the base's global
# block (everything above its first `[section]` header) and the variant's keys
# follow it; the dropped sections are listed in a comment. This also keeps the
# generated keys out of any section, where PROfit would silently ignore them.
# Pass subcommand options on the command line instead.
declare -a LABELS=()
if [ -n "$VARIANTS" ]; then
    [ "${#INIS[@]}" -eq 1 ] || { echo "error: --variants takes exactly one base config, got ${#INIS[@]}" >&2; exit 1; }
    VFILE="$(resolve "$VARIANTS")"
    BASE="$(resolve "${INIS[0]}")"
    BASE_NAME="$(basename "${INIS[0]}" .ini)"
    mkdir -p "$BATCH"
    INIS=()
    is_kw()  { case "$1" in exclude|inject-systs|inject-osc|syst-only) return 0 ;; *) return 1 ;; esac; }
    is_num() { [[ "$1" =~ ^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$ ]]; }
    vfail()  { echo "error: $VARIANTS: variant '$label': $*" >&2; exit 1; }
    ANY_SYSTONLY=0
    while IFS= read -r line || [ -n "$line" ]; do
        line="${line%%#*}"
        # shellcheck disable=SC2086
        set -- $line
        [ "$#" -eq 0 ] && continue
        label="$1"; shift
        [ "$#" -ge 1 ] || vfail "no sections (expected exclude / inject-systs / inject-osc / syst-only)"
        for t in "$@"; do
            [ "$t" = inject ] && vfail "bare 'inject' is ambiguous - use inject-systs (systematic pulls) or inject-osc (physics parameters)"
        done
        # Names cannot contain spaces: the line is split on whitespace and
        # quotes are not interpreted. Reject them rather than passing a token
        # like `"POT` through to PROfit.
        for t in "$@"; do
            case "$t" in *\"*|*\'*) vfail "'$t' contains a quote - names with spaces are not supported; list the member XML names instead" ;; esac
        done
        is_kw "$1" || vfail "'$1' is not a section keyword - start with exclude, inject-systs, inject-osc or syst-only (e.g. '$label exclude $*')"
        for l in ${LABELS[@]+"${LABELS[@]}"}; do
            [ "$l" = "$label" ] && { echo "error: duplicate variant label '$label' - outputs would clobber" >&2; exit 1; }
        done

        # Several exclude names are needed whenever one physical effect is
        # split across configs - the ICARUS detector systematics each exist as
        # a _Run2 and a _Run4 entry, and dropping only one answers nothing.
        excl=(); inj=(); osc=(); systonly=0; set_excl=0; set_inj=0; set_osc=0; sect=""
        while [ "$#" -gt 0 ]; do
            tok="$1"; shift
            case "$tok" in
                exclude)   [ "$set_excl" -eq 0 ] || vfail "'exclude' given twice"; sect=exclude; set_excl=1; continue ;;
                inject-systs) [ "$set_inj" -eq 0 ] || vfail "'inject-systs' given twice"; sect=inject-systs; set_inj=1; continue ;;
                inject-osc)   [ "$set_osc" -eq 0 ] || vfail "'inject-osc' given twice";   sect=inject-osc;   set_osc=1; continue ;;
                syst-only) sect=""; systonly=1; continue ;;
            esac
            case "$sect" in
                exclude) excl+=("$tok") ;;
                inject-systs|inject-osc)
                    # name/value pairs, checked here because PROfit's map parser
                    # would otherwise pair the wrong tokens without complaint.
                    { [ "$#" -ge 1 ] && ! is_kw "$1"; } || vfail "$sect '$tok' has no value"
                    is_num "$1" || vfail "$sect '$tok' value '$1' is not a number"
                    if [ "$sect" = inject-systs ]; then inj+=("$tok" "$1"); else osc+=("$tok" "$1"); fi
                    shift ;;
                *) vfail "'$tok' follows syst-only, which takes no arguments" ;;
            esac
        done
        [ "$set_excl" -eq 0 ] || [ "${#excl[@]}" -gt 0 ] || vfail "'exclude' with no names"
        [ "$set_inj"  -eq 0 ] || [ "${#inj[@]}"  -gt 0 ] || vfail "'inject-systs' with no name/sigma pairs"
        [ "$set_osc"  -eq 0 ] || [ "${#osc[@]}"  -gt 0 ] || vfail "'inject-osc' with no param/value pairs"
        [ "$systonly" -eq 1 ] && ANY_SYSTONLY=1

        # Strip from the base only the keys this variant overrides.
        strip="output|log"
        [ "$set_excl" -eq 1 ] && strip+="|exclude-systs"
        [ "$set_inj"  -eq 1 ] && strip+="|inject-systs"
        [ "$set_osc"  -eq 1 ] && strip+="|inject|i"
        [ "$systonly" -eq 1 ] && strip+="|syst-only"

        LABELS+=("$label")
        gen="$BATCH/${BASE_NAME}__${label}.ini"
        globals="$BATCH/.${label}.globals"
        {
            echo "# --- generated by submit_fanout.sh --variants $(basename "$VFILE") ---"
            echo "output = $label"
            [ "$set_excl" -eq 1 ] && echo "exclude-systs = ${excl[*]}"
            [ "$set_inj"  -eq 1 ] && echo "inject-systs = ${inj[*]}"
            [ "$set_osc"  -eq 1 ] && echo "inject = ${osc[*]}"
            [ "$systonly" -eq 1 ] && echo "syst-only = true"
            echo "log = ${label}.log"
        } > "$globals"
        # Keep the base's global block only; record which sections were dropped.
        sections="$( { grep -oE '^[[:space:]]*\[[^]]+\]' "$BASE" || true; } | tr -d ' \t' | tr '\n' ' ')"
        {
            sed -E "/^[[:space:]]*($strip)[[:space:]]*=/d" "$BASE" | awk '/^[[:space:]]*\[/ { exit } { print }'
            echo
            if [ -n "$sections" ]; then echo "# base sections dropped (would run as extra subcommands): $sections"; fi
            cat "$globals"
        } > "$gen"
        rm -f "$globals"
        INIS+=("$gen")
    done < "$VFILE"
    [ "${#INIS[@]}" -gt 0 ] || { echo "error: $VARIANTS contains no variants" >&2; exit 1; }

    # --syst-only pins the physics parameters at CV. PROfit's own surface code
    # notes it makes no sense there: the grid fixes physics point by point
    # anyway, so only the pre-fit global would change.
    if [ "$ANY_SYSTONLY" -eq 1 ] && [ "$SUBCOMMAND" = "surface" ]; then
        echo "warning: syst-only variants with 'surface' - PROfit only applies it to the pre-fit;" >&2
        echo "         it is meant for global/profile. Submitting anyway." >&2
    fi
fi

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
    path="$(resolve "$ini")"
    tag="$(tag_of "$path")"
    [ -n "$tag" ] || { echo "error: $ini declares no tag" >&2; exit 1; }
    if [ "$SHARED_CACHE" -eq 0 ]; then
        for t in ${TAGS[@]+"${TAGS[@]}"}; do
            [ "$t" = "$tag" ] && { echo "error: duplicate tag '$tag' - ranks would share a run directory. Use --shared-cache if that is intended." >&2; exit 1; }
        done
    fi
    TAGS+=("$tag")
done

# ---- shared-cache gate ------------------------------------------------------
if [ "$SHARED_CACHE" -eq 1 ]; then
    for tag in $(printf '%s\n' "${TAGS[@]}" | sort -u); do
        rd="$OUTPUT_ROOT/$tag"
        if [ ! -s "$rd/${tag}_prop.bin" ] || [ ! -s "$rd/${tag}_syst.bin" ]; then
            echo "error: --shared-cache requires an existing cache for tag '$tag'." >&2
            echo "       missing ${tag}_prop.bin / ${tag}_syst.bin in $rd" >&2
            echo "       Ranks sharing a tag with no cache would all process into it at once." >&2
            echo "       Build it first, e.g.:" >&2
            echo "         scripts/profit_run.sh <that config> $SUBCOMMAND" >&2
            exit 1
        fi
    done
    # Create the run directory and its symlinks once here, rather than letting N
    # ranks race on the same mkdir/ln -sfn at startup.
    for tag in $(printf '%s\n' "${TAGS[@]}" | sort -u); do
        mkdir -p "$OUTPUT_ROOT/$tag"
        ln -sfn "$REPO_ROOT/work" "$OUTPUT_ROOT/$tag/work"
        ln -sfn "$REPO_ROOT/xml"  "$OUTPUT_ROOT/$tag/xml"
    done
fi

mkdir -p "$BATCH"
TASKS="$BATCH/tasks.txt"
: > "$TASKS"
for ini in "${INIS[@]}"; do
    printf '%s %s %s %s\n' "$REPO_ROOT/scripts/profit_run.sh" "$(resolve "$ini")" "$SUBCOMMAND" "${EXTRA[*]}" >> "$TASKS"
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
echo "tasks:      $N run(s), subcommand '$SUBCOMMAND' ${EXTRA[*]:-}"
for i in "${!INIS[@]}"; do
    printf '  rank %-2s %-52s tag=%s\n' "$i" "$(basename "${INIS[$i]}")" "${TAGS[$i]}"
done
[ "$SHARED_CACHE" -eq 1 ] && echo "cache:      shared, verified present for $(printf '%s\n' "${TAGS[@]}" | sort -u | tr '\n' ' ')"
echo "submitting: $N node(s), walltime $WALLTIME, queue $QUEUE"
cd "$BATCH"
qsub "$JOB_SCRIPT"
