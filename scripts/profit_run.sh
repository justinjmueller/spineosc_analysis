#!/usr/bin/env bash
#
# Runs one PROfit subcommand for one job config.
#
# Usage: scripts/profit_run.sh <jobs/NAME.ini> <subcommand> [extra PROfit args...]
#
#   scripts/profit_run.sh jobs/varset1_mc.ini process
#   scripts/profit_run.sh jobs/varset1_mc.ini plot
#   scripts/profit_run.sh jobs/varset1_mc.ini plot --with-splines
#   scripts/profit_run.sh jobs/varset1_mc.ini process --dry-run
#
# --dry-run runs every check and prints the PROfit command without executing
# it. Worth using on a login node, where a stray `process` would start reading
# tens of GB of trees.
#
# Job settings live in the .ini (read by PROfit's own --config); machine paths
# live in input/site.conf. Arguments after the subcommand are passed straight
# through and override the .ini, since PROfit prefers the command line.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SITE_CONF="$REPO_ROOT/input/site.conf"

if [ "$#" -lt 2 ]; then
    echo "usage: scripts/profit_run.sh <jobs/NAME.ini> <subcommand> [extra PROfit args...]" >&2
    exit 2
fi

INI_ARG="$1"; shift
SUBCOMMAND="$1"; shift

DRY_RUN=0
PASSTHROUGH=()
for arg in "$@"; do
    if [ "$arg" = "--dry-run" ]; then
        DRY_RUN=1
    else
        PASSTHROUGH+=("$arg")
    fi
done
set -- ${PASSTHROUGH[@]+"${PASSTHROUGH[@]}"}

INI="$INI_ARG"
[ -f "$INI" ] || INI="$REPO_ROOT/$INI"
if [ ! -f "$INI" ]; then
    echo "error: job config not found: $INI_ARG" >&2
    exit 1
fi
INI="$(cd "$(dirname "$INI")" && pwd)/$(basename "$INI")"

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

PROFIT_BIN="$(get_setting profit_bin)"
OUTPUT_ROOT="$(get_setting output_root)"

for setting in PROFIT_BIN OUTPUT_ROOT; do
    if [ -z "${!setting}" ]; then
        echo "error: ${setting,,} is not set in $SITE_CONF" >&2
        exit 1
    fi
done

PROFIT="$PROFIT_BIN/PROfit"
if [ ! -x "$PROFIT" ]; then
    echo "error: PROfit not executable at $PROFIT (check profit_bin in $SITE_CONF)" >&2
    exit 1
fi

# PROfit segfaults rather than erroring cleanly when a work/<key>.root the XML
# references is absent, so fail here with something readable instead.
if [ ! -d "$REPO_ROOT/work" ] || [ -z "$(ls -A "$REPO_ROOT/work" 2>/dev/null)" ]; then
    echo "error: $REPO_ROOT/work is missing or empty." >&2
    echo "       run scripts/deploy.sh first (see README 'First-time setup')." >&2
    exit 1
fi

# The run directory shadows the repo: PROfit writes its outputs and .bin caches
# here, while work/ and xml/ symlinks let the repo-relative paths in the .ini
# resolve. That keeps GB-scale artifacts out of the repo without having to
# rewrite any paths.
#
# Keyed on the config's tag, not the .ini filename, because PROfit names its
# caches <tag>_prop.bin / <tag>_syst.bin / <tag>_data.bin. Two .ini files
# sharing a tag therefore share one cache, which is what lets an MC-only and a
# data/MC config reuse the same expensive processing run.
TAG="$(grep -oP '^\s*tag\s*=\s*\K\S+' "$INI" | tail -1)"
RUN_DIR="$OUTPUT_ROOT/${TAG:-$(basename "$INI" .ini)}"
mkdir -p "$RUN_DIR"
ln -sfn "$REPO_ROOT/work" "$RUN_DIR/work"
ln -sfn "$REPO_ROOT/xml" "$RUN_DIR/xml"

echo "config:     $INI"
echo "subcommand: $SUBCOMMAND"
echo "run dir:    $RUN_DIR"
echo "PROfit:     $PROFIT"
echo

cd "$RUN_DIR"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry run] all checks passed; would execute from $RUN_DIR:"
    printf '  '; printf '%q ' "$PROFIT" --config "$INI" "$SUBCOMMAND" "$@"; echo
    exit 0
fi

exec "$PROFIT" --config "$INI" "$SUBCOMMAND" "$@"
