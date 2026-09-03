#!/usr/bin/env bash
#
# Populates work/ with symlinks named after the logical file keys the XML
# configs under xml/ reference (work/<key>.root). Run this once per machine
# after editing input/sources.list to point at your local copies of the
# data/detsys files. Re-run any time input/sources.list changes.
#
# Usage: scripts/deploy.sh [path/to/sources.list]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCES_FILE="${1:-$REPO_ROOT/input/sources.list}"
WORK_DIR="$REPO_ROOT/work"

if [ ! -f "$SOURCES_FILE" ]; then
    echo "error: sources file not found: $SOURCES_FILE" >&2
    echo "       copy input/sources.list.example to input/sources.list and edit it for this machine." >&2
    exit 1
fi

mkdir -p "$WORK_DIR"

declare -a missing=()
declare -a linked=()

while IFS= read -r line || [ -n "$line" ]; do
    # strip comments and surrounding whitespace
    line="${line%%#*}"
    line="$(echo "$line" | xargs)"
    [ -z "$line" ] && continue

    if [[ "$line" != *"="* ]]; then
        echo "warning: skipping malformed line in $SOURCES_FILE: $line" >&2
        continue
    fi

    key="${line%%=*}"
    src="${line#*=}"

    if [ -z "$key" ] || [ -z "$src" ]; then
        echo "warning: skipping malformed line in $SOURCES_FILE: $line" >&2
        continue
    fi

    dest="$WORK_DIR/${key}.root"
    ln -sf "$src" "$dest"
    linked+=("$key")

    if [ ! -e "$src" ]; then
        missing+=("$key -> $src")
    fi
done < "$SOURCES_FILE"

echo "Linked ${#linked[@]} source(s) into $WORK_DIR:"
for key in "${linked[@]}"; do
    printf '  %-30s -> %s\n' "work/${key}.root" "$(readlink "$WORK_DIR/${key}.root")"
done

# Prune symlinks left behind by a since-renamed or since-removed key, so
# work/ never accumulates stale aliases from an earlier sources.list.
declare -a pruned=()
for existing in "$WORK_DIR"/*.root; do
    [ -e "$existing" ] || [ -L "$existing" ] || continue
    [ -L "$existing" ] || continue
    base="$(basename "$existing" .root)"
    found=0
    for key in "${linked[@]}"; do
        if [ "$key" = "$base" ]; then
            found=1
            break
        fi
    done
    if [ "$found" -eq 0 ]; then
        rm -f "$existing"
        pruned+=("$base")
    fi
done

if [ "${#pruned[@]}" -gt 0 ]; then
    echo
    echo "Pruned ${#pruned[@]} stale symlink(s) no longer in $SOURCES_FILE:"
    for key in "${pruned[@]}"; do
        echo "  work/${key}.root"
    done
fi

if [ "${#missing[@]}" -gt 0 ]; then
    echo
    echo "warning: ${#missing[@]} symlink(s) point at a file that does not exist:" >&2
    for entry in "${missing[@]}"; do
        echo "  $entry" >&2
    done
    exit 1
fi

echo
echo "All source files present. Run PROfit from $REPO_ROOT so the XML configs' relative work/ paths resolve."
