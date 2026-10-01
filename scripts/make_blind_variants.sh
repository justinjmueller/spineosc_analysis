#!/usr/bin/env bash
#
# Generates the blinded variants of a contour config from a single source XML.
#
# Usage: scripts/make_blind_variants.sh [xml/contours/contours_base.xml]
#
# Blinding is a 100% normalization systematic on one detector's subchannels: a
# `norm` allowlist entry with fraction 1.0 lets the fit rescale that detector
# freely, which discards its normalization information and leaves it
# contributing shape only. That gives three scenarios:
#
#   contours_blind.xml         SBND + ICARUS both shape-only
#   contours_sbndunblind.xml   ICARUS shape-only, SBND normalization used
#   <source>                   neither blinded - full sensitivity
#
# The variants are generated rather than committed because they differ from the
# source by one or two lines, and hand-maintained near-duplicates of a
# 2000-line XML drift.
#
# NOTE ON THE PATTERN: matching is an unanchored regex against subchannel
# fullnames (<mode>_<detector>_<channel>_<subchannel>). The detectors here are
# ICARUSRun2 and ICARUSRun4, so "nu_ICARUS" matches but "nu_ICARUS_" does NOT -
# the character after ICARUS is R, not an underscore. A non-matching norm
# pattern aborts this PROfit build (PROcreate.cxx: "matches NO subchannel
# fullname"), but older builds silently ignored it and ran effectively
# unblinded, so do not relax this without checking the match count.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${1:-$REPO_ROOT/xml/contours/contours_base.xml}"
OUT_DIR="$(dirname "$SRC")/generated"

[ -f "$SRC" ] || { echo "error: source XML not found: $SRC" >&2; exit 1; }

SBND_BLIND='    <allowlist type="norm" binning="var0" tag="Det" plotname="SBNDShapeBlind">nu_SBND_:1.0</allowlist>'
ICARUS_BLIND='    <allowlist type="norm" binning="var0" tag="Det" plotname="ICARUSShapeBlind">nu_ICARUS:1.0</allowlist>'

# Anchor on the first existing norm entry so the new ones land inside
# <variation_list> next to their own kind.
ANCHOR="$(grep -n '<allowlist type="norm"' "$SRC" | head -1 | cut -d: -f1)"
[ -n "$ANCHOR" ] || { echo "error: no <allowlist type=\"norm\"> anchor in $SRC" >&2; exit 1; }

mkdir -p "$OUT_DIR"

emit() {  # out_path, lines-to-insert...
    local out="$1"; shift
    awk -v anchor="$ANCHOR" -v ins="$*" 'NR==anchor{print ins} {print}' "$SRC" > "$out"
}

emit "$OUT_DIR/contours_blind.xml" "$SBND_BLIND"$'\n'"$ICARUS_BLIND"
emit "$OUT_DIR/contours_sbndunblind.xml" "$ICARUS_BLIND"

# Verify each variant carries the norm entries it should, and no others.
check() {
    local f="$1" want="$2"
    local got
    got="$(grep -c '<allowlist type="norm"' "$f")"
    if [ "$got" != "$want" ]; then
        echo "error: $(basename "$f") has $got norm entries, expected $want" >&2
        exit 1
    fi
    printf '  %-34s norm entries: %s\n' "$(basename "$f")" \
        "$(grep -oP '<allowlist type="norm"[^>]*>\K[^<]*' "$f" | tr '\n' ' ')"
}

BASE_NORMS="$(grep -c '<allowlist type="norm"' "$SRC")"
echo "generated in $OUT_DIR (source has $BASE_NORMS norm entries):"
check "$OUT_DIR/contours_blind.xml" "$((BASE_NORMS + 2))"
check "$OUT_DIR/contours_sbndunblind.xml" "$((BASE_NORMS + 1))"
echo
echo "source ($(basename "$SRC")) is the unblinded full-sensitivity config; use it directly."
