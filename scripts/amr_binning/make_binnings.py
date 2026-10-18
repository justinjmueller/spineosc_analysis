#!/usr/bin/env python3
"""
Generate binning variants of the AMR-comparison contour configurations.

Takes a base XML of ``runs/amr_compare/xml`` (the 1D configurations: ``contours_incl1D`` and ``contours_excl1D``) and
rewrites the edges of each channel's fit variable (the first ``<bins unit="RecoE" edges="...">`` of a channel that is not a
``plot="false"`` shadow), keeping the range of every channel. Everything else (samples, systematics, shadow variables) is
untouched, so the variants differ from the baseline in the binning only.

Variants (applied per channel, to that channel's baseline edges ``e_0 ... e_n``)::

    base         the baseline edges
    merge2       every second edge (pairs of bins merged), the last edge kept
    merge4       every fourth edge, the last edge kept
    uniform_W    uniform bins of width W GeV from e_0 to e_n (the last bin absorbs the remainder)
    fine_W       uniform bins of width W up to 1.5 GeV, the baseline edges above (the tail has too few events for fine bins)
    equalpop_N   N bins with equal expected populations, from a fine-binned histogram (``--population``)

Usage::

    python scripts/amr_binning/make_binnings.py --base runs/amr_compare/xml/contours_incl1D.xml --name incl1D \\
        --out runs/amr_binning/xml --variants base merge2 merge4 uniform_0.10 uniform_0.05
    python scripts/amr_binning/make_binnings.py --list-edges runs/amr_compare/xml/contours_excl1D.xml
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

#: Above this energy [GeV] the ``fine_W`` variant keeps the baseline's coarse edges (the high-energy tail has too few events).
FINE_LIMIT = 1.5

#: The fit variable's bins of a channel: not a ``plot="false"`` shadow, not a 2D binning.
FIT_BINS = re.compile(r'(<bins\s+unit="RecoE"\s+edges=")([^"]+)("\s*/>)')
CHANNEL = re.compile(r'<channel\s+name="([^"]+)"')


def channel_blocks(text: str) -> list[tuple[str, int, int]]:
    """(name, start, end) of each ``<channel ...>`` ... ``</channel>`` block."""
    out = []
    for match in CHANNEL.finditer(text):
        end = text.index("</channel>", match.end())
        out.append((match.group(1), match.start(), end))
    return out


def fit_edges(block: str) -> tuple[re.Match, np.ndarray]:
    """The fit-variable bins match of a channel block and its edges."""
    for match in FIT_BINS.finditer(block):
        if 'plot="false"' not in match.group(0):
            return match, np.array([float(x) for x in match.group(2).split()])
    raise ValueError("no fit-variable RecoE bins in channel block")


def make_edges(variant: str, edges: np.ndarray, population=None) -> np.ndarray:
    """New edges for one channel from its baseline edges."""
    lo, hi = edges[0], edges[-1]
    if variant == "base":
        return edges
    if variant.startswith("merge"):
        step = int(variant[5:])
        kept = list(edges[::step])
        if kept[-1] != hi:
            kept.append(hi)
        return np.array(kept)
    if variant.startswith("uniform_"):
        width = float(variant.split("_", 1)[1])
        n = max(1, int(round((hi - lo) / width)))
        return np.round(np.linspace(lo, hi, n + 1), 4)
    if variant.startswith("fine_"):
        width = float(variant.split("_", 1)[1])
        top = min(hi, FINE_LIMIT)
        low = lo + width * np.arange(int(round((top - lo) / width)) + 1)
        return np.round(np.r_[low, [e for e in edges if e > FINE_LIMIT + 1e-9]], 4)
    if variant.startswith("equalpop_"):
        n = int(variant.split("_", 1)[1])
        if population is None:
            raise ValueError("equalpop needs --population (centres, content) from a fine-binned histogram")
        centres, content = population
        inside = (centres > lo) & (centres < hi)
        cumulative = np.cumsum(content[inside]) / content[inside].sum()
        cut = np.interp(np.linspace(0, 1, n + 1)[1:-1], cumulative, centres[inside])
        return np.round(np.r_[lo, cut, hi], 4)
    raise ValueError(f"unknown variant {variant}")


def load_population(path: Path, channel_names: list[str]) -> dict:
    """
    Population per channel from a fine-binned PROfit ``*_PROplot_histograms.root``: the CV total summed over the detectors.

    ``channel_names`` gives the channel order of the XML (PROfit numbers channels in that order).
    """
    import uproot

    total = 4294967295
    with uproot.open(path) as handle:
        tree = {k: handle["hist1d"][k].array(library="np") for k in
                ("variable", "detector", "channel", "subchannel", "prefix", "bin_low_edge", "bin_high_edge", "bin_content")}
    out = {}
    for index, name in enumerate(channel_names):
        sel = (tree["variable"] == 0) & (tree["channel"] == index) & (tree["subchannel"] == total) & (tree["prefix"] == "CV")
        lows = np.unique(tree["bin_low_edge"][sel])
        centres, content = [], []
        for low in lows:
            hit = sel & (tree["bin_low_edge"] == low)
            first = np.where(hit)[0]
            per_detector = {}
            for i in first:  # the rows repeat per detector; sum one row per detector
                per_detector.setdefault(int(tree["detector"][i]), float(tree["bin_content"][i]))
            centres.append(0.5 * (low + float(tree["bin_high_edge"][first[0]])))
            content.append(sum(per_detector.values()))
        out[name] = (np.array(centres), np.array(content))
    return out


def rewrite(text: str, variant: str, population=None) -> tuple[str, dict]:
    """The XML text with every channel's fit edges replaced; also returns the new edges per channel."""
    if variant == "base":  # byte-identical, so the existing AMR results and caches are the baseline
        return text, {name: [float(x) for x in fit_edges(text[s:e])[1]] for name, s, e in channel_blocks(text)}
    result, cursor, summary = [], 0, {}
    for name, start, end in channel_blocks(text):
        block = text[start:end]
        match, edges = fit_edges(block)
        new = make_edges(variant, edges, None if population is None else population.get(name))
        result.append(text[cursor:start])
        result.append(block[:match.start(2)] + " ".join(f"{e:g}" for e in new) + block[match.end(2):])
        cursor = end
        summary[name] = [float(e) for e in new]
    result.append(text[cursor:])
    return "".join(result), summary


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", type=Path, help="baseline XML")
    parser.add_argument("--name", help="config name used in the output file names (e.g. incl1D)")
    parser.add_argument("--out", type=Path, default=Path("runs/amr_binning/xml"))
    parser.add_argument("--variants", nargs="+", default=["base", "merge2", "merge4", "uniform_0.10", "uniform_0.05"])
    parser.add_argument("--population", type=Path, help="fine-binned *_PROplot_histograms.root for the equalpop_N variants")
    parser.add_argument("--list-edges", type=Path, help="print the fit edges of each channel of an XML and exit")
    args = parser.parse_args()

    if args.list_edges:
        for name, start, end in channel_blocks(args.list_edges.read_text()):
            _, edges = fit_edges(args.list_edges.read_text()[start:end])
            print(f"{name}: {len(edges) - 1} bins  {' '.join(f'{e:g}' for e in edges)}")
        return
    text = args.base.read_text()
    args.out.mkdir(parents=True, exist_ok=True)
    population = None
    if args.population:
        population = load_population(args.population, [name for name, _, _ in channel_blocks(text)])
    for variant in args.variants:
        new, summary = rewrite(text, variant, population)
        path = args.out / f"{args.name}_{variant}.xml"
        path.write_text(new)
        counts = ", ".join(f"{k}: {len(v) - 1} bins" for k, v in summary.items())
        print(f"{path.name}  ({counts})")


if __name__ == "__main__":
    main()
