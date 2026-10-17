#!/usr/bin/env python3
"""
Ratios and fractional detector systematics for every smoothing variant of the closure study.

For one observable of the smoothing study (``muon_score_inclusive_s2`` or
``reco_energy_inclusive_s2``) this draws, for one role permutation of the three-way event split:

* ``ratios_<study>``: the varied / nominal ratio of each systematic under test as made from the
  ratio part (unsmoothed) and after each smoother, with the production full-sample ratio for the
  reco-energy study;
* ``frac_syst_<study>``: the standard PROfit "Detector Systematics" fractional-uncertainty panel
  for every variant (one line per systematic and their sum), read from the ``frac_syst`` tree PROfit
  writes next to every global fit.

Usage
-----
    python plotting/studies/detvar_smoothing_figures.py reco_energy_inclusive_s2 \
        --variants unsmoothed spline_gcv local3 poly_cv roll3 roll5 [--perm 0]

Outputs go to ``docs/detvar_smoothing/figures/``.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO, configure_style, save_figure  # noqa: E402

RUN = REPO / "runs" / "smoothing"
OUT = REPO / "docs" / "detvar_smoothing" / "figures"
SYSTS = ["ChargeScale15", "MattModYZ", "DENT", "WireModxThetaXW"]
COLORS = {"ChargeScale15": "#1f77b4", "MattModYZ": "#d62728", "DENT": "#2ca02c",
          "WireModxThetaXW": "#9467bd", "TriggerEmulation": "#8c564b"}
AXIS = {"muon_score_inclusive_s2": "Reconstructed muon score",
        "reco_energy_inclusive_s2": "Reconstructed visible energy [GeV]"}
PRODUCTION = REPO / "work" / "detsys_main.root"


def draw_ratios(study: str, variants: list[str], perm: int) -> None:
    """Ratio of every systematic under test, for every variant."""
    fig, axes = plt_subplots(len(SYSTS))
    for ax, syst in zip(axes.ravel(), SYSTS, strict=True):
        for index, variant in enumerate(variants):
            path = RUN / "ratios" / f"smm_{study}_{variant}_p{perm}.root"
            values, edges = uproot.open(path)[f"{syst}_sbnd_Run1"].to_numpy()
            centers = 0.5 * (edges[1:] + edges[:-1])
            if variant == "unsmoothed":
                ax.step(edges, np.append(values, values[-1]), where="post", color="0.55", lw=1.2,
                        label="Unsmoothed (one third)")
            else:
                ax.plot(centers, values, lw=1.8, label=variant, color=f"C{index}")
        if study.startswith("reco_energy"):
            values, edges = uproot.open(PRODUCTION)[f"{syst}_sbnd_Run1"].to_numpy()
            ax.step(edges, np.append(values, values[-1]), where="post", color="k", lw=1.0,
                    ls="--", label="Production (full sample)")
        ax.axhline(1.0, color="0.7", lw=0.8)
        ax.set_title(syst, fontsize=11)
        ax.set_xlabel(AXIS[study])
        ax.set_ylabel("Varied / nominal")
    axes[0, 0].legend(fontsize=8, ncol=1)
    save_figure(fig, OUT / f"ratios_{study}_p{perm}.png")


def plt_subplots(count: int):
    """A grid of panels (two columns) with the shared style applied."""
    import matplotlib.pyplot as plt

    configure_style()
    rows = (count + 1) // 2
    return plt.subplots(rows, 2, figsize=(11, 3.6 * rows), squeeze=False)


def fractional(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per-systematic fractional uncertainty (edges, values) of the Detector Systematics tag."""
    tree = uproot.open(path)["frac_syst"].arrays(library="np")
    keep = np.array([str(tag) == "Detector Systematics" for tag in tree["tag"]])
    out = {}
    for name in dict.fromkeys(str(n) for n in tree["systname"][keep]):
        rows = keep & np.array([str(n) == name for n in tree["systname"]])
        low, high = tree["bin_low_edge"][rows], tree["bin_high_edge"][rows]
        out[name.replace("_sbnd_Run1", "")] = (np.append(low, high[-1]), tree["bin_content"][rows])
    return out


def draw_fractional(study: str, variants: list[str], perm: int) -> None:
    """The Detector Systematics fractional panel for every variant."""
    fig, axes = plt_subplots(len(variants))
    for ax, variant in zip(axes.ravel(), variants, strict=False):
        path = RUN / f"smm_{study}_{variant}_p{perm}_null_nominal_PROplot_histograms.root"
        for name, (edges, values) in fractional(path).items():
            total = name == "SUM"
            ax.step(edges, np.append(values, values[-1]), where="post", lw=2.2 if total else 1.3,
                    color="k" if total else COLORS.get(name, "0.5"), label=name)
        ax.set_title(variant, fontsize=11)
        ax.set_xlabel(AXIS[study])
        ax.set_ylabel("Fractional uncertainty")
        ax.set_ylim(bottom=0)
    axes[0, 0].legend(fontsize=7)
    save_figure(fig, OUT / f"frac_syst_{study}_p{perm}.png")


def main() -> None:
    """Draw both figures."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("study", choices=sorted(AXIS))
    parser.add_argument("--variants", nargs="+", required=True)
    parser.add_argument("--perm", type=int, default=0)
    args = parser.parse_args()
    draw_ratios(args.study, args.variants, args.perm)
    draw_fractional(args.study, args.variants, args.perm)


if __name__ == "__main__":
    sys.exit(main())
