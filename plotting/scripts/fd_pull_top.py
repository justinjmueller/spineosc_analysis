#!/usr/bin/env python3
"""Top-10 nuisance pulls per fake dataset, in the style of the GUMP/GUMPLE slides: for the ten parameters with the largest pull in the null fit,
the CV closure fit (grey), the null fit (filled) and the oscillation-floating fit (open), error bars = MCMC posterior width, prior band +-1 sigma.

Combined analysis, data-XML route. Writes plotting/fakedata_2d_figures/pull_top10_<n>.pdf (six datasets per page) and tables/pull_top10.csv.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_pull_top.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fd2d as f2
from fd_figures2 import PANEL_STEMS

OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
TITLE = {"cpi_subthreshold_minus100": "sub-threshold CPI -100%", "cpi_subthreshold_plus100": "sub-threshold CPI +100%", "mec_bdt": "asymmetric MEC (BDT)",
         "jaesung_lowq2_pi_enhancement_prefsi": "low-Q2 pion, pre-FSI", "jaesung_lowq2_pi_enhancement_postfsi": "low-Q2 pion, post-FSI",
         "pi_fsi_ha2025": "pion FSI hA2025", "ub_cc1p0pi_dpt_dat_midW": "uB CC1p0pi mid-W", "ub_cc1p0pi_dpt_dat": "uB CC1p0pi all-W",
         "t2k_nc1pi_p_costh": "T2K NC1pi", "minerva_3dqelike_bnb": "MINERvA 3D QE-like", "qe_zexp_mva_to_lqcd": "LQCD"}
N = 10
ctrl = f2.read_fit("ds", "b", "fix", "ctrl")
rows, panels = [], []
for s in PANEL_STEMS:
    p = f2.pair("ds", "b", s)
    if not p:
        continue
    n, o = p
    names = sorted(f2.nuisance_names(n), key=lambda k: -abs(n["values"][k]))[:N]
    panels.append((s, names, n, o))
    for k in names:
        rows.append({"dataset": s, "parameter": k, "closure": ctrl["values"].get(k), "null pull": n["values"][k], "null width": n["sigma"].get(k),
                     "osc pull": o["values"].get(k), "osc width": o["sigma"].get(k)})
pd.DataFrame(rows).to_csv(OUT / "tables" / "pull_top10.csv", index=False)
for page, chunk in enumerate([panels[i:i + 6] for i in range(0, len(panels), 6)], start=1):
    fig, axes = plt.subplots(3, 2, figsize=(11, 11.5))
    for ax, (s, names, n, o) in zip(axes.flat, chunk):
        x = np.arange(len(names))
        ax.axhspan(-1, 1, color="0.92", zorder=0); ax.axhline(0, color="0.4", lw=0.6)
        ax.errorbar(x - 0.2, [ctrl["values"].get(k, np.nan) for k in names], yerr=[ctrl["sigma"].get(k, 0) for k in names], fmt="s", ms=4, color="0.35", lw=0.8, label="CV (closure)")
        ax.errorbar(x, [n["values"][k] for k in names], yerr=[n["sigma"].get(k, 0) for k in names], fmt="o", ms=5, color="C0", lw=0.9, label="fake data, oscillation fixed")
        ax.errorbar(x + 0.2, [o["values"].get(k, np.nan) for k in names], yerr=[o["sigma"].get(k, 0) for k in names], fmt="s", ms=5, mfc="none", color="C0", lw=0.9, label="fake data, oscillation floating")
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=40, ha="right", fontsize=7); ax.set_ylim(-1.6, 1.6)
        ax.set_title(f"{TITLE.get(s, s)}  ($\\chi^2_{{null}}$ = {n['chi2']:.1f})", fontsize=9); ax.set_ylabel(r"post-fit value [prior $\sigma$]", fontsize=8)
    for ax in list(axes.flat)[len(chunk):]:
        ax.axis("off")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(OUT / f"pull_top10_{page}.pdf", bbox_inches="tight"); plt.close(fig)
print("wrote", page, "pages;", len(panels), "datasets")
