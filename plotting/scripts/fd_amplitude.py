#!/usr/bin/env python3
"""Amplitude scan: null-fit chi2 and Delta chi2 against the injected amplitude a, where the fake weight is w' = max(0, 1 + a (w - 1)).

a = 0 is the closure (CV), a = 1 the plain dial, a = -1 the reflected dial. Stems: <dial>_<ampm1p0|amp0p5|amp1p5|amp2p0>; a = 1 is <dial>, a = 0 is ctrl.
Writes plotting/fakedata_2d_figures/amplitude_scan.pdf and tables/amplitude_scan.csv.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_amplitude.py
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

OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
DIALS = {"ub_cc1p0pi_dpt_dat_midW": "uB CC1p0pi mid-W", "mec_bdt": "asymmetric MEC (BDT)", "minerva_3dqelike_bnb": "MINERvA 3D QE-like",
         "jaesung_lowq2_pi_enhancement_prefsi": "low-Q2 pion, pre-FSI", "cpi_subthreshold_minus100": "sub-threshold CPI -100%"}
AMPS = {-1.0: "ampm1p0", 0.5: "amp0p5", 1.5: "amp1p5", 2.0: "amp2p0"}
CFGS = ["b", "1p", "Np"]


def stem(dial, a):
    return "ctrl" if a == 0 else (dial if a == 1 else f"{dial}_{AMPS[a]}")


rows = []
for d in DIALS:
    for a in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0]:
        for c in CFGS:
            p = f2.pair("ds", c, stem(d, a))
            if p:
                nm = f2.nuisance_names(p[0])
                j = max(nm, key=lambda k: abs(p[0]["values"][k]))
                rows.append({"dial": d, "amplitude": a, "config": c, "chi2 null": p[0]["chi2"], "chi2 osc": p[1]["chi2"],
                             "dchi2": p[0]["chi2"] - p[1]["chi2"], "max pull": p[0]["values"][j], "max-pull parameter": j})
df = pd.DataFrame(rows)
df.to_csv(OUT / "tables" / "amplitude_scan.csv", index=False)
fig, axes = plt.subplots(2, len(DIALS), figsize=(3.3 * len(DIALS), 6.2), sharex=True)
for k, (d, lab) in enumerate(DIALS.items()):
    for c in CFGS:
        s = df[(df.dial == d) & (df.config == c)].sort_values("amplitude")
        if len(s):
            axes[0, k].plot(s.amplitude, s["chi2 null"], "o-", ms=3, color=f2.CFG_COLOR[c], label=f2.CONFIGS[c])
            axes[1, k].plot(s.amplitude, s.dchi2, "o-", ms=3, color=f2.CFG_COLOR[c])
    axes[0, k].set_title(lab, fontsize=8); axes[0, k].set_yscale("symlog", linthresh=1)
    axes[1, k].axhline(1, color="0.7", lw=0.6, ls=":"); axes[1, k].axhline(4, color="0.7", lw=0.6, ls=":")
    axes[1, k].set_xlabel("injected amplitude a (1 = plain dial)", fontsize=8)
axes[1, len(DIALS) - 1].text(0.97, 0.05, "a > 1: weight clipped at 0\n(= a = 1)", transform=axes[1, len(DIALS) - 1].transAxes, ha="right", fontsize=6)
axes[0, 0].set_ylabel(r"null-fit $\chi^2$"); axes[1, 0].set_ylabel(r"$\Delta\chi^2$ (null - osc)"); axes[0, 0].legend(fontsize=6)
fig.tight_layout(); fig.savefig(OUT / "amplitude_scan.pdf", bbox_inches="tight")
print(df.pivot_table(index=["dial", "config"], columns="amplitude", values="chi2 null").round(1).to_string())
