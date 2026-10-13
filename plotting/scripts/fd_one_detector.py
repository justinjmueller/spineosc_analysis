#!/usr/bin/env python3
"""One-detector injections and all-at-once combinations.

one_detector.pdf / tables/one_detector.csv: Delta chi2 and null-fit chi2 when a dial is injected into both detectors (plain dial), SBND only, or ICARUS only.
tables/combinations.csv: the all-at-once datasets (theory5 = sub-threshold CPI -100%, asymmetric MEC, low-Q2 pion pre/post-FSI, hA2025; meas3 = uB CC1p0pi,
uB CC2p0pi, T2K NC1pi) against the sum of the individual datasets' chi2 and Delta chi2.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_one_detector.py
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
rows = []
for d in DIALS:
    for c in ("b", "1p"):
        for t, s in (("both detectors", d), ("SBND only", d + "_onlySBND"), ("ICARUS only", d + "_onlyICARUS")):
            p = f2.pair("ds", c, s)
            if p:
                rows.append({"dial": d, "config": c, "injected in": t, "chi2 null": p[0]["chi2"], "dchi2": p[0]["chi2"] - p[1]["chi2"]})
df = pd.DataFrame(rows)
df.to_csv(OUT / "tables" / "one_detector.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
w = 0.26
for ax, c in zip(axes, ("b", "1p")):
    for i, t in enumerate(("both detectors", "SBND only", "ICARUS only")):
        s = df[(df.config == c) & (df["injected in"] == t)].set_index("dial").reindex(list(DIALS))
        ax.bar(np.arange(len(DIALS)) + (i - 1) * w, s.dchi2.clip(lower=0.01), w, label=t)
    ax.axhline(4, color="0.5", ls=":", lw=0.8); ax.axhline(1, color="0.5", ls=":", lw=0.8); ax.set_yscale("log")
    ax.set_xticks(np.arange(len(DIALS))); ax.set_xticklabels(list(DIALS.values()), rotation=25, ha="right", fontsize=7); ax.set_title(f2.CONFIGS[c], fontsize=9)
axes[0].set_ylabel(r"$\Delta\chi^2$ (null - oscillation)"); axes[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT / "one_detector.pdf", bbox_inches="tight")
combos = {"theory5": ["cpi_subthreshold_minus100", "mec_bdt", "jaesung_lowq2_pi_enhancement_prefsi", "jaesung_lowq2_pi_enhancement_postfsi", "pi_fsi_ha2025"],
          "meas3": ["ub_cc1p0pi_dpt_dat", "ub_cc2p0pi_delta_pt", "t2k_nc1pi_p_costh"]}
rows = []
for name, ds in combos.items():
    for c in ("b", "1p", "Np"):
        p = f2.pair("ds", c, name)
        s = [f2.pair("ds", c, x) for x in ds]
        if p and all(s):
            rows.append({"combination": name, "config": c, "chi2 null (all at once)": p[0]["chi2"], "dchi2 (all at once)": p[0]["chi2"] - p[1]["chi2"],
                         "sum of single chi2 null": sum(x[0]["chi2"] for x in s), "sum of single dchi2": sum(x[0]["chi2"] - x[1]["chi2"] for x in s)})
pd.DataFrame(rows).to_csv(OUT / "tables" / "combinations.csv", index=False)
print(df.pivot_table(index=["dial", "config"], columns="injected in", values="dchi2").round(2).to_string())
