#!/usr/bin/env python3
"""Group attribution: increase of the null-fit chi2 when one systematic group (a PROfit tag) is removed from the fit.

Input: runs/fakedata_2d/jobs_grp.txt fits (null fits with `--exclude-systs <tag>`), parsed from logs/grp_*.out, against the full-model null fit.
Other groups can compensate, so each increase is a lower bound on how much the group matters. Writes
plotting/fakedata_2d_figures/group_attribution.pdf and tables/group_attribution.csv.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_group_attribution.py
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fd2d as f2

RUN = f2.RUN
OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
LABEL = {"cpi_subthreshold_minus100": "sub-threshold CPI -100%", "jaesung_lowq2_pi_enhancement_prefsi": "low-Q2 pion, pre-FSI", "mec_bdt": "asymmetric MEC (BDT)",
         "minerva_3dqelike_bnb": "MINERvA 3D QE-like", "ub_cc1p0pi_dpt_dat_midW": "uB CC1p0pi mid-W"}
GROUPS = {"MEC": "Cross-Section-MEC", "RES": "Cross-Section-RES", "FSI": "Cross-Section-FSI", "QE": "Cross-Section-QE", "Det": "Detector Systematics", "Flux": "Flux"}
rows = []
for line in open(RUN / "jobs_grp.txt"):
    cfg, stem, slug, _ = line.strip().split("|")
    out = RUN / "logs" / f"grp_d2{cfg}_{stem}_x{slug}.out"
    if not out.exists():
        continue
    txt = out.read_text(errors="ignore")
    m = re.findall(r"Global Best Fit chi\^2:\s*([0-9.eE+-]+)", txt)
    n = re.findall(r"resolves to (\d+) systematic", txt)
    if m:
        base = f2.read_fit("ds", cfg, "fix", stem)["chi2"]
        rows.append({"config": cfg, "dataset": stem, "group": slug, "n excluded": int(n[0]) if n else np.nan, "chi2 full model": base,
                     "chi2 group removed": float(m[-1]), "increase": float(m[-1]) - base, "relative increase [%]": 100 * (float(m[-1]) / base - 1)})
df = pd.DataFrame(rows)
df.to_csv(OUT / "tables" / "group_attribution.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
w = 0.13
for ax, cfg in zip(axes, ("b", "1p")):
    for i, g in enumerate(GROUPS):
        s = df[(df.config == cfg) & (df.group == g)].set_index("dataset").reindex(list(LABEL))
        ax.bar(np.arange(len(LABEL)) + (i - 2.5) * w, s["relative increase [%]"], w, label=f"{g} ({int(s['n excluded'].iloc[0])} params)")
    ax.set_xticks(np.arange(len(LABEL))); ax.set_xticklabels(list(LABEL.values()), rotation=25, ha="right", fontsize=7)
    ax.set_title(f2.CONFIGS[cfg], fontsize=9); ax.axhline(0, color="k", lw=0.6)
axes[0].set_ylabel("increase of the null-fit chi2 when the group is removed [%]"); axes[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT / "group_attribution.pdf", bbox_inches="tight")
print(df[df.config == "b"].pivot(index="dataset", columns="group", values="increase").round(1).to_string())
