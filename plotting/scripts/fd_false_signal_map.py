#!/usr/bin/env python3
"""False-signal map: where in (sin^2 2theta_mumu, dm^2_41) the oscillation fit lands for each non-oscillation fake dataset.

Background: the Asimov sensitivity surface of the full-exposure 2D analysis (runs/legacy_root/contours_full_v1_surf.root, made on
2026-09-14 with the same analysis XML family; a LEGACY production, not re-run on the current inputs) drawn at the Wilks 2-dof 90% and 95%
levels (not validated). Markers: best-fit oscillation point of the oscillation-floating fit, for fake datasets with Delta chi2 > 0.3 in the
combined analysis, sized by Delta chi2; stars: the three IceCube 2024 injection points. Left of a marker's size scale, small markers are
statistically insignificant. Writes plotting/fakedata_2d_figures/false_signal_map.pdf.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_false_signal_map.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).parent))
import fd2d as f2
from fd_figures2 import all_stems

OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
h = uproot.open("/Users/mueller/Projects/spineosc/runs/legacy_root/contours_full_v1_surf.root")["surf"]
v = h.values()
x = np.sqrt(h.axis(0).edges()[1:] * h.axis(0).edges()[:-1])
y = np.sqrt(h.axis(1).edges()[1:] * h.axis(1).edges()[:-1])
IC = {"dm2lo": (0.116, 2.07), "bestfit": (0.16, 3.5), "dm2hi": (0.187, 12.1)}
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
for ax, cfg in zip(axes, ("b", "1p")):
    ax.contour(x, y, v.T, levels=[4.61, 5.99], colors=["0.35", "0.55"], linestyles=["-", "--"])
    ax.contourf(x, y, v.T, levels=[5.99, 1e9], colors=["0.92"])
    n = 0
    for s in all_stems():
        p = f2.pair("ds", cfg, s)
        if not p:
            continue
        d = p[0]["chi2"] - p[1]["chi2"]
        if d < 0.3 or "osc_ic2024" in s:
            continue
        ax.scatter(p[1]["s22"], p[1]["dm2"], s=25 + 55 * d, alpha=0.6, edgecolor="k", lw=0.5, color="C3" if d > 4 else "C1")
        ax.annotate(s.replace("_", " ")[:22], (p[1]["s22"], p[1]["dm2"]), fontsize=5, xytext=(3, 3), textcoords="offset points"); n += 1
    for k, (s22, dm2) in IC.items():
        ax.plot(s22, dm2, "*", ms=13, color="C0", mec="k"); ax.annotate("IceCube " + k, (s22, dm2), fontsize=6, xytext=(5, -9), textcoords="offset points", color="C0")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.01, 1); ax.set_ylim(0.1, 100)
    ax.set_xlabel(r"$\sin^2 2\theta_{\mu\mu}$"); ax.set_title(f"{f2.CONFIGS[cfg]}: fits with $\\Delta\\chi^2$ > 0.3 ({n} datasets)", fontsize=9)
axes[0].set_ylabel(r"$\Delta m^2_{41}$ [eV$^2$]")
fig.suptitle("Best-fit oscillation point of non-oscillation fake data; contour: legacy full-exposure Asimov sensitivity (Wilks 90%, 95%)", fontsize=9)
fig.savefig(OUT / "false_signal_map.pdf", bbox_inches="tight")
