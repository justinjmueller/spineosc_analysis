#!/usr/bin/env python3
"""S1: baseline ratios with paired-bootstrap errors, comparison with the production histograms, selection flows. Writes docs/syst_study/baseline.json and figures."""
import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ana import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); FIG.mkdir(parents=True, exist_ok=True)
summ = json.loads((OUT / "pairs_summary.json").read_text()); res = {}
for name in pairs_available():
    s = summ[name]; df = load(name); r0, err, R, n0, v0 = boot_ratio(df, nb=300); fl = flows(df)
    prod = prod_hist(s["hist"]); rowratio = np.array(s["hist_v"]) / np.array(s["hist_n"])
    res[name] = dict(det=s["det"], label=s["label"], ratio=r0.tolist(), err=err.tolist(), prod=prod.tolist(), rowratio=rowratio.tolist(), n_nom=n0.tolist(), n_var=v0.tolist(),
                     lost=fl["lost"].tolist(), gained=fl["gained"].tolist(), migration=fl["migration"].tolist(), integrated=float(v0.sum() / n0.sum() - 1),
                     maxdiff_prod_vs_rows=float(np.abs(prod - rowratio).max()), maxdiff_prod_vs_unique=float(np.abs(prod - r0).max()),
                     lost_frac=float((df.in_n & ~df.in_v).sum() / df.in_n.sum()), gained_frac=float((~df.in_n & df.in_v).sum() / df.in_n.sum()), dup_frac=float(s["dup_n"] / s["n_n_rows"]))
    print(f"{name:34s} int {100*res[name]['integrated']:+6.2f}%  lost {100*res[name]['lost_frac']:5.2f}% gained {100*res[name]['gained_frac']:5.2f}%  |prod-rows| {res[name]['maxdiff_prod_vs_rows']:.1e} |prod-unique| {res[name]['maxdiff_prod_vs_unique']:.1e}  err(median) {100*np.median(err):.2f}%")
json.dump(res, open(FIG.parent / "syst_study" / "baseline.json", "w"))
for det in ("sbnd_run1", "icarus_run2", "icarus_run4"):
    names = [n for n in res if res[n]["det"] == det]
    if not names: continue
    nc = 3; nr = int(np.ceil(len(names) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(15, 3.6 * nr), squeeze=False)
    for ax, n in zip(axs.ravel(), names):
        r = res[n]; ax.errorbar(CTR, (np.array(r["ratio"]) - 1) * 100, yerr=np.array(r["err"]) * 100, fmt="o", color="k", ms=4)
        ax.plot(CTR, (np.array(r["prod"]) - 1) * 100, "-", color="#D55E00", lw=1.5); ax.axhline(0, color="0.5", lw=0.8)
        ax.text(0.04, 0.9, r["label"], transform=ax.transAxes, fontweight="bold", va="top"); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]")
    for ax in axs.ravel()[len(names):]: ax.axis("off")
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, nc - 1], det.replace("_", " ").upper(), alignment="right"); fig.tight_layout(); fig.savefig(FIG / f"s1_ratios_{det}.png", dpi=100); plt.close(fig)
