#!/usr/bin/env python3
"""S11: extra report figures: (a) ICARUS Run4 WireModxThetaXW duplicate artifact, (b) channel / class audit re-plotted from s6_audits.json (truth-less WireMod samples excluded from the signal curve)."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
from ana import CTR
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
d = json.loads((FIG / "icarus_run4_wiremodxthetaxw_dedup.json").read_text()); b = json.loads((FIG / "baseline.json").read_text())["icarus_run4__WireModxThetaXW"]
fig, ax = plt.subplots(figsize=(8.2, 5.2)); ax.plot(CTR, (np.array(d["prod"]) - 1) * 100, "-o", color="#D55E00", label="Production histogram (repeated events counted)"); ax.plot(CTR, (np.array(d["dedup"]) - 1) * 100, "-o", color="#0072B2", label="Exact duplicates removed")
ax.errorbar(CTR, (np.array(b["ratio"]) - 1) * 100, yerr=np.array(b["err"]) * 100, fmt="s", color="k", ms=4, label="Repeated keys dropped (pair table) with bootstrap error", alpha=0.7)
ax.axhline(0, color="0.5", lw=1); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.legend(fontsize=10, loc="upper left"); mark_axis(ax, r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(ax, "ICARUS Run 4 WireModxThetaXW", alignment="right"); fig.tight_layout(); fig.savefig(FIG / "s11_icarus_run4_duplicates.png", dpi=100); plt.close(fig)
res = json.loads((FIG / "s6_audits.json").read_text())
sel = {"sbnd_run1": ["ChargeScale15", "ChargeScale", "MattModYZ", "DENT", "SCE0x", "WireModxThetaXW"], "icarus_run2": ["ChargeScale", "SCE0x", "WireModxThetaXW"], "icarus_run4": ["ChargeScale", "SCE2x", "WireModxThetaXW"]}
for det, labs in sel.items():
    nc = 3; nr = int(np.ceil(len(labs) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(15, 3.9 * nr), squeeze=False)
    for ax, lab in zip(axs.ravel(), labs):
        o = res[f"{det}__{lab}"]; ax.plot(CTR, (np.array(o["combined"]["r"]) - 1) * 100, "k-o", ms=4, lw=1.8, label="Combined (as applied)")
        for k, c, t in (("1p", "#0072B2", "1p channel"), ("Mp", "#D55E00", "Mp channel"), ("signal", "#009E73", "Signal only")):
            if k == "signal" and det == "sbnd_run1" and lab in ("WireModxThetaXW", "WireModYZ"): continue
            ax.plot(CTR, (np.array(o[k]["r"]) - 1) * 100, "-", color=c, lw=1.4, label=t)
        if "filterEE" in o: ax.plot(CTR, (np.array(o["filterEE"]["r"]) - 1) * 100, "--", color="#CC79A7", lw=1.6, label="with filterEE")
        ax.axhline(0, color="0.6", lw=0.8); ax.text(0.04, 0.92, lab, transform=ax.transAxes, fontweight="bold", va="top"); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]")
    for ax in axs.ravel()[len(labs):]: ax.axis("off")
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=11)
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, nc - 1], det.replace("_", " ").upper(), alignment="right"); fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(FIG / f"s11_channel_audit_{det}.png", dpi=95); plt.close(fig)
print("ok")
