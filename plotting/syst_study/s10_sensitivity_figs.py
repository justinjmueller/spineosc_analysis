#!/usr/bin/env python3
"""S10: figures for the sensitivity discussion: (a) PROfit 'no X' scans (chi2 relative to the 'full' run in the 0.5<chi2<20 band) - note the common offset; (b) toy Asimov results: limit improvement from removing the SBND or ICARUS charge-scale nuisance versus the uncorrelated per-detector normalisation uncertainty; (c) event-weighted effect size of each active systematic."""
import json, numpy as np, uproot, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); P = REPO / "runs" / "legacy_root"
def L(n): return uproot.open(P / f"{n}_surf.root")["surf"].values()
full = L("contours_full_v1"); band = (full > 0.5) & (full < 20)
names = {"no SBND ChargeScale": "contours_full_nochargescale", "no ICARUS ChargeScale": "contours_full_noicchargescale", "no SBND DENT": "contours_full_nodent", "no SBND WireMod-YZ(*)": "contours_full_nowmyz", "no SBND WireMod x-theta": "contours_full_nowmxtxw",
         "no ICARUS WireMod YZ": "contours_full_noicwmyz", "no ICARUS WireMod x-theta": "contours_full_noicwmxtxw", "no ICARUS SCE": "contours_full_noicsce", "no SBND Trigger": "contours_full_notrigemu", "no ICARUS Trigger": "contours_full_noictrigemu"}
fig, axs = plt.subplots(1, 3, figsize=(19, 5.8))
ax = axs[0]; med = []; lo = []; hi = []
for k, n in names.items(): r = L(n)[band] / full[band]; med.append(np.median(r)); lo.append(np.percentile(r, 16)); hi.append(np.percentile(r, 84))
y = np.arange(len(names)); ax.errorbar(med, y, xerr=[np.array(med) - np.array(lo), np.array(hi) - np.array(med)], fmt="o", color="k", ms=6); ax.set_yticks(y); ax.set_yticklabels(list(names), fontsize=10); ax.invert_yaxis(); ax.axvline(1, color="0.5", lw=1)
ax.axvline(np.median(med), color="#D55E00", ls="--", lw=1.4); ax.set_xlabel(r"median $\chi^2$(no X) / $\chi^2$(full)"); ax.text(0.03, 0.04, "dashed: median over all variants", transform=ax.transAxes, color="#D55E00", fontsize=10)
# (b) toy
ax = axs[1]; res = {}
for tag, val in (("_norm0.003", 0.003), ("_norm0.01", 0.01), ("", 0.03), ("_norm0.10", 0.10)):
    d = json.loads((FIG / f"s7_toy_sensitivity{tag}.json").read_text()); base = np.array(d["all"]["limit"]); dm = np.array(d["dm2"]); sel = (dm > 0.5) & (dm < 20)
    res[val] = {k: 100 * np.nanmean((base[sel] - np.array(d[k]["limit"])[sel]) / base[sel]) for k in ("noSBNDCS", "noICARUSCS")}
xs = sorted(res); ax.plot([x * 100 for x in xs], [res[x]["noSBNDCS"] for x in xs], "-o", color="#D55E00", label="remove SBND charge scale"); ax.plot([x * 100 for x in xs], [res[x]["noICARUSCS"] for x in xs], "-o", color="#0072B2", label="remove ICARUS charge scale")
ax.set_xscale("log"); ax.set_xlabel("Per-detector normalisation uncertainty in toy [%]"); ax.set_ylabel("Limit improvement [%]"); ax.legend(fontsize=11, loc="upper right"); ax.text(0.04, 0.05, "Asimov toy (illustrative)", transform=ax.transAxes, fontweight="bold")
# (c) event-weighted effect size
ax = axs[2]; b = json.loads((FIG / "baseline.json").read_text()); items = []
for k, v in b.items():
    if k.split("__")[1] in ("GainHigh", "GainLow", "WireModYZ") and k.startswith("sbnd"): continue
    n = np.array(v["n_nom"]); r = np.array(v["ratio"]); items.append((k.replace("__", " ").replace("sbnd_run1", "SBND").replace("icarus_run", "ICARUS Run"), 100 * np.sqrt(np.sum(n * (r - 1) ** 2) / n.sum()), 100 * np.sqrt(np.mean((r - 1) ** 2)), 100 * v["integrated"]))
items.sort(key=lambda t: -t[1]); y = np.arange(len(items)); ax.barh(y - 0.2, [t[1] for t in items], 0.4, color="#009E73", label="event-weighted rms of ratio $-$ 1"); ax.barh(y + 0.2, [t[2] for t in items], 0.4, color="0.7", label="unweighted rms over bins")
ax.set_yticks(y); ax.set_yticklabels([t[0] for t in items], fontsize=9); ax.invert_yaxis(); ax.set_xlabel("[%]"); ax.legend(fontsize=10, loc="lower right")
mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[2], "Active detector systematics", alignment="right"); fig.tight_layout(w_pad=3.0); fig.savefig(FIG / "s10_sensitivity.png", dpi=95)
json.dump(dict(toy={str(k): v for k, v in res.items()}, profit_median_ratio=dict(zip(names, med))), open(FIG / "s10_sensitivity_summary.json", "w"))
print({k: (round(m, 3)) for k, m in zip(names, med)}); print(res)
