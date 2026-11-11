#!/usr/bin/env python3
"""Alternative views of the selection overlap between the nominal and varied samples (selected in both / nominal only / variation only), from docs/syst_deepdive/deepdive.json.
V1 stacked fractions of the union, V2 butterfly (nominal flows left, variation flows right), V3 churn-versus-net scatter, V4 heat map matrix, V5 integrated flow bars.
Run after syst_deepdive.py:  python plotting/dent_study/selection_overlap_views.py"""
import json, pathlib
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
OUT = pathlib.Path(__file__).resolve().parents[2] / "docs" / "syst_deepdive"; R = json.loads((OUT / "deepdive.json").read_text())
LAB = list(R); COL = {"both": "0.72", "nom": "#D55E00", "var": "#0072B2"}; VC = ["#0072B2", "#009E73", "#CC79A7", "#E69F00"]
def rows(label):
    out = []
    for ch in ("1p", "Mp"):
        b = R[label]["bins"][ch]; e = b["edges"]
        for i in range(len(b["total"])):
            lo = -b["lost"][i]; ga = b["gained"][i]; both = 1 - lo; out.append(dict(ch=ch, c=0.5 * (e[i] + e[i + 1]), lo=lo, ga=ga, both=both, uni=1 + ga, net=ga - lo))
    return out
def stamp(axs): mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[-1], "SBND Run 1", alignment="right")
# V1: stacked fractions of the union
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
for ax, l in zip(axs.ravel(), LAB):
    r = rows(l); x = np.arange(len(r)); uni = np.array([q["uni"] for q in r])
    ax.bar(x, [q["both"] / q["uni"] * 100 for q in r], color=COL["both"], label="Selected in both")
    ax.bar(x, [q["lo"] / q["uni"] * 100 for q in r], bottom=[q["both"] / q["uni"] * 100 for q in r], color=COL["nom"], label="Nominal only")
    ax.bar(x, [q["ga"] / q["uni"] * 100 for q in r], bottom=[(q["both"] + q["lo"]) / q["uni"] * 100 for q in r], color=COL["var"], label="Variation only")
    ax.set_ylim(60, 100); ax.set_ylabel("Share of union [%]"); ax.set_xticks(x); ax.set_xticklabels([f"{q['c']:.2f}" for q in r], rotation=60, fontsize=9)
    ax.set_xlabel("Reco visible energy [GeV]  (1p | Mp)"); ax.text(0.03, 0.06, l, transform=ax.transAxes, fontweight="bold", fontsize=13)
    ax.axvline(len(R[l]["bins"]["1p"]["total"]) - 0.5, color="k", lw=0.8, ls=":")
h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=11); stamp(axs.ravel()[[0, 1]]); fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(OUT / "views_v1_stacked_union.png", dpi=110); plt.close(fig)
# V2: butterfly: the nominal sample as 100 % = kept + lost (left), the variation sample = kept + gained (right)
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
for ax, l in zip(axs.ravel(), LAB):
    r = rows(l); y = np.arange(len(r))[::-1]
    ax.barh(y, [-q["both"] * 100 for q in r], color=COL["both"]); ax.barh(y, [q["both"] * 100 for q in r], color=COL["both"])
    ax.barh(y, [-q["lo"] * 100 for q in r], left=[-q["both"] * 100 for q in r], color=COL["nom"]); ax.barh(y, [q["ga"] * 100 for q in r], left=[q["both"] * 100 for q in r], color=COL["var"])
    ax.axvline(0, color="k", lw=0.8); ax.set_yticks(y); ax.set_yticklabels([f"{q['ch']} {q['c']:.2f}" for q in r], fontsize=8); ax.set_xlim(-125, 125)
    ax.set_xlabel("Nominal | varied [% of nominal]"); ax.text(0.03, 0.96, l, transform=ax.transAxes, fontweight="bold", fontsize=13, va="top")
    ax.set_xticks([-100, -50, 0, 50, 100]); ax.set_xticklabels(["100", "50", "0", "50", "100"])
stamp(axs.ravel()[[0, 1]]); fig.tight_layout(); fig.savefig(OUT / "views_v2_butterfly.png", dpi=110); plt.close(fig)
# V3: churn versus net
fig, axs = plt.subplots(1, 2, figsize=(13, 6.3))
ax = axs[0]
for l, c in zip(LAB, VC):
    r = rows(l); ax.scatter([q["lo"] * 100 for q in r], [q["ga"] * 100 for q in r], s=[40 if q["ch"] == "1p" else 40 for q in r], c=c, marker="o", edgecolor="k", label=l, zorder=3)
m = 24; ax.plot([0, m], [0, m], "k-", lw=1.2); ax.fill_between([0, m], [0, m], [m, m], color="0.93"); ax.text(2, 22, "variation gains net", fontsize=11); ax.text(14, 3, "variation loses net", fontsize=11)
ax.set_xlim(0, m); ax.set_ylim(0, m); ax.set_xlabel("Nominal only [% of nominal]"); ax.set_ylabel("Variation only [% of nominal]"); ax.legend(fontsize=10, loc="upper left", bbox_to_anchor=(0.0, 0.88)); ax.set_aspect("equal")
ax = axs[1]
for l, c in zip(LAB, VC):
    r = rows(l); ax.scatter([(q["lo"] + q["ga"]) * 100 for q in r], [q["net"] * 100 for q in r], c=c, edgecolor="k", s=40, label=l, zorder=3)
ax.axhline(0, color="k", lw=1); ax.set_xlabel("Churn (nominal only + variation only) [%]"); ax.set_ylabel("Net (var. $-$ nom.) [%]")
stamp(axs); fig.tight_layout(); fig.savefig(OUT / "views_v3_churn_vs_net.png", dpi=110); plt.close(fig)
# V4: heat maps
fig, axs = plt.subplots(3, 1, figsize=(12, 8.5), sharex=True)
rr = {l: rows(l) for l in LAB}; names = [f"{q['ch']} {q['c']:.2f}" for q in rr[LAB[0]]]
for ax, (key, title, cm, vm) in zip(axs, (("lo", "Nominal only [%]", "Oranges", None), ("ga", "Variation only [%]", "Blues", None), ("net", "Net (variation $-$ nominal) [%]", "RdBu_r", 8))):
    M = np.array([[q[key] * 100 for q in rr[l]] for l in LAB]); v = vm or M.max(); im = ax.imshow(M, aspect="auto", cmap=cm, vmin=-v if vm else 0, vmax=v)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]): ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=8, color="w" if abs(M[i, j]) > 0.6 * v else "k")
    ax.set_yticks(range(len(LAB))); ax.set_yticklabels(LAB, fontsize=9); ax.text(1.005, 0.5, title, transform=ax.transAxes, rotation=90, va="center", fontsize=10)
axs[-1].set_xticks(range(len(names))); axs[-1].set_xticklabels(names, rotation=60, fontsize=8); stamp(axs[:1].tolist() + [axs[0]]); fig.tight_layout(); fig.savefig(OUT / "views_v4_heatmaps.png", dpi=110); plt.close(fig)
# V5: integrated flow bars (whole selection)
fig, ax = plt.subplots(figsize=(11, 5.2))
for i, l in enumerate(LAB[::-1]):
    q = R[l]; lo, ga = q["lost_frac"] * 100, q["gained_frac"] * 100
    ax.barh(i + 0.18, 100 - lo, height=0.32, color=COL["both"]); ax.barh(i + 0.18, lo, left=100 - lo, height=0.32, color=COL["nom"])
    ax.barh(i - 0.18, 100 - lo, height=0.32, color=COL["both"]); ax.barh(i - 0.18, ga, left=100 - lo, height=0.32, color=COL["var"])
    ax.text(101, i + 0.18, f"nominal: {lo:.1f}% lost", va="center", fontsize=10); ax.text(101, i - 0.18, f"variation: {ga:.1f}% gained", va="center", fontsize=10)
ax.set_yticks(range(len(LAB))); ax.set_yticklabels(LAB[::-1]); ax.set_xlim(70, 135); ax.set_xlabel("Selected events, [% of nominal selected]"); ax.set_xticks([70, 80, 90, 100])
ax.plot([], [], "s", color=COL["both"], label="Selected in both"); ax.plot([], [], "s", color=COL["nom"], label="Nominal only"); ax.plot([], [], "s", color=COL["var"], label="Variation only"); ax.legend(fontsize=10, loc="upper right", bbox_to_anchor=(1.0, 1.0), frameon=True)
stamp([ax, ax]); fig.tight_layout(); fig.savefig(OUT / "views_v5_integrated_flows.png", dpi=110); plt.close(fig)
print("ok")
