#!/usr/bin/env python3
"""S8 (Q1/Q4): which flipped population carries the energy dependence?  Lost events (nominal only) and gained events (varied only) are classified by which leading candidate sits near its decision boundary:
muon candidate margin m_mu = 2*s_mu - max(s_pi, s_p) (selected muons have m_mu > 0 by construction, verified 100%), leading-proton margin m_p = s_p - max(2*s_mu, s_pi) (>0 for 100% of selected events).
Classes: 'proton-marginal' (m_p < T and smaller than m_mu), 'muon-marginal' (m_mu < T and smaller than m_p), 'neither'. Lost events use nominal-side features, gained events varied-side features.
Contribution to ratio-1 per coarse reconstructed-energy bin = (gained_k - lost_k)/N_nominal. T = 0.5 (scores range 0-1)."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ana import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
E6 = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0]); T = 0.5
def margins(df, s):
    g = lambda c: df[f"{s}_reco_leading_{c}"].to_numpy(float)
    return 2 * g("muon_muon_softmax") - np.maximum(g("muon_pion_softmax"), g("muon_proton_softmax")), g("proton_proton_softmax") - np.maximum(2 * g("proton_muon_softmax"), g("proton_pion_softmax"))
def klass(mm, mp):
    return np.where((mp < T) & (mp <= mm), 0, np.where(mm < T, 1, 2))      # 0 proton-marginal, 1 muon-marginal, 2 neither
res = {}
for name in pairs_available():
    df = load(name); inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); mmn, mpn = margins(df, "n"); mmv, mpv = margins(df, "v"); kn = klass(mmn, mpn); kv = klass(mmv, mpv)
    En = df.n_reco_visible_energy.to_numpy(float); Ev = df.v_reco_visible_energy.to_numpy(float)
    bn = np.digitize(En, E6) - 1; bv = np.digitize(Ev, E6) - 1; nb = len(E6) - 1; N = np.array([(inn & (bn == i)).sum() for i in range(nb)], float)
    out = dict(edges=E6.tolist(), N=N.tolist(), lost={}, gained={}, ratio=[])
    for c, lab in enumerate(("proton-marginal", "muon-marginal", "neither")):
        out["lost"][lab] = [float(((inn & ~inv) & (bn == i) & (kn == c)).sum() / N[i]) for i in range(nb)]
        out["gained"][lab] = [float(((~inn & inv) & (bv == i) & (kv == c)).sum() / N[i]) for i in range(nb)]
    out["both_in_class_frac"] = {lab: float(((inn & inv) & (kn == c)).sum() / (inn & inv).sum()) for c, lab in enumerate(("proton-marginal", "muon-marginal", "neither"))}
    out["net"] = [float(sum(out["gained"][l][i] - out["lost"][l][i] for l in out["lost"])) for i in range(nb)]
    res[name] = out
json.dump(res, open(FIG / "s8_populations.json", "w"))
sel = [("sbnd_run1__ChargeScale15", "SBND ChargeScale15"), ("sbnd_run1__ChargeScale", "SBND ChargeScale"), ("sbnd_run1__GainHigh", "SBND Gain +5%"), ("sbnd_run1__GainLow", "SBND Gain $-$5%"), ("icarus_run2__ChargeScale", "ICARUS Run2 ChargeScale"), ("icarus_run4__ChargeScale", "ICARUS Run4 ChargeScale")]
fig, axs = plt.subplots(2, 3, figsize=(16, 8.5)); cols = {"proton-marginal": "#D55E00", "muon-marginal": "#0072B2", "neither": "0.65"}
for ax, (n, lab) in zip(axs.ravel(), sel):
    o = res[n]; x = np.arange(len(E6) - 1); wd = 0.38
    for k, (grp, sgn) in enumerate((("lost", -1), ("gained", +1))):
        bot = np.zeros(len(x))
        for l in ("proton-marginal", "muon-marginal", "neither"):
            v = sgn * np.array(o[grp][l]) * 100; ax.bar(x + (k - 0.5) * wd, v, wd, bottom=bot, color=cols[l], edgecolor="k", lw=0.5, label=l if grp == "lost" else None, hatch="//" if grp == "gained" else None); bot = bot + v
    ax.plot(x, np.array(o["net"]) * 100, "k_", ms=22, mew=3, label="Net ratio $-$ 1"); ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"{a:.1f}-{b:.1f}" for a, b in zip(E6[:-1], E6[1:])], fontsize=9); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Fraction of nominal [%]"); ax.text(0.03, 0.95, lab, transform=ax.transAxes, fontweight="bold", va="top")
h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h + [plt.Rectangle((0, 0), 1, 1, fc="w", ec="k", hatch="//")], l + ["(hatched = gained, solid = lost)"], loc="lower center", ncol=5, frameon=False, fontsize=11)
mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 2], "Detector variation churn", alignment="right"); fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(FIG / "s8_populations.png", dpi=95)
def draw(selx, fname):
    nr = int(np.ceil(len(selx) / 3)); fig, axs = plt.subplots(nr, 3, figsize=(16, 4.4 * nr), squeeze=False)
    for ax, (n, lab) in zip(axs.ravel(), selx):
        o = res[n]; x = np.arange(len(E6) - 1); wd = 0.38
        for k, (grp, sgn) in enumerate((("lost", -1), ("gained", +1))):
            bot = np.zeros(len(x))
            for l in ("proton-marginal", "muon-marginal", "neither"):
                v = sgn * np.array(o[grp][l]) * 100; ax.bar(x + (k - 0.5) * wd, v, wd, bottom=bot, color=cols[l], edgecolor="k", lw=0.5, label=l if grp == "lost" else None, hatch="//" if grp == "gained" else None); bot = bot + v
        ax.plot(x, np.array(o["net"]) * 100, "k_", ms=22, mew=3, label="Net ratio $-$ 1"); ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels([f"{a:.1f}-{b:.1f}" for a, b in zip(E6[:-1], E6[1:])], fontsize=9); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Fraction of nominal [%]"); ax.text(0.03, 0.95, lab, transform=ax.transAxes, fontweight="bold", va="top")
    for ax in axs.ravel()[len(selx):]: ax.axis("off")
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h + [plt.Rectangle((0, 0), 1, 1, fc="w", ec="k", hatch="//")], l + ["(hatched = gained, solid = lost)"], loc="lower center", ncol=5, frameon=False, fontsize=11)
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 2], "Detector variation churn", alignment="right"); fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(FIG / fname, dpi=90); plt.close(fig)
draw([("sbnd_run1__MattModYZ", "SBND MattModYZ"), ("sbnd_run1__WireModYZ", "SBND WireModYZ"), ("sbnd_run1__WireModxThetaXW", "SBND WireModxThetaXW"), ("sbnd_run1__SCE0x", "SBND SCE0x"), ("sbnd_run1__SCE2x", "SBND SCE2x"), ("sbnd_run1__DENT", "SBND DENT")], "s8_populations_sbnd_other.png")
draw([("icarus_run2__SCE0x", "ICARUS Run2 SCE0x"), ("icarus_run2__SCE2x", "ICARUS Run2 SCE2x"), ("icarus_run2__WireModxThetaXW", "ICARUS Run2 WireModxThetaXW"), ("icarus_run2__WireModYZ", "ICARUS Run2 WireModYZ"), ("icarus_run4__SCE0x", "ICARUS Run4 SCE0x"), ("icarus_run4__WireModxThetaXW", "ICARUS Run4 WireModxThetaXW")], "s8_populations_icarus_other.png")
for n, lab in sel:
    o = res[n]; print(lab, "net%:", np.round(np.array(o["net"]) * 100, 2).tolist(), "| lost proton-marg %:", np.round(np.array(o["lost"]["proton-marginal"]) * 100, 2).tolist(), "| gained muon-marg %:", np.round(np.array(o["gained"]["muon-marginal"]) * 100, 2).tolist(), "| both-class", {k: round(v, 3) for k, v in o["both_in_class_frac"].items()})
