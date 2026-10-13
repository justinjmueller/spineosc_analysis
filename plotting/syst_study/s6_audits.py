#!/usr/bin/env python3
"""S6: audits of how the production ratio histogram is applied.
(a) ICARUS: the calculator selection omits reco_user_containment_cut_filterEE that the analysis XMLs apply -> ratio on the analysis selection.
(b) The single ratio histogram is applied to both proton-multiplicity channels (1p, Mp): channel-specific ratios (channel from each sample's own reco_proton_multiplicity).
(c) The ratio is built from signal+background mixture and applied to all subchannels: ratio for true target (signal) vs the rest."""
import json, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ana import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
def ratio_boot(bn, bv, nb=200, seed=2):
    n0, v0 = counts(bn, bv); r0 = v0 / n0; rng = np.random.default_rng(seed); R = np.empty((nb, NBIN))
    for i in range(nb):
        w = rng.poisson(1.0, len(bn)).astype(float); n, v = counts(bn, bv, w); R[i] = v / n
    return r0, R.std(0), n0
res = {}
for name in pairs_available():
    df = load(name); det = name.split("__")[0]; bn, bv = prep(df); out = {}
    r_all, e_all, n_all = ratio_boot(bn, bv)
    out["combined"] = dict(r=r_all.tolist(), e=e_all.tolist(), n=n_all.tolist())
    if det.startswith("icarus"):
        fn = df["n_reco_user_containment_cut_filterEE"].to_numpy(float) == 1; fv = df["v_reco_user_containment_cut_filterEE"].to_numpy(float) == 1
        r, e, n = ratio_boot(np.where(fn, bn, -1), np.where(fv, bv, -1)); out["filterEE"] = dict(r=r.tolist(), e=e.tolist(), n=n.tolist(), frac_pass_n=float((fn & df.in_n.to_numpy()).sum() / df.in_n.sum()))
    mn = df["n_reco_proton_multiplicity"].to_numpy(float); mv = df["v_reco_proton_multiplicity"].to_numpy(float)
    for ch, cn, cv in (("1p", mn == 1, mv == 1), ("Mp", mn > 1, mv > 1)):
        r, e, n = ratio_boot(np.where(cn, bn, -1), np.where(cv, bv, -1)); out[ch] = dict(r=r.tolist(), e=e.tolist(), n=n.tolist())
    def sig(side):
        cc = df[f"{side}_true_cc"].to_numpy(float) > 0; pdg = np.abs(df[f"{side}_true_pdg"].to_numpy(float)) == 14; cat = df[f"{side}_true_category"].to_numpy(float); return cc & pdg & np.isin(cat, [0, 2])
    sn = sig("n"); sv = sig("v"); sn = np.where(df.in_n.to_numpy(), sn, False); sv = np.where(df.in_v.to_numpy(), sv, False)
    for cl, cn, cv in (("signal", sn, sv), ("background", ~sn & df.in_n.to_numpy(), ~sv & df.in_v.to_numpy())):
        r, e, n = ratio_boot(np.where(cn, bn, -1), np.where(cv, bv, -1)); out[cl] = dict(r=r.tolist(), e=e.tolist(), n=n.tolist())
    # summary statistics: chi2 of (r_sub - r_all) with the quadrature of bootstrap errors (they are correlated, so the chi2 is an overestimate of significance; report also the mean shift)
    summ = {}
    for k in ("filterEE", "1p", "Mp", "signal", "background"):
        if k not in out: continue
        d = np.array(out[k]["r"]) - r_all; e = np.hypot(np.array(out[k]["e"]), e_all)
        summ[k] = dict(max_abs_dev=float(np.nanmax(np.abs(d))), mean_dev_weighted=float(np.nansum(d * n_all) / n_all.sum()), rms_dev=float(np.sqrt(np.nanmean(d ** 2))), max_dev_over_sigma=float(np.nanmax(np.abs(d) / e)))
    out["summary"] = summ; res[name] = out
    print(f"{name:30s} rms|ratio_comb-1|={100*np.sqrt(np.mean((r_all-1)**2)):.2f}%  " + "  ".join(f"{k}: rms dev {100*v['rms_dev']:.2f}% max {100*v['max_abs_dev']:.2f}%" for k, v in summ.items()), flush=True)
json.dump(res, open(FIG / "s6_audits.json", "w"))
# figure: for each detector, ratio curves of combined, 1p, Mp, signal for the main variations
sel = {"sbnd_run1": ["ChargeScale15", "MattModYZ", "WireModxThetaXW", "DENT", "SCE0x", "WireModYZ"], "icarus_run2": ["ChargeScale", "SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ"], "icarus_run4": ["ChargeScale", "SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ"]}
for det, labs in sel.items():
    nc = 3; nr = int(np.ceil(len(labs) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(15, 3.8 * nr), squeeze=False)
    for ax, lab in zip(axs.ravel(), labs):
        o = res[f"{det}__{lab}"]; ax.plot(CTR, (np.array(o["combined"]["r"]) - 1) * 100, "k-o", ms=4, lw=1.8, label="Combined (applied)")
        for k, c in (("1p", "#0072B2"), ("Mp", "#D55E00"), ("signal", "#009E73")): ax.plot(CTR, (np.array(o[k]["r"]) - 1) * 100, "-", color=c, lw=1.4, label={"1p": "1p channel", "Mp": "Mp channel", "signal": "Signal only"}[k])
        if "filterEE" in o: ax.plot(CTR, (np.array(o["filterEE"]["r"]) - 1) * 100, "--", color="#CC79A7", lw=1.6, label="with filterEE")
        ax.axhline(0, color="0.6", lw=0.8); ax.text(0.04, 0.92, lab, transform=ax.transAxes, fontweight="bold", va="top"); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]")
    for ax in axs.ravel()[len(labs):]: ax.axis("off")
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=5, frameon=False, fontsize=11)
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, nc - 1], det.replace("_", " ").upper(), alignment="right"); fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(FIG / f"s6_channel_audit_{det}.png", dpi=95); plt.close(fig)
