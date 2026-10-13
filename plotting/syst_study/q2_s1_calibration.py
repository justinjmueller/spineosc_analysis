#!/usr/bin/env python3
"""Q2 deliverable 1: calibration of the per-bin ratio uncertainty.
 (a) paired Poisson bootstrap (2000 hashed replicates) vs analytic delta-method variance (cross-classified table) vs the naive independent-samples error
 (b) split-half test: random event splits (hash) and random block splits; E[(rA-rB)^2]/E[sA^2+sB^2], chi2/ndf of (rA-rB)
 (c) block bootstrap / delete-one-block jackknife over (Run) [SBND] or (Run,Subrun) [ICARUS]; inflation factor vs the event-level error."""
from q2_lib import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
NSPLIT = 400
def groups(d, name, coarse=False):
    g = d["run"].astype(np.int64) * 100000 + (0 if det_of(name).startswith("sbnd") else d["sub"])
    if coarse: g = d["run"].astype(np.int64) // 25
    return np.unique(g, return_inverse=True)[1]
def half_stats(d, w):
    rA, CA, NA, VA = cov_analytic(d, w=w); rB, CB, NB_, VB = cov_analytic(d, w=1 - w); return rA, rB, np.diag(CA), np.diag(CB)
def split_test(d, mode, gid=None, nsplit=NSPLIT):
    dd, ss, zz, c2, vv = [], [], [], [], []
    key = d["key"] if gid is None else None; ug = None
    for s in range(nsplit):
        if gid is None: w = (uni(d["key"], 1000 + s) < 0.5).astype(float)
        else:
            ng = gid.max() + 1; ugrp = uni(np.arange(ng, dtype=np.uint64) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(7), 2000 + s); w = (ugrp[gid] < 0.5).astype(float)
        rA, rB, vA, vB = half_stats(d, w); dd.append(rA - rB); ss.append(vA + vB)
    dd = np.array(dd); ss = np.array(ss); z = dd / np.sqrt(ss)
    return dict(dsq=(dd ** 2).mean(0), pred=ss.mean(0), vr_bin=((dd ** 2).mean(0) / ss.mean(0)), chi2ndf=(z ** 2).mean(1), z=z)
def block_methods(d, gid, nb=1000):
    G = gid.max() + 1; bn, bv = d["bn"], d["bv"]; Ng = np.zeros((G, NB)); Vg = np.zeros((G, NB))
    for b in range(NB): Ng[:, b] = np.bincount(gid[bn == b], minlength=G); Vg[:, b] = np.bincount(gid[bv == b], minlength=G)
    N = Ng.sum(0); V = Vg.sum(0); r = V / N
    rj = (V - Vg) / (N - Ng); rbar = rj.mean(0); sj = np.sqrt((G - 1) / G * ((rj - rbar) ** 2).sum(0))   # delete-one-block jackknife
    gk = np.arange(G, dtype=np.uint64) * np.uint64(0xD6E8FEB86659FD93) + np.uint64(99); W = np.stack([pois_w(gk, k, 5) for k in range(nb)], 1)
    R = (W.T @ Vg) / (W.T @ Ng); return sj, R.std(0), G
res = {}
for n in NAMES:
    d = load(n); r, C, N, V = cov_analytic(d); san = np.sqrt(np.diag(C)); bz = np.load(RES / f"boot_{n}.npz"); R = bz["V"] / bz["N"]; sbs = R.std(0); sna = naive_sigma(N, V, r)
    out = dict(name=n, r=r.tolist(), N=N.tolist(), V=V.tolist(), sig_an=san.tolist(), sig_bs=sbs.tolist(), sig_naive=sna.tolist(), bs_bias_sigma=((R.mean(0) - r) / sbs).tolist(),
               bs_over_an_median=float(np.median(sbs / san)), bs_over_an_min=float((sbs / san).min()), bs_over_an_max=float((sbs / san).max()), naive_over_paired_median=float(np.median(sna / san)))
    Cbs = np.cov(R.T); out["corr_adjacent_an"] = float(np.mean([C[i, i + 1] / (san[i] * san[i + 1]) for i in range(NB - 1)])); out["corr_adjacent_bs"] = float(np.mean([Cbs[i, i + 1] / (sbs[i] * sbs[i + 1]) for i in range(NB - 1)]))
    ev = split_test(d, "event"); out["split_event"] = dict(vr_bin=ev["vr_bin"].tolist(), vr_all=float(ev["dsq"].sum() / ev["pred"].sum()), chi2ndf_mean=float(ev["chi2ndf"].mean()), chi2ndf_sd=float(ev["chi2ndf"].std()), z_mean=float(ev["z"].mean()), z_sd=float(ev["z"].std()))
    gid = groups(d, n); bm = split_test(d, "block", gid); out["split_block"] = dict(vr_bin=bm["vr_bin"].tolist(), vr_all=float(bm["dsq"].sum() / bm["pred"].sum()), chi2ndf_mean=float(bm["chi2ndf"].mean()), chi2ndf_sd=float(bm["chi2ndf"].std()), G=int(gid.max() + 1))
    sj, sbb, G = block_methods(d, gid); out["block"] = dict(G=int(G), jk_over_ev=(sj / san).tolist(), bb_over_ev=(sbb / san).tolist(), jk_med=float(np.median(sj / san)), bb_med=float(np.median(sbb / san)))
    if det_of(n).startswith("sbnd"):
        g2 = groups(d, n, True); sj2, sb2, G2 = block_methods(d, g2); out["block25"] = dict(G=int(G2), jk_med=float(np.median(sj2 / san)), bb_med=float(np.median(sb2 / san)), jk_over_ev=(sj2 / san).tolist())
    res[n] = out
    print(f"{n:30s} bs/an {out['bs_over_an_median']:.3f} [{out['bs_over_an_min']:.2f},{out['bs_over_an_max']:.2f}] naive/paired {out['naive_over_paired_median']:.1f} | split-ev VR {out['split_event']['vr_all']:.3f} chi2/ndf {out['split_event']['chi2ndf_mean']:.2f}+-{out['split_event']['chi2ndf_sd']:.2f} | "
          f"split-blk VR {out['split_block']['vr_all']:.3f} chi2/ndf {out['split_block']['chi2ndf_mean']:.2f} | blk JK/ev {out['block']['jk_med']:.3f} BB/ev {out['block']['bb_med']:.3f} G={G}" + (f" | Run//25 JK {out['block25']['jk_med']:.3f} G={out['block25']['G']}" if 'block25' in out else ''), flush=True)
json.dump(res, open(Q2 / "q2_calibration.json", "w"))
