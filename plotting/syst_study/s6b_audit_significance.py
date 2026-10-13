#!/usr/bin/env python3
"""S6b: significance of the audit deviations with a PAIRED bootstrap (same Poisson weights for the subset ratio and the combined ratio), full bootstrap covariance of the difference."""
import json, numpy as np
from ana import *
def rboot(bn, bv, W):
    out = np.empty((len(W), NBIN))
    for i, w in enumerate(W): n, v = counts(bn, bv, w); out[i] = v / n
    return out
res = {}; NB = 300
for name in pairs_available():
    df = load(name); det = name.split("__")[0]; bn, bv = prep(df); rng = np.random.default_rng(5); W = rng.poisson(1.0, (NB, len(bn))).astype(float)
    Ra = rboot(bn, bv, W); ra = counts(bn, bv)[1] / counts(bn, bv)[0]; subs = {}
    mn = df["n_reco_proton_multiplicity"].to_numpy(float); mv = df["v_reco_proton_multiplicity"].to_numpy(float)
    subs["1p"] = (mn == 1, mv == 1); subs["Mp"] = (mn > 1, mv > 1)
    if det.startswith("icarus"): subs["filterEE"] = (df["n_reco_user_containment_cut_filterEE"].to_numpy(float) == 1, df["v_reco_user_containment_cut_filterEE"].to_numpy(float) == 1)
    o = {}
    for k, (cn, cv) in subs.items():
        b1, b2 = np.where(cn, bn, -1), np.where(cv, bv, -1); n0, v0 = counts(b1, b2); r0 = v0 / n0; Rs = rboot(b1, b2, W); d0 = r0 - ra; D = Rs - Ra
        ok = np.isfinite(d0) & np.isfinite(D).all(0); C = np.cov(D[:, ok].T) + 1e-12 * np.eye(ok.sum()); chi2 = float(d0[ok] @ np.linalg.solve(C, d0[ok])); sd = D.std(0)
        o[k] = dict(d=d0.tolist(), sd=sd.tolist(), chi2=chi2, ndf=int(ok.sum()), max_sig=float(np.nanmax(np.abs(d0) / sd)))
    res[name] = o; print(f"{name:30s} " + "  ".join(f"{k}: chi2/ndf={v['chi2']:.0f}/{v['ndf']} max|d|/sd={v['max_sig']:.1f}" for k, v in o.items()), flush=True)
json.dump(res, open(FIG / "s6b_audit_significance.json", "w"))
