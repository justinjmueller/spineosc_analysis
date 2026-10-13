#!/usr/bin/env python3
"""S2 (Q1): why do ratios trend with energy?  For each pair: ratio decomposition (selection vs energy migration), mean fractional energy shift of events selected in both, the local spectral slope,
the first-order prediction r-1 = (alpha-1)*delta, and the share of the energy shift carried by the muon and leading-proton kinetic energies."""
import json, sys, numpy as np
from ana import *
res = {}
for name in pairs_available():
    df = load(name); both = (df.in_n & df.in_v).to_numpy(); bn, bv = prep(df)
    En = df["n_reco_visible_energy"].to_numpy(float); Ev = df["v_reco_visible_energy"].to_numpy(float)
    d = np.where(both, (Ev - En) / En, np.nan)
    n0 = np.bincount(bn[bn >= 0], minlength=NBIN).astype(float); w = np.diff(BINS)
    dens = n0 / w; lnE = np.log(CTR)
    alpha = np.full(NBIN, np.nan)
    for i in range(NBIN):
        lo, hi = max(0, i - 1), min(NBIN - 1, i + 1); alpha[i] = -(np.log(dens[hi]) - np.log(dens[lo])) / (lnE[hi] - lnE[lo])
    delta = np.array([np.nanmedian(d[(bn == i) & both]) if ((bn == i) & both).sum() > 20 else np.nan for i in range(NBIN)]); dmean = np.array([np.nanmean(d[(bn == i) & both]) for i in range(NBIN)])
    dsig = np.array([np.nanstd(d[(bn == i) & both]) for i in range(NBIN)])
    fl = flows(df); r0, err, R, n0_, v0 = boot_ratio(df, nb=200)
    # contribution of each particle's KE shift to the shift of E_vis (events selected in both, muon and leading proton)
    dmu = df["v_reco_leading_muon_ke"].to_numpy(float) - df["n_reco_leading_muon_ke"].to_numpy(float); dp = df["v_reco_leading_proton_ke"].to_numpy(float) - df["n_reco_leading_proton_ke"].to_numpy(float)
    dE = Ev - En
    sel = both & np.isfinite(dmu) & np.isfinite(dp) & np.isfinite(dE)
    X = np.c_[dmu[sel], dp[sel]]; coef = np.linalg.lstsq(X, dE[sel], rcond=None)[0]; resid = dE[sel] - X @ coef
    r2 = 1 - resid.var() / dE[sel].var()
    res[name] = dict(ratio=r0.tolist(), err=err.tolist(), sel=(fl["gained"] + fl["lost"]).tolist(), mig=fl["migration"].tolist(), delta_median=delta.tolist(), delta_mean=dmean.tolist(), delta_sigma=dsig.tolist(), alpha=alpha.tolist(),
                     pred_first_order=((alpha - 1) * dmean).tolist(), dE_vs_dmu_dp_coef=coef.tolist(), dE_r2=float(r2), mean_dmu=float(np.nanmean(dmu[sel])), mean_dp=float(np.nanmean(dp[sel])), mean_dE=float(np.nanmean(dE[sel])),
                     frac_E_in_mu=float(np.nanmean(df["n_reco_leading_muon_ke"][both] / df["n_reco_visible_energy"][both])))
    r = res[name]; print(f"{name:32s} <dE>={r['mean_dE']*1e3:+7.2f} MeV  <dmu_ke>={r['mean_dmu']*1e3:+7.2f}  <dp_ke>={r['mean_dp']*1e3:+7.2f}  coef(mu,p)={np.round(coef,2)} R2={r2:.2f}  delta(1.0-1.2 GeV) median={100*np.nanmedian(delta[11:14]):+.2f}%  alpha~{np.nanmean(alpha[11:14]):.1f}", flush=True)
json.dump(res, open(FIG / "s2_energy_trend.json", "w"))
