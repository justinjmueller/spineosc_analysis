#!/usr/bin/env python3
"""Q2 summary: classification of each variation, recommendation table, one summary JSON (docs/syst_study/q2/q2_summary.json)."""
from q2_lib import *
from scipy import stats
sn = json.load(open(Q2 / "q2_signal_noise.json")); cov = json.load(open(Q2 / "q2_covariance.json")); cal = json.load(open(Q2 / "q2_calibration.json")); pool = json.load(open(Q2 / "q2_pool.json")); val = json.load(open(Q2 / "q2_validation.json")); cons = json.load(open(Q2 / "q2_consistency.json"))
P01 = 0.01
def classify(o):
    tab = {t["deg"]: t for t in o["models"]}; bd = o["deg_bic"]
    if o["p_zero"] > P01 or bd == -1: cls = "iii"
    else:
        if bd <= 0: cls = "ii"
        else:
            k = tab[bd]["k"]; dchi = tab[0]["chi2"] - tab[bd]["chi2"]; pD = float(stats.chi2.sf(dchi, k - 1)); cls = "i" if pD < P01 else "ii"
    flags = []
    if o["resid_p_bic"] < 0.05: flags.append("excess non-smooth structure (p=%.3f)" % o["resid_p_bic"])
    if o["resid_p_bic"] > 0.995: flags.append("under-dispersed")
    return cls, flags
NAMEC = {"i": "(i) smooth physical trend", "ii": "(ii) constant offset", "iii": "(iii) consistent with zero/noise"}
summary = dict(calibration={}, per_variation={}, detector={}, consistency={}, pool={}, validation={})
rows = []
for n in ALL:
    o = sn[n]; c = cov["per_variation"][n]; cls, fl = classify(o); y = np.array(o["y"]); s = np.array(o["sigma"])
    exp_abs = float(np.mean(0.7979 * s)); obs_abs = float(np.mean(np.abs(y)))
    pv = dict(cls=cls, flags=fl, deg_bic=o["deg_bic"], deg_aic=o["deg_aic"], deg_loo=o["deg_loo"], chi2_zero=o["chi2_zero"], p_zero=o["p_zero"], chi2_const=o["chi2_const"], p_const=o["p_const"], resid_chi2_ndf=o["resid_chi2_ndf_bic"], resid_p=o["resid_p_bic"],
              n_within2sig=o["n_within2sig_zero"], rms_y=o["rms_y"], rms_noise=o["rms_noise"], rms_signal_fit=o["rms_signal_fit"], rms_signal_debiased=o["rms_signal_debiased"], snr=o["rms_signal_debiased"] / o["rms_noise"],
              mean_abs_y=obs_abs, expected_mean_abs_if_pure_noise=exp_abs, const=o["const_wmean"], const_err=o["const_err"], fnoise_d1=o["fnoise_d1"], fnoise_d2=o["fnoise_d2"], chi2_d1_ndf=o["chi2_d1_ndf"], chi2_d2_ndf=o["chi2_d2_ndf"],
              noise_trace_frac=c["noise_frac"], noise_trace_frac_bulk=c["noise_frac_bulk"], last_bin_share_of_noise=c["last_bin_share_of_noise"], trace_raw=c["trace_raw"], trace_fit=c["trace_fit"], trace_noise=c["trace_noise"], sign_flip_frac=c["sign_flip_frac"])
    summary["per_variation"][n] = pv; rows.append((n, pv))
for n in NAMES: 
    cc = cal[n]; summary["calibration"][n] = dict(bs_over_an_median=cc["bs_over_an_median"], naive_over_paired=cc["naive_over_paired_median"], split_vr=cc["split_event"]["vr_all"], split_chi2ndf=cc["split_event"]["chi2ndf_mean"], block_split_chi2ndf=cc["split_block"]["chi2ndf_mean"], jk_over_ev=cc["block"]["jk_med"], bb_over_ev=cc["block"]["bb_med"], G=cc["block"]["G"], sigma_range=[min(cc["sig_an"]), max(cc["sig_an"])], sigma_median=float(np.median(cc["sig_an"])))
summary["detector"] = cov["per_detector"]; summary["consistency"] = {t["label"]: {k: v for k, v in t.items() if not isinstance(v, list)} for t in cons["tests"]}; summary["consistency"]["icarus_pooled"] = cons["icarus_pooled"]
summary["pool"] = {n: dict(n_pool=pool[n]["n_pool"], common_chi2_ndf=(pool[n]["common"]["chi2"] / NB if "common" in pool[n] else None), common_p=(pool[n]["common"]["p"] if "common" in pool[n] else None), common_max_diff_over_sigma=(pool[n]["common"]["maxabs_diff_over_sigma_full"] if "common" in pool[n] else None),
              dup_frac=pool[n]["dup"]["frac_rows"], dup_rms_dr_over_sigma=pool[n]["dup"]["rms_dr_over_sigma"], flagged_runs=pool[n]["runscan"]["flagged"]) for n in NAMES}
summary["validation"] = dict(tau_hat=val["tau_hat"], tau_upper68=val["tau_interval"][1], methods={m: {tn: float(np.mean([val["res"][n][tn][m]["mse_ratio"] for n in NAMES])) for tn in val["res"][NAMES[0]]} for m in val["res"][NAMES[0]]["zero"]})
json.dump(summary, open(Q2 / "q2_summary.json", "w"), indent=1)
print(f"{'variation':34s} cls  degBIC  p_zero    chi2c/16  resid chi2/ndf  rms_y  rms_n  snr   noise-trace  const+-err    exp|y|  obs|y| within2s  flags")
for n, p in rows: print(f"{n:34s} {p['cls']:3s} {p['deg_bic']:3d}/{p['deg_aic']:2d}/{p['deg_loo']:2d} {p['p_zero']:9.2g} {p['chi2_const']/16:7.2f} {p['resid_chi2_ndf']:7.2f}(p{p['resid_p']:.2f}) {100*p['rms_y']:5.2f} {100*p['rms_noise']:5.2f} {p['snr']:5.1f} {p['noise_trace_frac']:6.2f} {100*p['const']:+6.2f}+-{100*p['const_err']:.2f}  {100*p['expected_mean_abs_if_pure_noise']:5.2f} {100*p['mean_abs_y']:5.2f}  {p['n_within2sig']:2d}  {'; '.join(p['flags'])}")
# ---- pooled statistics used in the report
from q2_fit import gls
pooled = {}
for dg in (2, 3, 4):
    c = sum(sn[n]["models"][dg + 1]["chi2"] for n in NAMES); nd = sum(17 - sn[n]["models"][dg + 1]["k"] for n in NAMES); cdf = np.array([stats.chi2.cdf(sn[n]["models"][dg + 1]["chi2"], 17 - sn[n]["models"][dg + 1]["k"]) for n in NAMES])
    pooled[f"fixed_deg{dg}"] = dict(chi2=float(c), ndf=int(nd), chi2_ndf=float(c / nd), p_pooled_two_sided=float(2 * min(stats.chi2.cdf(c, nd), stats.chi2.sf(c, nd))), ks_p=float(stats.kstest(cdf, "uniform").pvalue), n_cdf_lt_0p01=int((cdf < 0.01).sum()), n_cdf_gt_0p99=int((cdf > 0.99).sum()))
pooled["classes"] = {k: [n for n in NAMES if summary["per_variation"][n]["cls"] == k] for k in ("i", "ii", "iii")}
cm = [pool[n]["common"] for n in NAMES if "common" in pool[n]]; pooled["common_pool"] = dict(chi2=float(sum(c["chi2"] for c in cm)), ndf=NB * len(cm), p=float(stats.chi2.sf(sum(c["chi2"] for c in cm), NB * len(cm))), z_sd=float(np.std(np.concatenate([c["z"] for c in cm]))), max_abs_diff_over_sigma=float(max(c["maxabs_diff_over_sigma_full"] for c in cm)))
pooled["dup"] = dict(median_rms_dr_over_sigma_excl_run4wiremod=float(np.median([pool[n]["dup"]["rms_dr_over_sigma"] for n in NAMES if n != "icarus_run4__WireModxThetaXW"])), max_excl=float(max(pool[n]["dup"]["maxabs_dr_over_sigma"] for n in NAMES if n != "icarus_run4__WireModxThetaXW")))
summary["pooled"] = pooled; json.dump(summary, open(Q2 / "q2_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in pooled.items() if k != "classes"}, indent=0)[:1500]); print({k: len(v) for k, v in pooled["classes"].items()})
