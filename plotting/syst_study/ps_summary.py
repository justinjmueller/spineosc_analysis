#!/usr/bin/env python3
"""Collects core/*.json, pairs_compare.json, trigger.json into docs/syst_study/persyst/summary.json and prints the tables used in the report."""
from ps_common import *
CORE = PSFIG / "core"; C = {n: json.loads((CORE / f"{n}.json").read_text()) for n in PAIRS}
PC = json.loads((PSFIG / "pairs_compare.json").read_text()); TR = json.loads((PSFIG / "trigger.json").read_text())
A = lambda x: np.array([np.nan if v is None else v for v in x], float)
summ = {"pairs": {}}
def wfit(y, e):
    x = CTR; ok = np.isfinite(y) & np.isfinite(e) & (e > 0); w = 1 / e[ok] ** 2; X = np.vstack([np.ones(ok.sum()), x[ok]]).T
    cov = np.linalg.inv(X.T @ (X * w[:, None])); beta = cov @ (X.T @ (w * y[ok])); return float(beta[1]), float(np.sqrt(cov[1, 1])), float(np.sum(w * (y[ok] - X @ beta) ** 2)), int(ok.sum() - 2)
def excess(c, floor_bins):
    l = A(c["lost"]); n = A(c["n_nom"]); fl = np.nanmean(l[floor_bins]); ex = np.sum(n * np.maximum(0, l - fl)) / np.sum(n * l); return float(fl), float(ex)
for name, R in C.items():
    I = R["integrated"]; r = A(R["ratio"]) - 1; e = A(R["err"]); F = {k: A(v) for k, v in R["flows"].items()}
    w = 1 / e ** 2; mu = np.sum(w * r) / np.sum(w); var = np.sum(w * (r - mu) ** 2)
    share = {k: float(np.sum(w * (F[k] - np.sum(w * F[k]) / np.sum(w)) * (r - mu)) / var) for k in ("lost", "gained", "migration")}
    slope, eslope, chi2lin, dof = wfit(r, e)
    chi2_null = float(np.sum((r / e) ** 2))
    d = dict(integrated=I, share_var=share, trend_slope_per_GeV=slope, trend_slope_err=eslope, trend_chi2=chi2lin, trend_dof=dof, chi2_vs_unity=chi2_null,
             rms_dev=float(np.sqrt(np.mean(r ** 2))), low3=float(np.mean(r[:3])), high3=float(np.mean(r[-3:])), mid=float(np.mean(r[3:14])),
             churn_to_net=float(abs(I["net"]) / (I["lost_frac"] + I["gained_frac"])), E_unchanged=R["E_unchanged_frac"], sig_share=R["composition"]["nominal"][0],
             lost_sig_share=R["composition"]["lost"][0], truth_missing_gained=R["truth_missing_gained"], channel=R["channel"], mult_migration=R["mult_migration"],
             class_rates=R["class_rates"], side=R["side"], prod_match_maxdiff=R["prod_match_maxdiff"])
    d["shift_table"] = {k: (v["abs"]["mean"] if v and v.get("abs") else None) for k, v in R["shifts"].items() if k not in ("vertex_3d", "mu_end_3d", "mu_dir_angle_deg")}
    d["shift_rel"] = {k: (v["rel"]["mean"] if v.get("rel") else None) for k, v in R["shifts"].items() if k not in ("vertex_3d", "mu_end_3d", "mu_dir_angle_deg")}
    d["shift_rel_median"] = {k: (v["rel"]["median"] if v.get("rel") else None) for k, v in R["shifts"].items() if k not in ("vertex_3d", "mu_end_3d", "mu_dir_angle_deg")}
    d["vertex_shift_median_cm"] = R["shifts"]["vertex_3d"]["median"]; d["vertex_shift_p84_cm"] = R["shifts"]["vertex_3d"]["p84"]; d["dir_angle_median_deg"] = R["shifts"]["mu_dir_angle_deg"]["median"]
    d["frac_vertex_shift_gt5cm"] = None
    ex = {}
    for k, fb in (("mu_margin", slice(5, 10)), ("p_margin", slice(5, 10)), ("p_ke", slice(4, 9)), ("mu_ke", slice(3, 8))):
        if k in R["cond"]: ex[k] = dict(zip(("floor_lost_rate", "excess_share_of_lost"), excess(R["cond"][k], fb)))
    d["boundary"] = ex; d["smd"] = R["smd"]; d["thr"] = R["thr"]
    summ["pairs"][name] = d
summ["pairs_compare"] = {k: v for k, v in PC.items() if k in ("hamming", "joint_boot")}
summ["pairs_compare"]["SCE"] = {k: {a: b for a, b in v.items() if a not in ("odd", "even", "err", "lost_vs_drift")} for k, v in PC.items() if k.startswith("SCE_")}
summ["pairs_compare"]["YZ"] = {k: v for k, v in PC["YZ"].items() if k != "grid"}; summ["pairs_compare"]["YZ"]["grid_corr"] = PC["YZ"]["grid"]["corr_cells"]; summ["pairs_compare"]["YZ"]["grid_rms_ratio"] = PC["YZ"]["grid"]["ratio_rms"]
summ["trigger"] = {k: {a: b for a, b in v.items() if a not in ("pe_hist2d", "bins")} for k, v in TR.items() if k != "interplay"}
summ["trigger_interplay"] = {k: {a: b for a, b in v.items() if a not in ("r_all", "r_trig_own", "r_trig_nomPE")} for k, v in TR["interplay"].items()}
summ["fragility"] = json.loads((PSFIG / "fragile.json").read_text()); summ["fixdup"] = {k: v for k, v in json.loads((PSFIG / "fixdup.json").read_text()).items() if k.startswith("int") or k.startswith("lost") or k.startswith("gained")}
summ["duplicates"] = json.loads((PSFIG / "duplicates.json").read_text())
# charge-shift sensitivity of the energy trend: slope of ratio vs E_vis [%/GeV] and mean muon calo shift [%-pt]
sens = {}; bj = json.loads((REPO / "docs" / "syst_study" / "baseline.json").read_text()); sm = PC["YZ"]["shift_summary"]
for v in ("GainHigh", "ChargeScale", "ChargeScale15", "MattModYZ", "WireModYZ", "WireModxThetaXW"):
    b = bj[f"sbnd_run1__{v}"]; sl = wfit(np.array(b["ratio"], float) - 1, np.array(b["err"], float)); sens[f"sbnd_run1__{v}"] = dict(dmu_calo_pct=100 * sm[v]["dmu_calo"]["mean"], dmu_score_pct=100 * sm[v]["dmu_score"]["mean"], slope_pct_per_GeV=100 * sl[0], slope_err=100 * sl[1])
for n in PAIRS:
    if n.startswith("icarus"): sens[n] = dict(dmu_calo_pct=100 * summ["pairs"][n]["shift_table"]["mu_calo_csda"], dmu_score_pct=100 * summ["pairs"][n]["shift_table"]["mu_score"], slope_pct_per_GeV=100 * summ["pairs"][n]["trend_slope_per_GeV"], slope_err=100 * summ["pairs"][n]["trend_slope_err"])
summ["charge_sensitivity"] = sens
json.dump(summ, open(PS_JSON, "w"), indent=1)
# ------------------------- tables printed for the report
P = summ["pairs"]
print("TABLE1 net/churn")
print("pair | Nnom | lost% | gain% | net% | net/sig | |net|/churn | E_unch | rms r-1 | slope/GeV(err) | chi2 vs 1 | share lost/gain/mig")
for n, d in P.items():
    I = d["integrated"]; print(f"{n} | {I['n_nom']} | {100*I['lost_frac']:.2f} | {100*I['gained_frac']:.2f} | {100*I['net']:+.2f} | {I['net_sigma_churn']:+.1f} | {d['churn_to_net']:.2f} | {d['E_unchanged']:.2f} | {100*d['rms_dev']:.2f} | {100*d['trend_slope_per_GeV']:+.2f}({100*d['trend_slope_err']:.2f}) | {d['chi2_vs_unity']:.0f} | {d['share_var']['lost']:+.2f}/{d['share_var']['gained']:+.2f}/{d['share_var']['migration']:+.2f}")
print("TABLE2 shifts of kept events (mean): dE_rel% dmu_len_rel% dp_len_rel% | mu_score p_score | mu_calo p_calo | vtx_med_cm dir_med_deg")
for n, d in P.items():
    s = d["shift_rel"]; a = d["shift_table"]; print(f"{n} | {100*s['E_vis']['mean'] if isinstance(s['E_vis'],dict) else 100*s['E_vis']:+.3f} | {100*s['mu_len']:+.3f} | {100*s['p_len']:+.2f} | {a['mu_score']:+.4f} {a['p_score']:+.4f} | {a['mu_calo_csda']:+.4f} {a['p_calo_csda']:+.4f} | {d['vertex_shift_median_cm']:.2f} {d['dir_angle_median_deg']:.2f}")
print("TABLE3 boundary: floor_lost_rate / excess share (mu_margin, p_margin, p_ke)")
for n, d in P.items():
    b = d["boundary"]; print(n, {k: (round(100 * v["floor_lost_rate"], 1), round(100 * v["excess_share_of_lost"], 0)) for k, v in b.items()})
print("TABLE4 SMD lost/gained vs kept: p_margin, mu_margin, p_calo, mu_primary, p_primary, n_p, mu_len, p_ke")
for n, d in P.items():
    s = d["smd"]; print(n, {k: (round(s[k]["lost_vs_kept"], 2), round(s[k]["gained_vs_kept"], 2)) for k in ("p_margin", "mu_margin", "p_calo_csda", "mu_primary", "n_p", "mu_len", "p_ke") if k in s})
print("TABLE5 class rates (lost%, gained/nom%, net%) signal/otherCC/NC")
for n, d in P.items():
    cr = d["class_rates"]; print(n, "missG", round(d["truth_missing_gained"], 2), {k[:6]: (round(100 * v["lost"], 1), round(100 * v["gained_over_nom"], 1), round(100 * v["net"], 1)) for k, v in cr.items() if v["n_nom"] > 100})
print("TABLE6 channel: int all/1p/Mp, rms sys, rms d1p/dMp (noise-subtracted), chi2 d1p/dMp (n bins)")
for n, d in P.items():
    c = d["channel"]; print(f"{n} | {100*c['int_all']:+.2f} {100*c['int_1p']:+.2f} {100*c['int_Mp']:+.2f} | {100*c['rms_sys']:.2f} | {100*c['rms_d1p']:.2f}({100*c['rms_d1p_sub']:.2f}) {100*c['rms_dMp']:.2f}({100*c['rms_dMp_sub']:.2f}) | {c['chi2_d1p']:.0f}/{c['chi2_dMp']:.0f} ({c['nbin_ok']}) | mult mig {d['mult_migration']}")
