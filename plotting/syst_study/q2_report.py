#!/usr/bin/env python3
import pathlib
"""Writes docs/syst_study/Q2_statistical_content.md from the JSON results (tables) plus the analysis text."""
from q2_lib import *
S = json.load(open(Q2 / "q2_summary.json")); cal = json.load(open(Q2 / "q2_calibration.json")); cov = json.load(open(Q2 / "q2_covariance.json")); cons = json.load(open(Q2 / "q2_consistency.json")); pool = json.load(open(Q2 / "q2_pool.json")); val = json.load(open(Q2 / "q2_validation.json")); sc = json.load(open(Q2 / "q2_splitcv.json")); cl = json.load(open(Q2 / "q2_run4_wiremodx_clean.json")); pl = S["pooled"]
sh = lambda n: n.replace("sbnd_run1__", "SBND ").replace("icarus_run2__", "ICARUS-R2 ").replace("icarus_run4__", "ICARUS-R4 ").replace("clean", " (clean)")
CL = {"i": "(i) smooth trend", "ii": "(ii) constant", "iii": "(iii) ~zero/noise"}
REC = {"i": "smooth fit (poly-BIC), last bin per pre-test", "ii": "fit a constant (+ last-bin pre-test); quote offset +- err", "iii": "set to zero-model / drop; keep only as normalisation-type bound"}
T1 = "| pair | sigma bulk median [%] | sigma last bin [%] | bootstrap/analytic | naive/paired | split-half chi2/ndf (event / block) | jackknife/event | block-bootstrap/event |\n|---|---|---|---|---|---|---|---|\n"
for n in NAMES:
    c = cal[n]; s = np.array(c["sig_an"]); T1 += f"| {sh(n)} | {100*np.median(s[:-1]):.2f} | {100*s[-1]:.2f} | {c['bs_over_an_median']:.3f} | {c['naive_over_paired_median']:.1f} | {c['split_event']['chi2ndf_mean']:.2f} / {c['split_block']['chi2ndf_mean']:.2f} | {c['block']['jk_med']:.3f} | {c['block']['bb_med']:.3f} |\n"
T2 = "| variation | class | deg (BIC/AIC/LOO) | chi2(0)/17 | chi2(const)/16 | resid chi2/ndf (p) | #bins <2sigma | rms(r-1) [%] | rms noise [%] | signal/noise | frac. of d1/d2 variance that is noise | const +- err [%] |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
for n in ALL:
    p = S["per_variation"][n]; T2 += f"| {sh(n)} | {CL[p['cls']]}{' *' if p['flags'] else ''} | {p['deg_bic']}/{p['deg_aic']}/{p['deg_loo']} | {p['chi2_zero']/17:.1f} | {p['chi2_const']/16:.1f} | {p['resid_chi2_ndf']:.2f} ({p['resid_p']:.2f}) | {p['n_within2sig']} | {100*p['rms_y']:.2f} | {100*p['rms_noise']:.2f} | {p['snr']:.1f} | {p['fnoise_d1']:.2f} / {p['fnoise_d2']:.2f} | {100*p['const']:+.2f} +- {100*p['const_err']:.2f} |\n"
T3 = "| variation | trace raw [1e-4] | noise trace [1e-4] | noise/raw | fit trace [1e-4] | share of noise in last bin | sign-flip frac (bulk pairs) |\n|---|---|---|---|---|---|---|\n"
for n in ALL:
    c = cov["per_variation"][n]; T3 += f"| {sh(n)} | {1e4*c['trace_raw']:.1f} | {1e4*c['trace_noise']:.1f} | {c['noise_frac']:.2f} | {1e4*c['trace_fit']:.1f} | {c['last_bin_share_of_noise']:.2f} | {c['sign_flip_frac']:.2f} |\n"
T4 = "| detector / set | trace raw [1e-4] | trace fit | noise trace | noise/raw | eff. rank raw -> fit (MC of fit+noise) | frac. first eigenvalue raw -> fit | overlap e1(raw), e1(fit) (MC mean [16,84]%) | adjacent corr raw -> fit |\n|---|---|---|---|---|---|---|---|---|\n"
for g, d in cov["per_detector"].items():
    m = d["mc"]; T4 += f"| {g.replace(chr(124), ' / ')} | {1e4*d['trace_raw']:.1f} | {1e4*d['trace_fit']:.1f} | {1e4*d['trace_noise']:.1f} | {d['noise_frac']:.2f} | {d['eig_raw']['eff_rank']:.2f} -> {d['eig_fit']['eff_rank']:.2f} ({m['eff_rank_mean']:.2f}) | {d['eig_raw']['frac1']:.2f} -> {d['eig_fit']['frac1']:.2f} | {d['e1_overlap_raw_fit']:.3f} ({m['e1_overlap_mean']:.3f} [{m['e1_overlap_p16']:.3f},{m['e1_overlap_p84']:.3f}]) | {d['adj_corr_raw']:.2f} -> {d['adj_corr_fit']:.2f} |\n"
T5 = "| test | chi2/ndf | p | notes |\n|---|---|---|---|\n"
for t in cons["tests"]:
    if "chi2_equal" in t: T5 += f"| {t['label']}: equality | {t['chi2_equal']/NB:.2f} | {t['p_equal']:.2g} | z mean {t['z_mean']:+.2f}, z sd {t['z_sd']:.2f}; proportional model alpha={t['alpha_prop']:.2f}: chi2/16={t['chi2_prop']/16:.2f} (p={t['p_prop']:.3f}) |\n"
    if "chi2_antisym" in t: T5 += f"| {t['label']} | {t['chi2_antisym']/NB:.2f} | {t['p_antisym']:.3g} | rms linear {100*t['rms_lin']:.2f}% (sigma {100*t['rms_sigma_lin']:.2f}%), rms quadratic {100*t['rms_quad']:.2f}% (sigma {100*t['rms_sigma_quad']:.2f}%); quad=0 p={t['p_quad_zero']:.3g} |\n"
T6 = "| pair | pool [M ev] | rows in common pool | full-vs-common chi2/ndf (p) | max abs diff/sigma | dup. rows [%] | rms (all rows - unique)/sigma | flagged runs |\n|---|---|---|---|---|---|---|---|\n"
for n in NAMES:
    p = pool[n]; c = p.get("common"); T6 += f"| {sh(n)} | {p['n_pool']/1e6:.2f} | {('%.0f%%' % (100*c['frac_rows'])) if c else 'n/a'} | {('%.2f (%.2f)' % (c['chi2']/NB, c['p'])) if c else 'n/a'} | {('%.2f' % c['maxabs_diff_over_sigma_full']) if c else 'n/a'} | {100*p['dup']['frac_rows']:.2f} | {p['dup']['rms_dr_over_sigma']:.2f} | {p['runscan']['flagged'] or ''} |\n"
ns = NAMES; sb = lambda m, k: np.mean([sc[n][k][m] / sc[n]["floor_" + k] for n in ns])
T7 = "| estimator | bulk bins: T/floor (mean) | share of removable noise removed | last bin: T/floor (mean; #variations worse than raw) |\n|---|---|---|---|\n"
for m in sc[ns[0]]["bulk"]:
    ex = np.mean([(sc[n]["bulk"]["raw"] - sc[n]["bulk"][m]) / (sc[n]["bulk"]["raw"] - sc[n]["floor_bulk"]) for n in ns]); w = sum(sc[n]["last"][m] / sc[n]["floor_last"] > sc[n]["last"]["raw"] / sc[n]["floor_last"] + 1e-9 for n in ns)
    T7 += f"| {m} | {sb(m,'bulk'):.2f} | {ex:.2f} | {sb(m,'last'):.2f} ({w}) |\n"
ms = ["raw", "poly-BIC", "whit-CV", "poly-BIC, last raw", "whit-CV, last raw", "poly-BIC, last pretest(2sig)", "whit-CV, last pretest(2sig)"]; trs = list(val["res"][NAMES[0]])
T8 = "| estimator | " + " | ".join(trs) + " |\n|---|" + "---|" * len(trs) + "\n"
for m in ms: T8 += f"| {m} | " + " | ".join(f"{np.mean([val['res'][n][t][m]['mse_ratio'] for n in NAMES]):.2f}" for t in trs) + " |\n"
T9 = "| pair | class | recommended treatment | const +- err [%] | noise/raw trace |\n|---|---|---|---|---|\n"
for n in ALL:
    p = S["per_variation"][n]; extra = ""
    if n == "icarus_run4__WireModxThetaXW": extra = "; USE THE CLEAN VERSION (duplicate-event runs removed), the production histogram is biased by +0.5%"
    if p["cls"] == "i" and n.endswith("DENT"): extra = "; slope significant (p~1e-4) but amplitude (1%) < per-bin noise (1.5%): smoothing essential"
    if p["flags"]: extra += "; residual structure p<0.05, keep last-bin raw"
    T9 += f"| {sh(n)} | {CL[p['cls']]} | {REC[p['cls']]}{extra} | {100*p['const']:+.2f} +- {100*p['const_err']:.2f} | {p['noise_trace_frac']:.2f} |\n"
v = S["validation"]["methods"]; f2 = pl["fixed_deg2"]; f3 = pl["fixed_deg3"]; f4 = pl["fixed_deg4"]; cp = pl["common_pool"]; ic = cons["icarus_pooled"]
tx = open(pathlib.Path(__file__).parent / "q2_report_text.md").read()
for k, val_ in dict(T1=T1, T2=T2, T3=T3, T4=T4, T5=T5, T6=T6, T7=T7, T8=T8, T9=T9, TAU=f"{val['tau_hat']:.4f}", TAUUL=f"{val['tau_interval'][1]:.4f}", F3=f"{f3['chi2']:.1f}/{f3['ndf']} = {f3['chi2_ndf']:.2f}", F3KS=f"{f3['ks_p']:.2f}", F2=f"{f2['chi2_ndf']:.2f}", F4=f"{f4['chi2_ndf']:.2f}", CPC=f"{cp['chi2']:.0f}/{cp['ndf']} (p={cp['p']:.2f})", CPZ=f"{cp['z_sd']:.2f}",
                    ICC=f"{ic['chi2']:.0f}/{ic['ndf']} (p={ic['p']:.2f}), z sd {ic['z_sd']:.2f}", PRODCLEAN=f"{100*cl['const_wmean']:+.2f} +- {100*cl['const_err']:.2f}", PRODCHI=f"{cl['chi2']:.0f}").items(): tx = tx.replace("{{" + k + "}}", str(val_))
(FIG / "Q2_statistical_content.md").write_text(tx); print("report written", len(tx))
