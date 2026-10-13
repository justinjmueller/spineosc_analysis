#!/usr/bin/env python3
"""Q2 deliverable 2: signal vs noise in the 20 ratio shapes. GLS fits with the analytic (bootstrap-validated) full bin covariance; model selection by AIC/BIC/LOO-CV and by a split-half out-of-sample test."""
from q2_lib import *
from q2_fit import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
DEGS = list(range(-1, 7)); NSPL = 200
def second_diff_stats(y, S):
    D1 = np.zeros((NB - 1, NB)); D2 = np.zeros((NB - 2, NB))
    for i in range(NB - 1): D1[i, i] = -1; D1[i, i + 1] = 1
    for i in range(NB - 2): D2[i, i:i + 3] = [1, -2, 1]
    v1 = np.diag(D1 @ S @ D1.T); v2 = np.diag(D2 @ S @ D2.T); a1 = D1 @ y; a2 = D2 @ y
    return dict(fnoise_d1=float(min(1, v1.sum() / (a1 ** 2).sum())), chi2_d1_ndf=float((a1 ** 2 / v1).mean()), chi2_d2_ndf=float((a2 ** 2 / v2).mean()), fnoise_d2=float(min(1, v2.sum() / (a2 ** 2).sum())))
def splitcv(d, sfull, nsplit=NSPL):
    """fit on half A, score on half B (and v.v.). T = sum_b (yB - fA)^2/sfull_b^2 ; noise floor of B: sum sB^2/sfull^2 (= 2 per bin)."""
    w = 1 / sfull ** 2; res = {k: [] for k in ["raw", "zero"] + [f"deg{k}" for k in range(0, 7)] + ["bicsel", "whit", "floor"]}
    for s in range(nsplit):
        m = (uni(d["key"], 7000 + s) < 0.5).astype(float)
        for wa, wb in ((m, 1 - m), (1 - m, m)):
            rA, CA, _, _ = cov_analytic(d, w=wa); rB, CB, _, _ = cov_analytic(d, w=wb); yA, yB = rA - 1, rB - 1
            res["raw"].append(np.sum(w * (yB - yA) ** 2)); res["zero"].append(np.sum(w * yB ** 2)); res["floor"].append(np.sum(w * np.diag(CB)))
            bics = []
            for dg in range(0, 7):
                f, i = gls(yA, CA, dg); res[f"deg{dg}"].append(np.sum(w * (yB - f) ** 2)); bics.append((i["chi2"] + i["k"] * np.log(NB), f))
            f0, i0 = gls(yA, CA, -1); b0 = i0["chi2"]; best = min([(b0, f0)] + bics, key=lambda t: t[0]); res["bicsel"].append(np.sum(w * (yB - best[1]) ** 2))
            f, i = whit_cv(yA, CA); res["whit"].append(np.sum(w * (yB - f) ** 2))
    return {k: float(np.mean(v)) for k, v in res.items()}
out = {}
for n in ALL:
    d = load(n); r, C, N, V = cov_analytic(d); y = r - 1; s = np.sqrt(np.diag(C)); o = dict(name=n, y=y.tolist(), sigma=s.tolist(), N=N.tolist())
    fits = {}; tab = []
    for dg in DEGS:
        f, i = gls(y, C, dg); fits[dg] = f; tab.append(dict(deg=dg, k=i["k"], chi2=i["chi2"], aic=i["chi2"] + 2 * i["k"], bic=i["chi2"] + i["k"] * np.log(NB), loo=i["loo"], p=pval(i["chi2"], NB - i["k"])))
    fw, iw = whit_cv(y, C); o["whit"] = dict(f=fw.tolist(), lam=float(iw["lam"]), edf=iw["edf"], chi2=iw["chi2"], loo=iw["loo"])
    o["models"] = tab; bd = min(tab, key=lambda t: t["bic"]); ad = min(tab, key=lambda t: t["aic"]); cd = min(tab, key=lambda t: t["loo"])
    o["deg_bic"], o["deg_aic"], o["deg_loo"] = bd["deg"], ad["deg"], cd["deg"]
    chi0 = tab[0]; chic = tab[1]
    o["chi2_zero"], o["p_zero"] = chi0["chi2"], chi0["p"]; o["chi2_const"], o["p_const"] = chic["chi2"], chic["p"]
    o["n_within2sig_zero"] = int(np.sum(np.abs(y) < 2 * s)); o["n_within1sig_zero"] = int(np.sum(np.abs(y) < s))
    fb = fits[bd["deg"]] if bd["deg"] >= -1 else fits[-1]; o["fit_bic"] = fb.tolist(); kb = bd["k"]; o["resid_chi2_ndf_bic"] = bd["chi2"] / max(NB - kb, 1); o["resid_p_bic"] = bd["p"]
    # amplitudes
    o["rms_y"] = float(np.sqrt(np.mean(y ** 2))); o["rms_noise"] = float(np.sqrt(np.mean(s ** 2))); o["rms_signal_fit"] = float(np.sqrt(np.mean(fb ** 2))); o["rms_signal_whit"] = float(np.sqrt(np.mean(fw ** 2)))
    o["rms_signal_debiased"] = float(np.sqrt(max(0.0, np.mean(y ** 2) - np.mean(s ** 2)))); o["noise_trace_frac"] = float(np.sum(s ** 2) / np.sum(y ** 2))
    wc = 1 / s ** 2; o["const_wmean"] = float(np.sum(wc * y) / wc.sum()); o["const_err"] = float(1 / np.sqrt(wc.sum()))
    o["maxabs_over_sigma"] = float(np.max(np.abs(y) / s)); o["delta_chi2_const_to_best"] = float(chic["chi2"] - bd["chi2"]) if bd["deg"] >= 1 else 0.0
    o.update(second_diff_stats(y, C))
    # edge test: leave last bin out (and first bin), GLS deg 1..3 + whittaker, pull of the held-out bin
    edge = {}
    for tag, hold in (("last", NB - 1), ("first", 0)):
        keep = np.array([i for i in range(NB) if i != hold]); Sk = C[np.ix_(keep, keep)]; yk = y[keep]; xk = LOGX[keep]
        for dg in (1, 2, 3):
            f, i = gls(yk, Sk, dg, xk); beta = i["beta"]; Xh = basis(dg, LOGX[[hold]]); pred = float((Xh @ beta)[0])
            Xk = basis(dg, xk); Wk = whiten(Sk); cov_b = np.linalg.inv((Wk @ Xk).T @ (Wk @ Xk)); vp = float((Xh @ cov_b @ Xh.T)[0, 0])
            edge[f"{tag}_deg{dg}"] = dict(pred=pred, obs=float(y[hold]), pull=float((y[hold] - pred) / np.sqrt(C[hold, hold] + vp)), sd_pred=float(np.sqrt(vp)))
        fk, ik = whit_cv(yk, Sk, penalty(xk)); 
        # whittaker prediction by linear extrapolation of the last two smoothed points (2nd-diff penalty null space)
        if tag == "last": pred = fk[-1] + (fk[-1] - fk[-2]) * (LOGX[hold] - xk[-1]) / (xk[-1] - xk[-2])
        else: pred = fk[0] + (fk[1] - fk[0]) * (LOGX[hold] - xk[0]) / (xk[1] - xk[0])
        edge[f"{tag}_whit"] = dict(pred=float(pred), obs=float(y[hold]), pull=float((y[hold] - pred) / C[hold, hold] ** 0.5))
    o["edge"] = edge
    o["splitcv"] = splitcv(d, s); out[n] = o
    sc = o["splitcv"]
    print(f"{n:30s} rms y {100*o['rms_y']:.2f}% noise {100*o['rms_noise']:.2f}% | chi2 zero {chi0['chi2']:8.1f}/17 const {chic['chi2']:7.1f}/16 | deg BIC {bd['deg']:2d} AIC {ad['deg']:2d} LOO {cd['deg']:2d} resid chi2/ndf {o['resid_chi2_ndf_bic']:.2f} p {o['resid_p_bic']:.3f} | noise-trace frac {o['noise_trace_frac']:.3f} | within2sig {o['n_within2sig_zero']} | splitcv raw {sc['raw']/sc['floor']:.2f} bic {sc['bicsel']/sc['floor']:.2f} whit {sc['whit']/sc['floor']:.2f} deg1 {sc['deg1']/sc['floor']:.2f} (floor=1)", flush=True)
json.dump(out, open(Q2 / "q2_signal_noise.json", "w"))
