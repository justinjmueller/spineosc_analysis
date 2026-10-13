#!/usr/bin/env python3
"""Q2 deliverable 5: pool statistics. (i) ratio on the full per-variation pool vs on the COMMON pool of all (non-DENT) SBND variations / of all variations of that ICARUS run, with joint hashed bootstrap
(same event => same weight in both versions, so the covariance of the difference is right); (ii) effect of the ~0.6% duplicated keys (production = all rows vs pair table = unique keys);
(iii) anomalous-run scan (binomial test of the per-run lost/gained fractions) and the ratio without anomalous runs."""
from q2_lib import *
from scipy import stats
import pandas as pd, uproot
summ = json.loads((OUT / "pairs_summary.json").read_text()); sizes = json.load(open(Q2 / "q2_pool_sizes.json")); res = {}
NBS = 1000
for n in NAMES:
    d = load(n); det = det_of(n); lab = var_of(n); r, C, N, V = cov_analytic(d); o = dict(n_pool=summ[n]["n_pool"], n_common=sizes[det]["n_common"] if lab != "DENT" else None, n_rows_sel=int(len(d["key"])))
    # --- duplicates
    hn = np.array(summ[n]["hist_n"], float); hv = np.array(summ[n]["hist_v"], float); rr = hv / hn
    prod = uproot.open(REPO_PATH / "detvar" / "detsys-v1.1.0v2.root")[summ[n]["hist"]].values() if False else None
    dupn = hn - N; dupv = hv - V; sig = np.sqrt(np.diag(C)); dr = rr - r
    sig_dup = np.sqrt(np.clip(dupn * r ** 2 + dupv, 0, None)) / N     # conservative: dup rows treated as independent unpaired events
    o["dup"] = dict(frac_rows=float(summ[n]["dup_n"] / summ[n]["n_n_rows"]), dr=dr.tolist(), dr_over_sigma=(dr / sig).tolist(), maxabs_dr_over_sigma=float(np.abs(dr / sig).max()), rms_dr_over_sigma=float(np.sqrt(np.mean((dr / sig) ** 2))),
                    chi2_dr_over_stat=float(np.sum((dr / sig) ** 2)), dr_over_sigdup_rms=float(np.sqrt(np.mean(np.where(sig_dup > 0, dr / np.where(sig_dup > 0, sig_dup, 1), 0) ** 2))), rms_dr=float(np.sqrt(np.mean(dr ** 2))), rms_sigma=float(np.sqrt(np.mean(sig ** 2))))
    # --- common pool
    if lab != "DENT":
        flag = np.load(RES / f"common_{n}.npy"); rc = ratio(d, None, flag); nc, vc = counts(d, None, flag)
        Rf = np.empty((NBS, NB)); Rc = np.empty((NBS, NB))
        for k in range(NBS): w = pois_w(d["key"], k, 11); nn, vv = counts(d, w); Rf[k] = vv / nn; nn, vv = counts(d, w, flag); Rc[k] = vv / nn
        diff = r - rc; Cd = np.cov((Rf - Rc).T); sd = np.sqrt(np.diag(Cd)); z = diff / sd
        # analytic: independent halves A = common, B = removed; var(diff) = (NB/N)^2 (var rB + var rA) -> use bootstrap (joint); chi2 with bootstrap covariance
        chi2 = float(diff @ np.linalg.solve(Cd, diff)); sc = np.sqrt(np.diag(np.cov(Rc.T)))
        o["common"] = dict(frac_rows=float(flag.mean()), r_full=r.tolist(), r_common=rc.tolist(), diff=diff.tolist(), sd_diff=sd.tolist(), z=z.tolist(), chi2=chi2, ndf=NB, p=float(stats.chi2.sf(chi2, NB)), z_sd=float(z.std()), z_mean=float(z.mean()),
                           sigma_full=sig.tolist(), sigma_common_over_full=float(np.median(sc / np.sqrt(np.diag(np.cov(Rf.T))))), maxabs_diff_over_sigma_full=float(np.abs(diff / sig).max()), rms_diff=float(np.sqrt(np.mean(diff ** 2))))
    # --- anomalous run scan
    import pandas as pd
    x = pd.DataFrame(dict(run=d["run"], inn=d["in_n"], inv=d["in_v"])); x["g"] = x.run if det.startswith("icarus") else x.run // 10; x["lost"] = x.inn & ~x.inv; x["gain"] = ~x.inn & x.inv
    g = x.groupby("g").agg(n=("inn", "sum"), lost=("lost", "sum"), gain=("gain", "sum")); g = g[g.n >= 20]; p = g.lost.sum() / g.n.sum(); q = g.gain.sum() / g.n.sum()
    zl = (g.lost - p * g.n) / np.sqrt(g.n * p * (1 - p)); zg = (g.gain - q * g.n) / np.sqrt(g.n * q * (1 - q)); bad = list(zl.index[(zl.abs() > 5) | (zg.abs() > 5)])
    o["runscan"] = dict(n_blocks=int(len(g)), zl_sd=float(zl.std()), zg_sd=float(zg.std()), zl_maxabs=float(zl.abs().max()), zg_maxabs=float(zg.abs().max()), flagged=[int(b) for b in bad])
    if bad:
        runs = d["run"] if det.startswith("icarus") else d["run"] // 10; m = ~np.isin(runs, bad); rcl = ratio(d, None, m); nc, vc = counts(d, None, m)
        Cc = cov_analytic({k: (v[m] if hasattr(v, "__len__") and len(v) == len(m) else v) for k, v in d.items()})[1]
        o["runscan"].update(r_clean=rcl.tolist(), diff=(rcl - r).tolist(), z_vs_sigma=((rcl - r) / sig).tolist(), maxabs_diff_pct=float(100 * np.abs(rcl - r).max()), n_removed_nom=float(N.sum() - nc.sum()), frac_removed=float(1 - nc.sum() / N.sum()))
    res[n] = o
    cm = o.get("common"); print(f"{n:30s} dup: frac {100*o['dup']['frac_rows']:.2f}% max|dr|/sig {o['dup']['maxabs_dr_over_sigma']:.2f} rms {o['dup']['rms_dr_over_sigma']:.2f} | " + (f"common pool {100*cm['frac_rows']:.0f}% rows: chi2/ndf {cm['chi2']/NB:.2f} (p {cm['p']:.2f}) z-sd {cm['z_sd']:.2f} max|d|/sig {cm['maxabs_diff_over_sigma_full']:.2f}; sig_common/sig_full {cm['sigma_common_over_full']:.2f}" if cm else "no common pool (disjoint CV sample)") + (f" | FLAGGED runs {o['runscan']['flagged']}" if bad else ""), flush=True)
json.dump(res, open(Q2 / "q2_pool.json", "w"))
