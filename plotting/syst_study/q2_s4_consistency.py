#!/usr/bin/env python3
"""Q2 deliverable 4: consistency tests between related samples as noise probes. Covariances between ratio vectors of different pair tables use the exact analytic cross-covariance over shared events (SBND variations share the nominal sample and the pool events); ICARUS Run2/Run4 are independent samples (zero shared keys)."""
from q2_lib import *
from q2_fit import *
from scipy import stats
O = {}
D = {n: load(n) for n in NAMES}; Y = {}; S = {}
for n in NAMES: r, C, N, V = cov_analytic(D[n]); Y[n] = r - 1; S[n] = C
def chi2(delta, Sd):
    c = float(delta @ np.linalg.solve(Sd, delta)); return c
def prop_fit(ya, Sa, yb, Sb, X, sign=+1.0):
    """model ya = alpha*sign*yb ; scan alpha; Var(delta)=Sa + alpha^2 Sb - sign*alpha (X + X^T) with X=cov(ya,yb)."""
    best = None
    for a in np.linspace(0.0, 12.0, 2401):
        dl = ya - sign * a * yb; Sd = Sa + a * a * Sb - sign * a * (X + X.T); c = chi2(dl, Sd)
        if best is None or c < best[0]: best = (c, a)
    return best
def pair_test(a, b, label, kind):
    if kind == "indep": X = np.zeros((NB, NB))
    else: X = cross_cov(D[a], D[b])
    ya, yb = Y[a], Y[b]; Sa, Sb = S[a], S[b]; res = dict(label=label, a=a, b=b)
    if kind == "equal" or kind == "indep":
        dl = ya - yb; Sd = Sa + Sb - X - X.T; z = dl / np.sqrt(np.diag(Sd)); res.update(delta=dl.tolist(), sd_delta=np.sqrt(np.diag(Sd)).tolist(), chi2_equal=chi2(dl, Sd), z_mean=float(z.mean()), z_sd=float(z.std()), z=z.tolist())
        res["p_equal"] = float(stats.chi2.sf(res["chi2_equal"], NB))
        # smooth-vs-residual decomposition of the difference: GLS poly deg<=2 on dl with Sd
        for dg in (0, 1, 2, 3):
            f, i = gls(dl, Sd, dg); res[f"diff_deg{dg}_chi2"] = i["chi2"]
        res["rho_noise"] = float(np.mean([X[i, i] / np.sqrt(Sa[i, i] * Sb[i, i]) for i in range(NB)]))
        # proportionality (shape equal up to scale) for non-independent pairs
        c, al = prop_fit(ya, Sa, yb, Sb, X, +1.0); res.update(alpha_prop=float(al), chi2_prop=float(c), p_prop=float(stats.chi2.sf(c, NB - 1)))
    if kind == "antisym":
        dl = ya + yb; Sd = Sa + Sb + X + X.T; z = dl / np.sqrt(np.diag(Sd)); res.update(sum=dl.tolist(), sd_sum=np.sqrt(np.diag(Sd)).tolist(), chi2_antisym=chi2(dl, Sd), z_mean=float(z.mean()), z_sd=float(z.std()), z=z.tolist(), p_antisym=float(stats.chi2.sf(chi2(dl, Sd), NB)))
        c, al = prop_fit(ya, Sa, yb, Sb, X, -1.0); res.update(alpha_prop=float(al), chi2_prop=float(c), p_prop=float(stats.chi2.sf(c, NB - 1)))
        c2, al2 = prop_fit(ya, Sa, yb, Sb, X, +1.0); res.update(alpha_same=float(al2), chi2_same=float(c2))
        lin = 0.5 * (ya - yb); quad = 0.5 * (ya + yb); Sl = 0.25 * (Sa + Sb - X - X.T); Sq = 0.25 * (Sa + Sb + X + X.T)
        res.update(lin=lin.tolist(), quad=quad.tolist(), lin_sigma=np.sqrt(np.diag(Sl)).tolist(), quad_sigma=np.sqrt(np.diag(Sq)).tolist(), rms_lin=float(np.sqrt(np.mean(lin ** 2))), rms_quad=float(np.sqrt(np.mean(quad ** 2))),
                   rms_sigma_lin=float(np.sqrt(np.mean(np.diag(Sl)))), rms_sigma_quad=float(np.sqrt(np.mean(np.diag(Sq)))), rms_sigma_single=float(np.sqrt(np.mean(0.5 * (np.diag(Sa) + np.diag(Sb))))),
                   chi2_quad_zero=chi2(quad, Sq), p_quad_zero=float(stats.chi2.sf(chi2(quad, Sq), NB)), chi2_lin_zero=chi2(lin, Sl))
    return res
tests = []
for v in ICAR: tests.append(pair_test(f"icarus_run2__{v}", f"icarus_run4__{v}", f"ICARUS Run2 vs Run4: {v}", "indep"))
tests.append(pair_test("sbnd_run1__ChargeScale", "sbnd_run1__ChargeScale15", "SBND ChargeScale vs ChargeScale15", "equal"))
tests.append(pair_test("sbnd_run1__WireModYZ", "sbnd_run1__MattModYZ", "SBND WireModYZ vs MattModYZ", "equal"))
for det in ("sbnd_run1", "icarus_run2", "icarus_run4"): tests.append(pair_test(f"{det}__SCE0x", f"{det}__SCE2x", f"{det} SCE0x + SCE2x (antisymmetry)", "antisym"))
tests.append(pair_test("sbnd_run1__GainHigh", "sbnd_run1__GainLow", "SBND GainHigh + GainLow (antisymmetry)", "antisym"))
for t in tests:
    s = f"{t['label']:52s} "
    if "chi2_equal" in t: s += f"equal chi2/ndf {t['chi2_equal']/NB:7.2f} (p {t['p_equal']:.2g}) z-mean {t['z_mean']:+.2f} z-sd {t['z_sd']:.2f} rho_noise {t.get('rho_noise',0):+.2f} | diff smooth-deg: chi2 {t['diff_deg0_chi2']:.1f}/{t['diff_deg1_chi2']:.1f}/{t['diff_deg2_chi2']:.1f}"
    if "alpha_prop" in t and "chi2_prop" in t: s += f" | prop alpha {t['alpha_prop']:.3f} chi2/16 {t['chi2_prop']/16:.2f} (p {t['p_prop']:.2g})"
    if "chi2_same" in t: s += f" same-sign prop alpha {t['alpha_same']:.2f} chi2/16 {t['chi2_same']/16:.2f} |"
    if "chi2_antisym" in t: s += f"antisym chi2/ndf {t['chi2_antisym']/NB:6.2f} (p {t['p_antisym']:.2g}); rms lin {100*t['rms_lin']:.2f}% quad {100*t['rms_quad']:.2f}% sigma_lin {100*t['rms_sigma_lin']:.2f}% sigma_quad {100*t['rms_sigma_quad']:.2f}% single {100*t['rms_sigma_single']:.2f}%; quad=0: chi2 {t['chi2_quad_zero']:.1f} p {t['p_quad_zero']:.2g}"
    print(s)
ic = [t for t in tests if t["label"].startswith("ICARUS Run2")]; print("ICARUS Run2 vs Run4 pooled chi2/ndf", sum(t["chi2_equal"] for t in ic), 5 * NB, float(stats.chi2.sf(sum(t["chi2_equal"] for t in ic), 5 * NB)), "pooled z sd", float(np.std(np.concatenate([t["z"] for t in ic]))))
# noise correlation across SBND variations (joint hashed bootstrap): mean over bins of corr(y_a, y_b) over replicates
sb = [n for n in NAMES if n.startswith("sbnd")]; Rb = {n: (lambda z: z["V"] / z["N"])(np.load(RES / f"boot_{n}.npz")) for n in sb}
M = np.zeros((len(sb), len(sb)))
for i, a in enumerate(sb):
    for j, b in enumerate(sb): M[i, j] = np.mean([np.corrcoef(Rb[a][:, k], Rb[b][:, k])[0, 1] for k in range(0, NB)])
print("SBND noise correlation between variations (mean over bins), off-diagonal range:", np.round(M[np.triu_indices(len(sb), 1)].min(), 2), np.round(M[np.triu_indices(len(sb), 1)].max(), 2), "median", np.round(np.median(M[np.triu_indices(len(sb), 1)]), 2))
# check analytic cross-cov vs joint bootstrap for ChargeScale/ChargeScale15
a, b = "sbnd_run1__ChargeScale", "sbnd_run1__ChargeScale15"; X = cross_cov(D[a], D[b]); xb = np.array([np.cov(Rb[a][:, k], Rb[b][:, k])[0, 1] for k in range(NB)])
print("cross-cov diag analytic vs bootstrap (ratio):", np.round(np.diag(X) / xb, 2))
json.dump(dict(icarus_pooled=dict(chi2=float(sum(t["chi2_equal"] for t in ic)), ndf=5 * NB, p=float(stats.chi2.sf(sum(t["chi2_equal"] for t in ic), 5 * NB)), z_sd=float(np.std(np.concatenate([t["z"] for t in ic])))), tests=tests, sbnd_noise_corr=dict(names=sb, M=M.tolist()), crosscov_check=dict(analytic=np.diag(X).tolist(), bootstrap=xb.tolist())), open(Q2 / "q2_consistency.json", "w"))
