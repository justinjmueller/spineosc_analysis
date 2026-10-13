#!/usr/bin/env python3
"""Q2 deliverable 6 (validation): closure/simulation study of candidate smoothing algorithms. Noise draws = bootstrap deviations (R_k - mean R) of the real pair tables (2000 joint replicates),
so non-Gaussianity and bin correlations are those of the data. Truth scenarios: smooth (data-driven fit), smooth + last-bin edge offset, zero."""
from q2_lib import *
from q2_fit import *
from q2_smooth import *
import time
sn = json.load(open(Q2 / "q2_signal_noise.json"))
# --- empirical-Bayes estimate of the last-bin excess variance tau^2 (relative to a deg-2 extrapolation from bins 0..15)
E_, V_ = [], []
for n in NAMES:
    d = load(n); r, C, N, V = cov_analytic(d); y = r - 1; pred, vp, _ = edge_extrap(y[None, :], C, LOGX, 2); E_.append(y[-1] - pred[0]); V_.append(C[-1, -1] + vp)
E_ = np.array(E_); V_ = np.array(V_)
taus = np.linspace(0, 0.06, 601); ll = [np.sum(-0.5 * (np.log(V_ + t ** 2) + E_ ** 2 / (V_ + t ** 2))) for t in taus]; tau_hat = float(taus[int(np.argmax(ll))])
# profile-likelihood 68% interval
ll = np.array(ll); inside = taus[ll > ll.max() - 0.5]; tau_int = (float(inside.min()), float(inside.max()))
print("tau(last-bin excess) =", tau_hat, tau_int, " rms pull of deg2 extrapolation:", np.sqrt(np.mean(E_ ** 2 / V_)), flush=True)
TAU = max(tau_hat, tau_int[1]); TAU2 = TAU ** 2   # conservative: the 68% upper limit of the (zero-consistent) ML value
METHODS = {"raw": None, "poly-BIC": dict(mode="poly", last="smooth", dmax=4), "whit-CV": dict(mode="whit", last="smooth"), "poly-BIC, last raw": dict(mode="poly", last="raw", dmax=4), "whit-CV, last raw": dict(mode="whit", last="raw"),
           "poly-BIC, last shrunk": dict(mode="poly", last="shrink", tau2=TAU2, dmax=4), "whit-CV, last shrunk": dict(mode="whit", last="shrink", tau2=TAU2), "poly-BIC, last inflated": dict(mode="poly", last="inflate", tau2=TAU2, dmax=4), "whit-CV, last inflated": dict(mode="whit", last="inflate", tau2=TAU2), "poly-BIC, last pretest(2sig)": dict(mode="poly", last="pretest", tau2=TAU2, dmax=4), "whit-CV, last pretest(2sig)": dict(mode="whit", last="pretest", tau2=TAU2)}
res = {}
for n in NAMES:
    t0 = time.time(); bz = np.load(RES / f"boot_{n}.npz"); R = bz["V"] / bz["N"]; eps = R - R.mean(0); d = load(n); r, C, N, V = cov_analytic(d); sig = np.sqrt(np.diag(C)); y = r - 1
    # truths
    ft = np.array(sn[n]["fit_bic"]); fw = np.array(sn[n]["whit"]["f"])
    ft_cap = ft.copy()   # BIC poly fit may use deg up to 6 in S2; ok as a smooth truth
    def ed(k): e = np.zeros(NB); e[-1] = k * sig[-1]; return e
    truths = {"smooth(poly)": ft, "smooth(whit)": fw, "smooth+edge(+1sig)": ft + ed(1), "smooth+edge(-1sig)": ft - ed(1), "smooth+edge(+2sig)": ft + ed(2), "smooth+edge(-2sig)": ft - ed(2), "zero": np.zeros(NB), "zero+edge(2sig)": ed(2), "const(+0.5 rms-sigma)": 0.5 * np.sqrt(np.mean(sig ** 2)) * np.ones(NB), "const(+1 rms-sigma)": np.sqrt(np.mean(sig ** 2)) * np.ones(NB)}
    ests = {m: (Estimator(C, LOGX, **kw) if kw else None) for m, kw in METHODS.items()}; res[n] = {}
    for tn, tr in truths.items():
        Y = tr[None, :] + eps; res[n][tn] = {}
        for m, est in ests.items():
            F = Y if est is None else est(Y)[0]; err = F - tr; res[n][tn][m] = dict(mse_ratio=float((err ** 2).sum(1).mean() / (sig ** 2).sum()), last_mse_ratio=float((err[:, -1] ** 2).mean() / sig[-1] ** 2), bias_last_sig=float(err[:, -1].mean() / sig[-1]),
                 bias_rms_sig=float(np.sqrt(np.mean((err.mean(0) / sig) ** 2))), trace_ratio=float((F ** 2).sum(1).mean() / max((tr ** 2).sum(), 1e-12)), trace_over_noise=float((F ** 2).sum(1).mean() / (sig ** 2).sum()),
                 covbias=float(((F ** 2).sum(1).mean() - (tr ** 2).sum()) / (sig ** 2).sum()))
    print(n, round(time.time() - t0, 1), flush=True)
json.dump(dict(tau_hat=tau_hat, tau_interval=tau_int, tau_rms_pull=float(np.sqrt(np.mean(E_ ** 2 / V_))), edge_pulls=(E_ / np.sqrt(V_)).tolist(), res=res), open(Q2 / "q2_validation.json", "w"))
# aggregate
print("\nAGGREGATE (mean over the 20 variations): mse_ratio (raw=1) / last-bin mse_ratio / |bias|last [sigma] / covbias")
for tn in list(truths):
    print(tn)
    for m in [k for k in METHODS if 'shrunk' not in k and 'inflated' not in k]: a = [res[n][tn][m] for n in NAMES]; print(f"   {m:34s} {np.mean([x['mse_ratio'] for x in a]):.3f}  last {np.mean([x['last_mse_ratio'] for x in a]):.3f}  bias_last {np.mean([x['bias_last_sig'] for x in a]):+.2f}  rms-bias {np.mean([x['bias_rms_sig'] for x in a]):.2f} covbias/noise-trace {np.mean([x['covbias'] for x in a]):+.2f}")
