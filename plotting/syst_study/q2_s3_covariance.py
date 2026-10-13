#!/usr/bin/env python3
"""Q2 deliverable 3: noise contribution to the fractional covariance C = (r-1)(r-1)^T that PROfit builds from the ratio (per variation and summed per detector).
E[C_raw] = f f^T + Sigma_noise (cross terms vanish in expectation), so the noise adds exactly Sigma_noise (diag. dominated) to the trace. Smooth fit = recommended algorithm (poly-BIC, dmax 4, zero allowed, bins 0-15; last bin raw)."""
from q2_lib import *
from q2_fit import *
from q2_smooth import *
def rec_fit(y, C, last="pretest", mode="poly"):
    est = Estimator(C, LOGX, mode=mode, last=last, dmax=4, tau2=0.0095 ** 2); F, sel = est(y[None, :]); return F[0], sel[0]
def eig_stats(M):
    ev, U = np.linalg.eigh(M); ev = ev[::-1]; U = U[:, ::-1]; ev = np.clip(ev, 0, None); p = ev / ev.sum()
    return dict(ev=ev.tolist(), eff_rank=float(ev.sum() ** 2 / (ev ** 2).sum()), n90=int(np.searchsorted(np.cumsum(p), 0.9) + 1), frac1=float(p[0])), U[:, 0]
out = {}; per = {}; boots = {}
for n in ALL:
    d = load(n); r, C, N, V = cov_analytic(d); y = r - 1; s2 = np.diag(C); f, sel = rec_fit(y, C); fw, _ = rec_fit(y, C, mode="whit"); flr, _ = rec_fit(y, C, last="raw")
    o = dict(trace_raw=float((y ** 2).sum()), trace_noise=float(s2.sum()), trace_fit=float((f ** 2).sum()), trace_fit_lastraw=float((flr ** 2).sum()), trace_whit=float((fw ** 2).sum()), noise_frac=float(s2.sum() / (y ** 2).sum()),
             noise_frac_bulk=float(s2[:-1].sum() / (y[:-1] ** 2).sum()), last_bin_share_of_noise=float(s2[-1] / s2.sum()), last2_share_of_noise=float(s2[-2:].sum() / s2.sum()), last_bin_share_of_raw_trace=float(y[-1] ** 2 / (y ** 2).sum()),
             trace_debiased=float((y ** 2).sum() - s2.sum()), deg_sel=int(sel))
    # sign pattern agreement of off-diagonal entries raw vs fit (only where fit is non-zero)
    iu = np.triu_indices(NB - 1, 1); sr = np.sign(np.outer(y[:-1], y[:-1])[iu]); sf = np.sign(np.outer(f[:-1], f[:-1])[iu]); nz = sf != 0   # bulk bins only (last bin is kept raw)
    o["sign_flip_frac"] = float(np.mean(sr[nz] != sf[nz])) if nz.any() else float("nan"); o["frob_noise_frac"] = float(np.linalg.norm(np.outer(y, y) - np.outer(f, f)) / np.linalg.norm(np.outer(y, y)))
    o["rms_y"] = float(np.sqrt(np.mean(y ** 2))); o["rms_fit"] = float(np.sqrt(np.mean(f ** 2))); o["rms_noise"] = float(np.sqrt(np.mean(s2)))
    per[n] = dict(y=y, f=f, fw=fw, C=C, flr=flr); out[n] = o
# detector level
R4W = "icarus_run4__WireModxThetaXWclean"
det_groups = {"sbnd_run1": [n for n in NAMES if n.startswith("sbnd")], "icarus_run2": [n for n in NAMES if "run2" in n], "icarus_run4": [n for n in ALL if "run4" in n and n != "icarus_run4__WireModxThetaXW"]}
ACTIVE = ("ChargeScale", "SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ", "DENT")   # variations entering xml/standard.xml (SBND: no SCE; SCE only for ICARUS)
for det in list(det_groups): det_groups[det + "|xml-active"] = [n for n in det_groups[det] if var_of(n).replace("clean", "") in ACTIVE and not (det.startswith("sbnd") and var_of(n).startswith("SCE"))]
rng = np.random.default_rng(3); det_out = {}
for det, ns in det_groups.items():
    Craw = sum(np.outer(per[n]["y"], per[n]["y"]) for n in ns); Cfit = sum(np.outer(per[n]["f"], per[n]["f"]) for n in ns); Cwh = sum(np.outer(per[n]["fw"], per[n]["fw"]) for n in ns); Clr = sum(np.outer(per[n]["flr"], per[n]["flr"]) for n in ns); elr, _ = eig_stats(Clr); Cn = sum(per[n]["C"] for n in ns)
    er, e1r = eig_stats(Craw); ef, e1f = eig_stats(Cfit); ew, e1w = eig_stats(Cwh); en, _ = eig_stats(Cn)
    corr = lambda M: M / np.sqrt(np.outer(np.diag(M), np.diag(M)))
    adj = lambda M: float(np.mean([corr(M)[i, i + 1] for i in range(NB - 1)])); off = lambda M: float(np.mean(np.abs(corr(M)[np.triu_indices(NB, 1)])))
    # parametric/joint-bootstrap MC: truth = fit f_v, noise = bootstrap deviations of the SAME events (joint across variations)
    K = 1500; tr = np.zeros(K); ov = np.zeros(K); rk = np.zeros(K); fr1 = np.zeros(K); trN = np.zeros(K)
    eps = {}
    for n in ns:
        bz = np.load(RES / f"boot_{n}.npz"); R = (bz["V"] / bz["N"])[:K]; eps[n] = R - R.mean(0)
    for k in range(K):
        Ck = sum(np.outer(per[n]["f"] + eps[n][k], per[n]["f"] + eps[n][k]) for n in ns); ev = np.linalg.eigvalsh(Ck)[::-1]; w, U = np.linalg.eigh(Ck); tr[k] = np.trace(Ck); ov[k] = abs(U[:, -1] @ e1f); rk[k] = ev.sum() ** 2 / (ev ** 2).sum(); fr1[k] = ev[0] / ev.sum()
    dom = max(ns, key=lambda n: out[n]["trace_fit"])
    det_out[det] = dict(variations=ns, trace_raw=float(np.trace(Craw)), trace_fit=float(np.trace(Cfit)), trace_noise=float(np.trace(Cn)), noise_frac=float(np.trace(Cn) / np.trace(Craw)), eig_raw=er, eig_fit=ef, eig_whit=ew, eig_fit_lastraw=elr, trace_fit_lastraw=float(np.trace(Clr)), adj_corr_fit_lastraw=adj(Clr), eig_noise=en,
                        e1_overlap_raw_fit=float(abs(e1r @ e1f)), adj_corr_raw=adj(Craw), adj_corr_fit=adj(Cfit), meanabs_corr_raw=off(Craw), meanabs_corr_fit=off(Cfit), dominant=dom, dom_trace_fit=out[dom]["trace_fit"],
                        noise_over_dom=float(np.trace(Cn) / out[dom]["trace_fit"]), mc=dict(trace_mean=float(tr.mean()), trace_sd=float(tr.std()), trace_fit=float(np.trace(Cfit)), e1_overlap_mean=float(ov.mean()), e1_overlap_p16=float(np.percentile(ov, 16)), e1_overlap_p84=float(np.percentile(ov, 84)),
                        eff_rank_mean=float(rk.mean()), eff_rank_p16=float(np.percentile(rk, 16)), eff_rank_p84=float(np.percentile(rk, 84)), frac1_mean=float(fr1.mean()), eff_rank_true=ef["eff_rank"], frac1_true=ef["frac1"]),
                        e1_raw=e1r.tolist(), e1_fit=e1f.tolist())
    d_ = det_out[det]; print(f"{det}: trace raw {d_['trace_raw']*1e4:.2f}e-4 fit {d_['trace_fit']*1e4:.2f}e-4 noise {d_['trace_noise']*1e4:.2f}e-4 (frac {d_['noise_frac']:.3f}) | eff rank raw {er['eff_rank']:.2f} fit {ef['eff_rank']:.2f} MC {d_['mc']['eff_rank_mean']:.2f} | frac1 raw {er['frac1']:.3f} fit {ef['frac1']:.3f} noise-only {en['frac1']:.3f} | e1 overlap raw-fit {d_['e1_overlap_raw_fit']:.3f} MC {d_['mc']['e1_overlap_mean']:.3f} | adj corr raw {d_['adj_corr_raw']:.2f} fit {d_['adj_corr_fit']:.2f} | dominant {dom} noise/dom {d_['noise_over_dom']:.3f}")
for n, o in out.items(): print(f"{n:30s} trace raw {o['trace_raw']*1e4:8.2f} noise {o['trace_noise']*1e4:6.2f} fit {o['trace_fit']*1e4:8.2f} (e-4) noise/raw {o['noise_frac']:.2f} bulk {o['noise_frac_bulk']:.2f} last-bin share of noise {o['last_bin_share_of_noise']:.2f} last2 {o['last2_share_of_noise']:.2f} sign-flip {o['sign_flip_frac']:.2f} frob {o['frob_noise_frac']:.2f} deg {o['deg_sel']}")
json.dump(dict(per_variation=out, per_detector=det_out), open(Q2 / "q2_covariance.json", "w"))
