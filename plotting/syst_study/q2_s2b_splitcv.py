#!/usr/bin/env python3
"""Q2 deliverable 2 (out-of-sample): fit the smoothing algorithm on random event half A, score against the RAW ratio of the complementary half B (and v.v.).
T = sum_b (yB - fA)^2 / sigma_full,b^2 ; floor = sum_b sigma_B,b^2/sigma_full,b^2 (= irreducible noise of B ~ 2 per bin). Excess over the floor measures bias^2 + estimator variance (+ the noise of half A for the raw estimator)."""
from q2_lib import *
from q2_smooth import *
NSP = 100
METH = {"raw": None, "zero": "zero", "poly-BIC, last raw": dict(mode="poly", last="raw", dmax=4), "whit-CV, last raw": dict(mode="whit", last="raw"), "poly-BIC, last pretest(2sig)": dict(mode="poly", last="pretest", tau2=0.0095 ** 2, dmax=4), "poly-BIC, all bins": dict(mode="poly", last="smooth", dmax=4), "whit-CV, all bins": dict(mode="whit", last="smooth")}
out = {}
for n in ALL:
    d = load(n); r, C, N, V = cov_analytic(d); w = 1 / np.diag(C); acc = {m: [] for m in METH}; acc_last = {m: [] for m in METH}; floor = []; floor_last = []
    for s in range(NSP):
        mk = (uni(d["key"], 9100 + s) < 0.5).astype(float)
        for wa, wb in ((mk, 1 - mk), (1 - mk, mk)):
            rA, CA, _, _ = cov_analytic(d, w=wa); rB, CB, _, _ = cov_analytic(d, w=wb); yA, yB = rA - 1, rB - 1; floor.append(np.sum(w[:-1] * np.diag(CB)[:-1])); floor_last.append(w[-1] * CB[-1, -1])
            for m, kw in METH.items():
                f = yA if kw is None else (np.zeros(NB) if isinstance(kw, str) else Estimator(CA, LOGX, **kw)(yA[None, :])[0][0]); e = (yB - f) ** 2 * w; acc[m].append(e[:-1].sum()); acc_last[m].append(e[-1])
    fl = np.mean(floor); fll = np.mean(floor_last)
    out[n] = dict(floor_bulk=float(fl), floor_last=float(fll), bulk={m: float(np.mean(v)) for m, v in acc.items()}, last={m: float(np.mean(v)) for m, v in acc_last.items()})
    b = out[n]["bulk"]; print(f"{n:34s} bulk T/floor: " + " ".join(f"{m.split(',')[0][:9]}{'/L' if 'last raw' in m else ''}{'/A' if 'all' in m else ''}={b[m]/fl:5.2f}" for m in METH), flush=True)
json.dump(out, open(Q2 / "q2_splitcv.json", "w"))
