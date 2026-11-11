#!/usr/bin/env python3
"""Q3: 'direct' DetVar response (selected-attribute construction, observable-specific binning) versus the production hist1d treatment (per-event weight = ratio(E_vis) applied to the
CV events, then histogrammed in the observable) for several observables. All numpy on the pair tables. cvwgt not used."""
import sys, json, numpy as np, pandas as pd, uproot
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
Q3 = OUT / "q3"; rng = np.random.default_rng(7); Bt = 100
prodf = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root")
OBS = {"abs_muon_endx": ("reco_leading_muon_end_x", np.arange(0, 200.1, 10.0), np.abs, ("sbnd_run1",)),
       "muon_p": ("reco_leading_muon_p", np.arange(200, 1801, 100.0), None, None),
       "muon_softmax": ("reco_leading_muon_muon_softmax", np.array([0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0001]), None, None),
       "muon_length": ("reco_leading_muon_length", np.array([50, 75, 100, 150, 200, 300, 400, 600.0]), None, None)}
def bi(x, ed):
    i = np.digitize(x, ed) - 1; i[(i < 0) | (i >= len(ed) - 1) | ~np.isfinite(x)] = -1; return i
def rat(a, b): return np.where(b > 0, a / np.where(b > 0, b, 1), 1.0)
out = {}
for det, cfgs in configs().items():
    for lab, hname, *_ in cfgs:
        name = f"{det}__{lab}"; P = pd.read_pickle(PAIRS / f"{name}.pkl"); inn = P.in_n.to_numpy(); inv = P.in_v.to_numpy(); prod = prodf[hname].values()
        rE = np.where(inn, np.take(np.r_[prod, 1.0], np.where(inn, np.clip(np.digitize(np.nan_to_num(P.n_reco_visible_energy.to_numpy()), BINS) - 1, -1, len(prod)), -1) * 1), 1.0)
        ie = np.digitize(np.nan_to_num(P.n_reco_visible_energy.to_numpy()), BINS) - 1; ok = (ie >= 0) & (ie < len(prod)); w1d = np.where(inn & ok, prod[np.clip(ie, 0, len(prod) - 1)], np.where(inn, 1.0, 0.0))
        for on, (col, ed, fn, dets) in OBS.items():
            if dets and det not in dets: continue
            xn = P["n_" + col].to_numpy(float); xv = P["v_" + col].to_numpy(float)
            if fn is not None: xn, xv = fn(xn), fn(xv)
            bn = np.where(inn, bi(xn, ed), -1); bv = np.where(inv, bi(xv, ed), -1); nb = len(ed) - 1
            def H(b, w): m = b >= 0; return np.bincount(b[m], weights=w[m], minlength=nb)
            def both(wp):
                Hn = H(bn, wp); Hv = H(bv, wp); H1 = H(bn, wp * w1d); return rat(Hv, Hn), rat(H1, Hn), Hn
            direct, h1d, Hn = both(np.ones(len(P)))
            reps = [both(rng.poisson(1.0, len(P)).astype(float)) for _ in range(Bt)]; sd = np.std([r[0] for r in reps], 0, ddof=1); sh = np.std([r[1] for r in reps], 0, ddof=1); sdiff = np.std([r[0] - r[1] for r in reps], 0, ddof=1)
            d = direct - h1d; good = Hn >= 50
            out[f"{name}|{on}"] = dict(edges=ed.tolist(), direct=direct.tolist(), hist1d=h1d.tolist(), Hn=Hn.tolist(), sig_direct=sd.tolist(), sig_diff=sdiff.tolist(),
                                       max_abs=float(np.abs(d[good]).max()), rms=float(np.sqrt((d[good] ** 2).mean())), chi2=float(((d[good] / np.maximum(sdiff[good], 1e-9)) ** 2).sum()), ndf=int(good.sum()),
                                       max_resp=float(np.abs(direct[good] - 1).max()))
        print(name, {k.split("|")[1]: (round(v["max_abs"], 4), round(v["max_resp"], 3), round(v["chi2"], 1), v["ndf"]) for k, v in out.items() if k.startswith(name + "|")}, flush=True)
json.dump(out, open(Q3 / "observables.json", "w"))
