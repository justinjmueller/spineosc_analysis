#!/usr/bin/env python3
"""Q3: what the direct method gains when the analysis selection is stricter than the calculator selection: ratio in reco E_vis restricted to events with reco_proton_multiplicity >=1 (or ==1, >=2),
evaluated PER SAMPLE (so migration across the extra cut is included), versus the inclusive production ratio that the hist1d systematic applies to every subchannel/selection. cvwgt not used."""
import sys, json, numpy as np, pandas as pd, uproot
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
Q3 = OUT / "q3"; prodf = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root"); rng = np.random.default_rng(3); B = 150; NB = len(BINS) - 1; out = {}
def rat(a, b): return np.where(b > 0, a / np.where(b > 0, b, 1), 1.0)
for name, hname in (("sbnd_run1__ChargeScale15", "ChargeScale15_sbnd_Run1"), ("sbnd_run1__DENT", "DENT_sbnd_Run1"), ("sbnd_run1__GainHigh", "GainHigh_sbnd_Run1"), ("sbnd_run1__SCE2x", "SCE2x_sbnd_Run1"),
                    ("sbnd_run1__WireModYZ", "WireModYZ_sbnd_Run1"), ("icarus_run2__SCE2x", "SCE2x_icarus_Run2"), ("icarus_run2__ChargeScale", "ChargeScale_icarus_Run2")):
    P = pd.read_pickle(PAIRS / f"{name}.pkl"); prod = prodf[hname].values()
    for sel, fn in (("pm>=1", lambda m: m >= 1), ("pm==1", lambda m: m == 1), ("pm>=2", lambda m: m >= 2)):
        inn = P.in_n.to_numpy() & fn(P.n_reco_proton_multiplicity.fillna(-1).to_numpy()); inv = P.in_v.to_numpy() & fn(P.v_reco_proton_multiplicity.fillna(-1).to_numpy())
        bn = np.where(inn, np.clip(np.digitize(np.nan_to_num(P.n_reco_visible_energy.to_numpy()), BINS) - 1, -1, NB), -1); bv = np.where(inv, np.clip(np.digitize(np.nan_to_num(P.v_reco_visible_energy.to_numpy()), BINS) - 1, -1, NB), -1)
        bn[(bn < 0) | (bn >= NB)] = -1; bv[(bv < 0) | (bv >= NB)] = -1
        H = lambda b, w: np.bincount(b[b >= 0], weights=w[b >= 0], minlength=NB)
        f = lambda w: rat(H(bv, w), H(bn, w)); r = f(np.ones(len(P))); sg = np.std([f(rng.poisson(1.0, len(P)).astype(float)) for _ in range(B)], 0, ddof=1); Hn = H(bn, np.ones(len(P)))
        good = Hn >= 30; d = (r - prod)[good]
        out[f"{name}|{sel}"] = dict(ratio=r.tolist(), prod=prod.tolist(), sig=sg.tolist(), Hn=Hn.tolist(), max_abs=float(np.abs(d).max()), rms=float(np.sqrt((d ** 2).mean())), chi2=float(((d / np.maximum(sg[good], 1e-9)) ** 2).sum()), ndf=int(good.sum()))
    print(name, {k.split("|")[1]: (round(v["max_abs"], 4), round(v["chi2"], 1), v["ndf"]) for k, v in out.items() if k.startswith(name + "|")}, flush=True)
json.dump(out, open(Q3 / "selection.json", "w"))
