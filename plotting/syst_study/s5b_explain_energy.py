#!/usr/bin/env python3
"""S5b: do PID scores + proton multiplicity (no energy, no kinematics) explain the energy dependence of the lost / gained fractions?  Cross-validated GBDT; predicted vs observed fraction per energy bin."""
import json, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from ana import *
from s5_drivers import side  # noqa (re-runs nothing heavy: guarded below)
FEATS_NOE = ["mult", "mu_mu", "mu_pi", "mu_p", "mu_prim", "p_p", "p_mu", "p_pi", "p_prim", "m_mu", "m_p", "mu_calo", "p_calo"]
FEATS_KIN = FEATS_NOE + ["mu_ke", "p_ke", "mu_len", "p_len"]
res = {}
for name in ("sbnd_run1__ChargeScale15", "sbnd_run1__ChargeScale", "icarus_run2__ChargeScale", "sbnd_run1__GainHigh", "sbnd_run1__MattModYZ", "sbnd_run1__WireModYZ"):
    df = load(name); inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); Xn = side(df, "n"); Xv = side(df, "v"); out = {}
    for tag, mask, y, X in (("lost", inn, (inn & ~inv)[inn], {k: v[inn] for k, v in Xn.items()}), ("gained", inv, (~inn & inv)[inv], {k: v[inv] for k, v in Xv.items()})):
        o = {}
        for fs_name, cols in (("PID+mult", FEATS_NOE), ("PID+mult+kinematics", FEATS_KIN)):
            M = np.c_[[np.nan_to_num(X[c], nan=-1) for c in cols]].T; rng = np.random.default_rng(1); idx = rng.choice(len(y), min(250000, len(y)), replace=False)
            clf = HistGradientBoostingClassifier(max_iter=120, learning_rate=0.1, max_leaf_nodes=24, random_state=0)
            p = cross_val_predict(clf, M[idx], y[idx].astype(int), cv=StratifiedKFold(4, shuffle=True, random_state=1), method="predict_proba")[:, 1]
            E = X["E"][idx]; yy = y[idx]
            o[fs_name] = [(float(p[(E >= lo) & (E < hi)].mean()), float(yy[(E >= lo) & (E < hi)].mean()), int(((E >= lo) & (E < hi)).sum())) for lo, hi in zip([0.3, 0.5, 0.7, 0.9, 1.1, 1.5], [0.5, 0.7, 0.9, 1.1, 1.5, 2.0])]
        out[tag] = o
    res[name] = out
    for tag in ("lost", "gained"):
        for fs_name, v in out[tag].items(): print(f"{name:28s} {tag:6s} {fs_name:20s} predicted%: {[round(a*100,1) for a,_,_ in v]}  observed%: {[round(b*100,1) for _,b,_ in v]}", flush=True)
json.dump(res, open(FIG / "s5b_explain_energy.json", "w"))
