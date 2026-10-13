#!/usr/bin/env python3
"""S5 (Q1/Q4): what variable carries the energy dependence?  (a) ratio versus alternative variables (E_vis, muon KE, proton KE, muon/proton softmax, margins, proton multiplicity);
(b) cross-validated gradient-boosted models for 'lost' (selected in nominal only) and 'gained' (selected in the variation only) using nominal-side / varied-side features: AUC with feature sets
{E only, +multiplicity, PID scores, PID+kinematics, all}, and whether the model's predicted lost/gained fraction reproduces the observed energy dependence."""
import json, sys, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from ana import *
def side(df, s):
    g = lambda c: df[f"{s}_reco_leading_{c}"].to_numpy(float)
    X = dict(E=df[f"{s}_reco_visible_energy"].to_numpy(float), mult=df[f"{s}_reco_proton_multiplicity"].to_numpy(float), mu_ke=g("muon_ke"), p_ke=g("proton_ke"), mu_len=g("muon_length"), p_len=g("proton_length"),
             mu_mu=g("muon_muon_softmax"), mu_pi=g("muon_pion_softmax"), mu_p=g("muon_proton_softmax"), mu_prim=g("muon_primary_softmax"), p_p=g("proton_proton_softmax"), p_mu=g("proton_muon_softmax"), p_pi=g("proton_pion_softmax"),
             p_prim=g("proton_primary_softmax"), mu_calo=g("muon_calo_csda_diff"), p_calo=g("proton_calo_csda_diff"), cosz=g("muon_start_dir_z"), opening=df[f"{s}_reco_bivar_muon_proton_opening_angle"].to_numpy(float))
    X["m_mu"] = 2 * X["mu_mu"] - np.maximum(X["mu_pi"], X["mu_p"]); X["m_p"] = X["p_p"] - np.maximum(2 * X["p_mu"], X["p_pi"]); return X
SETS = {"E only": ["E"], "E+mult": ["E", "mult"], "PID scores": ["mu_mu", "mu_pi", "mu_p", "mu_prim", "p_p", "p_mu", "p_pi", "p_prim", "m_mu", "m_p"],
        "PID+calo": ["mu_mu", "mu_pi", "mu_p", "mu_prim", "p_p", "p_mu", "p_pi", "p_prim", "m_mu", "m_p", "mu_calo", "p_calo"],
        "kinematics": ["E", "mult", "mu_ke", "p_ke", "mu_len", "p_len", "cosz", "opening"], "all": None}
def auc_cv(X, y, cols, seed=0, nmax=250000):
    rng = np.random.default_rng(seed); idx = rng.choice(len(y), min(nmax, len(y)), replace=False); M = np.c_[[X[c][idx] for c in cols]].T; yy = y[idx]
    clf = HistGradientBoostingClassifier(max_iter=120, learning_rate=0.1, max_leaf_nodes=24, random_state=0)
    p = cross_val_predict(clf, M, yy, cv=StratifiedKFold(4, shuffle=True, random_state=seed), method="predict_proba")[:, 1]
    return roc_auc_score(yy, p), p, idx
if __name__ == '__main__':
    res = {}
    for name in pairs_available():
        df = load(name); inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); r = {}
        # (a) ratio versus alternative variables, paired bootstrap-free quick version (counts)
        Xn = side(df, "n"); Xv = side(df, "v"); ra = {}
        for var, edges in (("E", BINS), ("mu_ke", np.array([143, 200, 300, 400, 500, 650, 800, 1000, 2500])), ("p_ke", np.array([50, 70, 100, 130, 170, 220, 300, 450, 2000])),
                           ("mu_mu", np.array([0, .2, .4, .6, .8, .9, .95, .99, 1.0001])), ("p_p", np.array([0, .5, .8, .9, .95, .98, .995, 1.0001])), ("mult", np.array([1, 2, 3, 4, 8]))):
            a = np.histogram(Xn[var][inn & np.isfinite(Xn[var])], edges)[0].astype(float); b = np.histogram(Xv[var][inv & np.isfinite(Xv[var])], edges)[0].astype(float); ra[var] = dict(edges=edges.tolist(), ratio=(b / a).tolist(), n=a.tolist())
        r["ratio_vs"] = ra
        # (b) models: lost = selected in nominal not in varied; gained = in varied not in nominal
        out = {}
        for tag, mask, y, X in (("lost", inn, (inn & ~inv)[inn], {k: v[inn] for k, v in Xn.items()}), ("gained", inv, (~inn & inv)[inv], {k: v[inv] for k, v in Xv.items()})):
            o = {}
            for sname, cols in SETS.items():
                cc = cols or list(X)
                A, p, idx = auc_cv(X, y.astype(int), [c for c in cc if np.isfinite(X[c]).mean() > 0.9]); o[sname] = float(A)
                if sname == "all":
                    E = X["E"][idx]; yy = y[idx]; pe = np.array([(p[(E >= lo) & (E < hi)].mean(), yy[(E >= lo) & (E < hi)].mean(), ((E >= lo) & (E < hi)).sum()) for lo, hi in zip(BINS[:-1], BINS[1:])])
                    o["calibration_vs_E"] = pe.tolist()
            out[tag] = o
        r["models"] = out; res[name] = r
        print(f"{name:32s} lost AUC: " + " ".join(f"{k}={v:.3f}" for k, v in out["lost"].items() if k != "calibration_vs_E") + " | gained AUC: " + " ".join(f"{k}={v:.3f}" for k, v in out["gained"].items() if k != "calibration_vs_E"), flush=True)
    json.dump(res, open(FIG / "s5_drivers.json", "w"))
