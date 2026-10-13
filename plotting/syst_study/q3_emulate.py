#!/usr/bin/env python3
"""Q3 numpy emulation of the PROfit DetVar response for every detector/variation from the paired event tables.
Variants (all in the production by_vis_energy binning; cvwgt never used):
  prod        production ratio histogram (detsys-v1.1.0v2.root)
  pair        ratio of selected counts on the unique-key pool rows (= production minus the repeated-key rows)
  A_unmatched_pool  independent CV / variation selected samples RESTRICTED to the common pool (unmatched mode, pool-restricted trees) -> identical to 'pair'
  A_unmatched_full  independent whole selected trees (no pool restriction), raw counts, and POT-scaled
  B_matched   events selected in BOTH samples only (PROfit cv_variation_matching_vars, zero-weight events dropped)
  C_selattr   union rows on both sides with selected as weight (flag construction) -> identical to 'pair' by algebra (a check of the algebra, run natively elsewhere)
Errors: Poisson (per-event weights ~Poisson(1)) paired bootstrap over the union rows."""
import json, sys, numpy as np, pandas as pd, uproot
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
Q3 = OUT / "q3"; B = 300; rng = np.random.default_rng(20261010)
prodf = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root"); full = json.load(open(Q3 / "full_hists.json")); NB = len(BINS) - 1
def binidx(e): 
    i = np.digitize(e, BINS) - 1; i[(i < 0) | (i >= NB)] = -1; i[e == BINS[-1]] = NB - 1; return i
def ratio(num, den): return np.where(den > 0, num / np.where(den > 0, den, 1), 1.0)
def hist(i, w): m = i >= 0; return np.bincount(i[m], weights=w[m], minlength=NB)
res = {}
for det, cfgs in configs().items():
    for lab, hist_name, vtree, _, _ in cfgs:
        name = f"{det}__{lab}"; P = pd.read_pickle(PAIRS / f"{name}.pkl", compression=None) if False else pd.read_pickle(PAIRS / f"{name}.pkl")
        inn = P["in_n"].to_numpy(); inv = P["in_v"].to_numpy(); en = P["n_reco_visible_energy"].to_numpy(); ev = P["v_reco_visible_energy"].to_numpy()
        bn = np.where(inn, binidx(np.nan_to_num(en)), -1); bv = np.where(inv, binidx(np.nan_to_num(ev)), -1); both = inn & inv; N = len(P)
        prod = prodf[hist_name].values(); prod_nom = prodf[hist_name + "_nominal"].values(); prod_sc = prodf[hist_name + "_scaled"].values()
        one = np.ones(N)
        def methods(w):
            Hn = hist(bn, w); Hv = hist(bv, w); Hnb = hist(np.where(both, bn, -1), w); Hvb = hist(np.where(both, bv, -1), w)
            return dict(pair=ratio(Hv, Hn), matched=ratio(Hvb, Hnb)), Hn, Hv, Hnb, Hvb
        r0, Hn, Hv, Hnb, Hvb = methods(one)
        reps = {k: [] for k in r0}
        for _ in range(B):
            r, *_ = methods(rng.poisson(1.0, N).astype(float))
            for k in r: reps[k].append(r[k])
        sig = {k: np.std(v, axis=0, ddof=1) for k, v in reps.items()}
        fn, fv = full[name]["nom"], full[name]["var"]; Fn = np.array(fn["h"], float); Fv = np.array(fv["h"], float)
        r_full_raw = ratio(Fv, Fn); r_full_pot = ratio(Fv / fv["pot"], Fn / fn["pot"])
        sig_ind = r0["pair"] * np.sqrt(1 / np.maximum(Hn, 1) + 1 / np.maximum(Hv, 1))     # error if CV and variation were INDEPENDENT samples of this size
        sg = np.maximum(sig["pair"], 1e-12)
        def stat(r): d = r - prod; return dict(max_abs=float(np.abs(d).max()), rms=float(np.sqrt((d ** 2).mean())), chi2_paired=float(((d / sg) ** 2).sum()), chi2_indep=float(((d / np.maximum(sig_ind, 1e-12)) ** 2).sum()),
                                                max_pull_paired=float(np.abs(d / sg).max()), ndf=NB)
        res[name] = dict(prod=prod.tolist(), prod_nominal=prod_nom.tolist(), pair=r0["pair"].tolist(), matched=r0["matched"].tolist(), unmatched_full_raw=r_full_raw.tolist(), unmatched_full_pot=r_full_pot.tolist(),
                         sig_paired=sig["pair"].tolist(), sig_paired_matched=sig["matched"].tolist(), sig_indep=sig_ind.tolist(), Hn=Hn.tolist(), Hv=Hv.tolist(), Hn_both=Hnb.tolist(), Hv_both=Hvb.tolist(),
                         stats=dict(pair=stat(r0["pair"]), matched=stat(r0["matched"]), unmatched_full_raw=stat(r_full_raw), unmatched_full_pot=stat(r_full_pot)),
                         n=dict(rows=N, both=int(both.sum()), only_n=int((inn & ~inv).sum()), only_v=int((~inn & inv).sum()), full_nom_sel=fn["n_sel"], full_var_sel=fv["n_sel"], full_nom_events=fn["n_events"], full_var_events=fv["n_events"],
                                pot_nom=fn["pot"], pot_var=fv["pot"], prod_nom_sum=float(prod_nom.sum()), pair_nom_sum=float(Hn.sum())),
                         mean_sigma_ratio_ind_over_paired=float((sig_ind / sg).mean()))
        s = res[name]["stats"]; print(name, {k: (round(v["max_abs"], 4), round(v["chi2_paired"], 1)) for k, v in s.items()}, flush=True)
json.dump(res, open(Q3 / "emulation.json", "w"))
