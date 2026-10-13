#!/usr/bin/env python3
"""Extended paired tables (runs/syst_study/pairs2/<det>__<variation>.pkl): ALL rows of the production `selected` trees (before the per-detector post-cuts), restricted to the common pool,
paired by event key. Columns: tn/tv = row present in the selected tree of the nominal/varied sample; in_n/in_v = passes the post-cuts too (= the calculator's selection);
n_*/v_* features. Lets lost/gained events be attributed to (a) the production selection (tree membership) or (b) the post-cuts (fiducial, veto, crosser, end-x)."""
import json, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
P2 = OUT / "pairs2"; P2.mkdir(parents=True, exist_ok=True); only = sys.argv[1:]; k = KEYS + ["Ei"]
def load_all(fn, tree, det):
    t = uproot.open(fn)[f"events/{tree}/selected"]; have = set(t.keys()); cols = [c for c in KEYS + [EKEY] + FEATS if c in have]
    d = t.arrays(cols, library="pd"); d["Ei"] = np.round(d[EKEY].to_numpy() * 1e6).astype(np.int64); d["post"] = selection_mask(det, d).to_numpy(); return d
for det, cfgs in configs().items():
    cache = {}
    for lab, hist, vtree, (nfn, ntree), vfn in cfgs:
        name = f"{det}__{lab}"
        if only and not any(o in name for o in only): continue
        t0 = time.time()
        if (nfn, ntree) not in cache: cache[(nfn, ntree)] = (ekeys(nfn, ntree).drop_duplicates(), load_all(nfn, ntree, det))
        ek_n, dn = cache[(nfn, ntree)]; ek_v = ekeys(vfn, vtree).drop_duplicates(); dv = load_all(vfn, vtree, det)
        pool = ek_n.merge(ek_v, on=k, how="inner")
        sn = dn.merge(pool, on=k, how="inner"); sv = dv.merge(pool, on=k, how="inner")
        dn_ = sn.duplicated(k, keep=False); dv_ = sv.duplicated(k, keep=False); un = sn[~dn_]; uv = sv[~dv_]
        feats = [c for c in sn.columns if c not in k and c != EKEY]
        a = un.rename(columns={c: "n_" + c for c in feats}).drop(columns=[EKEY]); b = uv.rename(columns={c: "v_" + c for c in feats}).drop(columns=[EKEY])
        pair = a.merge(b, on=k, how="outer", indicator=True)
        pair["tn"] = pair["_merge"].isin(["both", "left_only"]); pair["tv"] = pair["_merge"].isin(["both", "right_only"]); pair = pair.drop(columns=["_merge"])
        pair["in_n"] = pair["tn"] & pair["n_post"].fillna(False).astype(bool); pair["in_v"] = pair["tv"] & pair["v_post"].fillna(False).astype(bool)
        pair = pair.copy()
        for c in pair.columns:
            if pair[c].dtype == np.float64: pair[c] = pair[c].astype(np.float32)
        pair.to_pickle(P2 / f"{name}.pkl")
        print(name, len(pair), dict(tree_both=int((pair.tn & pair.tv).sum()), tree_n_only=int((pair.tn & ~pair.tv).sum()), tree_v_only=int((~pair.tn & pair.tv).sum()), post_both=int((pair.in_n & pair.in_v).sum())), round(time.time() - t0, 1), flush=True)
