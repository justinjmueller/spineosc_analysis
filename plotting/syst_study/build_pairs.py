#!/usr/bin/env python3
"""Builds the paired event tables (runs/syst_study/pairs/<det>__<variation>.pkl) and a summary (runs/syst_study/pairs_summary.json).
Pool = events processed in BOTH samples (events-tree key intersection, as detsys/calculator). Selected rows restricted to the pool. Keys that repeat inside a selected tree
are excluded from the pairing (but counted in n_*_rows and in the production-style ratio)."""
import json, sys, time, pickle
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
PAIRS.mkdir(parents=True, exist_ok=True); summ = json.loads((OUT / "pairs_summary.json").read_text()) if (OUT / "pairs_summary.json").exists() else {}
only = sys.argv[1:]
for det, cfgs in configs().items():
    cache = {}
    for lab, hist, vtree, (nfn, ntree), vfn in cfgs:
        name = f"{det}__{lab}"
        if only and not any(o in name for o in only): continue
        t0 = time.time()
        if (nfn, ntree) not in cache:
            ek = ekeys(nfn, ntree); d, m = sel_tree(nfn, ntree, det); cache[(nfn, ntree)] = (ek.drop_duplicates(), d, m)
        ek_n, dn, mn = cache[(nfn, ntree)]
        ek_v = ekeys(vfn, vtree).drop_duplicates(); dv, mv = sel_tree(vfn, vtree, det)
        pool = ek_n.merge(ek_v, on=KEYS + ["Ei"], how="inner")
        sn = dn[mn].merge(pool, on=KEYS + ["Ei"], how="inner"); sv = dv[mv].merge(pool, on=KEYS + ["Ei"], how="inner")
        dupn = sn.duplicated(KEYS + ["Ei"], keep=False); dupv = sv.duplicated(KEYS + ["Ei"], keep=False)
        k = KEYS + ["Ei"]; un = sn[~dupn]; uv = sv[~dupv]
        feats = [c for c in sn.columns if c not in k and c != EKEY]
        a = un.rename(columns={c: "n_" + c for c in feats}).drop(columns=[EKEY]); b = uv.rename(columns={c: "v_" + c for c in feats}).drop(columns=[EKEY])
        pair = a.merge(b, on=k, how="outer", indicator=True)
        pair["in_n"] = pair["_merge"].isin(["both", "left_only"]); pair["in_v"] = pair["_merge"].isin(["both", "right_only"]); pair = pair.drop(columns=["_merge"])
        for c in pair.columns:
            if pair[c].dtype == np.float64: pair[c] = pair[c].astype(np.float32)
        pair.to_pickle(PAIRS / f"{name}.pkl")
        e = sn["reco_visible_energy"].to_numpy(); ev = sv["reco_visible_energy"].to_numpy()
        hn = np.histogram(e, BINS)[0]; hv = np.histogram(ev, BINS)[0]
        summ[name] = dict(det=det, label=lab, hist=hist, vtree=vtree, n_pool=int(len(pool)), n_n_rows=int(len(sn)), n_v_rows=int(len(sv)), dup_n=int(dupn.sum()), dup_v=int(dupv.sum()),
                          n_both=int((pair.in_n & pair.in_v).sum()), n_only_n=int((pair.in_n & ~pair.in_v).sum()), n_only_v=int((~pair.in_n & pair.in_v).sum()),
                          hist_n=hn.tolist(), hist_v=hv.tolist(), seconds=round(time.time() - t0, 1))
        print(name, {k_: v_ for k_, v_ in summ[name].items() if k_.startswith(("n_", "dup", "seconds"))}, flush=True)
        (OUT / "pairs_summary.json").write_text(json.dumps(summ))
