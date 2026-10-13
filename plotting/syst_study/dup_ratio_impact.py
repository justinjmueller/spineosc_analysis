#!/usr/bin/env python3
"""Impact of exact (full-row) duplicate rows on the 20 production ratio histograms: ratio with all rows (production) versus ratio after removing rows that are identical in every stored branch
(on both sides, after the pool restriction and the calculator cut)."""
import json, sys, numpy as np, pandas as pd, uproot
from common import *
k = KEYS + ["Ei"]; res = {}
def load(fn, tree):
    t = uproot.open(fn)[f"events/{tree}/selected"]; cols = [c for c in t.keys() if c.startswith(("Run", "Subrun", "Evt", "event_", "true_", "reco_"))]
    d = t.arrays(cols, library="pd"); d["Ei"] = np.round(d[EKEY].to_numpy(float) * 1e6).astype(np.int64); return d
for det, cfgs in configs().items():
    cache = {}
    for lab, hist, vtree, (nfn, ntree), vfn in cfgs:
        if (nfn, ntree) not in cache:
            dn = load(nfn, ntree); cache[(nfn, ntree)] = (ekeys(nfn, ntree).drop_duplicates(), dn)
        ek_n, dn = cache[(nfn, ntree)]; ek_v = ekeys(vfn, vtree).drop_duplicates(); dv = load(vfn, vtree); pool = ek_n.merge(ek_v, on=k, how="inner")
        out = {}
        for side, d in (("n", dn), ("v", dv)):
            m = selection_mask(det, d).to_numpy(); s = d[m].merge(pool, on=k, how="inner"); cols = [c for c in s.columns if c != "Ei"]
            full = s.duplicated(cols, keep="first"); sk = s.duplicated(k, keep="first")
            out[side] = dict(rows=len(s), K2_excess=int(sk.sum()), full_dup_excess=int(full.sum()), hist_all=np.histogram(s["reco_visible_energy"], BINS)[0], hist_dedup=np.histogram(s.loc[~full, "reco_visible_energy"], BINS)[0])
        r_all = out["v"]["hist_all"] / out["n"]["hist_all"]; r_dd = out["v"]["hist_dedup"] / out["n"]["hist_dedup"]
        res[f"{det}__{lab}"] = dict(n_rows=out["n"]["rows"], v_rows=out["v"]["rows"], n_K2_excess=out["n"]["K2_excess"], v_K2_excess=out["v"]["K2_excess"], n_full_dup=out["n"]["full_dup_excess"], v_full_dup=out["v"]["full_dup_excess"],
                                    integrated_all=float(out["v"]["hist_all"].sum() / out["n"]["hist_all"].sum() - 1), integrated_dedup=float(out["v"]["hist_dedup"].sum() / out["n"]["hist_dedup"].sum() - 1),
                                    max_bin_shift=float(np.abs(r_dd - r_all).max()), mean_bin_shift=float(np.mean(r_dd - r_all)), ratio_all=r_all.tolist(), ratio_dedup=r_dd.tolist())
        r = res[f"{det}__{lab}"]; print(f"{det}__{lab:16s} rows n/v {r['n_rows']}/{r['v_rows']}  K2-excess n/v {r['n_K2_excess']}/{r['v_K2_excess']}  full-row dups n/v {r['n_full_dup']}/{r['v_full_dup']}  integrated all {100*r['integrated_all']:+.2f}% dedup {100*r['integrated_dedup']:+.2f}%  max bin shift {100*r['max_bin_shift']:.2f}%", flush=True)
        json.dump(res, open(FIG / "duplicates_ratio_impact.json", "w"))
