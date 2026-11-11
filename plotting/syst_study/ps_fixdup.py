#!/usr/bin/env python3
"""Rebuilds the paired table of icarus_run4__WireModxThetaXW with duplicate (Run,Subrun,Evt,Ei) rows de-duplicated (keep first) instead of excluded.
Reason: the varied tree has 13,878 duplicated events (0.5%), only in this sample; excluding all copies from the varied side only (as build_pairs.py does) turns ~0.9k varied-selected events into spurious 'lost' events.
Output runs/syst_study/pairs_fixdup/icarus_run4__WireModxThetaXW.pkl and docs/syst_study/persyst/fixdup.json"""
from ps_common import *
from common import configs, ekeys, sel_tree
KEYK = KEYS + ["Ei"]; (OUT / "pairs_fixdup").mkdir(exist_ok=True)
cf = {c[0]: c for c in configs()["icarus_run4"]}["WireModxThetaXW"]; lab, hist, vtree, (nfn, ntree), vfn = cf
ek_n = ekeys(nfn, ntree).drop_duplicates(); ek_v = ekeys(vfn, vtree).drop_duplicates(); dn, mn = sel_tree(nfn, ntree, "icarus_run4"); dv, mv = sel_tree(vfn, vtree, "icarus_run4")
pool = ek_n.merge(ek_v, on=KEYK, how="inner")
sn = dn[mn].merge(pool, on=KEYK, how="inner").drop_duplicates(KEYK); sv = dv[mv].merge(pool, on=KEYK, how="inner").drop_duplicates(KEYK)
feats = [c for c in sn.columns if c not in KEYK and c != EKEY]
a = sn.rename(columns={c: "n_" + c for c in feats}).drop(columns=[EKEY]); b = sv.rename(columns={c: "v_" + c for c in feats}).drop(columns=[EKEY])
pair = a.merge(b, on=KEYK, how="outer", indicator=True); pair["in_n"] = pair["_merge"].isin(["both", "left_only"]); pair["in_v"] = pair["_merge"].isin(["both", "right_only"]); pair = pair.drop(columns=["_merge"])
for c in pair.columns:
    if pair[c].dtype == np.float64: pair[c] = pair[c].astype(np.float32)
pair.to_pickle(OUT / "pairs_fixdup" / "icarus_run4__WireModxThetaXW.pkl")
prod = prod_hist("WireModxThetaXW_icarus_Run4"); bn, bv = prep(pair); n0, v0 = counts(bn, bv); old = pd.read_pickle(OUT / "pairs" / "icarus_run4__WireModxThetaXW.pkl"); bo, bvo = prep(old); no, vo = counts(bo, bvo)
res = dict(n_dup_events_in_varied_events_tree=13878, ratio_dedup=(v0 / n0).tolist(), ratio_excl_groups=(vo / no).tolist(), ratio_prod=prod.tolist(), int_dedup=float(v0.sum() / n0.sum() - 1), int_excl=float(vo.sum() / no.sum() - 1),
           int_prod=float(prod @ n0 / n0.sum() - 1), lost_dedup=float((pair.in_n & ~pair.in_v).sum() / pair.in_n.sum()), gained_dedup=float((~pair.in_n & pair.in_v).sum() / pair.in_n.sum()),
           lost_excl=float((old.in_n & ~old.in_v).sum() / old.in_n.sum()), gained_excl=float((~old.in_n & old.in_v).sum() / old.in_n.sum()))
json.dump(res, open(PSFIG / "fixdup.json", "w")); print({k: (round(v, 4) if isinstance(v, float) else None) for k, v in res.items() if isinstance(v, float)}); print("mean per-bin prod-dedup %", 100 * np.mean(prod - v0 / n0), " prod-excl", 100 * np.mean(prod - vo / no))
