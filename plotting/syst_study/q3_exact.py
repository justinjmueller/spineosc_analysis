#!/usr/bin/env python3
"""Q3: exact reproduction of the production ratio including its repeated-key rows (all selected rows in the pool are kept in the trees).
Variants: A (unmatched, selected-only trees, equal POT) and Ce (matched, selected-attribute eps weight, unmodified PROfit). usage: q3_exact.py <det__variation>"""
import sys, numpy as np, pandas as pd, uproot, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import q3_native; from q3_native import *
name = sys.argv[1]; det, lab = name.split("__"); cfg = [c for c in configs()[det] if c[0] == lab][0]; _, hname, vtree, (nfn, ntree), vfn = cfg
dn, mn = sel_tree(nfn, ntree, det); dv, mv = sel_tree(vfn, vtree, det); k4 = KEYS + ["Ei"]
pool = ekeys(nfn, ntree).drop_duplicates().merge(ekeys(vfn, vtree).drop_duplicates(), on=k4)
sn = dn[mn].merge(pool, on=k4)[k4 + ["reco_visible_energy"]]; sv = dv[mv].merge(pool, on=k4)[k4 + ["reco_visible_energy"]]
d = Q3 / name / "exact"; d.mkdir(parents=True, exist_ok=True); fn = d / "trees.root"
keys_all = pd.concat([sn[k4], sv[k4]]).drop_duplicates()
def full_side(s):   # all selected rows + one placeholder (selected=0) row for union keys with no selected row on this side
    have = s[k4].drop_duplicates(); miss = keys_all.merge(have, on=k4, how="left", indicator=True); miss = miss[miss._merge == "left_only"][k4]
    a = s.assign(selected=1); b = miss.assign(reco_visible_energy=np.float32(0), selected=0); return pd.concat([a, b], ignore_index=True)
def write(f, tn, df):
    br = {c: df[c].to_numpy().astype(np.int32) for c in k4}; br["reco_visible_energy"] = np.nan_to_num(df["reco_visible_energy"].to_numpy()).astype(np.float32); br["selected"] = df["selected"].to_numpy().astype(np.int32)
    t = f.mktree(tn, {c: v.dtype for c, v in br.items()}); t.extend(br)
with uproot.recreate(fn) as f:
    write(f, "cvA", sn.assign(selected=1)); write(f, "varA", sv.assign(selected=1)); write(f, "cvU", full_side(sn)); write(f, "varU", full_side(sv))
q3_native.make_trees = lambda nm: fn
prod = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root")[hname].values(); out = {}
for variant in ("A", "Ce"):
    vd = Q3 / name / ("exact_" + variant); vd.mkdir(exist_ok=True); make_xml(name, variant, vd)
    for cmd, pre in ((["process"], []), (["plot"], ["-o", "nom"]), (["plot"], ["-o", "inj", "--inject-systs-cv", "Q3_native", "1.0"])):
        with open(vd / "steps.log", "ab") as fh: rc = subprocess.run([str(PROFIT), "--config", str(vd / "q3.ini")] + pre + cmd, cwd=vd, stdout=fh, stderr=subprocess.STDOUT).returncode
        assert rc == 0, (variant, cmd)
    q3_native.Q3 = Q3; r, nom, inj = response(name, "exact_" + variant); out[variant] = dict(max_abs=float(np.abs(r - prod).max()), rms=float(np.sqrt(((r - prod) ** 2).mean())), resp=r.tolist())
    print(name, variant, "max|native - production|", out[variant]["max_abs"], "rms", out[variant]["rms"], flush=True)
json.dump(dict(prod=prod.tolist(), **out), open(Q3 / name / "exact.json", "w"))
