#!/usr/bin/env python3
"""Q3: the construction the user describes -- an intermediate tree holding EVERY pool event (from the `events` trees) on both sides, with the `selected` attribute taken from the selected
trees. Tests scalability (6.3M rows/side) and correctness of the matched DetVar mode.   usage: q3_fullpool.py <det__variation> <Ce|Cp>"""
import sys, subprocess, time, resource, numpy as np, pandas as pd, uproot, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from q3_native import *
def make_full(name):
    det, lab = name.split("__"); cfg = [c for c in configs()[det] if c[0] == lab][0]; _, _, vtree, (nfn, ntree), vfn = cfg
    d = Q3 / name / "fullpool"; d.mkdir(parents=True, exist_ok=True); fn = d / "trees.root"
    if fn.exists(): return fn
    pool = ekeys(nfn, ntree).drop_duplicates().merge(ekeys(vfn, vtree).drop_duplicates(), on=KEYS + ["Ei"], how="inner")
    P = pd.read_pickle(PAIRS / f"{name}.pkl")[KEYS + ["Ei", "in_n", "in_v", "n_reco_visible_energy", "v_reco_visible_energy"]]
    pool = pool.merge(P, on=KEYS + ["Ei"], how="left"); inn = pool.in_n.fillna(False).to_numpy(bool); inv = pool.in_v.fillna(False).to_numpy(bool)
    k = {c: pool[c].to_numpy().astype(np.int32) for c in KEYS + ["Ei"]}
    with uproot.recreate(fn) as f:
        for tn, sel, e in (("cvU", inn, pool.n_reco_visible_energy), ("varU", inv, pool.v_reco_visible_energy)):
            br = dict(k); br["reco_visible_energy"] = np.where(sel, np.nan_to_num(e.to_numpy()), 0).astype(np.float32); br["selected"] = sel.astype(np.int32)
            t = f.mktree(tn, {c: v.dtype for c, v in br.items()}); t.extend(br)
    print("pool rows", len(pool), "selected cv/var", inn.sum(), inv.sum(), flush=True); return fn
if __name__ == "__main__":
    name, variant = sys.argv[1:3]; fn = make_full(name); d = Q3 / name / ("fullpool_" + variant); d.mkdir(exist_ok=True)
    # reuse the XML writer of q3_native by temporarily pointing make_trees at the full file
    import q3_native; q3_native.make_trees = lambda nm: fn
    q3_native.make_xml(name, variant, d)
    exe = PROFIT_PATCHED if variant == "Cp" else PROFIT
    for cmd, pre in ((["process"], []), (["plot"], ["-o", "nom"]), (["plot"], ["-o", "inj", "--inject-systs-cv", "Q3_native", "1.0"])):
        t0 = time.time(); r0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        with open(d / "steps.log", "ab") as fh: rc = subprocess.run([str(exe), "--config", str(d / "q3.ini")] + pre + cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
        print(name, variant, cmd[0], pre[:2], "rc", rc, f"{time.time()-t0:.0f}s maxrss_children_MB {resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1e6:.0f}", flush=True)
