#!/usr/bin/env python3
"""Q2 deliverable 5 (part 1): common-pool flags. For each detector the pool of variation v is P_v = events-tree keys in BOTH nominal and varied samples; the common pool is C = intersection of all P_v
(SBND: the nine variations sharing the main nominal; DENT has its own CV sample and is treated separately: C_DENT = P_DENT ∩ C). Writes per pair a boolean row flag 'in_common' and the pool sizes."""
from q2_lib import *
from common import configs, ekeys
import time
out = {}
for det, cfgs in configs().items():
    hs = {}; nom = {}
    for lab, hist, vtree, (nfn, ntree), vfn in cfgs:
        t0 = time.time()
        if (nfn, ntree) not in nom: e = ekeys(nfn, ntree).drop_duplicates(); nom[(nfn, ntree)] = np.unique(keyhash(*[e[c].to_numpy().astype(np.int64) for c in ("Run", "Subrun", "Evt", "Ei")]))
        e = ekeys(vfn, vtree).drop_duplicates(); hv = np.unique(keyhash(*[e[c].to_numpy().astype(np.int64) for c in ("Run", "Subrun", "Evt", "Ei")]))
        hs[lab] = np.intersect1d(nom[(nfn, ntree)], hv); print(det, lab, len(hs[lab]), round(time.time() - t0, 1), flush=True)
    main = [l for l in hs if l != "DENT"]; C = hs[main[0]]
    for l in main[1:]: C = np.intersect1d(C, hs[l])
    out[det] = dict(n_common=int(len(C)), pools={l: int(len(h)) for l, h in hs.items()})
    if "DENT" in hs: out[det]["n_common_dent_main"] = int(len(np.intersect1d(hs["DENT"], C)))
    for l in hs:
        d = np.load(SCR / f"{det}__{l}.npz")["key"]; flag = np.isin(d, C if l != "DENT" else np.intersect1d(hs["DENT"], C)); np.save(RES / f"common_{det}__{l}.npy", flag); out[det].setdefault("row_frac_in_common", {})[l] = float(flag.mean())
    print(det, out[det], flush=True)
json.dump(out, open(Q2 / "q2_pool_sizes.json", "w"))
