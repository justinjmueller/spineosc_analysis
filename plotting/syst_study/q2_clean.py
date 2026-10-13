#!/usr/bin/env python3
"""Q2: 'clean' ICARUS Run4 WireModxThetaXW: remove runs 12965-12970 where the varied sample contains duplicated events (events tree: 27756 rows duplicated; production counts duplicates twice, the
pair table drops the duplicated keys and so manufactures 'lost' events). Writes the clean cache table + bootstrap, and compares the production histogram with the clean ratio."""
from q2_lib import *
import uproot
n = "icarus_run4__WireModxThetaXW"; d = load(n); m = ~((d["run"] >= 12965) & (d["run"] <= 12970)); c = {k: v[m] for k, v in d.items()}
np.savez_compressed(SCR / f"{n}clean.npz", **{k: (v.astype(np.int8) if k in ("bn", "bv") else v) for k, v in c.items()})
N, V = boot_counts(c, 2000, 0); np.savez_compressed(RES / f"boot_{n}clean.npz", N=N, V=V)
r, C, Nn, Vv = cov_analytic(c); r0 = ratio(d); s = json.loads((OUT / "pairs_summary.json").read_text())[n]; prod = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root")[s["hist"]].values() if False else None
import pathlib; prod = uproot.open("/Users/mueller/Projects/spineosc/detvar/detsys-v1.1.0v2.root")[s["hist"]].values()
dl = prod - r; chi2 = float(dl @ np.linalg.solve(C, dl)); dl0 = prod - r0; _, C0, _, _ = cov_analytic(d)
print("removed events frac", 1 - len(c["key"]) / len(d["key"]), "| prod - clean (%):", np.round(100 * dl, 2), "chi2 vs clean", chi2, "of 17; const fit", float((dl / np.diag(C) ** 0.5 ** 1).mean()))
w = 1 / np.diag(C); print("weighted mean prod-clean %.3f%% +- %.3f%%" % (100 * (w * dl).sum() / w.sum(), 100 / np.sqrt(w.sum())))
print("pair - clean (%):", np.round(100 * (r0 - r), 2)); json.dump(dict(prod_minus_clean=dl.tolist(), chi2=chi2, pair_minus_clean=(r0 - r).tolist(), r_clean=r.tolist(), sigma_clean=np.sqrt(np.diag(C)).tolist(), const_wmean=float((w * dl).sum() / w.sum()), const_err=float(1 / np.sqrt(w.sum()))), open(Q2 / "q2_run4_wiremodx_clean.json", "w"))
