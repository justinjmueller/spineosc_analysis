#!/usr/bin/env python3
"""Q2: block-jackknife / block-bootstrap inflation for ICARUS Run4 WireModxThetaXW after removal of the duplicated-event runs (cf. q2_s1_calibration)."""
from q2_lib import *
def block(d, gid, nb=1000):
    G = gid.max() + 1; Ng = np.zeros((G, NB)); Vg = np.zeros((G, NB))
    for b in range(NB): Ng[:, b] = np.bincount(gid[d["bn"] == b], minlength=G); Vg[:, b] = np.bincount(gid[d["bv"] == b], minlength=G)
    N = Ng.sum(0); V = Vg.sum(0); rj = (V - Vg) / (N - Ng); sj = np.sqrt((G - 1) / G * ((rj - rj.mean(0)) ** 2).sum(0)); return sj
out = {}
for n in ["icarus_run4__WireModxThetaXW", "icarus_run4__WireModxThetaXWclean"]:
    d = load(n); gid = np.unique(d["run"].astype(np.int64) * 100000 + d["sub"], return_inverse=True)[1]; r, C, N, V = cov_analytic(d); sj = block(d, gid); out[n] = dict(jk_over_ev_median=float(np.median(sj / np.sqrt(np.diag(C)))), G=int(gid.max() + 1)); print(n, out[n])
json.dump(out, open(Q2 / "q2_block_clean.json", "w"))
