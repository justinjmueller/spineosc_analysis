#!/usr/bin/env python3
"""Q3: full-sample (NOT pool-restricted) selected histograms + POT + event counts for each detector/variation. This is what PROfit 'unmatched' mode sees if it is
given the whole CV and variation `selected` trees. cvwgt is not used."""
import json, sys, numpy as np, uproot
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
cols = ["reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x"]
def full(fn, tree, det):
    f = uproot.open(fn); t = f[f"events/{tree}/selected"]; have = set(t.keys()); d = t.arrays([c for c in cols if c in have], library="np")
    m = selection_mask(det, d); e = d["reco_visible_energy"][m]
    return dict(h=np.histogram(e, BINS)[0].tolist(), n_sel=int(m.sum()), n_events=int(f[f"events/{tree}/events"].num_entries), pot=float(f[f"events/{tree}/POT"].values()[0]))
out = {}; cache = {}
for det, cfgs in configs().items():
    for lab, hist, vtree, (nfn, ntree), vfn in cfgs:
        if (nfn, ntree) not in cache: cache[(nfn, ntree)] = full(nfn, ntree, det)
        out[f"{det}__{lab}"] = dict(nom=cache[(nfn, ntree)], var=full(vfn, vtree, det)); print(det, lab, out[f"{det}__{lab}"]["nom"]["n_sel"], out[f"{det}__{lab}"]["var"]["n_sel"], flush=True)
json.dump(out, open(OUT / "q3" / "full_hists.json", "w"))
