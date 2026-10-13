#!/usr/bin/env python3
"""Q3: diagnostics of the matching key used by PROfit (Run,Subrun,Evt,event_max_neutrino_energy; every value goes through float and is then cast with static_cast<int>(round(.)))."""
import sys, json, numpy as np, pandas as pd, uproot
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent)); from common import *
out = {}
for lab, fn, tr in (("sbnd_nominal_large", SBND_NOS, "sbnd_large"), ("sbnd_detvar_smear15", SBND_DET, "sbnd_detvar_smear15"), ("sbnd_newcv", SBND_DET, "sbnd_detvar_newcv"), ("icarus_run2_nominal", icarus_nos("run2"), "icarus_run2"),
                    ("icarus_run4_nominal", icarus_nos("run4"), "icarus_run4"), ("icarus_run4_wiremod_xthetaxw", icarus_det("run4"), "icarus_wiremod_xthetaxw_run4")):
    d = uproot.open(fn)[f"events/{tr}/events"].arrays(KEYS + [EKEY], library="np"); n = len(d["Run"]); E = d[EKEY]
    R, S, V = d["Run"].astype(np.int64), d["Subrun"].astype(np.int64), d["Evt"].astype(np.int64); Ei = np.round(E * 1e6).astype(np.int64)
    f32 = lambda a: a.astype(np.float32).astype(np.int64)
    def nuniq(*cols): return int(len(np.unique(np.stack(cols, 1), axis=0)))
    out[lab] = dict(n=n, max_run=int(R.max()), max_subrun=int(S.max()), max_evt=int(V.max()), max_Ei=int(Ei.max()), max_E_GeV=float(E.max()), neg_or_zero_E=int((E <= 0).sum()),
                    uniq_RSE=nuniq(R, S, V), uniq_RSE_Ei=nuniq(R, S, V, Ei), uniq_RSE_float32=nuniq(f32(R), f32(S), f32(V)), uniq_RSE_Ei_float32=nuniq(f32(R), f32(S), f32(V), f32(Ei)),
                    uniq_RSE_roundGeV=nuniq(R, S, V, np.round(E.astype(np.float32)).astype(np.int64)), float32_exact=bool(max(R.max(), S.max(), V.max(), Ei.max()) < 2 ** 24))
    print(lab, out[lab], flush=True)
json.dump(out, open(OUT / "q3" / "keys.json", "w"), indent=1)
