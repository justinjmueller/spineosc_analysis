#!/usr/bin/env python3
"""Q2 step 0: compact per-pair cache (key hash, run, subrun, bins of the nominal/varied candidates) in the scratchpad. cvwgt not used."""
import sys, numpy as np, pandas as pd
from ana import *
SCR = pathlib.Path("/private/tmp/claude-550/-Users-mueller-Projects-spineosc/32857cee-ac1d-4b8e-b98a-74f527609260/scratchpad/q2cache"); SCR.mkdir(exist_ok=True, parents=True)
def keyhash(run, sub, evt, ei):
    x = (run.astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)) ^ (sub.astype(np.uint64) * np.uint64(0xC2B2AE3D27D4EB4F)) ^ (evt.astype(np.uint64) * np.uint64(0x165667B19E3779F9)) ^ (ei.astype(np.uint64) * np.uint64(0x27D4EB2F165667C5))
    x ^= x >> np.uint64(33); x *= np.uint64(0xff51afd7ed558ccd); x ^= x >> np.uint64(33); x *= np.uint64(0xc4ceb9fe1a85ec53); x ^= x >> np.uint64(33); return x
if __name__ == "__main__":
    for name in pairs_available():
        df = load(name); bn, bv = prep(df)
        np.savez_compressed(SCR / f"{name}.npz", key=keyhash(*[df[c].to_numpy().astype(np.int64) for c in ("Run", "Subrun", "Evt", "Ei")]), run=df.Run.to_numpy().astype(np.int32), sub=df.Subrun.to_numpy().astype(np.int32),
                            bn=bn.astype(np.int8), bv=bv.astype(np.int8), in_n=df.in_n.to_numpy(), in_v=df.in_v.to_numpy(),
                            en=df.n_reco_visible_energy.to_numpy().astype(np.float32), ev=df.v_reco_visible_energy.to_numpy().astype(np.float32))
        print(name, len(df), flush=True)
