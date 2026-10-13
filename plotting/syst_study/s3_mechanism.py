#!/usr/bin/env python3
"""S3: where do the lost / gained events come from?  Uses pairs2 (all rows of the production `selected` trees, with post-cut flags).
Tier 1: the event is in the selected tree of only one sample (the production selection chain flipped: species assignment, thresholds, multiplicities, containment).
Tier 2: the event is in both selected trees but fails a post-cut (fiducial, veto, crosser, |x_end|) on one side.  Also: proton-multiplicity migration of events selected in both."""
import json, numpy as np
from ana import *
res = {}
for name in pairs_available(v2=True):
    df = load(name, v2=True); det = name.split("__")[0]; inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); tn = df.tn.to_numpy(); tv = df.tv.to_numpy()
    lost = inn & ~inv; gain = ~inn & inv
    def post_reason(side, m):
        p = "v_" if side == "n" else "n_"   # the side where the event FAILS is the other one
        fid = df[p + "reco_fiducialize_cathode"].to_numpy(float) != 1; veto = df[p + "reco_veto_sbnd_highy_highz"].to_numpy(float) != 1
        out = dict(fid=int((m & fid).sum()), veto=int((m & ~fid & veto).sum()))
        if det == "sbnd_run1":
            cr = df[p + "reco_cathode_crosser"].to_numpy(float) == 1; ex = np.abs(df[p + "reco_leading_muon_end_x"].to_numpy(float)) < 4.0
            out["crosser"] = int((m & ~fid & ~veto & cr).sum()); out["endx"] = int((m & ~fid & ~veto & ~cr & ex).sum())
        return out
    # lost: in tree on nominal; fails on varied either because absent from varied tree (tier 1) or fails post cut there (tier 2)
    lost_t1 = lost & ~tv; lost_t2 = lost & tv; gain_t1 = gain & ~tn; gain_t2 = gain & tn
    both = inn & inv; mn = df["n_reco_proton_multiplicity"].to_numpy(float); mv = df["v_reco_proton_multiplicity"].to_numpy(float)
    cat = lambda m: np.where(m == 1, 0, 1)
    mig = np.zeros((2, 2), int)
    for a in (0, 1):
        for b in (0, 1): mig[a, b] = int((both & (cat(mn) == a) & (cat(mv) == b)).sum())
    En = df["n_reco_visible_energy"].to_numpy(float)
    res[name] = dict(n_nom=int(inn.sum()), lost=int(lost.sum()), gained=int(gain.sum()), lost_tier1=int(lost_t1.sum()), lost_tier2=int(lost_t2.sum()), gained_tier1=int(gain_t1.sum()), gained_tier2=int(gain_t2.sum()),
                     lost_post=post_reason("n", lost_t2), gained_post=post_reason("v", gain_t2), mult_migration=mig.tolist(),
                     tier1_both_in_tree_frac=float((tn & tv).sum() / max(1, (tn | tv).sum())))
    r = res[name]; print(f"{name:32s} lost {100*r['lost']/r['n_nom']:5.2f}% (tier1 {100*r['lost_tier1']/r['lost']:4.1f}% / post-cut {100*r['lost_tier2']/r['lost']:4.1f}%: {r['lost_post']}) | gained {100*r['gained']/r['n_nom']:5.2f}% (tier1 {100*r['gained_tier1']/r['gained']:4.1f}%: post {r['gained_post']}) | 1p->Mp {mig[0,1]} Mp->1p {mig[1,0]}")
json.dump(res, open(FIG / "s3_mechanism.json", "w"))
