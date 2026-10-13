#!/usr/bin/env python3
"""Do the duplicated rows occur in streaks (wholly duplicated files)?  For each affected tree, in storage (entry) order:
 * map every repeated row (identical in all branches) to the entry of its first occurrence;
 * a 'streak' is a maximal run of consecutive repeated entries j, j+1, ... whose first occurrences i, i+1, ... are also consecutive (constant offset j - i): a contiguous block of the tree copied verbatim;
 * compare with the (Run, Subrun) segmentation of the tree (contiguous stretches of one (Run, Subrun) = presumed file/segment) and with event_is_first_in_subrun.
Output: docs/syst_study/duplicates_streaks.json + figure."""
import json, numpy as np, pandas as pd, uproot, pathlib
D = pathlib.Path("/Users/mueller/data/spineosc/v1.1.0"); OUT = pathlib.Path("/Users/mueller/Projects/spineosc/docs/syst_study")
CASES = [("ICARUS Run4 wiremod_xthetaxw (events)", "spineosc_numu_icarus_run4_v1.1.0_detvar_nosyst.root", "events/icarus_wiremod_xthetaxw_run4/events"),
         ("ICARUS Run4 dirt (events)", "spineosc_numu_icarus_run4_v1.1.0_dirt_nosyst.root", "events/icarus_dirt_run4/events"),
         ("ICARUS Run2 dirt (events)", "spineosc_numu_icarus_run2_v1.1.0_dirt_nosyst.root", "events/icarus_dirt_run2/events"),
         ("SBND large (events)", "spineosc_numu_sbnd_v1.1.0_nosyst.root", "events/sbnd_large/events"),
         ("ICARUS Run4 offbeam (events)", "spineosc_numu_icarus_run4_v1.1.0_nosyst.root", "events/icarus_offbeam_run4/events"),
         ("ICARUS Run4 data bnblight (events)", "spineosc_numu_icarus_run4_v1.1.0_bnblight.root", "events/icarus_bnblight_run4/events")]
res = {}
for lab, fn, tree in CASES:
    df = uproot.open(D / fn)[tree].arrays(library="pd"); cols = list(df.columns); n = len(df)
    h = pd.util.hash_pandas_object(df[cols], index=False).to_numpy(); first = {}; src = np.full(n, -1)
    for j, x in enumerate(h):
        if x in first: src[j] = first[x]
        else: first[x] = j
    rep = np.where(src >= 0)[0]; off = rep - src[rep]
    # streaks: consecutive repeated entries with constant offset
    streaks = []; s = 0
    for k in range(1, len(rep) + 1):
        if k == len(rep) or rep[k] != rep[k - 1] + 1 or off[k] != off[k - 1]:
            streaks.append((int(rep[s]), int(rep[k - 1]), int(off[s]), k - s)); s = k
    L = np.array([x[3] for x in streaks]) if streaks else np.array([0])
    # segmentation by (Run,Subrun) in entry order
    rs = (df.Run.astype("int64") * 100000 + df.Subrun.astype("int64")).to_numpy(); seg_start = np.r_[0, np.where(np.diff(rs) != 0)[0] + 1]; seg_end = np.r_[seg_start[1:], n] - 1
    seg_len = seg_end - seg_start + 1
    # does each streak coincide with a whole segment (start/end on segment boundaries)?
    starts = set(seg_start.tolist()); ends = set(seg_end.tolist()); big = [x for x in streaks if x[3] >= 10]
    aligned_both = sum(1 for a, b, d, l in big if a in starts and b in ends); aligned_start = sum(1 for a, b, d, l in big if a in starts); aligned_end = sum(1 for a, b, d, l in big if b in ends)
    # sorted-by-Evt within segment? and relation of the copy to the original (same segment? adjacent?)
    same_seg = np.mean([rs[j] == rs[src[j]] for j in rep]) if len(rep) else float("nan")
    offs = pd.Series(off).value_counts().head(5).to_dict() if len(rep) else {}
    res[lab] = dict(rows=n, repeated_rows=len(rep), n_streaks=len(streaks), streaks_ge10=len(big), rows_in_streaks_ge10=int(sum(x[3] for x in big)), frac_rows_in_streaks_ge10=float(sum(x[3] for x in big) / max(1, len(rep))), streak_len_median=float(np.median(L)), streak_len_max=int(L.max()), n_segments=len(seg_start), seg_len_median=float(np.median(seg_len)),
                    streaks_aligned_to_segment_both=aligned_both, aligned_start=aligned_start, aligned_end=aligned_end, copy_in_same_segment_frac=float(same_seg), top_offsets=offs, streaks=streaks[:200])
    print(f"{lab:40s} rows {n:8d} repeated {len(rep):7d} | streaks {len(streaks):5d} (>=10 rows: {len(big):4d}, covering {100*res[lab]['frac_rows_in_streaks_ge10']:.1f}% of repeats) median len {np.median(L):.0f} max {L.max()} | segments {len(seg_start)} (median len {np.median(seg_len):.0f}) streaks aligned start/end/both {aligned_start}/{aligned_end}/{aligned_both} of {len(big)} | copy in same (Run,Subrun) segment {100*same_seg:.1f}% | top offsets {list(offs.items())[:3]}", flush=True)
json.dump(res, open(OUT / "duplicates_streaks.json", "w"))
