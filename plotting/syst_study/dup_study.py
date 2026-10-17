#!/usr/bin/env python3
"""Duplicate-event study over every v1.1.0 MC sample (nosyst, detvar_nosyst, dirt, wsyst, wsyst_cvrw) plus index-only checks of the ICARUS data files and the SBND release detector-access file.
Per tree (events / selected / signal):
  K1 = (Run, Subrun, Evt)                      K2 = K1 + round(1e6*event_max_neutrino_energy) (the key the calculator uses)
  rows, unique K1, unique K2, 'excess' rows = rows - unique keys, rows sitting in duplicated groups,
  fraction of the K1 degeneracy broken by adding the energy, and for the groups that stay duplicated under K2: are they FULL-ROW duplicates (every branch identical) and which branches differ otherwise.
Output: docs/syst_study/duplicates_all.json (+ per-sample detail for the anomalous trees).  Index-only for data files (no physics columns read)."""
import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import DATA_DIR, REPO  # noqa: E402
import json, sys, time, pathlib, numpy as np, pandas as pd, uproot
D = DATA_DIR; OUT = REPO / "docs/syst_study"
K1 = ["Run", "Subrun", "Evt"]; EK = "event_max_neutrino_energy"
MC = sorted([p for p in D.glob("*.root") if any(t in p.name for t in ("_nosyst.root", "_wsyst.root", "_wsyst_cvrw.root", "_detvar_nosyst.root", "_detvar_wsyst.root", "_dirt_nosyst.root", "_dirt_wsyst.root")) and "DONOTLOOK" not in p.name])
DATA_IDX = [D / "spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root", D / "spineosc_numu_icarus_run2_v1.1.0_bnblight.root", D / "spineosc_numu_icarus_run4_v1.1.0_bnblight.root", D / "spineosc_numu_icarus_release-v1.1.0_detector_access_DONOTLOOK.root"]
def keystats(df, with_energy):
    k = K1 + (["Ei"] if with_energy else []); g = df.groupby(k, sort=False).size(); n = len(df)
    return dict(unique=int(len(g)), excess=int(n - len(g)), rows_in_dup_groups=int(g[g > 1].sum()), n_dup_groups=int((g > 1).sum()), max_group=int(g.max()))
def analyse(df):
    n = len(df); out = dict(rows=n)
    has_e = EK in df.columns
    if has_e: df["Ei"] = np.round(df[EK].to_numpy(float) * 1e6).astype(np.int64)
    out["K1"] = keystats(df, False)
    if has_e:
        out["K2"] = keystats(df, True); e1 = out["K1"]["excess"]; out["energy_breaks_frac"] = float((e1 - out["K2"]["excess"]) / e1) if e1 else None
        out["Ei_nonpositive_rows"] = int((df["Ei"] <= 0).sum())
        # full-row duplicates (all branches identical)
        cols = [c for c in df.columns if c != "Ei"]; h = pd.util.hash_pandas_object(df[cols], index=False).to_numpy()
        df["_h"] = h; kk = K1 + ["Ei"]; grp = df.groupby(kk, sort=False)
        sz = grp["_h"].transform("size").to_numpy(); nu = grp["_h"].transform("nunique").to_numpy(); dup = sz > 1
        out["K2_dup_rows"] = int(dup.sum()); out["K2_dup_rows_all_identical"] = int((dup & (nu == 1)).sum()); out["K2_dup_rows_not_identical"] = int((dup & (nu > 1)).sum())
        full = df.duplicated(cols, keep=False); out["full_row_dup_rows"] = int(full.sum()); out["full_row_excess"] = int(n - len(df.drop_duplicates(cols)))
        # which branches differ inside K2-duplicate groups that are not identical
        if (dup & (nu > 1)).any():
            sub = df[dup & (nu > 1)]; diffcols = {}
            for c in cols:
                if c in kk: continue
                v = sub.groupby(kk, sort=False)[c].nunique(dropna=False); diffcols[c] = int((v > 1).sum())
            out["differing_branches"] = {k: v for k, v in sorted(diffcols.items(), key=lambda kv: -kv[1]) if v > 0}
        # groups under K1 only: all-identical / differ only by energy / other
        df["_h1"] = h; g1 = df.groupby(K1, sort=False); sz1 = g1["_h1"].transform("size").to_numpy(); nu1 = g1["_h1"].transform("nunique").to_numpy(); d1 = sz1 > 1
        out["K1_dup_rows_all_identical"] = int((d1 & (nu1 == 1)).sum()); out["K1_dup_rows_not_identical"] = int((d1 & (nu1 > 1)).sum())
        # duplicated-run/subrun structure of the K2 duplicates
        if dup.any():
            rs = df[dup].groupby(["Run", "Subrun"], sort=False).size(); tot = df.groupby(["Run", "Subrun"], sort=False).size().reindex(rs.index)
            out["K2_dup_subruns"] = int(len(rs)); out["K2_dup_runs"] = sorted({int(r) for r, _ in rs.index})[:40]
            frac = (rs / tot); out["K2_dup_subrun_frac_rows_dup_median"] = float(frac.median()); out["K2_dup_subrun_frac_rows_dup_max"] = float(frac.max()); out["K2_subruns_gt50pct"] = int((frac > 0.5).sum())
    else:
        h = pd.util.hash_pandas_object(df.drop(columns=[c for c in df.columns if c == "Ei"]), index=False); df["_h"] = h; out["full_row_excess"] = int(n - len(df.drop_duplicates()))
    return out
res = {}; t0 = time.time()
for f in MC:
    k = uproot.open(f); trees = [x for x in k.keys(cycle=False) if x.count("/") == 2 and x.split("/")[2] in ("events", "selected", "signal")]
    for tpath in sorted(trees):
        try:
            t = k[tpath]; n = t.num_entries
            if n == 0: continue
            cols = [c for c in t.keys() if c.startswith(("Run", "Subrun", "Evt", "event_", "true_", "reco_")) or c in K1]
            if tpath.endswith("/events"): cols = t.keys()
            df = t.arrays(cols, library="pd"); 
            if "Evt" not in df.columns: continue
            r = analyse(df); r["file"] = f.name; r["tree"] = tpath; res[f"{f.name}::{tpath}"] = r
            e2 = r.get("K2", {}).get("excess"); print(f"{f.name[:48]:48s} {tpath:60s} rows {r['rows']:9d} K1excess {r['K1']['excess']:8d} K2excess {e2 if e2 is not None else -1:8d} fullrow_excess {r['full_row_excess']:8d}  [{time.time()-t0:.0f}s]", flush=True)
            del df
        except Exception as ex: print("FAIL", f.name, tpath, ex, flush=True)
    json.dump(res, open(OUT / "duplicates_all.json", "w"))
# index-only checks for data
for f in DATA_IDX:
    if not f.exists(): continue
    k = uproot.open(f)
    for tpath in sorted(x for x in k.keys(cycle=False) if x.count("/") == 2 and x.split("/")[2] in ("selected", "events")):
        t = k[tpath]
        if t.num_entries == 0 or not all(c in t.keys() for c in K1): continue
        df = t.arrays(K1, library="pd"); r = dict(file=f.name, tree=tpath, rows=len(df), K1=keystats(df, False), index_only=True); res[f"{f.name}::{tpath}"] = r
        print(f"[data index-only] {f.name[:50]:50s} {tpath:50s} rows {r['rows']:8d} K1excess {r['K1']['excess']}", flush=True)
json.dump(res, open(OUT / "duplicates_all.json", "w")); print("done", round(time.time() - t0))
