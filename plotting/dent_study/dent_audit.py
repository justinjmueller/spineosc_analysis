#!/usr/bin/env python3
"""Audit of the DENT narrative.  Tests: (T0) null test: detvar CV vs other nominal samples, same pairing, same bootstrap: is the pairing noise-free?
(T1) 2x2 factorial of the two cuts (crosser removal, |x_end|>=4 cm): DENT/CV - 1 in the coarse energy bins, integrated, with paired bootstrap errors.
(T2) the pile-up and crosser fractions by TPC side (sign of x_end) for data, production CV, detvar CV, DENT.
(T3) detvar CV against production CV (cvwgt) for the same fractions: is the detvar CV a fair stand-in?
Output: docs/dent_study/audit.json.  Run: python plotting/dent_study/dent_audit.py"""
import json, pathlib, sys
import numpy as np, pandas as pd, uproot
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import cathode_scan as cs
K = cs.K; NB = 300; RNG = np.random.default_rng(20261008)
COLS = cs.COLS + ["event_largest_flash_pe"]
DCOLS = [c for c in cs.COLS if c not in K and c != "reco_visible_energy"]
def sel(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="pd").drop_duplicates(K, keep=False).set_index(K)
    d = d[(d.reco_fiducialize_cathode == 1) & (d.reco_veto_sbnd_highy_highz == 1) & (0.647 * d.event_largest_flash_pe > 2000)]
    return d
SCEN = {"none": lambda x, c: np.ones_like(x, bool), "crosser": lambda x, c: c == 0, "endx4": lambda x, c: x >= 4, "both": lambda x, c: (c == 0) & (x >= 4)}
def pair(a, b):
    ev = a.index.union(b.index); ev = ev[ev.isin(cs.keys(cs.DET if True else None, "sbnd_detvar_newcv")) | True]
    return ev
def table(A, B, evset):
    """A, B: selected frames (before crosser/endx cuts) for the two samples; evset: common event index."""
    ev = A.index.union(B.index); ev = ev[ev.isin(evset)]; n = len(ev); out = {}
    ja, jb = A.reindex(ev), B.reindex(ev); boots = RNG.poisson(1.0, size=(NB + 1, n)).astype(np.float32); boots[0] = 1
    for scen, f in SCEN.items():
        res = {}
        for ch, ci in (("1p", 0), ("Mp", 1)):
            co = cs.COARSE[ch]
            def counts(j):
                E = j.reco_visible_energy.to_numpy(); pm = j.reco_proton_multiplicity.to_numpy(); x = np.abs(j.reco_leading_muon_end_x.to_numpy()); c = j.reco_cathode_crosser.to_numpy()
                ok = ~np.isnan(E) & f(np.nan_to_num(x), np.nan_to_num(c)) & ((pm == 1) if ci == 0 else (pm > 1))
                cols = [boots[:, ok & (E >= co[i]) & (E < co[i + 1])].sum(axis=1) for i in range(3)]; return np.stack(cols, axis=1)
            Na, Nb = counts(ja), counts(jb); r = Nb / Na - 1; tot = Nb.sum(1) / Na.sum(1) - 1
            res[ch] = dict(coarse=r[0].tolist(), coarse_err=r[1:].std(0).tolist(), integrated=float(tot[0]), integrated_err=float(tot[1:].std()))
        out[scen] = res
    return out
def fractions(d, w=None):
    """crosser fraction of all selected; pile-up (|x_end|<4 cm) fraction of non-crossers; by side of x_end."""
    w = np.ones(len(d)) if w is None else w; c = d.reco_cathode_crosser.to_numpy(); x = d.reco_leading_muon_end_x.to_numpy(); ax = np.abs(x)
    out = dict(n=float(w.sum()), n_eff=float(w.sum() ** 2 / np.sum(w ** 2)), crosser=float(np.sum(w * (c == 1)) / w.sum()))
    nx = c == 0; out["pileup"] = float(np.sum(w * (nx & (ax < 4))) / np.sum(w * nx)); out["pileup_err"] = float(np.sqrt(out["pileup"] * (1 - out["pileup"]) / (np.sum(w * nx) ** 2 / np.sum((w * nx) ** 2))))
    for side, m in (("east", x > 0), ("west", x < 0)):
        nxs = nx & m; out[f"pileup_{side}"] = float(np.sum(w * (nxs & (ax < 4))) / np.sum(w * nxs)); out[f"n_{side}"] = float(np.sum(w * nxs))
        out[f"crosser_{side}"] = float(np.sum(w * ((c == 1) & m)) / np.sum(w * m))
    return out
if __name__ == "__main__":
    D = cs.DET; res = {}
    print("loading", flush=True)
    s = {t: sel(D, t) for t in ("sbnd_detvar_newcv", "sbnd_detvar_newdent", "sbnd_detvar_cv")}
    s["sbnd_large_nosyst"] = sel(cs.NOS, "sbnd_large")
    k = {t: cs.keys(D, t) for t in ("sbnd_detvar_newcv", "sbnd_detvar_newdent", "sbnd_detvar_cv")}; k["sbnd_large_nosyst"] = cs.keys(cs.NOS, "sbnd_large")
    common = lambda a, b: k[a].intersection(k[b])
    print("T1 DENT/CV factorial", flush=True); res["T1_dent_over_cv"] = table(s["sbnd_detvar_newcv"], s["sbnd_detvar_newdent"], common("sbnd_detvar_newcv", "sbnd_detvar_newdent"))
    print("T0 nulls", flush=True)
    res["T0_null_newcv_vs_cv"] = table(s["sbnd_detvar_newcv"], s["sbnd_detvar_cv"], common("sbnd_detvar_newcv", "sbnd_detvar_cv"))
    res["T0_null_newcv_vs_production_nosyst"] = table(s["sbnd_detvar_newcv"], s["sbnd_large_nosyst"], common("sbnd_detvar_newcv", "sbnd_large_nosyst"))
    print("T2/T3 fractions", flush=True)
    data = uproot.open("/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root")["events/sbnd_bnblight_full/selected"].arrays(DCOLS, library="pd")
    data = data[(data.reco_fiducialize_cathode == 1) & (data.reco_veto_sbnd_highy_highz == 1) & (data.reco_proton_multiplicity >= 1)]
    prod = uproot.open(cs.REPO / "work" / "sbnd_run1_cvrw.root")["events/sbnd_large/selected"].arrays(DCOLS + ["event_largest_flash_pe", "cvwgt"], library="pd")
    prod = prod[(prod.reco_fiducialize_cathode == 1) & (prod.reco_veto_sbnd_highy_highz == 1) & (prod.reco_proton_multiplicity >= 1) & (0.647 * prod.event_largest_flash_pe > 2000)]
    pm = lambda d: d[d.reco_proton_multiplicity >= 1]
    res["T2_data"] = fractions(data); res["T3_production_cv_cvwgt"] = fractions(prod, prod.cvwgt.to_numpy()); res["T3_production_cv_noweight"] = fractions(prod)
    res["T2_detvar_cv"] = fractions(pm(s["sbnd_detvar_newcv"])); res["T2_detvar_dent"] = fractions(pm(s["sbnd_detvar_newdent"]))
    pathlib.Path(cs.REPO / "docs" / "dent_study" / "audit.json").write_text(json.dumps(res, indent=1))
    for t in ("T1_dent_over_cv", "T0_null_newcv_vs_cv", "T0_null_newcv_vs_production_nosyst"):
        print(t)
        for sc, r in res[t].items(): print(f"  {sc:8s}", {ch: ([round(100 * v, 1) for v in r[ch]["coarse"]], round(100 * r[ch]["integrated"], 2), round(100 * r[ch]["integrated_err"], 2), [round(100 * v, 1) for v in r[ch]["coarse_err"]]) for ch in ("1p", "Mp")})
    for t in ("T2_data", "T3_production_cv_cvwgt", "T3_production_cv_noweight", "T2_detvar_cv", "T2_detvar_dent"): print(t, {k_: round(v, 4) for k_, v in res[t].items()})
