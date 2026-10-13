#!/usr/bin/env python3
"""Cut-value scan for the SBND muon end-point cut (see docs/cathode_motivation/PLAN.md for the pre-registered criterion).
Writes docs/cathode_motivation/scan.json. Run with the mutools environment:  python plotting/dent_study/cathode_scan.py"""
import json, pathlib, sys
import numpy as np, pandas as pd, uproot
REPO = pathlib.Path(__file__).resolve().parents[2]
D = pathlib.Path("/Users/mueller/data/spineosc/v1.1.0")
DET = D / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"; NOS = D / "spineosc_numu_sbnd_v1.1.0_nosyst.root"
K = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
XS = [0, 1, 2, 3, 4, 5, 6, 8, 10]
FINE = {"1p": np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5]),
        "Mp": np.array([0.4, 0.5, 0.59, 0.67, 0.75, 0.83, 0.91, 0.99, 1.07, 1.16, 1.26, 1.37, 1.5, 1.7, 2.0])}
COARSE = {"1p": [0.3, 0.7, 1.0, 1.5], "Mp": [0.4, 0.8, 1.2, 2.0]}
NBOOT = 500; RNG = np.random.default_rng(20261007)
COLS = K + ["reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", "reco_proton_multiplicity"]

def keys(fn, tree):
    return pd.MultiIndex.from_frame(uproot.open(fn)[f"events/{tree}/events"].arrays(K, library="pd").drop_duplicates(K))
def sel(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="pd").drop_duplicates(K, keep=False).set_index(K)
    return d[(d.reco_fiducialize_cathode == 1) & (d.reco_veto_sbnd_highy_highz == 1)]   # crossers kept; each mode filters on the flag below

def fine_index(E, pm):
    out = np.full(len(E), -1)
    for ci, (ch, m) in enumerate((("1p", pm == 1), ("Mp", pm > 1))):
        e = FINE[ch]; b = np.digitize(E, e) - 1; ok = m & (b >= 0) & (b < len(e) - 1); out[ok] = ci * 100 + b[ok]
    return out

def main():
    print("loading detvar CV and DENT ...", flush=True)
    common = keys(DET, "sbnd_detvar_newcv").intersection(keys(DET, "sbnd_detvar_newdent"))
    cv = sel(DET, "sbnd_detvar_newcv"); dt = sel(DET, "sbnd_detvar_newdent")
    cv = cv[cv.index.isin(common)]; dt = dt[dt.index.isin(common)]
    ev = cv.index.union(dt.index); n = len(ev); print("events selected in either:", n)
    def cols(d):
        j = d.reindex(ev); return j.reco_visible_energy.to_numpy(), j.reco_proton_multiplicity.to_numpy(), np.abs(j.reco_leading_muon_end_x.to_numpy()), j.reco_cathode_crosser.to_numpy()
    Ec, Pc, Xc, Cc = cols(cv); Ed, Pd, Xd, Cd = cols(dt)
    in_c = ~np.isnan(Ec); in_d = ~np.isnan(Ed)
    nx_c = in_c & (Cc == 0); nx_d = in_d & (Cd == 0)      # non-crossers
    res = {"X": [-1] + XS, "channels": {}, "band": {}}   # X = -1: no crosser cut
    fc = fine_index(np.nan_to_num(Ec), np.nan_to_num(Pc)); fd = fine_index(np.nan_to_num(Ed), np.nan_to_num(Pd))
    boots = RNG.poisson(1.0, size=(NBOOT + 1, n)).astype(np.float32); boots[0] = 1.0     # row 0 = nominal
    for ch, ci in (("1p", 0), ("Mp", 1)):
        e = FINE[ch]; nb = len(e) - 1; cen = 0.5 * (e[1:] + e[:-1])
        co = COARSE[ch]; cmap = [(np.where((e[:-1] >= co[i] - 1e-9) & (e[1:] <= co[i + 1] + 1e-9))[0]) for i in range(3)]
        out = []
        for X in [-1] + XS:
            passc = (in_c if X < 0 else nx_c & (Xc >= X)); passd = (in_d if X < 0 else nx_d & (Xd >= X))
            ic = np.where(passc & (fc // 100 == ci), fc % 100, -1); idd = np.where(passd & (fd // 100 == ci), fd % 100, -1)
            # counts per bootstrap replica and fine bin
            Nc = np.zeros((NBOOT + 1, nb)); Nd = np.zeros((NBOOT + 1, nb))
            for b in range(nb):
                Nc[:, b] = boots[:, ic == b].sum(axis=1); Nd[:, b] = boots[:, idd == b].sum(axis=1)
            coarse_c = np.stack([Nc[:, m].sum(axis=1) for m in cmap], axis=1); coarse_d = np.stack([Nd[:, m].sum(axis=1) for m in cmap], axis=1)
            r = coarse_d / coarse_c - 1; tot = Nd.sum(axis=1) / Nc.sum(axis=1) - 1
            cc = np.array([0.5 * (co[i] + co[i + 1]) for i in range(3)])
            slope = np.array([np.polyfit(cc, r[k], 1)[0] for k in range(NBOOT + 1)])
            out.append(dict(X=X, coarse=r[0].tolist(), coarse_err=r[1:].std(axis=0).tolist(), integrated=float(tot[0]), integrated_err=float(tot[1:].std()),
                            slope=float(slope[0]), slope_err=float(slope[1:].std()), kept_cv=float(Nc[0].sum() / (Nc[0].sum() if X == 0 else 1)),
                            fine=(Nd[0] / np.maximum(Nc[0], 1) - 1).tolist(), fine_err=(Nd[1:] / np.maximum(Nc[1:], 1)).std(axis=0).tolist(), n_cv=float(Nc[0].sum())))
        n0 = out[1]["n_cv"]      # crosser-cut selection (X = 0)
        for o in out: o["kept_cv"] = o["n_cv"] / n0; del o["n_cv"]
        res["channels"][ch] = out
    # band ratio: DENT/CV of events with |x_end| in [X, X+4) among non-crossers (all channels)
    for X in XS:
        cb = nx_c & (Xc >= X) & (Xc < X + 4); db = nx_d & (Xd >= X) & (Xd < X + 4)
        r = np.array([boots[k, db].sum() / boots[k, cb].sum() for k in range(NBOOT + 1)])
        res["band"][str(X)] = dict(ratio=float(r[0]), err=float(r[1:].std()))
    # efficiency cost, reco: selected true-signal events (categories 0, 2) with crosser==0, flash threshold on MC
    s = uproot.open(NOS)["events/sbnd_large/selected"].arrays(["reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", "true_category", "event_largest_flash_pe", "reco_proton_multiplicity"], library="pd")
    s = s[(s.reco_fiducialize_cathode == 1) & (s.reco_veto_sbnd_highy_highz == 1) & (s.reco_cathode_crosser == 0) & (0.647 * s.event_largest_flash_pe > 2000) & s.true_category.isin([0, 2]) & (s.reco_proton_multiplicity >= 1)]
    ax = np.abs(s.reco_leading_muon_end_x.to_numpy()); res["cost_reco"] = {str(X): float((ax < X).mean()) for X in XS}
    sg = uproot.open(NOS)["events/sbnd_large/signal"].arrays(["true_fiducialize_cathode", "true_veto_sbnd_highy_highz", "true_cathode_crosser", "true_leading_muon_end_x"], library="pd")
    sg = sg[(sg.true_fiducialize_cathode == 1) & (sg.true_veto_sbnd_highy_highz == 1) & (sg.true_cathode_crosser == 0)]
    at = np.abs(sg.true_leading_muon_end_x.to_numpy()); res["cost_true"] = {str(X): float((at < X).mean()) for X in XS}
    # data / CV / DENT fractions of selected events that are non-crossers with |x_end| < X (energy-free; MC with flash threshold)
    sys.path.insert(0, str(REPO / "plotting" / "dent_validation")); import dentval as dv
    d = dv.load(); fr = {}
    for k, df in d.items():
        nc = (df.reco_cathode_crosser == 0).to_numpy(); ae = df.abs_end_x.to_numpy(); N = len(df)
        fr[k] = {str(X): [float((nc & (ae < X)).sum() / N), float(np.sqrt(max((nc & (ae < X)).sum(), 1)) / N)] for X in XS}
        fr[k]["crosser_fraction"] = [float((~nc).sum() / N), float(np.sqrt((~nc).sum()) / N)]
    res["fractions"] = fr
    out = REPO / "docs" / "cathode_motivation" / "scan.json"; out.write_text(json.dumps(res, indent=1)); print("wrote", out)

if __name__ == "__main__": main()
