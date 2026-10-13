#!/usr/bin/env python3
"""v1.1.0 binding-energy ratio histograms in reco_visible_energy, all three detectors plus the pooled BindingEnergy_SBN.

Differences from the v0 production of these histograms (make_be_syst.py run on the v1.0.10 samples):
  * inputs are the v1.1.0 CV-reweighted beam MC (work/*_cvrw.root) and every event is weighted by cvwgt, because the CV
    reweighting changes the interaction-mode mix and the shift is N_struck x 25 MeV;
  * the selection mirrors xml/systematics/systematics_suite.xml: >= 1 proton, cathode fiducialization, high-y/z veto, the SBND
    cathode-crosser cut (SBND only), the ICARUS Run 2 EE-filter containment, and the flash-PE trigger emulation per detector.
The shift itself (n_struck, deterministic half-sample, 25 MeV) is unchanged and imported from make_be_syst.py.
Writes a ROOT file with TH1D objects named like the entries in detvar/detsys-v1.1.0v1.root; inject_be.C copies them in."""
import argparse, pathlib, sys
import numpy as np, uproot
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from make_be_syst import shifted_mask, n_struck, D_BE, AXES

EDGES = AXES["by_vis_energy"][1]
ROOT = pathlib.Path("/Users/mueller/Projects/spineosc/work")
DETS = {  # label, file, tree
    "sbnd_Run1":   ("sbnd_run1_cvrw.root",   "events/sbnd_large/selected"),
    "icarus_Run2": ("icarus_run2_cvrw.root", "events/icarus_run2/selected"),
    "icarus_Run4": ("icarus_run4_cvrw.root", "events/icarus_run4/selected"),
}
TRIGGER = {"sbnd_Run1": (0.647, 2000.0), "icarus_Run2": (0.615, 6000.0), "icarus_Run4": (0.361, 3000.0)}
BR = ["reco_visible_energy", "reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser",
      "reco_user_containment_cut_filterEE", "reco_leading_muon_end_x", "event_largest_flash_pe", "true_interaction_mode", "true_neutrino_energy", "cvwgt"]

def selection(d, det, sbnd_endx_cut=0.0):
    sel = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1)
    if det == "sbnd_Run1":
        sel &= d["reco_cathode_crosser"] == 0
        if sbnd_endx_cut > 0: sel &= np.abs(d["reco_leading_muon_end_x"]) >= sbnd_endx_cut   # scenario: drop muons ending near the cathode
    if det == "icarus_Run2": sel &= d["reco_user_containment_cut_filterEE"] == 1
    scale, thr = TRIGGER[det]
    return sel & (scale * d["event_largest_flash_pe"] > thr)

def run(det, weighted=True, sbnd_endx_cut=0.0):
    fn, tree = DETS[det]
    d = uproot.open(ROOT / fn)[tree].arrays(BR, library="np")
    sel = selection(d, det, sbnd_endx_cut)
    d = {k: v[sel] for k, v in d.items()}
    w = d["cvwgt"] if weighted else np.ones(len(d["cvwgt"]))
    shift = shifted_mask(d["true_neutrino_energy"])
    e_cv = d["reco_visible_energy"]; e_var = e_cv + n_struck(d["true_interaction_mode"]) * np.where(shift, D_BE, 0.0) / 1000.0
    ok = np.isfinite(e_cv) & np.isfinite(e_var) & np.isfinite(w)
    nom, _ = np.histogram(e_cv[ok], EDGES, weights=w[ok]); sc, _ = np.histogram(e_var[ok], EDGES, weights=w[ok])
    print(f"{det}: {sel.sum()} selected events, weighted={weighted}, cvwgt mean {w[ok].mean():.3f}")
    return nom, sc

def write(path, results):
    with uproot.recreate(path) as f:
        pooled_n = sum(r[0] for r in results.values()); pooled_s = sum(r[1] for r in results.values())
        items = dict(results); items["SBN"] = (pooled_n, pooled_s)
        for det, (n, s) in items.items():
            ratio = np.where(n > 0, s / np.maximum(n, 1e-12), 1.0)
            for suffix, v in (("", ratio), ("_nominal", n), ("_scaled", s)):
                f[f"BindingEnergy_{det}{suffix}"] = (v, EDGES)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--output", required=True); ap.add_argument("--unweighted", action="store_true"); ap.add_argument("--sbnd-endx-cut", type=float, default=0.0, help="SBND only: require |reco_leading_muon_end_x| >= this [cm]")
    a = ap.parse_args()
    res = {det: run(det, not a.unweighted, a.sbnd_endx_cut) for det in DETS}
    write(a.output, res)
    for det, (n, s) in {**res, "SBN": (sum(r[0] for r in res.values()), sum(r[1] for r in res.values()))}.items():
        print(det, np.round(s / n, 3))
