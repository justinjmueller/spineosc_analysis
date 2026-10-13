"""Shared definitions for the detector-systematic study (overnight 2026-10-10).
Paired event tables: for every detector/variation, all events common to the nominal (CV) and varied samples, with a selected flag on each side and the
reconstructed features of the candidate on the side(s) where it is selected. cvwgt is NOT used anywhere in this study."""
import pathlib, numpy as np, pandas as pd, uproot
REPO = pathlib.Path(__file__).resolve().parents[2]
D = pathlib.Path("/Users/mueller/data/spineosc/v1.1.0")
OUT = REPO / "runs" / "syst_study"; PAIRS = OUT / "pairs"; FIG = REPO / "docs" / "syst_study"
KEYS = ["Run", "Subrun", "Evt"]; EKEY = "event_max_neutrino_energy"
BINS = np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5, 2.0])   # production ratio binning (by_vis_energy)
SBND_DET = D / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"; SBND_NOS = D / "spineosc_numu_sbnd_v1.1.0_nosyst.root"
def icarus_det(run): return D / f"spineosc_numu_icarus_{run}_v1.1.0_detvar_nosyst.root"
def icarus_nos(run): return D / f"spineosc_numu_icarus_{run}_v1.1.0_nosyst.root"
# (label, production ratio histogram name, varied tree, nominal (file, tree), varied file)
def configs():
    out = {}
    sb = lambda lab, hist, vt: (lab, hist, vt, (SBND_NOS, "sbnd_large"), SBND_DET)
    out["sbnd_run1"] = [sb("ChargeScale", "ChargeScale_sbnd_Run1", "sbnd_detvar_charge_scale"), sb("ChargeScale15", "ChargeScale15_sbnd_Run1", "sbnd_detvar_smear15"),
        sb("SCE0x", "SCE0x_sbnd_Run1", "sbnd_detvar_sce0x"), sb("SCE2x", "SCE2x_sbnd_Run1", "sbnd_detvar_sce2x"),
        sb("WireModxThetaXW", "WireModxThetaXW_sbnd_Run1", "sbnd_detvar_wiremod_xthetaxw"), sb("WireModYZ", "WireModYZ_sbnd_Run1", "sbnd_detvar_wiremod_yz"),
        sb("MattModYZ", "MattModYZ_sbnd_Run1", "sbnd_detvar_mattmod_yz"), sb("GainHigh", "GainHigh_sbnd_Run1", "sbnd_detvar_p5_gain"), sb("GainLow", "GainLow_sbnd_Run1", "sbnd_detvar_m5_gain"),
        ("DENT", "DENT_sbnd_Run1", "sbnd_detvar_newdent", (SBND_DET, "sbnd_detvar_newcv"), SBND_DET)]
    for run, R in (("run2", "Run2"), ("run4", "Run4")):
        ic = lambda lab, hist, vt: (lab, f"{hist}_icarus_{R}", f"{vt}_{run}", (icarus_nos(run), f"icarus_{run}"), icarus_det(run))
        out[f"icarus_{run}"] = [ic("ChargeScale", "ChargeScale", "icarus_detvar_charge_scale"), ic("SCE0x", "SCE0x", "icarus_detvar_sce0x"), ic("SCE2x", "SCE2x", "icarus_detvar_sce2x"),
                                ic("WireModxThetaXW", "WireModxThetaXW", "icarus_wiremod_xthetaxw"), ic("WireModYZ", "WireModYZ", "icarus_wiremod_yz")]
    return out
def selection_mask(det, d):
    """Exactly the calculator's by_vis_energy cut (detsys/calculator_endx4)."""
    m = (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1)
    if det == "sbnd_run1": m &= (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4.0)
    return m
FEATS = ["reco_visible_energy", "reco_energy_qel", "reco_proton_multiplicity", "reco_muon_multiplicity", "reco_pion_multiplicity", "reco_photon_multiplicity", "reco_electron_multiplicity",
         "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_user_containment_cut_filterEE", "event_largest_flash_pe",
         "reco_vertex_x", "reco_vertex_y", "reco_vertex_z", "reco_pn", "reco_dpT"]
for p in ("muon", "proton"):
    FEATS += [f"reco_leading_{p}_{v}" for v in ("ke", "p", "length", "muon_softmax", "pion_softmax", "proton_softmax", "primary_softmax", "calo_csda_diff", "end_x", "end_y", "end_z", "start_dir_x", "start_dir_y", "start_dir_z", "polar_angle")]
FEATS += ["reco_leading_muon_mcs_csda_diff", "reco_leading_proton_mip_softmax", "reco_leading_proton_hadron_softmax", "reco_bivar_muon_proton_opening_angle",
          "true_category", "true_cc", "true_pdg", "true_neutrino_energy", "true_leading_muon_ke", "true_leading_proton_ke", "true_leading_muon_end_x", "true_interaction_mode", "true_detector"]
def ekeys(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/events"].arrays(KEYS + [EKEY], library="pd")
    d["Ei"] = np.round(d[EKEY].to_numpy() * 1e6).astype(np.int64); return d[KEYS + ["Ei"]]
def sel_tree(fn, tree, det):
    t = uproot.open(fn)[f"events/{tree}/selected"]; have = set(t.keys()); cols = [c for c in KEYS + [EKEY] + FEATS if c in have]
    d = t.arrays(cols, library="pd"); d["Ei"] = np.round(d[EKEY].to_numpy() * 1e6).astype(np.int64)
    return d, selection_mask(det, d).to_numpy()
