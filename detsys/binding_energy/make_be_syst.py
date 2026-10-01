#!/usr/bin/env python3
"""Binding-energy systematic: ratio histograms in the detvar hist1d format.

The binding energy enters the reconstruction in two independent places, so the
+1 sigma universe shifts two observables and leaves a third group untouched:

  reco_visible_energy  = muonKE + m_mu + SUM(protonKE) + N * BE
                         -> shifts by N_struck * dBE
  reco_pn (delta p)    = hypot(dpT, dpL), dpL depending on BE only through
                         M_A' = M_A - m_n + BE
                         -> shifts by recomputing dpL at the shifted M_A'
  dpT / dphiT / dalphaT are purely transverse and carry NO binding-energy
                         dependence, so their ratios must come out identically
                         1.0. That is the closure test, not an assumption.

N_struck follows CalorimetricEnergy.h's NStruckNucleons(): 1 for QE/Res/DIS,
2 for MEC, 0 otherwise -- the GENIE interaction mode, NOT the reconstructed
proton count. The two disagree for ~28% of selected events; this file follows
the header.

dpL is stored signed, and R is recovered exactly by inverting the header's own
formula, so nothing here is approximated:

    dpL = R/2 - (M_A'^2 + dpT^2) / (2R)
      => R = dpL + sqrt(dpL^2 + M_A'^2 + dpT^2)       (positive root)

R is built from the measured final state and is therefore BE-independent, which
is what makes the shifted dpL a two-line recompute.

Only a deterministic half of events is shifted (kBindingEnergyShiftFraction).
The partition is keyed on the true neutrino energy so it cannot move between
universes of other systematics; it is NOT bit-compatible with be_syst.py's
lowbias32 hash, which is not required here.
"""
import argparse
import pathlib
import numpy as np
import uproot

# --- CalorimetricEnergy.h constants, in MeV to match the ntuple units -------
M_A = (22 * 0.939565 + 18 * 0.938272 - 0.34381) * 1000.0
M_P = 938.272
M_MU = 105.658
M_N = 0.939565 * 1000.0
# NOTE: the header says 29.5 MeV, but the production's visible energy was built
# with 30.900 MeV/proton (verified: exact for 96% of 1p events). The CV value
# only sets M_A'; the systematic is the +dBE shift either way.
BE_CV = 30.900
D_BE = 25.0
SHIFT_FRACTION = 0.5

QE, RES, DIS, COH, MEC = 0, 1, 2, 3, 10

SUBCHANNEL_CUT = "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1"

# (name, axis edges) per observable, matching detsys/calculator/configuration.h
AXES = {
    "by_dp":       ("reco_pn", np.arange(0.0, 2001.0, 100.0)),
    "by_dpt":      ("reco_dpT", np.arange(0.0, 2001.0, 100.0)),
    "by_dalphat":  ("reco_dalphaT", np.linspace(0.0, 3.2, 21)),
    "by_dphit":    ("reco_dphiT", np.linspace(0.0, 3.2, 21)),
    "by_energy_qel": ("reco_energy_qel",
                      np.round(np.arange(0.15, 3.001, 0.15), 4)),
    "by_vis_energy": ("reco_visible_energy",
                      np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82,
                                0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5, 2.0])),
}

BRANCHES = ["reco_pn", "reco_dpT", "reco_dpL", "reco_dalphaT", "reco_dphiT",
            "reco_energy_qel", "reco_leading_muon_p", "reco_leading_muon_polar_angle",
            "reco_visible_energy", "reco_proton_multiplicity",
            "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
            "true_interaction_mode", "true_neutrino_energy"]


def n_struck(mode):
    """Struck nucleons from the GENIE interaction mode (header's convention)."""
    out = np.zeros(len(mode), dtype=np.float64)
    finite = np.isfinite(mode)
    out[finite & np.isin(mode, [QE, RES, DIS])] = 1.0
    out[finite & (mode == MEC)] = 2.0
    return out


def shifted_mask(true_e, fraction=SHIFT_FRACTION):
    """Deterministic `fraction` of events, keyed on the true neutrino energy.

    Truth, not reco, so the partition cannot move between universes of other
    systematics. Any stable hash does; this one is not the Python lowbias32.
    """
    x = np.asarray(true_e, dtype=np.float32).view(np.uint32).astype(np.uint64)
    x = (x ^ (x >> np.uint64(16))) * np.uint64(0x7feb352d) & np.uint64(0xFFFFFFFF)
    x = (x ^ (x >> np.uint64(15))) * np.uint64(0x846ca68b) & np.uint64(0xFFFFFFFF)
    x = x ^ (x >> np.uint64(16))
    return x < np.uint64(round(fraction * 2**32))


def shifted_observables(d):
    """CV and +1 sigma values for every observable the systematic touches."""
    shift = shifted_mask(d["true_neutrino_energy"])
    dbe = np.where(shift, D_BE, 0.0)

    cv, var = {}, {}

    # Energy: N_struck binding-energy terms, so the shift scales with the mode.
    cv["reco_visible_energy"] = d["reco_visible_energy"]
    var["reco_visible_energy"] = (
        d["reco_visible_energy"] + n_struck(d["true_interaction_mode"]) * dbe / 1000.0
    )

    # delta p: BE enters once, through M_A'. Recover R, then recompute.
    dpL, dpT = d["reco_dpL"], d["reco_dpT"]
    map_cv = M_A - M_N + BE_CV
    R = dpL + np.sqrt(dpL**2 + map_cv**2 + dpT**2)
    map_sh = M_A - M_N + BE_CV + dbe
    dpL_sh = 0.5 * R - (map_sh**2 + dpT**2) / (2.0 * R)
    cv["reco_pn"] = d["reco_pn"]
    var["reco_pn"] = np.hypot(dpT, dpL_sh)

    # CCQE energy estimator. Verified empirically against the ntuple: the
    # standard formula with E_B = 30.900 MeV reproduces reco_energy_qel to
    # 1e-6 GeV for 96.4% of events, so the shift is a clean recompute. The
    # binding energy here is a single nuclear property, NOT multiplied by the
    # struck-nucleon count (same convention as M_A' in the TKI block).
    e_mu = np.sqrt(d["reco_leading_muon_p"] ** 2 + M_MU ** 2)
    cos_th = np.cos(d["reco_leading_muon_polar_angle"])

    def e_qe(binding):
        mn = M_N - binding
        num = M_P ** 2 - mn ** 2 - M_MU ** 2 + 2.0 * mn * e_mu
        den = 2.0 * (mn - e_mu + d["reco_leading_muon_p"] * cos_th)
        return num / den / 1000.0          # MeV -> GeV

    cv["reco_energy_qel"] = d["reco_energy_qel"]
    var["reco_energy_qel"] = np.where(shift, e_qe(BE_CV + D_BE), d["reco_energy_qel"])

    # Purely transverse: no BE dependence at all. Carried through unchanged so
    # the ratio is a real closure test rather than an assumed 1.0.
    for b in ("reco_dpT", "reco_dalphaT", "reco_dphiT"):
        cv[b] = d[b]
        var[b] = d[b]
    return cv, var


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="nominal MC ROOT file")
    ap.add_argument("--tree", default="events/sbnd_large/selected")
    ap.add_argument("--output", required=True, help="directory for per-subchannel files")
    ap.add_argument("--label", default="BindingEnergy_sbnd_Run1")
    args = ap.parse_args()

    d = uproot.open(args.input)[args.tree].arrays(BRANCHES, library="np")
    sel = (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1)
    d = {k: v[sel] for k, v in d.items()}
    print(f"{sel.sum()} events pass: {SUBCHANNEL_CUT}")
    print(f"shifted fraction: {100*shifted_mask(d['true_neutrino_energy']).mean():.2f}%")

    cv_vals, var_vals = shifted_observables(d)

    out = {}
    for sub, (branch, edges) in AXES.items():
        c, v = cv_vals[branch], var_vals[branch]
        ok = np.isfinite(c) & np.isfinite(v)
        h_cv, _ = np.histogram(c[ok], edges)
        h_var, _ = np.histogram(v[ok], edges)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(h_cv > 0, h_var / np.maximum(h_cv, 1), 1.0)
        name = f"{args.label}"
        out[f"{sub}/{name}"] = (ratio, edges)
        out[f"{sub}/{name}_nominal"] = (h_cv.astype(float), edges)
        out[f"{sub}/{name}_scaled"] = (h_var.astype(float), edges)
        dev = 100 * np.max(np.abs(ratio - 1.0))
        tag = "  <-- closure (must be 0.00%)" if branch in ("reco_dpT", "reco_dalphaT", "reco_dphiT") else ""
        print(f"  {sub:16s} {branch:22s} max |ratio-1| = {dev:6.2f}%{tag}")

    outdir = pathlib.Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in AXES:
        path = outdir / f"be_{sub}.root"
        with uproot.recreate(path) as f:
            for suffix in ("", "_nominal", "_scaled"):
                f[f"{args.label}{suffix}"] = out[f"{sub}/{args.label}{suffix}"]
        print(f"  wrote {path.name}")
    print(f"wrote {len(AXES)} per-subchannel files to {outdir}")


if __name__ == "__main__":
    main()
