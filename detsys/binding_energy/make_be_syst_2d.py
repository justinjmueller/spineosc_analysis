#!/usr/bin/env python3
"""Binding-energy systematic as a 2D ratio in (visible energy, delta p).

A single +25 MeV shift moves BOTH observables on the same event: visible
energy by N_struck * dBE, and delta p through M_A' in the longitudinal
component. Either 1D projection therefore throws away half the response, and
the two projections disagree in sign across delta-phi_T -- the symptom of a
joint effect forced through one axis. This writes the joint ratio instead.

Output is a TH2 named <label>, consumed by a type="hist2d" allowlist entry
with xvar = the visible-energy variable and yvar = the delta-p variable.
"""
import argparse
import pathlib
import numpy as np
import uproot

M_A = (22 * 0.939565 + 18 * 0.938272 - 0.34381) * 1000.0
M_N = 0.939565 * 1000.0
BE_CV = 30.900
D_BE = 25.0
SHIFT_FRACTION = 0.5
QE, RES, DIS, MEC = 0, 1, 2, 10

# Coarse enough that every cell is populated: verified min 46 events, no empty
# cells, 100% of events inside the axes. An empty CV cell would silently give a
# ratio of 1 rather than erroring.
XEDGES = np.array([0.30, 0.50, 0.65, 0.80, 0.95, 1.15, 1.50, 3.00])   # visible energy [GeV]
YEDGES = np.array([0.0, 120.0, 240.0, 360.0, 500.0, 700.0, 1000.0, 3000.0])  # delta p [MeV/c]

BRANCHES = ["reco_pn", "reco_dpT", "reco_dpL", "reco_visible_energy",
            "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
            "true_interaction_mode", "true_neutrino_energy"]


def n_struck(mode):
    out = np.zeros(len(mode))
    fin = np.isfinite(mode)
    out[fin & np.isin(mode, [QE, RES, DIS])] = 1.0
    out[fin & (mode == MEC)] = 2.0
    return out


def shifted_mask(true_e, fraction=SHIFT_FRACTION):
    x = np.asarray(true_e, dtype=np.float32).view(np.uint32).astype(np.uint64)
    x = (x ^ (x >> np.uint64(16))) * np.uint64(0x7feb352d) & np.uint64(0xFFFFFFFF)
    x = (x ^ (x >> np.uint64(15))) * np.uint64(0x846ca68b) & np.uint64(0xFFFFFFFF)
    x = x ^ (x >> np.uint64(16))
    return x < np.uint64(round(fraction * 2 ** 32))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--tree", default="events/sbnd_large/selected")
    ap.add_argument("--output", required=True)
    ap.add_argument("--label", default="BindingEnergy2D_sbnd_Run1")
    a = ap.parse_args()

    d = uproot.open(a.input)[a.tree].arrays(BRANCHES, library="np")
    sel = (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1)
    d = {k: v[sel] for k, v in d.items()}
    shift = shifted_mask(d["true_neutrino_energy"])
    dbe = np.where(shift, D_BE, 0.0)
    print(f"{sel.sum()} events; shifted fraction {100*shift.mean():.2f}%")

    # visible energy: N_struck binding-energy terms
    ve_cv = d["reco_visible_energy"]
    ve_var = ve_cv + n_struck(d["true_interaction_mode"]) * dbe / 1000.0

    # delta p: BE enters once, through M_A'. R is built from the measured final
    # state and is BE-independent, so recover it and recompute.
    dpL, dpT = d["reco_dpL"], d["reco_dpT"]
    map_cv = M_A - M_N + BE_CV
    R = dpL + np.sqrt(dpL ** 2 + map_cv ** 2 + dpT ** 2)
    map_sh = M_A - M_N + BE_CV + dbe
    dp_var = np.hypot(dpT, 0.5 * R - (map_sh ** 2 + dpT ** 2) / (2.0 * R))
    dp_cv = d["reco_pn"]

    good = np.isfinite(ve_cv) & np.isfinite(dp_cv) & np.isfinite(ve_var) & np.isfinite(dp_var)
    nom, _, _ = np.histogram2d(ve_cv[good], dp_cv[good], [XEDGES, YEDGES])
    sca, _, _ = np.histogram2d(ve_var[good], dp_var[good], [XEDGES, YEDGES])
    ratio = np.where(nom > 0, sca / np.maximum(nom, 1), 1.0)

    pathlib.Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(a.output) as f:
        f[a.label] = (ratio, XEDGES, YEDGES)
        f[f"{a.label}_nominal"] = (nom, XEDGES, YEDGES)
        f[f"{a.label}_scaled"] = (sca, XEDGES, YEDGES)
    print(f"wrote {a.output}: {ratio.shape[0]}x{ratio.shape[1]}, "
          f"max |ratio-1| = {100*np.max(np.abs(ratio-1)):.2f}%, empty CV cells {int((nom==0).sum())}")
    print("\nratio (rows = visible energy, cols = delta p):")
    print("        " + "".join(f"{y:>8.0f}" for y in YEDGES[1:]))
    for i in range(ratio.shape[0]):
        print(f"{XEDGES[i]:5.2f}  " + "".join(f"{v:8.3f}" for v in ratio[i]))


if __name__ == "__main__":
    main()
