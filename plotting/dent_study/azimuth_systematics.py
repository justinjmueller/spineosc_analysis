#!/usr/bin/env python3
"""Expected azimuth dependence of the SBND detector variations: mean-score shift versus muon azimuth for each detvar sample against its nominal, with the
same selection and slices as azimuth_scan.py, overlaid on the observed data - MC offset (shape of the offset after removing its azimuth-averaged value).
Run after azimuth_scan.py:  python plotting/dent_study/azimuth_systematics.py"""
import json, pathlib, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import azimuth_scan as az, cathode_scan as cs
from mutools.plotting.helpers import mark_axis
DET = cs.DET; NOS = cs.NOS; K = cs.K
VARS = [("Charge scale", "sbnd_detvar_smear15", "sbnd_large"), ("WireMod", "sbnd_detvar_wiremod_xthetaxw", "sbnd_large"),
        ("Non-uniformity", "sbnd_detvar_mattmod_yz", "sbnd_large"), ("DENT", "sbnd_detvar_newdent", "newcv")]   # SBND detector variations active in the analysis (no Gain, no SCE)
COLS = K + az.BASE + ["reco_cathode_crosser", "reco_leading_muon_end_x", "event_largest_flash_pe", "reco_leading_proton_azimuthal_angle_dent"]
def sel_load(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="np")
    ok = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1) & (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4) & (0.647 * d["event_largest_flash_pe"] > 2000)
    return {k: v[ok] for k, v in d.items()}
def keyset(d): return set(zip(d["Run"].astype(int), d["Subrun"].astype(int), d["Evt"].astype(int), np.round(d["event_max_neutrino_energy"], 6)))
def profile(d, nphi, col, mask=None, proton=False):
    phi = d["reco_leading_proton_azimuthal_angle_dent"] if proton else np.arctan2(d["reco_leading_muon_start_dir_x"], d["reco_leading_muon_start_dir_y"]); e = np.linspace(-np.pi, np.pi, nphi + 1); j = np.clip(np.digitize(phi, e) - 1, 0, nphi - 1)
    x = d[col]; m = np.ones(len(x), bool) if mask is None else mask; m &= np.isfinite(x)
    mean = np.array([x[(j == i) & m].mean() for i in range(nphi)]); err = np.array([x[(j == i) & m].std(ddof=1) / np.sqrt(((j == i) & m).sum()) for i in range(nphi)]); return mean, err
if __name__ == "__main__":
    nphi = 24; res = json.loads((az.OUT / "azimuth_scan.json").read_text())["sbnd"]; edges = np.array(res["edges"]); ctr = 0.5 * (edges[1:] + edges[:-1]) * 180 / np.pi
    nom = {"sbnd_large": sel_load(NOS, "sbnd_large"), "newcv": sel_load(DET, "sbnd_detvar_newcv")}
    out = {}
    for label, tree, nomkey in VARS:
        v = sel_load(DET, tree); n = nom[nomkey]; common = keyset(v) & keyset(n)
        mv = np.array([k in common for k in zip(v["Run"].astype(int), v["Subrun"].astype(int), v["Evt"].astype(int), np.round(v["event_max_neutrino_energy"], 6))])
        mn = np.array([k in common for k in zip(n["Run"].astype(int), n["Subrun"].astype(int), n["Evt"].astype(int), np.round(n["event_max_neutrino_energy"], 6))])
        out[label] = {}
        for nm, col in az.SCORES.items():
            pr = az.PHI_OF[nm] == "phi_p"; a, ea = profile(v, nphi, col, mv.copy(), pr); b, eb = profile(n, nphi, col, mn.copy(), pr); out[label][nm] = ((a - b).tolist(), np.hypot(ea, eb).tolist())
        print(label, "common events", len(common))
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    colors = ["#0072B2", "#009E73", "#CC79A7", "#E69F00"]
    for ax, (nm, r) in zip(axs.ravel(), res["scores"].items()):
        d = np.array(r["delta"]) - r["const"]; de = np.array(r["delta_err"])
        shifts = [np.array(out[label][nm][0]) for label, _, _ in VARS]; shifts = [s - s.mean() for s in shifts]
        env = np.sqrt(np.sum([s ** 2 for s in shifts], axis=0)) * 100                     # quadrature sum of the active variations about the CV
        ax.fill_between(ctr, -env, env, color="0.82", alpha=0.7, lw=0, label="Systematic envelope", zorder=1)
        out_ = np.abs(d * 100) > env + de * 100                                          # data not covered by the envelope within its statistical error
        ax.errorbar(ctr[~out_], d[~out_] * 100, yerr=de[~out_] * 100, fmt="o", color="k", ms=4.5, label="Data $-$ CV", zorder=5)
        ax.errorbar(ctr[out_], d[out_] * 100, yerr=de[out_] * 100, fmt="o", color="#D55E00", ms=6, mec="k", label="Not covered", zorder=6)
        for (label, _, _), s, c in zip(VARS, shifts, colors): ax.plot(ctr, s * 100, color=c, lw=1.6, label=label, zorder=3)
        ax.axhline(0, color="k", lw=2.2, zorder=4, label="CV"); ax.set_ylabel(r"Shift relative to CV [%]"); ax.set_xlim(-180, 180); ax.set_xlabel(az.XLAB[nm])
        ax.text(0.03, 0.80, az.SHORT[nm], transform=ax.transAxes, fontsize=13, va="top", fontweight="bold")
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.35 * (hi - lo))
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, fontsize=11, frameon=False)
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 1], "SBND Run 1", alignment="right"); fig.tight_layout(rect=(0, 0.07, 1, 1))
    print({nm: int(np.sum(np.abs(np.array(res["scores"][nm]["delta"]) - res["scores"][nm]["const"]) * 100 > 0)) for nm in az.SCORES})
    fig.savefig(az.OUT / "azimuth_mean_sbnd_detvar.pdf"); fig.savefig(az.OUT / "azimuth_mean_sbnd_detvar.png", dpi=120)
    (az.OUT / "azimuth_detvar.json").write_text(json.dumps(out)); print("done")
    for nm in az.SCORES:
        obs = np.array(res["scores"][nm]["delta"]) - res["scores"][nm]["const"]
        print(nm, "observed peak-to-peak %.2f%%" % (100 * np.ptp(obs)), "| detvar peak-to-peak %:", {l: round(100 * float(np.ptp(np.array(out[l][nm][0]))), 2) for l, _, _ in VARS})
