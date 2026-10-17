#!/usr/bin/env python3
"""Azimuthal scan of the PID and primary scores with the full detector-access statistics (SBND Run 1 and ICARUS Run 2).
Data: the release detector-access files (energy-free columns only). MC: the CV-reweighted production samples (cvwgt), flash/trigger emulation on the MC only.
Per azimuthal slice and score, the mean score of data and MC (shape information, area-normalized by construction) and the pull map of the score
distribution against the MC shape after removing the global (azimuth-independent) data/MC ratio in each score bin. Statistical uncertainties only for the data/MC
differences; the expected azimuth-dependent size of the detector systematic variations is taken from the detvar samples (SBND).
Run:  python plotting/dent_study/azimuth_scan.py          Outputs: docs/azimuth_scan/*"""
import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import DATA_DIR, REPO  # noqa: E402
import json, pathlib, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
OUT = REPO / "docs" / "azimuth_scan"; OUT.mkdir(parents=True, exist_ok=True)
D = DATA_DIR
SCORES = {"muon_pid": "reco_leading_muon_muon_softmax", "muon_primary": "reco_leading_muon_primary_softmax",
          "proton_pid": "reco_leading_proton_proton_softmax", "proton_primary": "reco_leading_proton_primary_softmax"}
LABEL = {"muon_pid": "Muon PID score", "muon_primary": "Muon primary score", "proton_pid": "Leading-proton PID score", "proton_primary": "Leading-proton primary score"}
BASE = ["reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_leading_muon_start_dir_x", "reco_leading_muon_start_dir_y", *SCORES.values()]
PHI_OF = {"muon_pid": "phi_mu", "muon_primary": "phi_mu", "proton_pid": "phi_p", "proton_primary": "phi_p"}   # muon scores vs muon azimuth, proton scores vs proton azimuth
CFG = {
  "sbnd": dict(label="SBND Run 1", data=(D / "spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root", "events/sbnd_bnblight_full/selected"),
               mc=(REPO / "work" / "sbnd_run1_cvrw.root", "events/sbnd_large/selected"), extra=["reco_cathode_crosser", "reco_leading_muon_end_x"], flash=(0.647, 2000.0)),
  "icarus": dict(label="ICARUS Run 2", data=(D / "spineosc_numu_icarus_release-v1.1.0_detector_access_DONOTLOOK.root", "events/icarus_run2_bnblight_full/selected"),
               mc=(REPO / "work" / "icarus_run2_cvrw.root", "events/icarus_run2/selected"), extra=["reco_user_containment_cut_filterEE"], flash=(0.615, 6000.0)),
}

def load(det, which):
    c = CFG[det]; fn, tree = c[which]; cols = BASE + c["extra"] + (["cvwgt", "event_largest_flash_pe", "reco_leading_proton_azimuthal_angle_dent"] if which == "mc" else ["reco_leading_proton_start_dir_x", "reco_leading_proton_start_dir_y"])
    d = uproot.open(fn)[tree].arrays(cols, library="np")
    sel = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1)
    if det == "sbnd": sel &= (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4)
    else: sel &= d["reco_user_containment_cut_filterEE"] == 1
    if which == "mc": sel &= c["flash"][0] * d["event_largest_flash_pe"] > c["flash"][1]
    out = {k: v[sel].astype(float) for k, v in d.items()}
    out["phi_mu"] = np.arctan2(out["reco_leading_muon_start_dir_x"], out["reco_leading_muon_start_dir_y"])   # same convention as the *_azimuthal_angle_dent branches
    out["phi_p"] = out["reco_leading_proton_azimuthal_angle_dent"] if which == "mc" else np.arctan2(out["reco_leading_proton_start_dir_x"], out["reco_leading_proton_start_dir_y"])
    out["w"] = out["cvwgt"] if which == "mc" else np.ones(sel.sum())
    return out

def analyse(det, nphi):
    data, mc = load(det, "data"), load(det, "mc"); edges = np.linspace(-np.pi, np.pi, nphi + 1); res = {"n_data": int(len(data["phi_mu"])), "edges": edges.tolist(), "scores": {}}
    for nm, col in SCORES.items():
        jd = np.clip(np.digitize(data[PHI_OF[nm]], edges) - 1, 0, nphi - 1); jm = np.clip(np.digitize(mc[PHI_OF[nm]], edges) - 1, 0, nphi - 1)
        xd, xm, wm = data[col], mc[col], mc["w"]; ok_d = np.isfinite(xd); ok_m = np.isfinite(xm)
        # score bins: 10 equal-population bins of the MC distribution
        q = np.quantile(xm[ok_m], np.linspace(0, 1, 11)); q[0] = -np.inf; q[-1] = np.inf; q = np.unique(q); K = len(q) - 1
        kd = np.clip(np.digitize(xd, q) - 1, 0, K - 1); km = np.clip(np.digitize(xm, q) - 1, 0, K - 1)
        mean_d = np.zeros(nphi); err_d = np.zeros(nphi); mean_m = np.zeros(nphi); err_m = np.zeros(nphi); nd = np.zeros(nphi)
        Dk = np.zeros((nphi, K)); Mk = np.zeros((nphi, K)); M2k = np.zeros((nphi, K))
        for j in range(nphi):
            a = (jd == j) & ok_d; b = (jm == j) & ok_m; nd[j] = a.sum()
            mean_d[j] = xd[a].mean(); err_d[j] = xd[a].std(ddof=1) / np.sqrt(a.sum())
            wb = wm[b]; mean_m[j] = np.average(xm[b], weights=wb); neff = wb.sum() ** 2 / (wb ** 2).sum()
            err_m[j] = np.sqrt(np.average((xm[b] - mean_m[j]) ** 2, weights=wb) / neff)
            Dk[j] = np.bincount(kd[a], minlength=K); Mk[j] = np.bincount(km[b], weights=wb, minlength=K); M2k[j] = np.bincount(km[b], weights=wb ** 2, minlength=K)
        delta = mean_d - mean_m; derr = np.hypot(err_d, err_m); cst = np.sum(delta / derr ** 2) / np.sum(1 / derr ** 2)
        chi2_const = float(np.sum(((delta - cst) / derr) ** 2))
        # expected counts: MC shape in each slice x the global data/MC ratio in each score bin, normalized to the slice's data total
        Rg = Dk.sum(axis=0) / np.maximum(Mk.sum(axis=0), 1e-12); E = Mk * Rg[None, :]; E *= (Dk.sum(axis=1) / E.sum(axis=1))[:, None]
        var = E + (E / np.maximum(Mk, 1e-12)) ** 2 * M2k; pull = (Dk - E) / np.sqrt(np.maximum(var, 1e-12))
        res["scores"][nm] = dict(delta=delta.tolist(), delta_err=derr.tolist(), const=float(cst), chi2_const=chi2_const, ndf_const=nphi - 1,
                                 pull=pull.tolist(), chi2_pull=float(np.sum(pull ** 2)), ndf_pull=int((nphi - 1) * (K - 1)), n_data=nd.tolist(), qedges=[float(v) for v in q[1:-1]])
    return res

SHORT = {"muon_pid": "Muon PID score", "muon_primary": "Muon primary score", "proton_pid": "Proton PID score", "proton_primary": "Proton primary score"}
XLAB = {"muon_pid": r"Muon azimuth $\phi_\mu$ [deg]", "muon_primary": r"Muon azimuth $\phi_\mu$ [deg]", "proton_pid": r"Proton azimuth $\phi_p$ [deg]", "proton_primary": r"Proton azimuth $\phi_p$ [deg]"}

def plots(det, res):
    e = np.array(res["edges"]); ctr = 0.5 * (e[1:] + e[:-1]) * 180 / np.pi; w = np.diff(e) * 90 / np.pi
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.5))
    for ax, (nm, r) in zip(axs.ravel(), res["scores"].items()):
        ax.errorbar(ctr, np.array(r["delta"]) * 100, xerr=w, yerr=np.array(r["delta_err"]) * 100, fmt="o", color="k", ms=4.5, elinewidth=1.2)
        ax.axhline(r["const"] * 100, color="#0072B2", lw=1.8); ax.axhline(0, color="0.5", lw=1)
        ax.text(0.03, 0.80, SHORT[nm], transform=ax.transAxes, fontsize=13, va="top", fontweight="bold")
        ax.text(0.97, 0.06, rf"$\chi^2$/ndf = {r['chi2_const']:.0f}/{r['ndf_const']}", transform=ax.transAxes, fontsize=12, ha="right")
        ax.set_xlim(-180, 180); ax.set_xlabel(XLAB[nm]); ax.set_ylabel(r"Data $-$ MC [%]")
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.35 * (hi - lo))
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 1], CFG[det]["label"], alignment="right"); fig.tight_layout()
    fig.savefig(OUT / f"azimuth_mean_{det}.pdf"); fig.savefig(OUT / f"azimuth_mean_{det}.png", dpi=120); plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.5))
    for ax, (nm, r) in zip(axs.ravel(), res["scores"].items()):
        p = np.array(r["pull"]); im = ax.imshow(p.T, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-4, vmax=4, extent=[-180, 180, 0, p.shape[1]])
        ax.text(0.03, 0.97, SHORT[nm], transform=ax.transAxes, fontsize=12, va="top", fontweight="bold"); ax.text(0.97, 0.97, rf"$\chi^2$/ndf = {r['chi2_pull']:.0f}/{r['ndf_pull']}", transform=ax.transAxes, fontsize=11, va="top", ha="right")
        ax.set_xlabel(XLAB[nm]); ax.set_ylabel("Score bin")
    fig.tight_layout(rect=(0, 0, 0.9, 0.94)); fig.colorbar(im, ax=axs.ravel().tolist(), label="Pull", shrink=0.8)
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 1], CFG[det]["label"], alignment="right")
    fig.savefig(OUT / f"azimuth_pulls_{det}.pdf"); fig.savefig(OUT / f"azimuth_pulls_{det}.png", dpi=120); plt.close(fig)

if __name__ == "__main__":
    allres = {}
    for det, nphi in (("sbnd", int(sys.argv[1]) if len(sys.argv) > 1 else 24), ("icarus", int(sys.argv[2]) if len(sys.argv) > 2 else 16)):
        r = analyse(det, nphi); allres[det] = r; plots(det, r)
        print(det, "data events", r["n_data"], "| per slice", int(min(r["scores"]["muon_pid"]["n_data"])), "to", int(max(r["scores"]["muon_pid"]["n_data"])))
        for nm, s in r["scores"].items(): print(f"   {nm:15s} mean-offset const {100*s['const']:+.2f}%  chi2_const {s['chi2_const']:.1f}/{s['ndf_const']}  pull-map chi2 {s['chi2_pull']:.0f}/{s['ndf_pull']}")
    (OUT / "azimuth_scan.json").write_text(json.dumps(allres, indent=1))
