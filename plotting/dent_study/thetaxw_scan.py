#!/usr/bin/env python3
"""Does the WireMod X-thetaXW variation improve data/MC agreement as a function of theta_xW?  SBND Run 1 (detector-access data), per wire plane.
theta_xW = atan2(d . w_perp, |d_x|): the track angle from the drift axis in the plane spanned by x and the in-plane direction perpendicular to the wires.
SBND wires: collection vertical (along y); the two induction planes are at +-60 deg from the vertical. Mean score versus theta_xW: data, MC (CV), and MC + the
WireMod / ChargeScale shift taken from the detvar samples (full strength, a fixed point, no strength fit). Muon scores use the muon direction, proton scores the proton direction.
Run: python plotting/dent_study/thetaxw_scan.py         Outputs: docs/azimuth_scan/thetaxw_*"""
import json, pathlib, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import azimuth_scan as az, cathode_scan as cs
from mutools.plotting.helpers import mark_axis
OUT = az.OUT; NB = 18; EDGES = np.linspace(-90, 90, NB + 1); CTR = 0.5 * (EDGES[1:] + EDGES[:-1])
PART = {"muon_pid": "muon", "muon_primary": "muon", "proton_pid": "proton", "proton_primary": "proton"}
DIRS = [f"reco_leading_{p}_start_dir_{a}" for p in ("muon", "proton") for a in "xyz"]
PLANES = {"Collection": 0.0, "Induction 1 (+60$^\\circ$)": +60.0, "Induction 2 ($-$60$^\\circ$)": -60.0}   # wire angle from the vertical; collection is vertical
def theta(d, part, ang):
    dx, dy, dz = (d[f"reco_leading_{part}_start_dir_{a}"] for a in "xyz"); a = np.radians(ang)
    u = np.array([np.cos(a), np.sin(a)])                       # wire direction in (y, z)
    perp = -u[1] * dy + u[0] * dz                              # component along the in-plane perpendicular to the wires
    return np.degrees(np.arctan2(perp, np.abs(dx)))
def sel(d, mc):
    ok = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1) & (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4)
    if mc: ok &= 0.647 * d["event_largest_flash_pe"] > 2000
    return {k: v[ok].astype(float) for k, v in d.items()}
COLS = ["reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", *az.SCORES.values(), *DIRS]
PCOLS = [c for c in COLS if "proton_start_dir" not in c]
def sel_dirs(fn, tree, extra, mc):
    """MC and detvar samples carry the proton direction as polar angle and phi_dent (x = sin th sin phi, y = sin th cos phi, z = cos th; checked on the muon); data carries the components."""
    t = uproot.open(fn)[tree]
    if mc:
        d = t.arrays(PCOLS + ["reco_leading_proton_polar_angle", "reco_leading_proton_azimuthal_angle_dent"] + list(extra), library="np")
        th, ph = d["reco_leading_proton_polar_angle"], d["reco_leading_proton_azimuthal_angle_dent"]
        d["reco_leading_proton_start_dir_x"] = np.sin(th) * np.sin(ph); d["reco_leading_proton_start_dir_y"] = np.sin(th) * np.cos(ph); d["reco_leading_proton_start_dir_z"] = np.cos(th)
    else: d = t.arrays(COLS, library="np")
    return sel(d, mc)
def load(fn, tree, extra=()): return sel_dirs(fn, tree, extra, bool(extra))
def prof(x, j, w, m):
    mean = np.full(NB, np.nan); err = np.full(NB, np.nan)
    for i in range(NB):
        k = (j == i) & m
        if k.sum() > 5: ww = w[k]; mu = np.sum(ww * x[k]) / ww.sum(); neff = ww.sum() ** 2 / np.sum(ww ** 2); mean[i] = mu; err[i] = np.sqrt(np.sum(ww * (x[k] - mu) ** 2) / ww.sum() / neff)
    return mean, err
def chi2_flat(r, e):
    ok = np.isfinite(r); w = 1 / e[ok] ** 2; c = np.sum(w * r[ok]) / w.sum(); return float(np.sum(w * (r[ok] - c) ** 2)), int(ok.sum() - 1)
if __name__ == "__main__":
    c = az.CFG["sbnd"]; data = load(*c["data"]); mc = load(*c["mc"], extra=["cvwgt", "event_largest_flash_pe"])
    nos = cs.NOS; det = cs.DET
    VAR = {"WireMod X-$\\theta_{xW}$": "sbnd_detvar_wiremod_xthetaxw", "Charge scale": "sbnd_detvar_smear15"}
    KEYS = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
    def lk(fn, tree): return load(fn, f"events/{tree}/selected", extra=["event_largest_flash_pe"])
    def lkk(fn, tree):
        return sel_dirs(fn, f"events/{tree}/selected", KEYS + ["event_largest_flash_pe"], True)
    n0 = lkk(nos, "sbnd_large"); vs = {l: lkk(det, t) for l, t in VAR.items()}
    key = lambda d: list(zip(d["Run"].astype(int), d["Subrun"].astype(int), d["Evt"].astype(int), np.round(d["event_max_neutrino_energy"], 6)))
    res = {}
    fig, axs = plt.subplots(4, 3, figsize=(15, 16))
    for ci, (pl, ang) in enumerate(PLANES.items()):
        res[pl] = {}
        for ri, (nm, col) in enumerate(az.SCORES.items()):
            part = PART[nm]; jd = np.clip(np.digitize(theta(data, part, ang), EDGES) - 1, 0, NB - 1); jm = np.clip(np.digitize(theta(mc, part, ang), EDGES) - 1, 0, NB - 1)
            md, ed = prof(data[col], jd, np.ones(len(jd)), np.isfinite(data[col])); mm, em = prof(mc[col], jm, mc["cvwgt"], np.isfinite(mc[col]))
            delta = (md - mm) / np.nanmean(mm); derr = np.hypot(ed, em) / np.nanmean(mm)
            shifts = {}
            kn = set(key(n0))
            for l, v in vs.items():
                common = kn & set(key(v)); mv = np.array([k in common for k in key(v)]); mn = np.array([k in common for k in key(n0)])
                jv = np.clip(np.digitize(theta(v, part, ang), EDGES) - 1, 0, NB - 1); jn = np.clip(np.digitize(theta(n0, part, ang), EDGES) - 1, 0, NB - 1)
                a, _ = prof(v[col], jv, np.ones(len(jv)), mv & np.isfinite(v[col])); b, _ = prof(n0[col], jn, np.ones(len(jn)), mn & np.isfinite(n0[col])); shifts[l] = (a - b) / np.nanmean(mm)
            wm = shifts["WireMod X-$\\theta_{xW}$"]
            c0, nd0 = chi2_flat(delta, derr); c1, nd1 = chi2_flat(delta - wm, derr)
            res[pl][nm] = dict(delta=delta.tolist(), err=derr.tolist(), wm=wm.tolist(), cs=shifts["Charge scale"].tolist(), chi2_cv=c0, chi2_wm=c1, ndf=nd0)
            ax = axs[ri, ci]; ax.errorbar(CTR, delta * 100, xerr=5, yerr=derr * 100, fmt="o", color="k", ms=4, label="Data $-$ MC")
            ax.plot(CTR, (wm - np.nanmean(wm) + np.nanmean(delta)) * 100, color="#009E73", lw=2, label="WireMod shift")
            ax.plot(CTR, (shifts["Charge scale"] - np.nanmean(shifts["Charge scale"]) + np.nanmean(delta)) * 100, color="#0072B2", lw=1.6, label="Charge-scale shift")
            ax.text(0.03, 0.80, az.SHORT[nm], transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")
            ax.text(0.97, 0.06, rf"$\chi^2$/ndf: MC {c0:.0f}, MC+WM {c1:.0f} / {nd0}", transform=ax.transAxes, fontsize=10, ha="right")
            ax.set_xlim(-90, 90); ax.set_ylabel(r"Data $-$ MC [%]"); lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.35 * (hi - lo))
            ax.set_xlabel(rf"$\theta_{{xW}}$ ({part}) [deg]")
            if ri == 0: ax.set_title(pl, fontsize=13)
            print(pl, nm, f"chi2 flat (offset removed): MC {c0:.0f}  MC+WireMod {c1:.0f}  / {nd0}")
    axs[0, 0].legend(fontsize=9, loc="upper right"); mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 2], "SBND Run 1", alignment="right")
    fig.tight_layout(); fig.savefig(OUT / "thetaxw_sbnd.pdf"); fig.savefig(OUT / "thetaxw_sbnd.png", dpi=110)
    (OUT / "thetaxw_scan.json").write_text(json.dumps(res))
