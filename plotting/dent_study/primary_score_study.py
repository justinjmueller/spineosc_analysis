#!/usr/bin/env python3
"""What drives the SBND muon primary-score mismodeling? Data vs CV-reweighted MC (production), selection as in azimuth_scan.py.
1. the score distribution and where the data/MC ratio departs from 1; 2. the high-tail fraction (score above a threshold) and the mean score in bins of
event features (vertex/end-point positions, direction, proton multiplicity, other scores); 3. the shift produced by each detector variation (paired, detvar vs nominal).
Outputs docs/azimuth_scan/primary_*.png/.json. Run in the mutools environment."""
import json, pathlib, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import azimuth_scan as az, azimuth_systematics as sy, cathode_scan as cs
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
OUT = az.OUT; T = 0.975; COL = "reco_leading_muon_primary_softmax"
FEATS = {"vertex x [cm]": ("reco_vertex_x", None), "vertex y [cm]": ("reco_vertex_y", None), "vertex z [cm]": ("reco_vertex_z", None),
         "muon end x [cm]": ("reco_leading_muon_end_x", None), "muon end y [cm]": ("reco_leading_muon_end_y", None), "muon end z [cm]": ("reco_leading_muon_end_z", None),
         "muon direction z": ("reco_leading_muon_start_dir_z", None), "proton multiplicity": ("reco_proton_multiplicity", None),
         "muon PID score": ("reco_leading_muon_muon_softmax", None), "leading-proton primary score": ("reco_leading_proton_primary_softmax", None)}
EXTRA = ["reco_vertex_x", "reco_vertex_y", "reco_vertex_z", "reco_leading_muon_end_y", "reco_leading_muon_end_z", "reco_leading_muon_start_dir_z"]
def load(which):
    c = az.CFG["sbnd"]; fn, tree = c[which]; cols = az.BASE + c["extra"] + EXTRA + (["cvwgt", "event_largest_flash_pe"] if which == "mc" else [])
    d = uproot.open(fn)[tree].arrays(cols, library="np")
    sel = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1) & (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4)
    if which == "mc": sel &= 0.647 * d["event_largest_flash_pe"] > 2000
    o = {k: v[sel].astype(float) for k, v in d.items()}; o["w"] = o["cvwgt"] if which == "mc" else np.ones(sel.sum()); return o
data, mc = load("data"), load("mc")
xd, xm, wm = data[COL], mc[COL], mc["w"]
# 1. distribution
e = np.concatenate([np.linspace(0.8, 0.95, 4)[:-1], np.linspace(0.95, 1.0, 11)]); hd = np.histogram(xd, e)[0].astype(float); hm = np.histogram(xm, e, weights=wm)[0]; hm2 = np.histogram(xm, e, weights=wm ** 2)[0]
sc = hd.sum() / hm.sum(); ratio = hd / (hm * sc); rerr = np.sqrt(hd + (hd ** 2 / np.maximum(hm, 1e-12) ** 2) * hm2 / np.maximum(hm, 1e-12) ** 0) / (hm * sc) if False else np.sqrt(hd) / (hm * sc)
print("score bin: data/MC ratio"); [print(f"  {e[i]:.4f}-{e[i+1]:.4f}: {ratio[i]:.3f} +- {rerr[i]:.3f}") for i in range(len(e) - 1)]
# 2. features
def tail_stats(x, w, m): 
    n = w[m].sum(); f = w[m & (x > T)].sum() / n; neff = n ** 2 / (w[m] ** 2).sum(); return f, np.sqrt(f * (1 - f) / neff), np.average(x[m], weights=w[m])
res = {"ratio": ratio.tolist(), "ratio_err": rerr.tolist(), "edges": e.tolist(), "features": {}}
fig, axs = plt.subplots(2, 5, figsize=(22, 8)); 
fd0, _, md0 = tail_stats(xd, np.ones(len(xd)), np.ones(len(xd), bool)); fm0, _, mm0 = tail_stats(xm, wm, np.ones(len(xm), bool)); print(f"global high-tail fraction data {fd0:.4f} MC {fm0:.4f} ratio {fd0/fm0:.3f}; mean {md0:.4f} {mm0:.4f}")
for ax, (lab, (col, _)) in zip(axs.ravel(), FEATS.items()):
    q = np.unique(np.quantile(mc[col][np.isfinite(mc[col])], np.linspace(0, 1, 11) if col != "reco_proton_multiplicity" else [0, .45, .7, .85, .95, 1])) if col != "reco_proton_multiplicity" else np.array([0.5, 1.5, 2.5, 3.5, 4.5, 20])
    cd, cm = [], []; ctr = 0.5 * (q[1:] + q[:-1]); chi = 0.0
    for i in range(len(q) - 1):
        a = np.isfinite(xd) & (data[col] >= q[i]) & (data[col] < q[i + 1] + (1e-9 if i == len(q) - 2 else 0)); b = np.isfinite(xm) & (mc[col] >= q[i]) & (mc[col] < q[i + 1] + (1e-9 if i == len(q) - 2 else 0))
        fd, ed, mdn = tail_stats(xd, np.ones(len(xd)), a); fm, em, mmn = tail_stats(xm, wm, b); r = fd / fm; er = r * np.hypot(ed / fd, em / fm); cd.append((r, er)); chi += ((fd - fm) / np.hypot(ed, em)) ** 2
    r = np.array([c[0] for c in cd]); er = np.array([c[1] for c in cd]); mean_r = np.sum(r / er ** 2) / np.sum(1 / er ** 2); c2 = float(np.sum(((r - mean_r) / er) ** 2))
    ax.errorbar(ctr, r, yerr=er, fmt="o", color="k"); ax.axhline(mean_r, color="#0072B2"); ax.axhline(1, color="0.5", lw=1); ax.set_xlabel(lab); ax.set_ylabel(f"data/MC high-tail fraction (>{T})"); ax.text(0.04, 0.06, rf"const. $\chi^2$/ndf={c2:.1f}/{len(r)-1}", transform=ax.transAxes, fontsize=11)
    res["features"][lab] = dict(centers=ctr.tolist(), ratio=r.tolist(), err=er.tolist(), chi2_const=c2, ndf=len(r) - 1)
mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 4], "SBND Run 1", alignment="right"); fig.tight_layout(); fig.savefig(OUT / "primary_features.png", dpi=100); fig.savefig(OUT / "primary_features.pdf"); plt.close(fig)
# 3. detector variations: change of the high-tail fraction and mean (paired, same selection)
nom = {"sbnd_large": sy.sel_load(cs.NOS, "sbnd_large"), "newcv": sy.sel_load(cs.DET, "sbnd_detvar_newcv")}; sysres = {}
for label, tree, nk in sy.VARS:
    v = sy.sel_load(cs.DET, tree); n = nom[nk]; common = sy.keyset(v) & sy.keyset(n)
    kv = np.array([k in common for k in zip(v["Run"].astype(int), v["Subrun"].astype(int), v["Evt"].astype(int), np.round(v["event_max_neutrino_energy"], 6))]); kn = np.array([k in common for k in zip(n["Run"].astype(int), n["Subrun"].astype(int), n["Evt"].astype(int), np.round(n["event_max_neutrino_energy"], 6))])
    fv = (v[COL][kv] > T).mean(); fn = (n[COL][kn] > T).mean(); sysres[label] = dict(dtail_rel=float(fv / fn - 1), dmean=float(v[COL][kv].mean() - n[COL][kn].mean()), n=int(kv.sum()))
    print(f"  {label:16s} high-tail fraction {100*(fv/fn-1):+.2f}%  mean shift {100*(v[COL][kv].mean()-n[COL][kn].mean()):+.3f}%")
res["detvar"] = sysres; res["global"] = dict(tail_data=fd0, tail_mc=fm0, mean_data=md0, mean_mc=mm0)
(OUT / "primary_study.json").write_text(json.dumps(res, indent=1))
