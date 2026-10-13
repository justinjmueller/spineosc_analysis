#!/usr/bin/env python3
"""Statistical resolution of the SBND ChargeScale15 detector ratio at high reconstructed energy (selection: fiducial, veto, crosser == 0, |end x| >= 4 cm).
Compares the ratio in the calculator binning (single 1.5-2.0 GeV bin) with finer bins and a smooth fit, with paired Poisson-bootstrap errors.
Writes docs/cathode_motivation/chargescale_highE.json and .pdf. Run in the mutools environment."""
import json, pathlib, sys
import numpy as np, pandas as pd, uproot
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import cathode_scan as cs
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
REPO = cs.REPO; OUT = REPO / "docs" / "cathode_motivation"
cols = cs.K + ["reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x"]
def load(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(cols, library="pd").drop_duplicates(cs.K, keep=False).set_index(cs.K)
    return d[(d.reco_fiducialize_cathode == 1) & (d.reco_veto_sbnd_highy_highz == 1) & (d.reco_cathode_crosser == 0) & (d.reco_leading_muon_end_x.abs() >= 4)]
common = cs.keys(cs.NOS, "sbnd_large").intersection(cs.keys(cs.DET, "sbnd_detvar_smear15"))
cv = load(cs.NOS, "sbnd_large"); sm = load(cs.DET, "sbnd_detvar_smear15"); cv = cv[cv.index.isin(common)]; sm = sm[sm.index.isin(common)]
ev = cv.index.union(sm.index); n = len(ev); Ec = cv.reindex(ev).reco_visible_energy.to_numpy(); Es = sm.reindex(ev).reco_visible_energy.to_numpy()
rng = np.random.default_rng(20261007); NB = 400; W = rng.poisson(1.0, size=(NB + 1, n)).astype(np.float32); W[0] = 1.0
def ratio(edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        c = W[:, (Ec >= lo) & (Ec < hi)].sum(axis=1); s = W[:, (Es >= lo) & (Es < hi)].sum(axis=1); r = s / c
        out.append((r[0], r[1:].std(), c[0]))
    return out
res = {}
res["calculator_bins"] = ratio([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5, 2.0])
res["split_bins"] = ratio([1.2, 1.3, 1.5, 1.7, 2.0])
# smooth: weighted linear fit of the calculator-bin ratio above 0.7 GeV, bootstrap for the band
e = np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5, 2.0]); ctr = 0.5 * (e[1:] + e[:-1]); sel = ctr > 0.7
def fit(k):
    c = np.array([W[k, (Ec >= lo) & (Ec < hi)].sum() for lo, hi in zip(e[:-1], e[1:])]); s = np.array([W[k, (Es >= lo) & (Es < hi)].sum() for lo, hi in zip(e[:-1], e[1:])])
    r = s / c; w = c[sel]; return np.polyfit(ctr[sel], r[sel], 1, w=np.sqrt(w))
coef = [fit(k) for k in range(NB + 1)]; xs = np.linspace(0.7, 2.0, 60)
band = np.array([np.polyval(c, xs) for c in coef[1:]]); res["fit"] = dict(coef=coef[0].tolist(), x=xs.tolist(), y=np.polyval(coef[0], xs).tolist(), err=band.std(axis=0).tolist())
(OUT / "chargescale_highE.json").write_text(json.dumps(res, indent=1, default=float))
fig, ax = plt.subplots(figsize=(8.5, 6))
cb = res["calculator_bins"]; ax.stairs([b[0] for b in cb], e, color="#0072B2", lw=2, label="Ratio as in the detector systematic", baseline=None)
ax.errorbar(ctr, [b[0] for b in cb], yerr=[b[1] for b in cb], fmt="none", ecolor="#0072B2", elinewidth=1, alpha=0.6)
sp = res["split_bins"]; se = np.array([1.2, 1.3, 1.5, 1.7, 2.0]); ax.stairs([b[0] for b in sp], se, color="#D55E00", lw=2, ls="--", label="Split at 1.7 GeV", baseline=None)
ax.errorbar(0.5 * (se[1:] + se[:-1]), [b[0] for b in sp], yerr=[b[1] for b in sp], fmt="none", ecolor="#D55E00", elinewidth=1, alpha=0.6)
f = res["fit"]; ax.plot(f["x"], f["y"], color="k", lw=1.6, label="Linear fit above 0.7 GeV"); ax.fill_between(f["x"], np.array(f["y"]) - np.array(f["err"]), np.array(f["y"]) + np.array(f["err"]), color="0.8", zorder=0)
ax.axhline(1, color="0.4", lw=1); ax.set_ylim(0.78, 1.06); ax.set_xlim(0.3, 2.0); ax.set_xlabel("Reconstructed Visible Energy [GeV]"); ax.set_ylabel("Charge scale / smearing: variation / CV"); ax.legend(loc="lower left", fontsize=12)
mark_axis(ax, r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(ax, "SBND Simulation", alignment="right")
fig.savefig(OUT / "chargescale_highE.pdf"); fig.savefig(OUT / "chargescale_highE.png", dpi=130)
print("calculator last bins:", [(round(b[0], 3), round(b[1], 3), int(b[2])) for b in cb[-3:]]); print("split:", [(round(b[0], 3), round(b[1], 3), int(b[2])) for b in sp]); print("fit coef (slope, intercept):", np.round(coef[0], 3))
