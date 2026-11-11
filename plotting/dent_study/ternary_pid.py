#!/usr/bin/env python3
"""Score-simplex (muon, pion, proton) views of the selected SBND muon candidates, using the ternary tooling of /Users/mueller/Projects/spineosc_pid.
(1) Scatter of the nominal selected muon candidates, signal and background separately, with the argmax (alpha=1) and 2x-boosted (alpha=2) decision boundaries.
(2) Paired migration field of the candidates between the nominal and each *active* SBND detector variation (charge scale, WireMod X-thetaxW, non-uniformity, DENT).
(3) Fractions of candidates in the boost band (muon score wins only with the 2x boost), nominal vs variations, and the rate in the band by truth class.
Run: python plotting/dent_study/ternary_pid.py        Outputs: docs/ternary_pid/*"""
import json, pathlib, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/mueller/Projects/spineosc_pid"); sys.path.insert(0, str(pathlib.Path(__file__).parent))
from pidtools import ternary as T
import cathode_scan as cs
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
OUT = pathlib.Path(__file__).resolve().parents[2] / "docs" / "ternary_pid"; OUT.mkdir(parents=True, exist_ok=True)
S = ["reco_leading_muon_muon_softmax", "reco_leading_muon_pion_softmax", "reco_leading_muon_proton_softmax"]
KEYS = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
COLS = KEYS + S + ["true_cc", "true_pdg", "true_category", "reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", "event_largest_flash_pe"]
VARS = [("Charge scale", "sbnd_detvar_smear15", "sbnd_large"), ("WireMod X-$\\theta_{xW}$", "sbnd_detvar_wiremod_xthetaxw", "sbnd_large"),
        ("Non-uniformity", "sbnd_detvar_mattmod_yz", "sbnd_large"), ("DENT", "sbnd_detvar_newdent", "newcv")]
def load(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="np")
    ok = (d["reco_proton_multiplicity"] >= 1) & (d["reco_fiducialize_cathode"] == 1) & (d["reco_veto_sbnd_highy_highz"] == 1) & (d["reco_cathode_crosser"] == 0) & (np.abs(d["reco_leading_muon_end_x"]) >= 4) & (0.647 * d["event_largest_flash_pe"] > 2000)
    d = {k: v[ok].astype(float) for k, v in d.items()}
    d["sig"] = (d["true_cc"] > 0) & (np.abs(d["true_pdg"]) == 14) & np.isin(d["true_category"], [0, 2])
    d["key"] = list(zip(d["Run"].astype(int), d["Subrun"].astype(int), d["Evt"].astype(int), np.round(d["event_max_neutrino_energy"], 6)))
    return d
def scores(d, m=None):
    a = [d[s] if m is None else d[s][m] for s in S]; return a
def band(d, lo=1.0, hi=2.0):
    m, p, q = d[S[0]], d[S[1]], d[S[2]]; tot = m + p + q
    return (lo * m <= np.maximum(p, q)) & (hi * m > np.maximum(p, q))   # muon wins only with the boost
def scatter(ax, sc, n=30000, seed=3):
    x, y, good = T.to_ternary(*sc); rng = np.random.default_rng(seed); i = rng.choice(len(x), min(n, len(x)), replace=False)
    gx, gy, inside = T.triangle_grid(130); dens = T.density(x, y, gx, gy, inside, 0.018); v = T.density_at(x[i], y[i], dens, (gx, gy, inside)) if True else None
    o = np.argsort(v); ax.scatter(x[i][o], y[i][o], c=v[o], s=0.8, cmap="viridis", rasterized=True, linewidths=0)
def frame(ax, boundaries=(1.0, 2.0)):
    T.draw_frame(ax)
    for a, ls in zip(boundaries, ("-", "--")):
        bx, by = T.boundary(a); ax.plot(bx, by, color="#D55E00", lw=1.8, ls=ls)
if __name__ == "__main__":
    nom = load(cs.NOS, "sbnd_large"); ncv = load(cs.DET, "sbnd_detvar_newcv"); res = {}
    fig, axs = plt.subplots(1, 2, figsize=(12, 6.3))
    for ax, (nm, m) in zip(axs, (("Signal", nom["sig"]), ("Background", ~nom["sig"]))):
        scatter(ax, scores(nom, m)); frame(ax); ax.text(0.02, 0.97, nm, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")
    axs[0].plot([], [], color="#D55E00", lw=1.8, label=r"$\alpha=1$ (argmax)"); axs[0].plot([], [], color="#D55E00", ls="--", lw=1.8, label=r"$\alpha=2$ (muon score $\times$2)"); axs[1].legend(loc="upper right", fontsize=11, frameon=False)
    mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[1], "SBND Run 1", alignment="right"); fig.tight_layout()
    fig.savefig(OUT / "ternary_scatter_sbnd.png", dpi=130); fig.savefig(OUT / "ternary_scatter_sbnd.pdf"); plt.close(fig)
    b = band(nom); res["nominal"] = dict(band_frac_signal=float(b[nom["sig"]].mean()), band_frac_background=float(b[~nom["sig"]].mean()),
                                         signal_share_in_band=float(nom["sig"][b].mean()), signal_share_outside=float(nom["sig"][~b].mean()), n=int(len(b)))
    fig, axs = plt.subplots(2, 2, figsize=(12, 11))
    for ax, (label, tree, nk) in zip(axs.ravel(), VARS):
        v = load(cs.DET, tree); n = nom if nk == "sbnd_large" else ncv
        kv = {k: i for i, k in enumerate(v["key"])}; pairs = [(i, kv[k]) for i, k in enumerate(n["key"]) if k in kv]
        ia = np.array([p[0] for p in pairs]); ib = np.array([p[1] for p in pairs]); sig = n["sig"][ia]
        scatter(ax, scores(n, ia), n=12000); frame(ax)
        tx, ty, dx, dy, w = T.displacement_field(scores(n, ia), scores(v, ib), bins=9, minimum=40)
        ax.quiver(tx, ty, dx, dy, angles="xy", scale_units="xy", scale=1.0, width=0.006, color="k", zorder=5)
        ax.text(0.02, 0.97, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")
        bn = band(n)[ia]; bv = band(v)[ib]; mus_n = (n[S[0]][ia] * 2 > np.maximum(n[S[1]][ia], n[S[2]][ia])); mus_v = (v[S[0]][ib] * 2 > np.maximum(v[S[1]][ib], v[S[2]][ib]))
        res[label] = dict(n_pairs=int(len(ia)), boost_band_nom=float(bn.mean()), boost_band_var=float(bv.mean()), muon_region_boost_nom=float(mus_n.mean()), muon_region_boost_var=float(mus_v.mean()),
                          muon_region_argmax_nom=float((n[S[0]][ia] > np.maximum(n[S[1]][ia], n[S[2]][ia])).mean()), muon_region_argmax_var=float((v[S[0]][ib] > np.maximum(v[S[1]][ib], v[S[2]][ib])).mean()),
                          mean_mu_nom=float(n[S[0]][ia].mean()), mean_mu_var=float(v[S[0]][ib].mean()))
        print(label, json.dumps({k: round(x, 4) if isinstance(x, float) else x for k, x in res[label].items()}))
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 1], "SBND Run 1", alignment="right"); fig.tight_layout()
    fig.savefig(OUT / "ternary_migration_sbnd.png", dpi=120); fig.savefig(OUT / "ternary_migration_sbnd.pdf"); plt.close(fig)
    (OUT / "ternary_pid.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res["nominal"]))
