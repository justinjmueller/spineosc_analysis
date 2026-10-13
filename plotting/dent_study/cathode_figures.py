#!/usr/bin/env python3
"""Polished figures motivating the SBND cathode-crosser and 4 cm muon end-point cuts. No track length or proxy is shown; only crosser fractions and
the distance of the end point from the cathode, |x_end|. Area-normalized quantities only. Run (mutools environment):
    python plotting/dent_study/cathode_figures.py [name ...]      # default: all
Inputs: docs/cathode_motivation/scan.json (cathode_scan.py), the dv_new CSV dumps (plotting/dent_validation/dump_samples.C), detsys files in detvar/.
Outputs: docs/cathode_motivation/*.pdf"""
import json, pathlib, sys
import numpy as np, pandas as pd, uproot
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mutools.plotting import use_style
from mutools.plotting.helpers import mark_axis
REPO = pathlib.Path(__file__).resolve().parents[2]; OUT = REPO / "docs" / "cathode_motivation"; OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "plotting" / "dent_validation")); sys.path.insert(0, str(pathlib.Path(__file__).parent))
use_style("rootlike")
CUT = 4.0
COL = {"data": "k", "cv": "#0072B2", "dent": "#E69F00"}; LAB = {"data": "Data", "cv": "Simulation (CV)", "dent": "Simulation (DENT)"}
STAMP = "medulla v1.1.0"

def stamp(ax, left=r"$\bf{SBN}$ Analysis-in-Progress", right="SBND Simulation", size=8):
    mark_axis(ax, left)
    if right: mark_axis(ax, right, alignment="right")

def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf"); fig.savefig(OUT / f"{name}.png", dpi=130); plt.close(fig); print("wrote", name)

def ratio_axes(figsize=(8, 7)):
    fig, (a, r) = plt.subplots(2, 1, figsize=figsize, sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05)); return fig, a, r

def samples():
    import dentval as dv
    return dv.load()

def hist_fraction(df, col, edges, mask):
    """fraction of ALL selected events of the sample per unit of `col` (per cm), with Poisson errors"""
    N = len(df); c = np.histogram(df.loc[mask, col].to_numpy(), edges)[0].astype(float); w = np.diff(edges)
    return c / N / w, np.sqrt(np.maximum(c, 1)) / N / w

# ------------------------------------------------------------------------------------------------ stage 2: data/MC near the cathode
def fig_endx_distance():
    d = samples(); e = np.arange(0, 42, 2.0); ctr = 0.5 * (e[1:] + e[:-1])
    fig, a, r = ratio_axes()
    P = {k: hist_fraction(d[k], "abs_end_x", e, (d[k].reco_cathode_crosser == 0).to_numpy()) for k in d}
    for k in ("cv", "dent"): a.stairs(P[k][0], e, color=COL[k], lw=2, label=LAB[k])
    a.errorbar(ctr, P["data"][0], yerr=P["data"][1], fmt="o", color="k", ms=4.5, label=LAB["data"], zorder=5)
    base = np.where(P["cv"][0] > 0, P["cv"][0], np.nan)
    r.stairs(P["dent"][0] / base, e, color=COL["dent"], lw=2); r.axhline(1, color=COL["cv"], lw=2)
    r.errorbar(ctr, P["data"][0] / base, yerr=P["data"][1] / base, fmt="o", color="k", ms=4.5, zorder=5)
    for ax in (a, r): ax.axvline(CUT, color="0.45", ls="--", lw=1.2)
    a.text(CUT + 0.6, 0.92 * a.get_ylim()[1], f"{CUT:.0f} cm", color="0.35", fontsize=11)
    a.set_ylabel("Fraction of selected events / cm"); r.set_ylabel("Ratio to CV"); r.set_ylim(0.4, 2.4); r.set_xlabel(r"Muon end point distance from the cathode $|x_{end}|$ [cm]")
    a.set_xlim(0, 40); a.set_ylim(0, None); a.legend(loc="upper right", fontsize=12); stamp(a, right="SBND, non-crossers"); a.text(0.98, 0.64, STAMP, transform=a.transAxes, ha="right", fontsize=8, color="0.45")
    save(fig, "cathode_endx_distance")

def fig_fractions():
    J = json.load(open(OUT / "scan.json"))["fractions"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (key, title) in zip(axs, (("crosser_fraction", "Flagged as cathode crossers"), ("4", r"Non-crossers with $|x_{end}| < 4$ cm"))):
        for i, k in enumerate(("data", "cv", "dent")):
            v, err = J[k][key]; ax.errorbar([i], [100 * v], yerr=[100 * err], fmt="o" if k == "data" else "s", color=COL[k], ms=9, capsize=3, lw=1.8)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["Data", "CV", "DENT"]); ax.set_xlim(-0.6, 2.6); ax.set_ylabel("Fraction of selected events [%]"); ax.text(0.04, 0.92, title, transform=ax.transAxes, fontsize=13, va="top")
        ax.set_ylim(0, 1.25 * max(100 * J[k][key][0] for k in J if isinstance(J[k], dict) and key in J[k]))
    stamp(axs[0], right=None); mark_axis(axs[1], "SBND", alignment="right")
    fig.tight_layout(); save(fig, "cathode_fractions")

def fig_after_cuts():
    d = samples(); e = np.arange(-40, 42, 2.0); ctr = 0.5 * (e[1:] + e[:-1])
    fig, axs = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0.05, wspace=0.22))
    for j, (col, xl) in enumerate((("reco_leading_muon_end_x", r"Muon end point $x$ [cm]"), ("reco_vertex_x", r"Interaction vertex $x$ [cm]"))):
        a, r = axs[0, j], axs[1, j]; P = {}
        for k, df in d.items():
            keep = (df.reco_cathode_crosser == 0).to_numpy() & (df.abs_end_x >= CUT).to_numpy(); sub = df[keep]
            c = np.histogram(sub[col].to_numpy(), e)[0].astype(float); n = c.sum(); w = np.diff(e); P[k] = (c / n / w, np.sqrt(np.maximum(c, 1)) / n / w)
        for k in ("cv", "dent"): a.stairs(P[k][0], e, color=COL[k], lw=2, label=LAB[k] if j == 0 else None)
        a.errorbar(ctr, P["data"][0], yerr=P["data"][1], fmt="o", color="k", ms=4, label=LAB["data"] if j == 0 else None, zorder=5)
        base = np.where(P["cv"][0] > 0, P["cv"][0], np.nan)
        r.stairs(P["dent"][0] / base, e, color=COL["dent"], lw=2); r.axhline(1, color=COL["cv"], lw=2); r.errorbar(ctr, P["data"][0] / base, yerr=P["data"][1] / base, fmt="o", color="k", ms=4, zorder=5)
        a.set_ylabel("Fraction of selected events / cm"); r.set_ylabel("Ratio to CV"); r.set_ylim(0.6, 1.4); r.set_xlabel(xl); a.set_ylim(0, None)
    h, l = axs[0, 0].get_legend_handles_labels(); axs[0, 1].legend(h, l, loc="lower left", fontsize=12); stamp(axs[0, 0], right=None); mark_axis(axs[0, 1], r"SBND, crossers and $|x_{end}|<4$ cm removed", alignment="right")
    save(fig, "cathode_after_cuts")

# ------------------------------------------------------------------------------------------------ stage 1: DENT systematic
def fig_dent_ratio():
    J = json.load(open(OUT / "scan.json")); e = {"1p": np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13, 1.2, 1.3, 1.5]),
         "Mp": np.array([0.4, 0.5, 0.59, 0.67, 0.75, 0.83, 0.91, 0.99, 1.07, 1.16, 1.26, 1.37, 1.5, 1.7, 2.0])}
    sel = ((-1, "No crosser cut", "#009E73"), (0, "Crosser cut", "#0072B2"), (4, r"Crosser cut + $|x_{end}| \geq 4$ cm", "#D55E00"))
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.4), sharey=True)
    for ax, ch, title in zip(axs, ("1p", "Mp"), ("Single-proton", "Multi-proton")):
        ctr = 0.5 * (e[ch][1:] + e[ch][:-1])
        for X, lab, c in sel:
            o = next(o for o in J["channels"][ch] if o["X"] == X); y = 100 * np.array(o["fine"]); err = 100 * np.array(o["fine_err"])
            ax.stairs(y, e[ch], color=c, lw=2, label=lab); ax.errorbar(ctr, y, yerr=err, fmt="none", ecolor=c, elinewidth=1, alpha=0.55)
        ax.axhline(0, color="0.4", lw=1); ax.set_xlabel("Reconstructed Visible Energy [GeV]"); ax.set_ylim(-22, 16); stamp(ax, right=f"SBND Simulation, {title}")
    axs[0].set_ylabel("DENT / CV $-$ 1 [%]"); axs[0].legend(loc="lower left", fontsize=12); save(fig, "dent_ratio_vs_energy")

def fig_dent_uncertainty():
    files = (("detsys-v1.1.0v0.root", "No crosser cut", "#009E73"), ("detsys-v1.1.0v1.root", "Crosser cut", "#0072B2"), ("detsys-v1.1.0v1_endx4.root", r"Crosser cut + $|x_{end}| \geq 4$ cm", "#D55E00"))
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for fn, lab, c in files:
        h = uproot.open(REPO / "detvar" / fn)["DENT_sbnd_Run1"]; ed = h.axis().edges(); v = np.abs(h.values() - 1)
        ax.stairs(v, ed, color=c, lw=2, label=lab)
    ax.set_xlabel("Reconstructed Visible Energy [GeV]"); ax.set_ylabel("DENT fractional uncertainty"); ax.set_xlim(0.3, 2.0); ax.set_ylim(0, None); ax.legend(loc="upper left", fontsize=12)
    stamp(ax, right="SBND Simulation"); save(fig, "dent_fractional_uncertainty")

def pid_impact():
    """matched detvar CV/DENT events (fiducial + veto), DENT/CV - 1 per PID-score bin for crossers and non-crossers, paired Poisson bootstrap"""
    import cathode_scan as cs
    D = cs.DET; cols = cs.K + ["reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_muon_softmax", "reco_leading_muon_primary_softmax", "reco_leading_proton_proton_softmax", "reco_leading_proton_primary_softmax"]
    def load(tree):
        d = uproot.open(D)[f"events/{tree}/selected"].arrays(cols, library="pd").drop_duplicates(cs.K, keep=False).set_index(cs.K)
        return d[(d.reco_fiducialize_cathode == 1) & (d.reco_veto_sbnd_highy_highz == 1)]
    common = cs.keys(D, "sbnd_detvar_newcv").intersection(cs.keys(D, "sbnd_detvar_newdent"))
    cv = load("sbnd_detvar_newcv"); dt = load("sbnd_detvar_newdent"); cv = cv[cv.index.isin(common)]; dt = dt[dt.index.isin(common)]
    ev = cv.index.union(dt.index); c = cv.reindex(ev); t = dt.reindex(ev); n = len(ev)
    rng = np.random.default_rng(20261007); NB = 300; W = rng.poisson(1.0, size=(NB + 1, n)).astype(np.float32); W[0] = 1.0
    edges = np.linspace(0, 1, 21); out = {}
    for var in ("reco_leading_muon_muon_softmax", "reco_leading_muon_primary_softmax", "reco_leading_proton_proton_softmax", "reco_leading_proton_primary_softmax"):
        out[var] = {}
        for nm, fc in (("nocut", None), ("crossers", 1), ("noncrossers", 0)):
            vals = {}
            for s, d in (("cv", c), ("dent", t)):
                x = d[var].to_numpy(); ok = ~np.isnan(x); ok &= (np.ones(n, bool) if fc is None else (d.reco_cathode_crosser.to_numpy() == fc))
                b = np.digitize(x, edges) - 1; b[(b < 0) | (b >= 20)] = -1; b[~ok] = -1
                vals[s] = np.stack([W[:, b == i].sum(axis=1) for i in range(20)], axis=1)
            ratio = vals["dent"] / np.maximum(vals["cv"], 1) - 1
            out[var][nm] = (ratio[0], ratio[1:].std(axis=0), vals["cv"][0])
    return edges, out

def fig_pid_impact():
    edges, P = pid_impact(); ctr = 0.5 * (edges[1:] + edges[:-1])
    names = {"reco_leading_muon_muon_softmax": "Muon PID score", "reco_leading_muon_primary_softmax": "Muon primary score",
             "reco_leading_proton_proton_softmax": "Leading-proton PID score", "reco_leading_proton_primary_softmax": "Leading-proton primary score"}
    fig, axs = plt.subplots(2, 2, figsize=(13, 9), sharey=False)
    sel = (("nocut", "All events", "k"), ("crossers", "Cathode crossers", "#D55E00"), ("noncrossers", "Non-crossers", "#0072B2"))
    for ax, (var, title) in zip(axs.ravel(), names.items()):
        for k, lab, c in sel:
            y, err, ncv = P[var][k]; ok = ncv > 30       # drop bins too empty to carry a ratio
            ax.stairs(np.where(ok, 100 * y, np.nan), edges, color=c, lw=2, label=lab); ax.errorbar(ctr[ok], 100 * y[ok], yerr=100 * err[ok], fmt="none", ecolor=c, elinewidth=1, alpha=0.5)
        ax.axhline(0, color="0.4", lw=1); ax.set_xlabel(title); ax.set_ylabel("DENT / CV $-$ 1 [%]"); pop = np.where(P[var]["nocut"][2] > 30)[0]; ax.set_xlim(edges[max(pop[0] - 1, 0)], 1)
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=3, fontsize=13, frameon=False, bbox_to_anchor=(0.5, 1.0)); stamp(axs[0, 0], right=None); mark_axis(axs[0, 1], "SBND Simulation", alignment="right")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); save(fig, "dent_pid_impact")

def fig_scan():
    J = json.load(open(OUT / "scan.json")); X = np.array([x for x in J["X"] if x >= 0])
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.2))
    a = axs[0, 0]
    for ch, c in (("1p", "#0072B2"), ("Mp", "#D55E00")):
        o = [o for o in J["channels"][ch] if o["X"] >= 0]; a.errorbar(X, [100 * q["integrated"] for q in o], yerr=[100 * q["integrated_err"] for q in o], fmt="o-", color=c, label={"1p": "Single-proton", "Mp": "Multi-proton"}[ch], lw=1.8)
    a.axhspan(-1.5, 1.5, color="0.85", zorder=0); a.set_ylabel("Integrated DENT / CV $-$ 1 [%]"); a.legend(fontsize=11, loc="upper right")
    b = axs[0, 1]; bx = [x for x in X]; b.errorbar(bx, [J["band"][str(x)]["ratio"] for x in bx], yerr=[J["band"][str(x)]["err"] for x in bx], fmt="o-", color="#E69F00", lw=1.8)
    b.axhspan(0.9, 1.1, color="0.85", zorder=0); b.axhline(1, color="0.4", lw=1); b.set_ylabel("DENT / CV in the next 4 cm band")
    c_ = axs[1, 0]
    for k in ("data", "cv", "dent"):
        v = [100 * J["fractions"][k][str(x)][0] for x in X]; er = [100 * J["fractions"][k][str(x)][1] for x in X]; c_.errorbar(X, v, yerr=er, fmt="o-" if k == "data" else "s-", color=COL[k], label=LAB[k], lw=1.6)
    c_.set_ylabel(r"Non-crossers, $|x_{end}| < X$ [%]"); c_.legend(fontsize=11, loc="upper left")
    d_ = axs[1, 1]; d_.plot(X, [100 * J["cost_reco"][str(x)] for x in X], "o-", color="#0072B2", label="Reconstructed", lw=1.8); d_.plot(X, [100 * J["cost_true"][str(x)] for x in X], "s--", color="#009E73", label="True", lw=1.8)
    d_.axhline(3, color="0.4", ls=":", lw=1.2); d_.set_ylabel("Signal removed [%]"); d_.legend(fontsize=11, loc="upper left")
    for ax in axs.ravel(): ax.axvline(CUT, color="0.45", ls="--", lw=1.2); ax.set_xlabel(r"End-point cut value $X$ [cm]")
    stamp(axs[0, 0], right=None); mark_axis(axs[0, 1], "SBND Simulation", alignment="right"); fig.tight_layout(); save(fig, "cathode_cut_scan")

ALL = dict(endx_distance=fig_endx_distance, fractions=fig_fractions, after_cuts=fig_after_cuts, dent_ratio=fig_dent_ratio, dent_uncertainty=fig_dent_uncertainty, pid_impact=fig_pid_impact, scan=fig_scan)
if __name__ == "__main__":
    for k in (sys.argv[1:] or ALL): ALL[k]()
