#!/usr/bin/env python3
"""Q3: collect all results into docs/syst_study/q3/q3_summary.json and make the figures."""
import sys, json, numpy as np, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from q3_native import response, Q3, BINS, REPO, OUT
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); FIG = REPO / "docs" / "syst_study" / "q3"; FIG.mkdir(parents=True, exist_ok=True)
E_ = json.load(open(Q3 / "emulation.json")); O_ = json.load(open(Q3 / "obs_repro.json")); S_ = json.load(open(Q3 / "selection.json")); K_ = json.load(open(Q3 / "keys.json")); M_ = json.load(open(Q3 / "sbnd_run1__SCE_pm" / "result.json"))
cen = 0.5 * (BINS[1:] + BINS[:-1]); C = dict(prod="k", nat="#d62728", matched="#1f77b4", ufull="#2ca02c", pair="#ff7f0e")
summary = dict(meta=dict(binning=BINS.tolist(), note="all deviations are native/emulated minus production ratio in by_vis_energy bins unless stated; chi2 uses the paired Poisson-bootstrap sigma of the ratio (B=300), ndf=17"))
# 1. emulation table
summary["emulation"] = {k: dict(stats=v["stats"], n=v["n"], mean_sigma_indep_over_paired=v["mean_sigma_ratio_ind_over_paired"], median_sigma_paired=float(np.median(v["sig_paired"])),
                               bias_unmatched_full_pot_mean=float(np.mean(np.array(v["unmatched_full_pot"]) - np.array(v["prod"]))), migration_fraction=(v["n"]["only_n"] + v["n"]["only_v"]) / (v["n"]["both"] + v["n"]["only_n"] + v["n"]["only_v"])) for k, v in E_.items()}
# 2. native E_vis results
nat = {}
def nstat(name, variant, ref):
    r, nom, inj = response(name, variant); p = np.array(E_[name][ref]); sg = np.maximum(np.array(E_[name]["sig_paired"]), 1e-12); d = r - p
    return dict(max_abs=float(np.abs(d).max()), rms=float(np.sqrt((d ** 2).mean())), chi2=float(((d / sg) ** 2).sum()), resp=r.tolist())
for name in ("sbnd_run1__ChargeScale15", "sbnd_run1__DENT", "icarus_run2__SCE2x", "icarus_run4__WireModxThetaXW"):
    nat[name] = {}
    for var, ref in (("A", "pair"), ("Bm", "matched"), ("Ce", "pair"), ("U_equal", "unmatched_full_raw"), ("U_real", "unmatched_full_pot")):
        try: nat[name][var] = dict(vs_emulation=nstat(name, var, ref)); nat[name][var]["vs_production"] = nstat(name, var, "prod")
        except Exception as e: nat[name][var] = dict(error=str(e))
for name, vs in (("sbnd_run1__ChargeScale15", ("Cz", "Cp", "fullpool_Cp", "fullpool_Ce")), ("sbnd_run1__DENT", ("Cp",)), ("icarus_run2__SCE2x", ("Cp",))):
    for var in vs:
        try: nat[name][var] = dict(vs_emulation=nstat(name, var, "pair"), vs_production=nstat(name, var, "prod"))
        except Exception as e: nat[name][var] = dict(error=str(e))
for name in nat:
    p = Q3 / name / "exact.json"
    if p.exists(): nat[name]["exact_with_repeated_keys"] = {k: v for k, v in json.load(open(p)).items() if k != "prod"}
def strip(d): return {k: strip(v) for k, v in d.items() if k != "resp"} if isinstance(d, dict) else d
summary["native_E_vis"] = strip(nat)
summary["native_runtime"] = dict(note="fullpool = 6,274,248 events on each side (ChargeScale15), matched mode", process_seconds=dict(Cp=36, Ce=38), plot_seconds_each=19, max_rss_GB=dict(Cp=3.7, Ce=3.9), pair_union_trees="2-4 s process, <0.6 GB", full_selected_trees_unmatched="4 s, 0.5 GB")
summary["multi_universe_SCE"] = {k: v for k, v in M_.items() if not isinstance(v, list)} | dict(max_dev_native_vs_numpy_commonpool=dict(m1=float(np.abs(np.array(M_["m1"]) - np.array(M_["np_m1"])).max()), p1=float(np.abs(np.array(M_["p1"]) - np.array(M_["np_p1"])).max())),
    max_dev_native_vs_production=dict(m1=float(np.abs(np.array(M_["m1"]) - np.array(M_["prod_m1"])).max()), p1=float(np.abs(np.array(M_["p1"]) - np.array(M_["prod_p1"])).max())))
summary["per_observable"] = {k: dict(stats=v["stats"], n=v["n"]) for k, v in O_.items()}
summary["selection_gain_proton_multiplicity"] = {k: {a: b for a, b in v.items() if a in ("max_abs", "rms", "chi2", "ndf")} for k, v in S_.items()}
summary["keys"] = K_
json.dump(summary, open(FIG / "q3_summary.json", "w"), indent=1)
# ---- Figure 1: ratio comparisons for 4 cases
cases = [("sbnd_run1__ChargeScale15", "SBND ChargeScale15"), ("sbnd_run1__DENT", "SBND DENT"), ("icarus_run2__SCE2x", "ICARUS Run2 SCE2x"), ("icarus_run4__WireModxThetaXW", "ICARUS Run4 WireMod x#thetaXW")]
fig, axs = plt.subplots(2, 4, figsize=(20, 8), gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.06), sharex=True)
for j, (name, lab) in enumerate(cases):
    e = E_[name]; p = np.array(e["prod"]); sg = np.array(e["sig_paired"]); ax, bx = axs[0, j], axs[1, j]
    ax.fill_between(BINS, np.r_[p - sg, (p - sg)[-1]], np.r_[p + sg, (p + sg)[-1]], step="post", color="0.85", label="paired bootstrap #pm1#sigma")
    ax.step(BINS, np.r_[p, p[-1]], where="post", color="k", lw=2, label="production ratio")
    ax.step(BINS, np.r_[e["matched"], e["matched"][-1]], where="post", color=C["matched"], lw=1.5, label="matched mode (both selected only)")
    ax.step(BINS, np.r_[e["unmatched_full_pot"], e["unmatched_full_pot"][-1]], where="post", color=C["ufull"], lw=1.5, label="unmatched, full trees, POT-scaled")
    ra, _, _ = response(name, "Ce"); ax.step(BINS, np.r_[ra, ra[-1]], where="post", color=C["nat"], lw=1.5, ls="--", label="PROfit native: selected-attribute (matched)")
    ax.set_ylabel("varied / nominal"); mark_axis(ax, lab)
    if j == 0: ax.legend(fontsize=9, loc="lower left")
    for key, col, ls, arr in (("matched", C["matched"], "-", np.array(e["matched"])), ("ufull", C["ufull"], "-", np.array(e["unmatched_full_pot"])), ("nat", C["nat"], "--", ra)):
        bx.step(BINS, np.r_[(arr - p) / sg, ((arr - p) / sg)[-1]], where="post", color=col, ls=ls, lw=1.3)
    bx.axhline(0, color="k", lw=0.8); bx.set_ylim(-6, 6) if name != "sbnd_run1__ChargeScale15" else bx.set_ylim(-60, 60); bx.set_xlabel("Reco visible energy [GeV]"); bx.set_ylabel("(method - prod) / #sigma")
fig.savefig(FIG / "q3_fig1_ratio_comparison.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# ---- Figure 2: summary over all 20 cases
names = list(E_); fig, axs = plt.subplots(1, 2, figsize=(17, 6.5))
x = np.arange(len(names)); w = 0.27
for ax, key, yl in ((axs[0], "max_abs", "max |method - production| over bins"), (axs[1], "chi2_paired", "#chi^{2} vs production (ndf = 17), paired #sigma")):
    for i, (m, c, l) in enumerate((("pair", C["pair"], "selected-attribute / pool-restricted unmatched"), ("matched", C["matched"], "matched (both selected only)"), ("unmatched_full_pot", C["ufull"], "unmatched, full trees, POT-scaled"))):
        ax.bar(x + (i - 1) * w, [E_[n]["stats"][m][key] + (1e-4 if key == "max_abs" else 0.02) for n in names], w, color=c, label=l)
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([n.replace("sbnd_run1__", "SBND ").replace("icarus_run", "ICARUS R") .replace("__", " ") for n in names], rotation=70, ha="right", fontsize=9); ax.set_ylabel(yl)
    if key == "chi2_paired": ax.axhline(17, color="k", ls=":", lw=1)
axs[0].legend(fontsize=9, loc="upper right"); mark_axis(axs[0], "Direct DetVar variants vs production, 20 detector/variation pairs")
fig.savefig(FIG / "q3_fig2_summary_all.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# ---- Figure 3: matched-mode failure vs migration fraction
fig, ax = plt.subplots(figsize=(7.5, 6)); mf = [summary["emulation"][n]["migration_fraction"] for n in names]; md = [E_[n]["stats"]["matched"]["max_abs"] for n in names]
ax.scatter(mf, md, c=["#1f77b4" if n.startswith("sbnd") else "#9467bd" for n in names], s=40)
for n, a, b in zip(names, mf, md):
    if b > 0.09 or n.endswith("DENT"): ax.annotate(n.split("__")[1], (a, b), fontsize=9, xytext=(4, 3), textcoords="offset points")
ax.set_xlabel("selection-migration fraction (CV-only + var-only)/(union rows)"); ax.set_ylabel("matched mode: max |dev| from production"); mark_axis(ax, "Matched mode loses the migration (blue SBND, purple ICARUS)")
fig.savefig(FIG / "q3_fig3_matched_failure.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# ---- Figure 4: per-observable reproduction
fig, axs = plt.subplots(2, 4, figsize=(20, 8), gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.06))
sel = [("sbnd_run1__ChargeScale15|muon_score", "SBND ChargeScale15, muon score"), ("sbnd_run1__DENT|muon_score", "SBND DENT, muon score"), ("sbnd_run1__DENT|muon_endx_sbnd", "SBND DENT, muon end-x"), ("icarus_run2__ChargeScale|muon_score", "ICARUS Run2 ChargeScale, muon score")]
for j, (k, lab) in enumerate(sel):
    d = O_[k]; ed = np.array(d["edges"]); ax, bx = axs[0, j], axs[1, j]; p = np.array(d["prod"]); sg = np.maximum(np.array(d["sig"]), 1e-9); good = np.array(d["Hn"]) > 0
    st = lambda a: np.r_[a, a[-1]]
    ax.fill_between(ed, st(p - sg), st(p + sg), step="post", color="0.85"); ax.step(ed, st(p), where="post", color="k", lw=2, label="production (per-observable file)")
    ax.step(ed, st(np.array(d["matched"])), where="post", color=C["matched"], label="matched"); ax.step(ed, st(np.array(d["unmatched_full_pot"])), where="post", color=C["ufull"], label="unmatched full, POT")
    ax.step(ed, st(np.array(d["native"]["Ce"])), where="post", color=C["nat"], ls="--", lw=1.5, label="PROfit native selected-attribute")
    lo, hi = np.percentile(p[good], [0, 100]); ax.set_ylim(min(lo, 1) - 0.15, hi + 0.15); ax.set_ylabel("varied / nominal"); mark_axis(ax, lab)
    for arr, col, ls in ((np.array(d["matched"]), C["matched"], "-"), (np.array(d["unmatched_full_pot"]), C["ufull"], "-"), (np.array(d["native"]["Ce"]), C["nat"], "--")):
        bx.step(ed, st(np.where(good, (arr - p) / sg, 0)), where="post", color=col, ls=ls)
    bx.axhline(0, color="k", lw=0.8); bx.set_ylim(-8, 8); bx.set_ylabel("(method - prod)/#sigma"); bx.set_xlabel("observable")
    if j == 0: ax.legend(fontsize=8, loc="best")
fig.savefig(FIG / "q3_fig4_per_observable.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# ---- Figure 5: gain from a stricter analysis selection
fig, axs = plt.subplots(1, 3, figsize=(18, 5.5), sharey=False)
for ax, name, lab in zip(axs, ("sbnd_run1__ChargeScale15", "sbnd_run1__GainHigh", "icarus_run2__ChargeScale"), ("SBND ChargeScale15", "SBND GainHigh", "ICARUS Run2 ChargeScale")):
    d = S_[f"{name}|pm==1"]; sg = np.array(d["sig"]); r = np.array(d["ratio"]); p = np.array(d["prod"]); st = lambda a: np.r_[a, a[-1]]
    ax.fill_between(BINS, st(r - sg), st(r + sg), step="post", color="#d62728", alpha=0.25); ax.step(BINS, st(r), where="post", color="#d62728", lw=1.8, label="direct response, reco_proton_multiplicity == 1")
    ax.step(BINS, st(p), where="post", color="k", lw=2, label="production inclusive ratio (what hist1d applies)"); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("varied / nominal"); mark_axis(ax, lab)
    if name == "sbnd_run1__ChargeScale15": ax.legend(fontsize=8, loc="lower left")
fig.savefig(FIG / "q3_fig5_selection_gain.png", dpi=130, bbox_inches="tight"); plt.close(fig)
# ---- Figure 6: independent-sample noise vs paired noise
fig, ax = plt.subplots(figsize=(8, 5.5)); ax.bar(x, [E_[n]["mean_sigma_ratio_ind_over_paired"] for n in names], color="#7f7f7f"); ax.set_xticks(x); ax.set_xticklabels([n.replace("sbnd_run1__", "SBND ").replace("icarus_run", "ICARUS R").replace("__", " ") for n in names], rotation=70, ha="right", fontsize=9)
ax.set_ylabel("#sigma_{independent}/#sigma_{paired} (bin average)"); mark_axis(ax, "Statistical noise if CV and variation were independent samples")
fig.savefig(FIG / "q3_fig6_noise_ratio.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("done")
