#!/usr/bin/env python3
"""Q2 figures (docs/syst_study/q2/*.png)."""
import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from spineosc_paths import REPO  # noqa: E402
from q2_lib import *
from q2_fit import *
from q2_smooth import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
WM = r"$\bf{SBN}$ Analysis-in-Progress"; CC = {"i": "#0072B2", "ii": "#E69F00", "iii": "#999999"}; DC = {"sbnd_run1": "#0072B2", "icarus_run2": "#D55E00", "icarus_run4": "#009E73"}
summ = json.load(open(Q2 / "q2_summary.json")); cal = json.load(open(Q2 / "q2_calibration.json")); cov = json.load(open(Q2 / "q2_covariance.json")); cons = json.load(open(Q2 / "q2_consistency.json")); pool = json.load(open(Q2 / "q2_pool.json")); val = json.load(open(Q2 / "q2_validation.json"))
sn = json.load(open(Q2 / "q2_signal_noise.json")); short = lambda n: n.replace("sbnd_run1__", "SBND ").replace("icarus_run2__", "IC2 ").replace("icarus_run4__", "IC4 ").replace("clean", " (clean)")
DATA = {}
for n in ALL:
    d = load(n); r, C, N, V = cov_analytic(d); y = r - 1; f, sel = Estimator(C, LOGX, mode="poly", last="raw", dmax=4)(y[None, :]); fw, _ = Estimator(C, LOGX, mode="whit", last="raw")(y[None, :]); DATA[n] = dict(y=y, C=C, s=np.sqrt(np.diag(C)), f=f[0], fw=fw[0], N=N)
# ---- Fig 1: calibration
fig, axs = plt.subplots(2, 2, figsize=(13, 10))
ax = axs[0, 0]
for det in DC:
    ns = [n for n in NAMES if det_of(n) == det]; a = np.concatenate([cal[n]["sig_an"] for n in ns]); b = np.concatenate([cal[n]["sig_bs"] for n in ns]); ax.plot(100 * a, 100 * b, ".", color=DC[det], label=det.replace("_", " "), ms=5)
ax.plot([0.1, 5], [0.1, 5], "k-", lw=1); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"analytic $\sigma(r_b)$ [%]"); ax.set_ylabel(r"paired Poisson bootstrap $\sigma(r_b)$ [%]"); ax.legend(loc="upper left"); mark_axis(ax, WM)
ax = axs[0, 1]; x = np.arange(len(NAMES)); ax.axhline(1, color="k", lw=1)
ax.errorbar(x - 0.15, [cal[n]["split_event"]["chi2ndf_mean"] for n in NAMES], yerr=[0.35 / np.sqrt(8) for n in NAMES], fmt="o", color="#0072B2", label="split by event hash"); ax.plot(x + 0.15, [cal[n]["split_block"]["chi2ndf_mean"] for n in NAMES], "s", color="#D55E00", label="split by Run (SBND) / (Run,Subrun) (ICARUS)")
ax.set_xticks(x); ax.set_xticklabels([short(n) for n in NAMES], rotation=90, fontsize=8); ax.set_ylabel(r"split-half $\chi^2$/ndf of $(r_A-r_B)$"); ax.set_ylim(0.6, 1.5); ax.legend(loc="upper left", fontsize=9)
ax = axs[1, 0]; ax.axhline(1, color="k", lw=1); ax.plot(x, [cal[n]["block"]["jk_med"] for n in NAMES], "o", color="#0072B2", label="delete-one-block jackknife / event-level"); ax.plot(x, [cal[n]["block"]["bb_med"] for n in NAMES], "s", color="#D55E00", label="block bootstrap / event-level")
bc = json.load(open(Q2 / "q2_block_clean.json")); ax.plot([NAMES.index("icarus_run4__WireModxThetaXW")], [bc["icarus_run4__WireModxThetaXWclean"]["jk_over_ev_median"]], "*", color="#009E73", ms=14, label="Run4 WireModxThetaXW without dup. runs")
ax.set_xticks(x); ax.set_xticklabels([short(n) for n in NAMES], rotation=90, fontsize=8); ax.set_ylabel(r"median$_b$ $\sigma_{\rm block}/\sigma_{\rm event}$"); ax.set_ylim(0.9, 1.25); ax.legend(loc="upper left", fontsize=8)
ax = axs[1, 1]; ax.bar(x, [cal[n]["naive_over_paired_median"] for n in NAMES], color=[DC[det_of(n)] for n in NAMES]); ax.set_xticks(x); ax.set_xticklabels([short(n) for n in NAMES], rotation=90, fontsize=8); ax.set_ylabel(r"median $\sigma_{\rm independent\ samples}/\sigma_{\rm paired}$")
fig.tight_layout(); fig.savefig(Q2 / "q2_fig1_calibration.png", dpi=110); plt.close(fig)
# ---- Fig 2: shapes
def shapes(names, fn, nc=3, title=""):
    nr = int(np.ceil(len(names) / nc)); fig, axs = plt.subplots(nr, nc, figsize=(5.2 * nc, 3.9 * nr), squeeze=False)
    for ax, n in zip(axs.ravel(), names):
        D = DATA[n]; p = summ["per_variation"][n]
        ax.axhline(0, color="0.5", lw=0.8); ax.errorbar(CTR, 100 * D["y"], yerr=100 * D["s"], fmt="o", color="k", ms=4, lw=1, label=r"$r-1$ $\pm\,\sigma$ (analytic)")
        ax.plot(CTR[:-1], 100 * D["f"][:-1], "-", color="#D55E00", lw=2, label="smooth fit (poly-BIC, bins 1-16)"); ax.plot(CTR[:-1], 100 * D["fw"][:-1], "--", color="#0072B2", lw=1.3, label="Whittaker-CV (alt.)"); ax.plot([CTR[-1]], [100 * D["f"][-1]], "D", color="#D55E00", ms=6, mfc="none", label="recommended last-bin value")
        ax.text(0.03, 0.97, short(n) + "\n" + {"i": "(i) smooth trend", "ii": "(ii) constant", "iii": "(iii) zero/noise"}[p["cls"]] + rf"; $\chi^2_{{0}}$/17={p['chi2_zero']/17:.1f}", transform=ax.transAxes, va="top", fontsize=9, fontweight="bold")
        lim = np.max(np.abs(100 * D["y"][:-1])) * 1.4 + 0.5; ax.set_ylim(-lim, lim); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.text(0.97, 0.04, f"last bin: {100*D['y'][-1]:+.1f} $\\pm$ {100*D['s'][-1]:.1f} %", transform=ax.transAxes, ha="right", fontsize=8)
    h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="upper right", fontsize=9, ncol=2); mark_axis(axs[0, 0], WM)
    for ax in axs.ravel()[len(names):]: ax.axis("off")
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(Q2 / fn, dpi=100); plt.close(fig)
shapes([n for n in NAMES if n.startswith("sbnd")], "q2_fig2_shapes_sbnd.png"); shapes([n for n in ALL if "run2" in n], "q2_fig2_shapes_icarus_run2.png"); shapes([n for n in ALL if "run4" in n], "q2_fig2_shapes_icarus_run4.png")
# ---- Fig 3: signal vs noise
fig, axs = plt.subplots(2, 2, figsize=(13, 10)); ax = axs[0, 0]; x = np.arange(len(ALL))
ax.bar(x, [summ["per_variation"][n]["noise_trace_frac"] for n in ALL], color=[CC[summ["per_variation"][n]["cls"]] for n in ALL]); ax.axhline(1, color="k", lw=1, ls=":"); ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([short(n) for n in ALL], rotation=90, fontsize=8)
ax.set_ylabel("noise trace / raw trace"); mark_axis(ax, WM)
ax = axs[0, 1]; e = np.array([summ["per_variation"][n]["expected_mean_abs_if_pure_noise"] for n in ALL]) * 100; o = np.array([summ["per_variation"][n]["mean_abs_y"] for n in ALL]) * 100
ax.scatter(e, o, c=[CC[summ["per_variation"][n]["cls"]] for n in ALL], s=40); ax.plot([0.1, 8], [0.1, 8], "k-", lw=1); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"noise-only expectation $\langle|r-1|\rangle=0.80\,\langle\sigma_b\rangle$ [%]"); ax.set_ylabel(r"observed $\langle|r-1|\rangle$ [%]"); [ax.scatter([], [], c=CC[k], label=NAMEC) for k, NAMEC in (("i", "(i) smooth trend"), ("ii", "(ii) constant offset"), ("iii", "(iii) consistent with zero"))]; ax.legend(loc="lower right", fontsize=9)
ax = axs[1, 0]; pulls = np.concatenate([((DATA[n]["y"] - DATA[n]["f"]) / DATA[n]["s"])[:-1] for n in ALL if summ["per_variation"][n]["deg_bic"] >= -1]); h = np.linspace(-4, 4, 33); ax.hist(pulls, h, color="#0072B2", alpha=0.7, density=True, label=f"bulk-bin residual pulls (N={len(pulls)}, sd={pulls.std():.2f})"); xx = np.linspace(-4, 4, 200); ax.plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), "k-", label="N(0,1)"); ax.set_xlabel(r"$(r_b-1-\hat f_b)/\sigma_b$"); ax.legend(fontsize=8)
ax = axs[1, 1]; sc = json.load(open(Q2 / "q2_splitcv.json")) if (Q2 / "q2_splitcv.json").exists() else None
if sc:
    ms = ["raw", "poly-BIC, last raw", "whit-CV, last raw", "poly-BIC, all bins"]; xx = np.arange(len(ALL))
    for k, m in enumerate(ms): ax.plot(xx + 0.15 * (k - 1.5), [sc[n]["bulk"][m] / sc[n]["floor_bulk"] for n in ALL], "o", ms=4, label=m)
    ax.axhline(1, color="k", lw=1); ax.set_ylim(0.5, 2.4); ax.set_xticks(xx); ax.set_xticklabels([short(n) for n in ALL], rotation=90, fontsize=7); ax.set_ylabel("held-out half: $T$ / noise floor (bulk bins)"); ax.legend(fontsize=8, ncol=2, loc="upper right")
fig.tight_layout(); fig.savefig(Q2 / "q2_fig3_signal_noise.png", dpi=110); plt.close(fig)
# ---- Fig 4: covariance
fig, axs = plt.subplots(2, 3, figsize=(16, 9.5)); gl = [("sbnd_run1|xml-active", "SBND (xml-active set)"), ("icarus_run2|xml-active", "ICARUS Run2"), ("icarus_run4|xml-active", "ICARUS Run4")]
for ax, (g, t) in zip(axs[0], gl):
    d_ = cov["per_detector"][g]; k = np.arange(1, 8); er = np.array(d_["eig_raw"]["ev"][:7]); ef = np.array(d_["eig_fit"]["ev"][:7]); en = np.array(d_["eig_noise"]["ev"][:7])
    ax.bar(k - 0.27, er / sum(d_["eig_raw"]["ev"]), 0.27, label="raw $\\sum_v (r-1)(r-1)^T$", color="#D55E00"); ax.bar(k, ef / sum(d_["eig_fit"]["ev"]), 0.27, label="smooth fit", color="#0072B2"); ax.bar(k + 0.27, en / sum(d_["eig_noise"]["ev"]), 0.27, label="noise only ($\\sum_v\\Sigma_v$)", color="#999999")
    ax.set_xlabel("eigenvalue index"); ax.set_ylabel("fraction of trace"); ax.set_title(f"{t}: eff. rank raw {d_['eig_raw']['eff_rank']:.2f} / fit {d_['eig_fit']['eff_rank']:.2f}", fontsize=10)
    if g.startswith("sbnd"): ax.legend(fontsize=8)
    if g.startswith("sbnd"): mark_axis(ax, WM)
for ax, (g, t) in zip(axs[1], gl):
    ns = cov["per_detector"][g]["variations"]; Craw = sum(np.outer(DATA[n]["y"], DATA[n]["y"]) for n in ns); Cfit = sum(np.outer(DATA[n]["f"], DATA[n]["f"]) for n in ns)
    cr = lambda M: M / np.sqrt(np.outer(np.diag(M), np.diag(M))); A = cr(Craw); B = cr(Cfit); M = np.triu(A, 1) + np.tril(B, -1) + np.eye(NB) * 1.0
    im = ax.imshow(M, origin="lower", cmap="RdBu_r", norm=TwoSlopeNorm(0, -1, 1)); ax.set_title(f"{t}: upper = raw, lower = smooth fit\nadjacent corr. raw {cov['per_detector'][g]['adj_corr_raw']:.2f} / fit {cov['per_detector'][g]['adj_corr_fit']:.2f}", fontsize=10); ax.set_xlabel("bin"); ax.set_ylabel("bin")
fig.colorbar(im, ax=axs[1].tolist(), fraction=0.02); fig.savefig(Q2 / "q2_fig4_covariance.png", dpi=100); plt.close(fig)
# ---- Fig 5: consistency
T = {t["label"]: t for t in cons["tests"]}; fig, axs = plt.subplots(2, 5, figsize=(22, 8.5))
for ax, v in zip(axs[0], ICAR):
    t = T[f"ICARUS Run2 vs Run4: {v}"]; a, b = f"icarus_run2__{v}", f"icarus_run4__{v}"
    ax.errorbar(CTR - 0.004, 100 * DATA[a]["y"], yerr=100 * DATA[a]["s"], fmt="o", color="#D55E00", ms=3, label="Run2"); ax.errorbar(CTR + 0.004, 100 * DATA[b]["y"], yerr=100 * DATA[b]["s"], fmt="s", color="#009E73", ms=3, label="Run4"); ax.axhline(0, color="0.5", lw=0.8)
    ax.set_title(f"{v}: $\\chi^2$/ndf(Run2$-$Run4) = {t['chi2_equal']/NB:.2f}", fontsize=10); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.set_ylim(-4, 4)
    if v == "ChargeScale": ax.legend(fontsize=8); mark_axis(ax, WM)
ax = axs[1, 0]; sbn = lambda v: DATA[f"sbnd_run1__{v}"]
t = T["SBND ChargeScale vs ChargeScale15"]; ax.errorbar(CTR, 100 * sbn("ChargeScale15")["y"], yerr=100 * sbn("ChargeScale15")["s"], fmt="o", color="#0072B2", ms=3, label="ChargeScale15"); ax.errorbar(CTR, 100 * sbn("ChargeScale")["y"] / t["alpha_prop"], yerr=100 * sbn("ChargeScale")["s"] / t["alpha_prop"], fmt="s", color="#D55E00", ms=3, label=f"ChargeScale / {t['alpha_prop']:.2f}")
ax.set_title(f"proportional fit: $\\chi^2$/16 = {t['chi2_prop']/16:.2f} (p={t['p_prop']:.2f})", fontsize=10); ax.legend(fontsize=8); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]")
ax = axs[1, 1]; t = T["SBND WireModYZ vs MattModYZ"]; ax.errorbar(CTR, 100 * sbn("WireModYZ")["y"], yerr=100 * sbn("WireModYZ")["s"], fmt="o", color="#0072B2", ms=3, label="WireModYZ"); ax.errorbar(CTR, 100 * sbn("MattModYZ")["y"] * t["alpha_prop"], yerr=100 * sbn("MattModYZ")["s"] * t["alpha_prop"], fmt="s", color="#D55E00", ms=3, label=f"MattModYZ $\\times$ {t['alpha_prop']:.1f}")
ax.set_title(f"proportional fit: $\\chi^2$/16 = {t['chi2_prop']/16:.2f} (p={t['p_prop']:.3f})", fontsize=10); ax.legend(fontsize=8); ax.set_xlabel("Reco visible energy [GeV]")
for ax, det in zip(axs[1, 2:], ("sbnd_run1", "icarus_run2", "icarus_run4")[:3]):
    pass
for ax, det in zip(axs[1, 2:5], ("sbnd_run1", "icarus_run2", "icarus_run4")):
    t = T[f"{det} SCE0x + SCE2x (antisymmetry)"]; a = DATA[f"{det}__SCE0x"]; b = DATA[f"{det}__SCE2x"]
    ax.errorbar(CTR - 0.004, 100 * a["y"], yerr=100 * a["s"], fmt="o", color="#0072B2", ms=3, label="SCE0x"); ax.errorbar(CTR + 0.004, -100 * b["y"], yerr=100 * b["s"], fmt="s", color="#D55E00", ms=3, label="$-$SCE2x"); ax.axhline(0, color="0.5", lw=0.8)
    ax.set_title(f"{det}: $\\chi^2$/ndf(0x+2x) = {t['chi2_antisym']/NB:.2f}; quad: {100*t['rms_quad']:.2f}% vs $\\sigma$ {100*t['rms_sigma_quad']:.2f}%", fontsize=9); ax.set_xlabel("Reco visible energy [GeV]"); ax.legend(fontsize=8); ax.set_ylim(-5, 5)
fig.tight_layout(); fig.savefig(Q2 / "q2_fig5_consistency.png", dpi=100); plt.close(fig)
# ---- Fig 6: pool
fig, axs = plt.subplots(2, 2, figsize=(13, 10)); x = np.arange(len(NAMES))
ax = axs[0, 0]; ax.bar(x, [pool[n]["n_pool"] / 1e6 for n in NAMES], color=[DC[det_of(n)] for n in NAMES]); ax.set_xticks(x); ax.set_xticklabels([short(n) for n in NAMES], rotation=90, fontsize=8); ax.set_ylabel("common-event pool [M events]"); mark_axis(ax, WM)
ax = axs[0, 1]; nn = [n for n in NAMES if "common" in pool[n]]; ax.plot(range(len(nn)), [pool[n]["common"]["chi2"] / NB for n in nn], "o", color="#0072B2"); ax.axhline(1, color="k", lw=1); ax.set_xticks(range(len(nn))); ax.set_xticklabels([short(n) for n in nn], rotation=90, fontsize=8); ax.set_ylabel(r"$\chi^2$/ndf (full pool $-$ common pool), joint bootstrap")
ax = axs[1, 0]; ax.bar(x, [pool[n]["dup"]["rms_dr_over_sigma"] for n in NAMES], color=[DC[det_of(n)] for n in NAMES]); ax.set_xticks(x); ax.set_xticklabels([short(n) for n in NAMES], rotation=90, fontsize=8); ax.set_ylabel(r"rms$_b$ [(all rows $-$ unique keys) ratio]/$\sigma_b$"); ax.set_yscale("log")
ax = axs[1, 1]; cl = json.load(open(Q2 / "q2_run4_wiremodx_clean.json")); import uproot; prod = uproot.open(str(REPO / "detvar/detsys-v1.1.0v2.root"))["WireModxThetaXW_icarus_Run4"].values()
ax.plot(CTR, 100 * (prod - 1), "o-", color="#D55E00", label="production histogram (all rows, dup. events counted twice)"); ax.plot(CTR, 100 * DATA["icarus_run4__WireModxThetaXW"]["y"], "s-", color="#0072B2", label="pair table (duplicated keys dropped)"); ax.errorbar(CTR, 100 * DATA["icarus_run4__WireModxThetaXWclean"]["y"], yerr=100 * DATA["icarus_run4__WireModxThetaXWclean"]["s"], fmt="^-", color="#009E73", label="runs 12965-12970 removed")
ax.axhline(0, color="0.5"); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.set_title("ICARUS Run4 WireModxThetaXW", fontsize=10); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(Q2 / "q2_fig6_pool.png", dpi=110); plt.close(fig)
# ---- Fig 7: validation heat table
ms = [m for m in val["res"][NAMES[0]]["zero"] if "shrunk" not in m]; trs = list(val["res"][NAMES[0]]); M = np.array([[np.mean([val["res"][n][tn][m]["mse_ratio"] for n in NAMES]) for tn in trs] for m in ms])
fig, ax = plt.subplots(figsize=(14, 6.5)); im = ax.imshow(M, cmap="viridis_r", vmin=0, vmax=1.4, aspect="auto"); ax.set_xticks(range(len(trs))); ax.set_xticklabels(trs, rotation=40, ha="right", fontsize=8); ax.set_yticks(range(len(ms))); ax.set_yticklabels(ms, fontsize=9)
for i in range(len(ms)):
    for j in range(len(trs)): ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", color="w" if M[i, j] > 0.6 else "k", fontsize=8)
mark_axis(ax, WM); fig.tight_layout(); fig.colorbar(im, label="MSE vs truth / sum sigma_b^2 (raw ~ 1); mean over 20 variations"); fig.savefig(Q2 / "q2_fig7_validation.png", dpi=110); plt.close(fig)
print("figs done")
