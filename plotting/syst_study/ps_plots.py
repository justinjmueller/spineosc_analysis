#!/usr/bin/env python3
"""Figures for the per-systematic report -> docs/syst_study/persyst/*.png  (reads core/*.json, pairs_compare.json, trigger.json)"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ps_common import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
plt.rcParams.update({"axes.labelsize": 10, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 7.5, "axes.titlesize": 9})
CORE = PSFIG / "core"; C = {n: json.loads((CORE / f"{n}.json").read_text()) for n in PAIRS}
PC = json.loads((PSFIG / "pairs_compare.json").read_text()); TR = json.loads((PSFIG / "trigger.json").read_text())
A = lambda x: np.array([np.nan if v is None else v for v in x], float)
COL = {"SCE0x": "#0072B2", "SCE2x": "#D55E00", "WireModxThetaXW": "#009E73", "WireModYZ": "#CC79A7", "MattModYZ": "#E69F00", "DENT": "#56B4E9"}
def _tl(fig):
    try: _tl(fig)
    except RuntimeError: pass
def lab(n): return n.split("__")[1]
def tag(n): return n.split("__")[0].replace("_run", " Run").replace("sbnd", "SBND").replace("icarus", "ICARUS")
def mk(ax, txt=None, right=None):
    mark_axis(ax, r"$\bf{SBN}$ Analysis-in-Progress")
    if right: mark_axis(ax, right, alignment="right")
def hline(ax): ax.axhline(0, color="0.5", lw=0.8)

def fig_ratio_flows(pairs, fname, title):
    nr = len(pairs); fig, axs = plt.subplots(nr, 3, figsize=(15, 3.6 * nr), squeeze=False)
    for i, n in enumerate(pairs):
        R = C[n]; c = COL.get(lab(n), "k"); r = 100 * (A(R["ratio"]) - 1); e = 100 * A(R["err"]); F = {k: 100 * A(v) for k, v in R["flows"].items()}
        ax = axs[i, 0]; ax.errorbar(CTR, r, e, fmt="o", color="k", ms=4, label="all (applied)")
        ax.errorbar(CTR - 0.004, 100 * (A(R["r1p"]) - 1), 100 * A(R["e1p"]), fmt="s", color="#0072B2", ms=3, alpha=0.8, label="1p")
        ax.errorbar(CTR + 0.004, 100 * (A(R["rMp"]) - 1), 100 * A(R["eMp"]), fmt="^", color="#D55E00", ms=3, alpha=0.8, label="Mp")
        hline(ax); ax.set_ylabel("Ratio $-$ 1 [%]"); lo, hi = np.nanmin(r) - 1.5, np.nanmax(r) + 1.5; ax.set_ylim(max(lo, -15), min(hi, 15)); ax.legend(fontsize=8, ncol=3, loc="upper left")
        ax.set_title(f"{tag(n)} {lab(n)}", fontsize=10)
        ax = axs[i, 1]; ax.plot(CTR, F["lost"], "-o", ms=3, color="#D55E00", label="lost (nom only)"); ax.plot(CTR, F["gained"], "-o", ms=3, color="#0072B2", label="gained (var only)"); ax.plot(CTR, F["migration"], "-o", ms=3, color="#009E73", label="migration")
        ax.plot(CTR, r, "k-", lw=2, label="sum = ratio-1"); hline(ax); ax.set_ylabel("fraction of nominal [%]"); ax.legend(fontsize=7, ncol=2, loc="lower center" if True else None)
        ax = axs[i, 2]; ch = R["channel"]; d1 = 100 * A(ch["d1p"]); dM = 100 * A(ch["dMp"]); ax.errorbar(CTR - 0.004, d1, 100 * A(ch["e_d1p"]), fmt="s", color="#0072B2", ms=3, label=r"$r_{1p}-r_{all}$")
        ax.errorbar(CTR + 0.004, dM, 100 * A(ch["e_dMp"]), fmt="^", color="#D55E00", ms=3, label=r"$r_{Mp}-r_{all}$"); ax.fill_between(CTR, -np.abs(r), np.abs(r), color="0.85", label=r"$\pm|r_{all}-1|$"); hline(ax)
        ax.set_ylabel("channel $-$ combined [%]"); ax.set_ylim(-6, 6); ax.legend(fontsize=8, ncol=3, loc="upper left")
        for a in axs[i]: a.set_xlabel("Reco visible energy [GeV]")
    mk(axs[0, 0]); fig.suptitle(title, y=1.0, fontsize=11); _tl(fig); savefig(fig, fname)

def cond_panel(ax, pairs, key, xlab, logx=False, extra="cond", xmid=None, ylim=None):
    for n in pairs:
        R = C[n]; d = R[extra][key]; ed = A(d["edges"]); x = xmid(ed) if xmid else 0.5 * (ed[1:] + ed[:-1]); c = COL.get(lab(n), "k")
        ax.errorbar(x, 100 * A(d["lost"]), 100 * A(d["lost_err"]), fmt="-o", ms=3, color=c, label=f"{lab(n)} lost | nominal bin")
        ax.errorbar(x, 100 * A(d["gained"]), 100 * A(d["gained_err"]), fmt="--s", ms=3, color=c, alpha=0.7, label=f"{lab(n)} new | varied bin")
    ax.set_xlabel(xlab); ax.set_ylabel("fraction [%]")
    if logx: ax.set_xscale("log")
    if ylim: ax.set_ylim(*ylim)

def fig_mech(pairs, fname, title, theta=True):
    fig, axs = plt.subplots(2, 4, figsize=(17, 7.4))
    cond_panel(axs[0, 0], pairs, "p_margin", "proton margin $s_p-\\max(2s_\\mu,s_\\pi)$"); cond_panel(axs[0, 1], pairs, "mu_margin", "muon margin $2s_\\mu-\\max(s_\\pi,s_p)$")
    cond_panel(axs[0, 2], pairs, "p_ke", "leading-proton KE [MeV]", logx=True); cond_panel(axs[0, 3], pairs, "n_p", "proton multiplicity (nominal)")
    axs[0, 0].legend(fontsize=7)
    for n in pairs:
        R = C[n]; c = COL.get(lab(n), "k"); P = R["pos"]["drift"]; ed = A(P["edges"]); x = 0.5 * (ed[1:] + ed[:-1]); rows = P["rows"]
        axs[1, 0].plot(x, 100 * A([r["dL"] for r in rows]), "-o", ms=3, color=c, label=lab(n)); axs[1, 1].plot(x, 100 * A(P["r_pos"]), "-o", ms=3, color=c, label=lab(n))
        d = R["cond"]["drift"]; e2 = A(d["edges"]); x2 = 0.5 * (e2[1:] + e2[:-1]); axs[1, 2].errorbar(x2, 100 * A(d["lost"]), 100 * A(d["lost_err"]), fmt="-o", ms=3, color=c); axs[1, 2].errorbar(x2, 100 * A(d["gained"]), 100 * A(d["gained_err"]), fmt="--s", ms=3, color=c, alpha=0.7)
        if theta:
            T = R["pos"]["theta_xz"]; et = A(T["edges"]); xt = 0.5 * (et[1:] + et[:-1]); axs[1, 3].plot(xt, 100 * A([r["dcalo"] for r in T["rows"]]), "-o", ms=3, color=c, label=lab(n) + " $\\Delta$calo(mu) [%-pt]")
            axs[1, 3].plot(xt, 100 * A([r["dscore"] for r in T["rows"]]), "--s", ms=3, color=c, alpha=0.7, label=lab(n) + " $\\Delta s_\\mu$ [%-pt]")
    axs[1, 0].set_ylabel("mean $\\Delta L_\\mu/L_\\mu$ (kept) [%]"); axs[1, 1].set_ylabel("net: varied/nominal $-$ 1 [%]"); axs[1, 2].set_ylabel("lost (solid) / new (dashed) [%]"); axs[1, 3].set_ylabel("shift of kept events"); axs[1, 3].set_xlabel("muon $\\theta_{xz}$ [deg] (90 = along drift)")
    for a in axs[1, :3]: a.set_xlabel("distance to cathode [cm]")
    for a in (axs[1, 0], axs[1, 1], axs[1, 3]): hline(a)
    axs[1, 0].legend(fontsize=8); axs[1, 3].legend(fontsize=7); mk(axs[0, 0]); fig.suptitle(title, y=1.0, fontsize=11); _tl(fig); savefig(fig, fname)

def fig_sce_linear(fname):
    fig, axs = plt.subplots(2, 3, figsize=(16, 7.4)); J = PC["joint_boot"]
    for j, det in enumerate(("sbnd_run1", "icarus_run2", "icarus_run4")):
        ax = axs[0, j]; r0 = 100 * (A(C[f"{det}__SCE0x"]["ratio"]) - 1); r2 = 100 * (A(C[f"{det}__SCE2x"]["ratio"]) - 1); e0 = 100 * A(C[f"{det}__SCE0x"]["err"]); e2 = 100 * A(C[f"{det}__SCE2x"]["err"])
        ax.errorbar(CTR, r0, e0, fmt="o", color=COL["SCE0x"], ms=4, label="SCE0x"); ax.errorbar(CTR, r2, e2, fmt="s", color=COL["SCE2x"], ms=4, label="SCE2x"); ax.plot(CTR, -r0, "--", color=COL["SCE0x"], alpha=0.7, label="$-($SCE0x$-1)$")
        hline(ax); ax.set_title(tag(det + "__x"), fontsize=10); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.set_xlabel("Reco visible energy [GeV]"); ax.legend(fontsize=8, ncol=3, loc="upper center")
        ax = axs[1, j]; ax.errorbar(CTR - 0.004, 100 * A(J[det]["odd"]), 100 * A(J[det]["e_odd"]), fmt="o", color="k", ms=4, label=f"odd $=(r_{{2x}}-r_{{0x}})/2$  [$\\chi^2$={J[det]['chi2_odd']:.0f}/{NBIN}]")
        ax.errorbar(CTR + 0.004, 100 * A(J[det]["even"]), 100 * A(J[det]["e_even"]), fmt="s", color="#CC79A7", ms=4, label=f"even $=(r_{{2x}}+r_{{0x}})/2$  [$\\chi^2$={J[det]['chi2_even']:.0f}/{NBIN}]")
        hline(ax); ax.set_ylabel("odd / even part of ratio $-$ 1 [%]"); ax.set_xlabel("Reco visible energy [GeV]"); ax.legend(fontsize=8, loc="upper left")
    mk(axs[0, 0]); _tl(fig); savefig(fig, fname)

def fig_sce_event(fname):
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6)); H = PC["hamming"]; dets = ["sbnd_run1", "icarus_run2", "icarus_run4"]; x = np.arange(3); w = 0.25
    for k, (key, nm, c) in enumerate((("d_n0_over_N", "nominal vs 0x", COL["SCE0x"]), ("d_n2_over_N", "nominal vs 2x", COL["SCE2x"]), ("d_02_over_N", "0x vs 2x", "k"))): axs[0].bar(x + (k - 1) * w, [100 * H[d][key] for d in dets], w, color=c, label=nm)
    axs[0].set_xticks(x); axs[0].set_xticklabels(["SBND", "ICARUS R2", "ICARUS R4"]); axs[0].set_ylabel("|selection set difference| / N$_{nom}$ [%]"); axs[0].legend(fontsize=8); axs[0].set_title("Hamming distance between selected sets", fontsize=10)
    S = PC; names = ["dE", "dmu_len", "dp_len", "dmu_score", "dp_score", "dmu_calo", "dvx", "dvz", "dmu_endx"]; ax = axs[1]
    for k, det in enumerate(dets):
        ec = PC[f"SCE_{det}"]["event_corr"]; ax.bar(np.arange(len(names)) + (k - 1) * w, [ec[n]["pearson"] for n in names], w, label=f"{tag(det + '__x')} Pearson", color=["#0072B2", "#009E73", "#E69F00"][k])
        ax.plot(np.arange(len(names)) + (k - 1) * w, [ec[n]["spearman"] for n in names], "o", color="k", mfc="w", ms=5)
    ax.axhline(0, color="0.5", lw=0.8); ax.set_xticks(np.arange(len(names))); ax.set_xticklabels(["$\\Delta E$", "$\\Delta L_\\mu$", "$\\Delta L_p$", "$\\Delta s_\\mu$", "$\\Delta s_p$", "$\\Delta$calo$_\\mu$", "$\\Delta v_x$", "$\\Delta v_z$", "$\\Delta x^{end}_\\mu$"], rotation=30)
    ax.set_ylabel("corr(shift 0x, shift 2x)"); ax.legend(fontsize=7); ax.set_title("bars Pearson, ticks Spearman (odd response: negative)", fontsize=10)
    ax = axs[2]
    for k, det in enumerate(dets):
        st = PC[f"SCE_{det}"]["states"]; ax.bar(k - 0.2, 100 * st["P_lostB_given_lostA"], 0.2, color="#D55E00"); ax.bar(k, 100 * st["P_lostB_base"], 0.2, color="0.6"); ax.bar(k + 0.2, 100 * st["P_gainedB_given_gainedA_lowerbound"], 0.2, color="#0072B2")
    from matplotlib.patches import Patch; ax.legend(handles=[Patch(color="#D55E00", label="P(lost in 2x | lost in 0x)"), Patch(color="0.6", label="P(lost in 2x) overall"), Patch(color="#0072B2", label="P(new in 2x | new in 0x), lower bound")], fontsize=7, loc="upper right")
    ax.set_xticks(range(3)); ax.set_xticklabels(["SBND", "ICARUS R2", "ICARUS R4"]); ax.set_ylabel("[%]"); ax.set_title("churn overlap between 0x and 2x", fontsize=10)
    mk(axs[0]); _tl(fig); savefig(fig, fname)

def fig_yz(fname):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8)); g = PC["YZ"]["grid"]; ey = np.array(g["edges_y"]); ez = np.array(g["edges_z"])
    for j, (k, nm) in enumerate((("WireModYZ", "WireModYZ"), ("MattModYZ", "MattModYZ"))):
        ax = axs[0, j]; im = ax.pcolormesh(ez, ey, 100 * A(np.array(g[k]).ravel()).reshape(len(ey) - 1, len(ez) - 1), cmap="RdBu_r", vmin=-8, vmax=8); ax.set_xlabel("vertex z [cm]"); ax.set_ylabel("vertex y [cm]"); ax.set_title(f"{nm}: varied/nominal $-$ 1 [%] per cell (>300 evts)", fontsize=10); fig.colorbar(im, ax=ax)
    ax = axs[0, 2]; sm = PC["YZ"]["shift_summary"]; names = ["dmu_calo", "dp_calo", "dmu_score", "dp_score", "dE"]
    for k, (v, c) in enumerate((("WireModYZ", COL["WireModYZ"]), ("MattModYZ", COL["MattModYZ"]), ("ChargeScale", "0.4"), ("GainHigh", "0.7"))):
        ax.bar(np.arange(len(names)) + (k - 1.5) * 0.2, [100 * sm[v][n]["mean"] for n in names], 0.2, color=c, label=v)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(["$\\Delta$calo$_\\mu$", "$\\Delta$calo$_p$", "$\\Delta s_\\mu$", "$\\Delta s_p$", "$\\Delta E/E$"]); ax.axhline(0, color="0.5", lw=0.8); ax.set_ylabel("mean shift [%-pt or %]"); ax.legend(fontsize=8)
    # charge-scale proportionality of the energy slope (SBND)
    b = json.loads((REPO / "docs" / "syst_study" / "baseline.json").read_text()); pts = {}
    def slope(r, e):
        ok = np.isfinite(r) & (e > 0); w = 1 / e[ok] ** 2; X = np.vstack([np.ones(ok.sum()), CTR[ok]]).T; cov = np.linalg.inv(X.T @ (X * w[:, None])); be = cov @ (X.T @ (w * r[ok])); return 100 * be[1], 100 * np.sqrt(cov[1, 1])
    for v in ("GainHigh", "GainLow", "ChargeScale", "ChargeScale15", "MattModYZ", "WireModYZ", "WireModxThetaXW"):
        sh = sm[v]["dmu_calo"]["mean"] * 100 if v in sm else None; bb = b[f"sbnd_run1__{v}"]; s, es = slope(np.array(bb["ratio"], float) - 1, np.array(bb["err"], float)); pts[v] = (sh, s, es)
    ax = axs[1, 0]
    for v, (sh, s, es) in pts.items():
        if sh is None: continue
        ax.errorbar(sh, s, es, fmt="o", ms=6, color=COL.get(v, "k")); ax.annotate(v, (sh, s), fontsize=7, xytext=(4, 4), textcoords="offset points")
    xs = np.linspace(-3, 5, 10); ax.plot(xs, 3.5 * xs, "k:", lw=1); ax.axhline(0, color="0.5", lw=0.8); ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("mean $\\Delta$calo$_\\mu$ of kept events [%-pt]"); ax.set_ylabel("slope of ratio vs $E_{vis}$ [%/GeV]"); ax.set_title("energy trend vs charge shift (SBND); dotted: 3.5%/GeV per %-pt", fontsize=10)
    ax = axs[1, 1]
    for n in ("sbnd_run1__WireModYZ", "sbnd_run1__MattModYZ"):
        R = C[n]; ax.plot(CTR, 100 * A(R["flows"]["lost"]), "-", color=COL[lab(n)], label=lab(n) + " lost"); ax.plot(CTR, 100 * A(R["flows"]["gained"]), "--", color=COL[lab(n)], label=lab(n) + " gained")
    hline(ax); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("fraction of nominal [%]"); ax.legend(fontsize=8)
    ax = axs[1, 2]; st = PC["YZ"]; names = ["dE", "dmu_len", "dmu_score", "dp_score", "dmu_calo", "dp_calo", "dvx"]; co = st["WireModYZ_vs_MattModYZ"]["corr"]
    ax.bar(np.arange(len(names)), [co[n]["pearson"] for n in names], 0.4, color="0.4", label="WireModYZ vs MattModYZ"); co2 = st["MattModYZ_vs_ChargeScale"]["corr"]; ax.bar(np.arange(len(names)) + 0.4, [co2[n]["pearson"] for n in names], 0.4, color="0.75", label="MattModYZ vs ChargeScale")
    ax.set_xticks(np.arange(len(names)) + 0.2); ax.set_xticklabels(["$\\Delta E$", "$\\Delta L_\\mu$", "$\\Delta s_\\mu$", "$\\Delta s_p$", "$\\Delta$calo$_\\mu$", "$\\Delta$calo$_p$", "$\\Delta v_x$"], rotation=30); ax.axhline(0, color="0.5", lw=0.8); ax.set_ylabel("Pearson corr. of event shifts"); ax.legend(fontsize=8)
    mk(axs[0, 0]); _tl(fig); savefig(fig, fname)

def fig_trigger(fname):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8)); dets = ["sbnd_run1", "icarus_run2", "icarus_run4"]
    for j, det in enumerate(dets):
        T = TR[det]; ax = axs[0, j]; H = np.array(T["pe_hist2d"]); ye = np.linspace(2.5, 5.5, 61); im = ax.pcolormesh(BINS, ye, np.log10(H.T + 1), cmap="viridis"); ax.axhline(np.log10(T["pe_threshold"]), color="w", lw=1.5, ls="--")
        ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("log$_{10}$ largest-flash PE"); ax.set_title(f"{tag(det + '__x')}: nominal threshold {T['pe_threshold']:.0f} PE ({T['cfg']['scale']}$\\times$PE>{T['cfg']['thr']})", fontsize=9)
        ax = axs[1, j]; ax.errorbar(CTR, 100 * (A(T["M1"]["ratio"]) - 1), fmt="o-", ms=3, color="#0072B2", label=f"M1 (scale $-$ {T['cfg']['sigma']})"); ax.errorbar(CTR, 100 * (A(T["P1"]["ratio"]) - 1), fmt="s-", ms=3, color="#D55E00", label=f"P1 (scale $+$ {T['cfg']['sigma']})")
        ax2 = ax.twinx(); ax2.plot(CTR, 100 * np.array([b["fail"] for b in T["bins"]]), "k:", label="selected evts failing nominal trigger [%]"); ax2.set_ylabel("failing nominal trigger [%]"); ax2.set_ylim(0, max(25, 100 * T["bins"][0]["fail"] * 1.1)); ax2.legend(fontsize=7, loc="lower right")
        hline(ax); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.legend(fontsize=8, loc="upper right" if j == 0 else "center right")
    mk(axs[0, 0]); _tl(fig); savefig(fig, fname)

def fig_dent(fname):
    n = "sbnd_run1__DENT"; R = C[n]; fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axs[0]; sd = R["side"]
    for k, (s, c) in enumerate((("x<0", "#0072B2"), ("x>0", "#D55E00"))): ax.bar(np.arange(3) + (k - 0.5) * 0.35, [100 * sd[s]["lost"], 100 * sd[s]["gained"], 100 * sd[s]["net"]], 0.35, color=c, label=f"vertex {s}")
    ax.set_xticks(range(3)); ax.set_xticklabels(["lost", "new (/nominal)", "net"]); ax.set_ylabel("[% of nominal events in the half]"); ax.axhline(0, color="0.5", lw=0.8); ax.legend(fontsize=8); ax.set_title("DENT churn is the same in both TPC halves", fontsize=10)
    ax = axs[1]; P = R["pos"]["drift"]; ed = A(P["edges"]); x = 0.5 * (ed[1:] + ed[:-1]); d = R["cond"]["drift"]; e2 = A(d["edges"]); x2 = 0.5 * (e2[1:] + e2[:-1])
    ax.errorbar(x2, 100 * A(d["lost"]), 100 * A(d["lost_err"]), fmt="-o", ms=3, color=COL["DENT"], label="DENT lost"); ax.errorbar(x2, 100 * A(d["gained"]), 100 * A(d["gained_err"]), fmt="--s", ms=3, color=COL["DENT"], alpha=0.7, label="DENT new")
    for v in ("SCE0x", "MattModYZ"):
        dd = C[f"sbnd_run1__{v}"]["cond"]["drift"]; ee = A(dd["edges"]); ax.errorbar(0.5 * (ee[1:] + ee[:-1]), 100 * A(dd["lost"]), 100 * A(dd["lost_err"]), fmt="-o", ms=3, color=COL[v], alpha=0.6, label=f"{v} lost")
    ax.set_xlabel("distance to cathode [cm]"); ax.set_ylabel("fraction [%]"); ax.legend(fontsize=7)
    ax = axs[2]; pn = ["mu_margin", "p_margin"]
    cond_panel(ax, [n], "p_margin", "leading-proton margin"); ax.set_title("DENT: species-boundary dependence", fontsize=10); ax.legend(fontsize=7)
    mk(axs[0]); _tl(fig); savefig(fig, fname)

def fig_overview(fname):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.4)); order = list(PAIRS)
    ax = axs[0, 0]; S = json.loads(PS_JSON.read_text())["pairs"]
    for n in order:
        I = S[n]["integrated"]; ax.scatter(100 * (I["lost_frac"] + I["gained_frac"]), 100 * I["net"], color=COL.get(lab(n), "k"), marker={"sbnd_run1": "o", "icarus_run2": "s", "icarus_run4": "^"}[n.split("__")[0]], s=50)
    ax.axhline(0, color="0.5", lw=0.8); ax.set_xlabel("churn: lost + new [% of nominal]"); ax.set_ylabel("net change [%]"); ax.set_title("net is a few % of the churn (marker: o SBND, s ICARUS R2, ^ R4)", fontsize=9)
    for v, c in COL.items(): ax.scatter([], [], color=c, label=v)
    ax.legend(fontsize=7, ncol=2)
    ax = axs[0, 1]; keys = ["p_margin", "mu_margin", "p_calo_csda", "mu_primary", "mu_len", "p_ke", "n_p"]
    M = np.array([[S[n]["smd"][k]["lost_vs_kept"] for k in keys] for n in order]); im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto"); ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys, rotation=40, ha="right", fontsize=8); ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{tag(n)} {lab(n)}" for n in order], fontsize=7); fig.colorbar(im, ax=ax); ax.set_title("standardized mean difference, lost vs kept (nominal features)", fontsize=9)
    ax = axs[0, 2]
    for n in order:
        b = S[n]["boundary"]; ax.scatter(100 * b["p_margin"]["floor_lost_rate"], 100 * b["p_margin"]["excess_share_of_lost"], color=COL.get(lab(n), "k"), marker={"sbnd_run1": "o", "icarus_run2": "s", "icarus_run4": "^"}[n.split("__")[0]], s=50)
    ax.set_xlabel("lost rate of confident events (floor) [%]"); ax.set_ylabel("share of lost above floor, by leading-proton margin [%]"); ax.set_title("how much churn sits at the species boundary", fontsize=9)
    ax = axs[1, 0]
    for n in order:
        c = S[n]["channel"]; ax.scatter(100 * c["rms_sys"], 100 * max(c["rms_d1p"], c["rms_dMp"]), color=COL.get(lab(n), "k"), marker={"sbnd_run1": "o", "icarus_run2": "s", "icarus_run4": "^"}[n.split("__")[0]], s=50)
    xs = np.linspace(0, 3.5, 3); ax.plot(xs, xs, "k:"); ax.set_xlabel("rms of applied ratio $-$ 1 [%]"); ax.set_ylabel("larger of rms(1p$-$all), rms(Mp$-$all) [%]"); ax.set_title("channel mismatch vs systematic size (above dotted: mismatch > effect)", fontsize=9)
    ax = axs[1, 1]
    for n in order:
        if n.startswith("sbnd"):
            R = C[n]; ax.plot(CTR, 100 * (A(R["ratio"]) - 1), "-o", ms=3, color=COL[lab(n)], label=lab(n))
    hline(ax); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.set_title("SBND", fontsize=10); ax.legend(fontsize=7, ncol=2)
    ax = axs[1, 2]
    for n in order:
        if n.startswith("icarus_run2") or n.startswith("icarus_run4"):
            R = C[n]; ax.plot(CTR, 100 * (A(R["ratio"]) - 1), "-" if "run2" in n else "--", marker="o", ms=2, color=COL[lab(n)], alpha=0.8, label=("R2 " if "run2" in n else "R4 ") + lab(n))
    hline(ax); ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel("Ratio $-$ 1 [%]"); ax.set_title("ICARUS", fontsize=10); ax.legend(fontsize=6, ncol=2)
    mk(axs[0, 0]); _tl(fig); savefig(fig, fname)

if __name__ == "__main__":
    fig_ratio_flows(["sbnd_run1__SCE0x", "sbnd_run1__SCE2x"], "sbnd_sce_ratio.png", "SBND SCE 0x / 2x")
    fig_mech(["sbnd_run1__SCE0x", "sbnd_run1__SCE2x"], "sbnd_sce_mech.png", "SBND SCE: where events flip and what moves")
    fig_sce_linear("sce_linear.png"); fig_sce_event("sce_event_level.png")
    fig_ratio_flows(["sbnd_run1__WireModxThetaXW"], "sbnd_wmxt_ratio.png", "SBND WireModxThetaXW"); fig_mech(["sbnd_run1__WireModxThetaXW"], "sbnd_wmxt_mech.png", "SBND WireModxThetaXW")
    fig_ratio_flows(["sbnd_run1__WireModYZ", "sbnd_run1__MattModYZ"], "sbnd_yz_ratio.png", "SBND WireModYZ (computed) vs MattModYZ (used)"); fig_yz("sbnd_yz_compare.png")
    fig_mech(["sbnd_run1__WireModYZ", "sbnd_run1__MattModYZ"], "sbnd_yz_mech.png", "SBND WireModYZ vs MattModYZ")
    fig_trigger("trigger.png"); fig_ratio_flows(["sbnd_run1__DENT"], "sbnd_dent_ratio.png", "SBND DENT"); fig_dent("sbnd_dent_mech.png")
    fig_ratio_flows(["icarus_run2__SCE0x", "icarus_run2__SCE2x", "icarus_run4__SCE0x", "icarus_run4__SCE2x"], "icarus_sce_ratio.png", "ICARUS SCE 0x / 2x (Run 2, Run 4)")
    fig_mech(["icarus_run2__SCE0x", "icarus_run2__SCE2x"], "icarus_sce_mech.png", "ICARUS Run 2 SCE: where events flip")
    fig_ratio_flows(["icarus_run2__WireModxThetaXW", "icarus_run4__WireModxThetaXW", "icarus_run2__WireModYZ", "icarus_run4__WireModYZ"], "icarus_wm_ratio.png", "ICARUS WireMod x-theta_xw and YZ")
    fig_mech(["icarus_run2__WireModxThetaXW", "icarus_run2__WireModYZ"], "icarus_wm_mech.png", "ICARUS Run 2 WireMod")
    fig_overview("overview.png"); print("done")
