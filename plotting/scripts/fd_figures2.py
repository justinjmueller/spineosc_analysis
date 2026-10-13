#!/usr/bin/env python3
"""Summary figures for the 2D fake-data campaign, built from the existing fits (no new PROfit runs):

  stoplight.pdf         datasets x analyses, coloured by Delta chi2 (green < 1, amber 1-4, red > 4), pull > 1 sigma flagged
  chi2_decomp.pdf/.csv  null-fit chi2 split into a prior part (sum of squared pulls) and the rest (data + MC statistics)
  absorbed_fraction.*   1 - (data part of the null chi2) / (stat-only chi2 of the injection), per dataset
  spectra/<stem>.pdf    injected / best-fit / residual maps (1p: 2D RecoE x pn; Np: 1D) per detector, null fit

Datasets are the data-XML ("ds") route: the injected dial is never a nuisance parameter.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_figures2.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).parent))
import fd2d as f2

OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
TAB = OUT / "tables"
(OUT / "spectra").mkdir(exist_ok=True, parents=True)
CFGS = ["b", "1p", "Np", "m"]
DET = {"SBND": "SBND", "ICARUSRun2": "ICARUS Run 2", "ICARUSRun4": "ICARUS Run 4"}
PANEL_STEMS = ["cpi_subthreshold_minus100", "cpi_subthreshold_plus100", "mec_bdt", "jaesung_lowq2_pi_enhancement_prefsi",
               "jaesung_lowq2_pi_enhancement_postfsi", "pi_fsi_ha2025", "ub_cc1p0pi_dpt_dat_midW", "ub_cc1p0pi_dpt_dat", "t2k_nc1pi_p_costh",
               "minerva_3dqelike_bnb", "qe_zexp_mva_to_lqcd"]


def all_stems():
    out = []
    for m in sorted((f2.RUN / "done").glob("ds_b_*")):
        s = m.name[len("ds_b_"):]
        if s != "ctrl" and not s.endswith(("_ratenorm",)) and "_amp" not in s and f2.available("ds", "b", "fix", s):
            out.append(s)
    return out


def global_files(cfg, stem, variant="fix"):
    tag = f2._paths("ds", cfg, variant, stem)[0]
    return f2.RUN / f"{tag}_PROglobal.root"


def counts(h):
    v = h.values()
    ex = h.axis(0).edges()
    if v.ndim == 2:
        ey = h.axis(1).edges()
        return v * np.outer(np.diff(ex), np.diff(ey)), ex, ey
    return v * np.diff(ex), ex, None


def spectra(cfg, stem, variant="fix"):
    """{det: {'1p': (data, cv, bf, ex, ey), 'Np': (...)}} in counts, or None."""
    p = global_files(cfg, stem, variant)
    if not p.exists():
        return None
    out = {}
    with uproot.open(p) as f:
        for det in DET:
            o = {}
            for ch, (nd, nc, nb) in (("1p", ("data2d", "cv2d", "bestfit2d")), ("Np", ("data", "cv1d", "bestfit"))):
                name = "ccqelike" if ch == "1p" else "nproton"
                try:
                    d, ex, ey = counts(f[f"nu_{det}_{name}_{nd};1"])
                    c, _, _ = counts(f[f"nu_{det}_{name}_{nc};1"])
                    b, _, _ = counts(f[f"nu_{det}_{name}_{nb};1"])
                except KeyError:
                    continue
                o[ch] = (d, c, b, ex, ey)
            out[det] = o
    return out


def stat_chi2(sp):
    tot = 0.0
    for det in sp.values():
        for d, c, b, ex, ey in det.values():
            tot += float(np.sum((d - c) ** 2 / np.maximum(d, 1e-9)))
    return tot


def decomposition():
    rows = []
    for s in all_stems():
        for c in CFGS:
            if not f2.available("ds", c, "fix", s):
                continue
            n = f2.read_fit("ds", c, "fix", s)
            pull2 = sum(n["values"][k] ** 2 for k in f2.nuisance_names(n))
            sp = spectra(c, s)
            S = stat_chi2(sp) if sp else np.nan
            rows.append({"dataset": s, "config": c, "chi2 null": n["chi2"], "sum pull^2": pull2, "chi2 data+MCstat": n["chi2"] - pull2,
                         "stat-only chi2 of injection": S, "absorbed fraction": 1 - (n["chi2"] - pull2) / S if S > 0 else np.nan})
    return pd.DataFrame(rows)


def stoplight(stems):
    cmap = ListedColormap(["#8fd19e", "#ffd966", "#e06666"])
    norm = BoundaryNorm([0, 1, 4, 1e9], cmap.N)
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(stems) + 1.6))
    for i, s in enumerate(stems):
        for j, c in enumerate(CFGS):
            p = f2.pair("ds", c, s)
            if not p:
                continue
            d = p[0]["chi2"] - p[1]["chi2"]
            nm = f2.nuisance_names(p[0])
            mp = max(abs(p[0]["values"][k]) for k in nm)
            ax.add_patch(plt.Rectangle((j, i), 1, 1, color=cmap(norm(d)), ec="w"))
            ax.text(j + 0.5, i + 0.5, f"{p[0]['chi2']:.1f} | {d:.2f}" + (" *" if mp > 1 else ""), ha="center", va="center", fontsize=7)
    ax.set_xlim(0, len(CFGS)); ax.set_ylim(len(stems), 0)
    ax.set_xticks(np.arange(len(CFGS)) + 0.5); ax.set_xticklabels([f2.CONFIGS[c] for c in CFGS], fontsize=8)
    ax.xaxis.tick_top(); ax.set_yticks(np.arange(len(stems)) + 0.5); ax.set_yticklabels(stems, fontsize=7)
    ax.set_title(r"null-fit $\chi^2$ | $\Delta\chi^2$ (green < 1, amber 1-4, red > 4); * = a pull above 1 $\sigma$", fontsize=9, pad=28)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    fig.savefig(OUT / "stoplight.pdf", bbox_inches="tight"); plt.close(fig)


def heat(ax, v, ex, ey, vmax, cmap, title):
    im = ax.pcolormesh(np.arange(len(ex)), np.arange(len(ey)), v.T, cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(0, len(ex), 3)); ax.set_xticklabels([f"{x:.2f}" for x in ex[::3]], fontsize=6)
    ax.set_yticks(np.arange(len(ey))); ax.set_yticklabels([f"{y:.2f}" for y in ey], fontsize=6)
    ax.set_title(title, fontsize=7)
    return im


def spectra_page(pdf, cfg, stem):
    sp = spectra(cfg, stem)
    if not sp:
        return
    fit = f2.read_fit("ds", cfg, "fix", stem)
    fig = plt.figure(figsize=(11, 10.5))
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 1.15], hspace=0.7, wspace=0.35)
    axs_ratio, axs_res, im1, im3 = [], [], None, None
    vmax = min(0.6, max(0.1, round(max((float(np.nanmax(np.abs(v["1p"][0] / v["1p"][1] - 1))) for v in sp.values() if "1p" in v), default=0.1), 1)))
    for r, det in enumerate(("SBND", "ICARUSRun4", "ICARUSRun2")):
        if "1p" not in sp.get(det, {}):
            continue
        d, c, b, ex, ey = sp[det]["1p"]
        res = (d - b) / np.sqrt(np.maximum(b, 1e-9))
        a0, a1, a2 = (fig.add_subplot(gs[r, k]) for k in range(3))
        im1 = heat(a0, d / c - 1, ex, ey, vmax, "RdBu_r", f"{DET[det]} 1p: injected / CV - 1")
        heat(a1, b / c - 1, ex, ey, vmax, "RdBu_r", "best fit / CV - 1")
        im3 = heat(a2, res, ex, ey, 3, "PuOr", r"residual $(D-BF)/\sqrt{BF}$")
        axs_ratio += [a0, a1]; axs_res.append(a2)
        a0.set_ylabel(r"$p_n$ [GeV/c]", fontsize=6)
        if r == 2:
            for a in (a0, a1, a2):
                a.set_xlabel("RecoE [GeV]", fontsize=6)
    if im1 is not None:
        fig.colorbar(im1, ax=axs_ratio, fraction=0.025, pad=0.02, shrink=0.8, label="fractional change")
        fig.colorbar(im3, ax=axs_res, fraction=0.05, pad=0.02, shrink=0.8, label="sigma")
    for k, det in enumerate(("SBND", "ICARUSRun4")):
        if "Np" not in sp.get(det, {}):
            continue
        d, c, b, ex, ey = sp[det]["Np"]
        ax = fig.add_subplot(gs[3, 0:2] if k == 0 else gs[3, 2])
        x = 0.5 * (ex[1:] + ex[:-1])
        ax.step(x, d / c, where="mid", color="k", label="injected / CV"); ax.step(x, b / c, where="mid", color="C3", label="best fit / CV")
        ax.axhline(1, color="0.6", lw=0.6); ax.set_title(f"{DET[det]} Np (1D RecoE)", fontsize=7); ax.tick_params(labelsize=6)
        if k == 0:
            ax.legend(fontsize=6); ax.set_xlabel("RecoE [GeV]", fontsize=6)
    p = fit["chi2"]
    fig.suptitle(f"{f2.CONFIGS[cfg]} | {stem} | null fit (oscillation fixed): chi2 = {p:.2f}", fontsize=9)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def proj(cfg, stem, variant):
    """{(det, ch): dict(x, w, data, cv, bf, cv_err, bf_err)} from the PROfit global-fit histogram file (1D RecoE projections, events per GeV)."""
    tag = f2._paths("ds", cfg, variant, stem)[0]
    p = f2.RUN / f"{tag}_PROglobal_histograms.root"
    if not p.exists():
        return {}
    with uproot.open(p) as f:
        h = f["hist1d"].arrays(library="np")
        e = f["errorband"].arrays(library="np")
    pre = np.asarray(h["prefix"]).astype(str)
    epre = np.asarray(e["prefix"]).astype(str)
    out = {}
    for d in range(3):
        for c in range(2):
            def get(P):
                m = (pre == P) & (h["detector"] == d) & (h["channel"] == c)
                if not m.any():
                    return None
                _, o = np.unique(h["bin_index"][m], return_index=True)
                return h["bin_content"][m][o], h["bin_low_edge"][m][o], h["bin_high_edge"][m][o]
            def band(P):
                m = (epre == P) & (e["detector"] == d) & (e["channel"] == c)
                if not m.any():
                    return None
                o = np.argsort(e["point_index"][m])
                return e["error_y_high"][m][o], e["y_value"][m][o]
            dat, bf, cv = get("DATA"), get("BF"), None
            m = (pre == "CV") & (h["detector"] == d) & (h["channel"] == c) & (h["subchannel"] == 4294967295)
            cvb = band("CV"); bfb = band("BF")
            if dat is None or bf is None:
                continue
            cvv = cvb[1] if cvb else None
            out[(d, c)] = dict(lo=dat[1], hi=dat[2], data=dat[0], bf=bf[0], cv=cvv, cv_err=cvb[0] if cvb else None, bf_err=bfb[0] if bfb else None)
    return out


POT = ["1.0e20 POT", "2.0e20 POT", "3.0e20 POT"]      # nominal exposures of contours_base.xml (SBND, ICARUS Run 2, Run 4)


POT = ["1.0e20 POT", "2.0e20 POT", "3.0e20 POT"]      # nominal exposures of contours_base.xml (SBND, ICARUS Run 2, Run 4)


def overlay_pdf(cfg_stem, fit_kind):
    """Data/MC-style pages, ONE fit per plot (mutools.plotting.fit_histograms): injected data, stacked CV prediction by true subchannel,
    prefit band and the best fit of the null ('fix') or the oscillation-floating ('free') fit; one page per detector and channel."""
    import mutools.plotting as mp
    from mutools.plotting import fit_histograms, use_style
    from mutools.plotting.profit import ProfitPlotData
    stem = cfg_stem
    mp.saver.configure(fmt="pdf", dpi=300, rasterized=True, ytick_precision=2)
    use_style("rootlike")
    mp.saver.configure(mark_hadj=0.055)
    label = {"fix": "Best fit (null)", "free": "Best fit (oscillation)"}[fit_kind]
    with PdfPages(OUT / "spectra" / f"{stem}_{'null' if fit_kind == 'fix' else 'osc'}.pdf") as pdf:
        for cfg in ("b", "1p", "Np"):
            if not f2.available("ds", cfg, fit_kind, stem):
                continue
            rd = ProfitPlotData(str(f2.RUN / f"{f2._paths('ds', cfg, fit_kind, stem)[0]}_PROglobal_histograms.root"), scale_by_width=True)
            chans = {"b": [(0, "Single-Proton"), (1, "Multi-Proton")], "1p": [(0, "Single-Proton")], "Np": [(0, "Multi-Proton")]}[cfg]
            figs = fit_histograms(rd, pot=[stem] * 3, channels=chans, fit_label=label, code_version="v3.0.4", selection_version="v1.1.0")   # right watermark: the exact fake dataset name
            for (d, ch), fig in sorted(figs.items(), key=lambda kv: (kv[0][1], kv[0][0])):
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    stems = all_stems()
    stoplight(stems)
    dec = decomposition()
    dec.to_csv(TAB / "chi2_decomposition.csv", index=False)
    b = dec[dec.config == "b"].sort_values("chi2 null", ascending=False).head(25)
    fig, ax = plt.subplots(figsize=(8, 7))
    y = np.arange(len(b))
    ax.barh(y, b["chi2 data+MCstat"], color="C0", label="data + MC-stat part")
    ax.barh(y, b["sum pull^2"], left=b["chi2 data+MCstat"], color="C1", label=r"prior part ($\sum$ pull$^2$)")
    ax.set_yticks(y); ax.set_yticklabels(b.dataset, fontsize=7); ax.invert_yaxis(); ax.set_xscale("log"); ax.legend(fontsize=8)
    ax.set_xlabel(r"null-fit $\chi^2$ (combined analysis)"); fig.savefig(OUT / "chi2_decomp.pdf", bbox_inches="tight"); plt.close(fig)
    a = dec[dec.config == "b"].dropna(subset=["absorbed fraction"]).sort_values("absorbed fraction")
    fig, ax = plt.subplots(figsize=(8, 0.22 * len(a) + 1.5))
    ax.hlines(np.arange(len(a)), 0, a["absorbed fraction"] * 100, color="C0"); ax.plot(a["absorbed fraction"] * 100, np.arange(len(a)), "o", color="C0")
    ax.set_yticks(np.arange(len(a))); ax.set_yticklabels(a.dataset, fontsize=6.5); ax.set_xlabel("fraction of the injection's stat-only chi2 absorbed by the null fit [%]")
    ax.axvline(100, color="0.6", lw=0.6); fig.savefig(OUT / "absorbed_fraction.pdf", bbox_inches="tight"); plt.close(fig)
    for s in PANEL_STEMS:
        with PdfPages(OUT / "spectra" / f"{s}.pdf") as pdf:
            for c in ("b", "1p", "Np"):
                if f2.available("ds", c, "fix", s):
                    spectra_page(pdf, c, s)
    for s_ in PANEL_STEMS:
        for kind in ("fix", "free"):
            overlay_pdf(s_, kind)
    print("wrote", [p.name for p in sorted(OUT.glob("stoplight.pdf"))], len(list((OUT / 'spectra').glob('*.pdf'))), "spectra files")
