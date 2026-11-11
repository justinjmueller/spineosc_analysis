"""Analysis helpers for the DENT data/MC study (see docs/dent_study/PREREGISTRATION.md).

Reads the per-study result JSON written by scripts/dent_study/run_study.py and the PROplot_histograms files, and renders
PROfit-style panels with mutools exactly as plotting/studies/data_access.py does.
"""
import glob, json, pathlib, shutil
import numpy as np, pandas as pd
from scipy.stats import chi2 as chi2_dist

REPO = pathlib.Path(__file__).resolve().parents[2]
RUN = REPO / "runs" / "dent_study"
OBS = {o["name"]: o for o in json.loads((REPO / "scripts/dent_study/observables.json").read_text())}
PRIMARY = json.loads((REPO / "scripts/dent_study/primary_selection.json").read_text())["primary"]
VARIANTS_B = {"v1": "CV, DENT in covariance (standard)", "exDENT": "CV, DENT removed from covariance",
              "inj0p5": "CV + 0.5 x DENT", "inj1p0": "CV + 1.0 x DENT"}

def p_value(chi2, ndf):
    return float(chi2_dist.sf(chi2, ndf))

def load_results():
    """One row per (study, variant, detector, channel)."""
    rows = []
    for f in sorted(glob.glob(str(RUN / "results" / "*.json"))):
        r = json.loads(pathlib.Path(f).read_text())
        tag = r["tag"]
        if tag.startswith("dentA_"):
            tier, name, which = "A", tag.split("_dd_")[1].rsplit("_", 1)[0], tag.rsplit("_", 1)[1]
        elif "_dd_" in tag:
            tier, name, which = "B", tag.split("_dd_")[1], ""
        else:
            tier, name, which = "baseline", r["study"], ""
        for label, lst in r["chi2"].items():
            for x in lst:
                rows.append(dict(tag=tag, tier=tier, det=r["det"], study=name, which=which, variant=label, detector=x["detector"],
                                 channel=x["channel"], chi2=x["chi2"], nbins=x["nbins"], ndf=x["ndf"], p=p_value(x["chi2"], x["ndf"]),
                                 shapeonly_live=r.get("shapeonly_live")))
    return pd.DataFrame(rows)

def load_tierC(key="tierC_pe"):
    """key: "tierC_pe" (pseudo-experiments, the calibrated version) or "tierC" (Poisson-only).
    Paired fake-data throws: one row per (study, seed, truth) with chi2 at test strength 0 and 1 (channel 0 summed over channels)."""
    out = []
    for f in sorted(glob.glob(str(RUN / "results" / "*.json"))):
        r = json.loads(pathlib.Path(f).read_text())
        if not r.get(key): continue
        df = pd.DataFrame(r[key]); df["chi2"] = df["chi2"].apply(sum)
        piv = df.pivot_table(index=["seed", "truth"], columns="test", values="chi2").reset_index()
        piv.columns = ["seed", "truth", "chi2_0", "chi2_1"]
        piv["dchi2"] = piv["chi2_0"] - piv["chi2_1"]; piv["study"] = r["study"]
        out.append(piv)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()

def two_sided(x, obs):
    """Two-sided empirical p of obs under the sample x, with the +1 correction (floor 2/(n+1))."""
    n = len(x); lo = (np.sum(x <= obs) + 1) / (n + 1); hi = (np.sum(x >= obs) + 1) / (n + 1)
    return float(min(1.0, 2 * min(lo, hi)))

def tierC_pvalues(tc, data_dchi2):
    """For one study: the data's dchi2 placed in the fake-data distributions for truth CV (theta=0), DENT x0.5 and DENT x1.0.
    p_* are TWO-sided empirical p-values (floor 2/(n+1)); the one-sided p_cv/p_dent of the plan are kept for reference."""
    g = {t: tc[tc.truth == t]["dchi2"].to_numpy() for t in ("cv", "half", "dent")}
    a, b = g["cv"], g["dent"]; out = dict(n=len(a))
    for t, x in g.items():
        if len(x) > 1: out.update({f"mean_{t}": float(x.mean()), f"sd_{t}": float(x.std(ddof=1)), f"p2_{t}": two_sided(x, data_dchi2)})
    out.update(p_cv=(np.sum(a >= data_dchi2) + 1) / (len(a) + 1), p_dent=(np.sum(b <= data_dchi2) + 1) / (len(b) + 1))
    return out

def holm(pvals):
    """Holm-Bonferroni adjusted p-values for a list; returned in input order."""
    p = np.asarray(pvals, float); order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i]); adj[i] = min(1.0, run)
    return adj.tolist()

# ----------------------------------------------------------------------------------------------- mutools panels
SUBCH = ['Target', r'Other $\nu$ CC', r'$\nu_e$ CC', r'$\nu$ NC', 'Cosmics', 'In-time Cosmics']   # reversed XML order, as in the notebook
COLORS = ['C0', 'C1', 'C2', 'C4', 'C3', 'C5']

def peak_height(path, with_data=True):
    from mutools.plotting.profit import ProfitPlotData, TraceType
    d = ProfitPlotData(path, True); key = '0:0:0:0:total'
    cv = d.get_trace(f'{key}:CV', TraceType.HIST_CONTENTS); band = d.get_trace(f'{key}:CV', TraceType.HIST_ERROR_BAND)
    stack = [cv[:, 3], band[:, 1] + band[:, 3]]
    if with_data:
        dat = d.get_trace(f'{key}:DATA', TraceType.HIST_CONTENTS); stack.append(dat[:, 3] + (dat[:, 4] if dat.shape[1] > 4 else np.sqrt(np.maximum(dat[:, 3], 0))))
    return float(np.nanmax(np.concatenate(stack)))

def render_panel(path, outdir, name, xlabel, xlim=None, chi2_label=r'$\chi^2_{shape}$', legend_loc='upper right'):
    """One PROfit-style data/MC panel (mutools), saved as <outdir>/<name>.png. y-limit = 1.55 x peak, as in the notebook."""
    import matplotlib; matplotlib.use("Agg")
    import mutools.plotting as mp
    from mutools.plotting import run, use_style
    use_style('rootlike'); outdir = pathlib.Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    mp.saver.configure(fmt='png', dpi=200, rasterized=True, ytick_precision=2, mark_hadj=0.055)
    for old in outdir.glob('hist_*_v3.0.4.png'): old.unlink()
    run({'general': {'input': str(path), 'output': str(outdir), 'savefig': True, 'code_version': 'v3.0.4', 'selection_version': 'v1.1.0',
                     'subchannels': SUBCH, 'detectors': ['SBND'], 'channels': ['Inclusive'], 'scale_by_width': 'enabled', 'counter_index': 4},
         'plot': [{'type': 'histogram', 'variable': 0, 'channel': 0, 'detectors': [0], 'colors': COLORS, 'xlabel': xlabel, 'ylabel': 'Area Normalized',
                   'xlim': xlim, 'ylim': [0.0, 1.55 * peak_height(path)], 'rlim': [0.0, 2.0], 'ratio': 'data', 'scale_by_width': 'enabled',
                   'disable_systematics': False, 'band_label': 'Total Uncertainty', 'counter_fmt': '.1%', 'watermark': r'$\bf{SBN}$ Analysis-in-Progress',
                   'show_chi2': True, 'chi2_label': chi2_label, 'show_pvalue': True, 'legend_ncol': 2, 'legend_loc': legend_loc}]})
    staged = sorted(outdir.glob('hist_*_v3.0.4.png'), key=lambda q: q.stat().st_mtime)[-1]
    target = outdir / f'{name}.png'; shutil.move(staged, target); return target

# ----------------------------------------------------------------------------------------------- summary figures
def style():
    import matplotlib; matplotlib.use("Agg")
    from mutools.plotting import use_style
    use_style('rootlike')

def plot_pvalues(df, observables, labels, variants=('v1', 'exDENT', 'inj0p5', 'inj1p0'), colors=('C0', 'C3', 'C2', 'C1'), title='Data Access',
                 savepath=None, ylabel='Shape-only p-value'):
    """p-value per observable for several prediction variants, notebook style (log y, p = 0.05 line)."""
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis
    style()
    fig, ax = plt.subplots(figsize=(11, 5)); x = np.arange(len(observables))
    names = {'v1': 'CV (DENT in covariance)', 'exDENT': 'CV (DENT removed)', 'inj0p5': 'CV + 0.5 DENT', 'inj1p0': 'CV + 1.0 DENT'}
    for v, c, dx in zip(variants, colors, np.linspace(-0.18, 0.18, len(variants))):
        sub = df[df.variant == v].set_index('study').reindex(observables)
        ax.scatter(x + dx, sub['p'].clip(lower=1e-12), s=40, color=c, label=names.get(v, v), zorder=3)
    ax.axhline(0.05, color='#eb6834', lw=1.2, ls='--', zorder=2, label='p = 0.05')
    ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha='right'); ax.set_ylabel(ylabel)
    pmin = float(df[df.variant.isin(variants)].p.clip(lower=1e-300).min()); ax.set_xlim(-0.6, len(observables) - 0.4); ax.set_ylim(min(pmin / 5, 1e-2), 3.0)
    ax.legend(loc='lower left', ncol=2, fontsize=9)
    mark_axis(ax, rf'$\bf{{SBND}}$ {title}', vadj=0.05)
    if savepath: fig.savefig(savepath, bbox_inches='tight', dpi=200)
    return fig

def plot_dchi2(table, savepath=None):
    """table: DataFrame with study, label, dchi2_data, mean_cv, sd_cv, mean_dent, sd_dent. Data point vs the Tier C +-1 sd bands."""
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis
    style(); fig, ax = plt.subplots(figsize=(11, 5)); x = np.arange(len(table))
    ax.errorbar(x - 0.12, table.mean_cv, yerr=table.sd_cv, fmt='s', color='C0', ms=6, capsize=3, label='fake data from CV (mean $\\pm$ 1 sd)', zorder=3)
    ax.errorbar(x + 0.12, table.mean_dent, yerr=table.sd_dent, fmt='D', color='C1', ms=6, capsize=3, label='fake data from CV + DENT (mean $\\pm$ 1 sd)', zorder=3)
    ax.scatter(x, table.dchi2_data, marker='o', s=70, color='k', label='data', zorder=4)
    ax.axhline(0, color='0.5', lw=1); ax.set_yscale('symlog', linthresh=10)
    ax.set_xticks(x); ax.set_xticklabels(table.label, rotation=35, ha='right'); ax.set_ylabel(r'$\Delta\chi^2 = \chi^2(\theta{=}0) - \chi^2(\theta{=}1)$')
    ax.legend(loc='upper right', fontsize=9); mark_axis(ax, r'$\bf{SBND}$ Data Access', vadj=0.05)
    if savepath: fig.savefig(savepath, bbox_inches='tight', dpi=200)
    return fig

def plot_tierC_hist(tc, data_dchi2, title, ax=None):
    import matplotlib.pyplot as plt
    style(); ax = ax or plt.subplots(figsize=(5.5, 4))[1]
    a = tc[tc.truth == 'cv'].dchi2; b = tc[tc.truth == 'dent'].dchi2
    lo, hi = min(a.min(), b.min(), data_dchi2), max(a.max(), b.max(), data_dchi2); pad = 0.05 * (hi - lo + 1e-9)
    bins = np.linspace(lo - pad, hi + pad, 21)
    ax.hist(a, bins=bins, histtype='step', lw=1.8, color='C0', label='truth: CV'); ax.hist(b, bins=bins, histtype='step', lw=1.8, color='C1', label='truth: CV + DENT')
    ax.axvline(data_dchi2, color='k', lw=2, label='data'); ax.set_xlabel(r'$\Delta\chi^2$'); ax.set_ylabel('Fake experiments'); ax.set_title(title, loc='left', fontsize=10)
    return ax

# ----------------------------------------------------------------------------------------------- display helpers (notebook or static build)
_SINK = None            # set to a list by build_report.py; in a notebook it stays None and output is displayed inline

def show_table(df, caption='', floatfmt=None):
    if _SINK is not None:
        _SINK.append(('table', caption, df.to_html(index=False, border=0, float_format=floatfmt or (lambda v: f'{v:.4g}'), classes='tbl'))); return
    from IPython.display import display, HTML
    if caption: display(HTML(f'<b>{caption}</b>'))
    display(df)

def show_figure(path, caption=''):
    if _SINK is not None:
        _SINK.append(('figure', caption, str(path))); return
    from IPython.display import display, Image, HTML
    if caption: display(HTML(f'<i>{caption}</i>'))
    display(Image(filename=str(path)))

def plot_tierA(tab, savepath=None, title='Tier A'):
    """tab: study, label, chi2_cv, chi2_dent, ndf.  Horizontal bars of chi2_cv - chi2_dent (positive: DENT sample closer to data)."""
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis
    style(); tab = tab.sort_values('dchi2'); fig, ax = plt.subplots(figsize=(9, 0.42 * len(tab) + 1.8)); y = np.arange(len(tab))
    ax.barh(y, tab.dchi2, color=['C1' if v > 0 else 'C0' for v in tab.dchi2], height=0.62)
    ax.set_xscale('symlog', linthresh=10); ax.axvline(0, color='0.3', lw=1); ax.set_yticks(y); ax.set_yticklabels(tab.label)
    ax.set_xlabel(r'$\chi^2_{shape}$(CV) $-$ $\chi^2_{shape}$(DENT)' + '\n[positive: DENT sample closer to data]')
    for yi, (v, a, b, n) in enumerate(zip(tab.dchi2, tab.chi2_cv, tab.chi2_dent, tab.ndf)):
        ax.text(-1.5 if v > 0 else 1.5, yi, f'{a:.0f} → {b:.0f}  (ndf {n})', va='center', ha='right' if v > 0 else 'left', fontsize=8)
    mark_axis(ax, rf'$\bf{{SBND}}$ {title}', vadj=0.05)
    if savepath: fig.savefig(savepath, bbox_inches='tight', dpi=200)
    return fig
