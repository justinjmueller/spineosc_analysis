#!/usr/bin/env python3
"""Build (and execute) plotting/fakedata_channel_studies.ipynb.

Technote-style fake-data fits (null vs oscillation fit, Delta chi2, pulls) for every spinefd dial and the native
IceCube injections, repeated for three analysis configurations: both channels, single-proton only, multi-proton only.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/build_fakedata_channel_notebook.py [--no-exec]
"""

import argparse
import sys
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_fakedata_injection_nb import execute  # noqa: E402  (the shared cell executor)

NB = HERE.parents[0] / 'fakedata_channel_studies.ipynb'
CELLS = []


def md(t):
    CELLS.append(('md', t.strip('\n')))


def code(t):
    CELLS.append(('code', t.strip('\n')))


md(r"""
# Fake-data fits and the single-channel question

This notebook follows the fake-data section of the technote (`spineosc_technote/sections/fake_data_studies.tex`) and
extends it. What the technote does: each fake dataset is fit twice, a **null fit** (nuisance parameters only) and an
**oscillation fit** ($\Delta m^2_{41}$ and $\sin^2 2\theta_{\mu\mu}$ also float); the test statistic is
$\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm osc}$ with Wilks $p$-values for two degrees of freedom (to be replaced by
Feldman-Cousins); post-fit nuisance pulls and post-fit spectra are shown for both fits. Its studies are a Poisson
fluctuation of the CV, $\pm1\sigma$ shifts of the QE interference and spectral-function models, and the three IceCube
2024 oscillation points. (Its numbers predate the v1.1.0 selection; `xml/MANIFEST.toml` notes those reruns are pending.)

Here the same two fits are run on every fake dataset made with the reBruce dials (`spinefd`) and on PROfit's native
IceCube injections, for **three analysis configurations**: both channels, **single-proton only** (`ccqelike`) and
**multi-proton only** (`nproton`). The single-channel analyses are separate XMLs (the other channel removed from the
`xml/systematics/systematics_suite.xml` structure), processed and fit exactly like the two-channel one, with the **full
systematics suite** on the CV. Fake data is the CV MC reweighted by the dial at +1 sigma (`--use-fake-data
--inject-systs`), or by PROfit's own oscillation injection, or a Poisson throw of the CV; the injected dial is excluded
from the fit. No real beam data is used anywhere.
""")

code(r'''
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, '/Users/mueller/Projects/spineosc/plotting/scripts')
import fd_fits as fd
from mutools.plotting import use_style

use_style('rootlike')
REPO = Path('/Users/mueller/Projects/spineosc')
OUT = REPO / 'plotting' / 'fakedata_channel_figures'
OUT.mkdir(parents=True, exist_ok=True)

CFG = fd.CONFIGS                      # b / 1p / Np
CFG_COLOR = {'b': 'k', '1p': 'C0', 'Np': 'C3'}

# every injection with a label; stems are the PROfit output stems used by fitjob.sh
LABELS = {
    'ctrl': 'No injection (CV)',
    'phys_lo': 'IceCube low $\\Delta m^2$ (PROfit --inject)', 'phys_best': 'IceCube best fit (PROfit --inject)',
    'phys_hi': 'IceCube high $\\Delta m^2$ (PROfit --inject)',
    'osc_ic2024_pebl_dm2lo': 'IceCube low $\\Delta m^2$ (dial, per-event L)',
    'osc_ic2024_pebl_bestfit': 'IceCube best fit (dial, per-event L)',
    'osc_ic2024_pebl_dm2hi': 'IceCube high $\\Delta m^2$ (dial, per-event L)',
    'osc_ic2024_dm2lo': 'IceCube low $\\Delta m^2$ (dial, L = 110 m)', 'osc_ic2024_bestfit': 'IceCube best fit (dial, L = 110 m)',
    'osc_ic2024_dm2hi': 'IceCube high $\\Delta m^2$ (dial, L = 110 m)',
    'osc_ic2024_unf_dm2lo': 'IceCube low $\\Delta m^2$ (unfolded dial)', 'osc_ic2024_unf_bestfit': 'IceCube best fit (unfolded dial)',
    'osc_ic2024_unf_dm2hi': 'IceCube high $\\Delta m^2$ (unfolded dial)',
    'mec_bdt': 'MEC BDT', 'qe_zexp_mva_to_lqcd': 'QE axial FF MINERvA$\\to$LQCD', 'ccqe_xsec_corr': 'CCQE bugfix',
    'pi_fsi_ha2025': 'Pion FSI hA2025', 'jaesung_lowq2_pi_enhancement_postfsi': 'Low-$Q^2$ pion enh. (post-FSI)',
    'jaesung_lowq2_pi_enhancement_prefsi': 'Low-$Q^2$ pion enh. (pre-FSI)',
    'flux_horn_nightmare': 'Horn nightmare', 'flux_horn_171p5kA': 'Horn 171.5 kA', 'flux_hadprod_sw_piplus_u387': 'Hadron production u387',
    'minerva_3dqelike_bnb': 'MINERvA 3D QE-like', 'minerva_3dqelike_bnb_pzmarg': 'MINERvA QE-like (pz marg.)',
}
for k in range(1, 11):
    LABELS[f'pois{k:02d}'] = f'Poisson throw {k}'


def label(stem):
    if stem in LABELS:
        return LABELS[stem]
    return stem.replace('_', ' ')


JOBS = [l.split()[2] for l in (fd.RUN / 'jobs_fits.txt').read_text().splitlines()]
STEMS = list(dict.fromkeys(JOBS))                   # injection order of the queue
POIS = [s for s in STEMS if s.startswith('pois')]
INJ = [s for s in STEMS if not s.startswith('pois')]
done = {c: sum(fd.available(c, v, s) for s in STEMS for v in ('free', 'fix')) for c in CFG}
print('completed fits per configuration (of', 2 * len(STEMS), 'each):', done)
''')

md(r"""
## 1. Delta chi2 of every injection, per analysis configuration

$\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm osc}$ and its Wilks $p$-value (2 dof), as in the technote tables. The
oscillation fit also reports the fitted point $(\Delta m^2_{41},\sin^2 2\theta_{\mu\mu})$.
""")

code(r'''
def table(cfg, stems):
    rows = []
    for s in stems:
        if not (fd.available(cfg, 'free', s) and fd.available(cfg, 'fix', s)):
            continue
        n, o = fd.read_fit(cfg, 'fix', s), fd.read_fit(cfg, 'free', s)
        d = n['chi2'] - o['chi2']
        rows.append({'injection': s, 'chi2 null': n['chi2'], 'chi2 osc': o['chi2'], 'dchi2': d,
                     'p (Wilks, 2 dof)': stats.chi2.sf(max(d, 0), 2), 'fit dm2 [eV2]': o['dm2'], 'fit sin2 2th': o['s22']})
    return pd.DataFrame(rows)


tabs = {c: table(c, INJ) for c in CFG}
with pd.option_context('display.width', 220, 'display.max_rows', 200, 'display.float_format', lambda v: f'{v:.4g}'):
    for c in CFG:
        print(f'=== {CFG[c]} ===')
        print(tabs[c].to_string(index=False))
        print()
''')

md(r"""
### Cross-check against the technote's IceCube numbers

The technote quotes $\Delta\chi^2$ = 26.85, 21.57, 33.14 for the three IceCube points (older selection, different
production). Ours, for PROfit's native injection with both channels, are below; the floating fit recovers the injected
point exactly in every configuration where it is sensitive.
""")

code(r'''
tech = {'phys_lo': 26.85, 'phys_best': 21.57, 'phys_hi': 33.14}
rows = []
for s, t in tech.items():
    r = {'point': label(s), 'technote': t}
    for c in CFG:
        r[CFG[c]] = fd.delta_chi2(c, s)
    rows.append(r)
print(pd.DataFrame(rows).round(2).to_string(index=False))
''')

md(r"""
## 2. Poisson fluctuations: the null distribution of Delta chi2

Ten Poisson throws of the CV per configuration. If Wilks' theorem held, $\Delta\chi^2$ would follow a $\chi^2$
distribution with 2 degrees of freedom (CDF in black). Ten throws only give a coarse check, enough to see whether
the Wilks $p$-values in the technote are reasonable, and what the single-channel analyses do to the distribution.
""")

code(r'''
fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
rows = []
for c in CFG:
    d = np.array([fd.delta_chi2(c, s) for s in POIS])
    d = d[np.isfinite(d)]
    n0 = np.array([fd.read_fit(c, 'fix', s)['chi2'] for s in POIS if fd.available(c, 'fix', s)])
    if len(d) == 0:
        continue
    xs = np.sort(d)
    axes[0].step(xs, np.arange(1, len(xs) + 1) / len(xs), where='post', color=CFG_COLOR[c], label=f'{CFG[c]} (n = {len(xs)})')
    ks = stats.kstest(d, stats.chi2(2).cdf)
    rows.append({'config': CFG[c], 'throws': len(d), 'mean dchi2': d.mean(), 'max dchi2': d.max(),
                 'mean chi2 null': n0.mean(), 'KS p vs chi2(2)': ks.pvalue})
    axes[1].hist(n0, bins=8, histtype='step', color=CFG_COLOR[c], label=CFG[c])
xx = np.linspace(0, 12, 200)
axes[0].plot(xx, stats.chi2(2).cdf(xx), 'k:', label=r'$\chi^2$, 2 dof')
axes[0].set_xlabel(r'$\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm osc}$'); axes[0].set_ylabel('cumulative fraction'); axes[0].legend(fontsize=9)
axes[1].set_xlabel(r'null-fit $\chi^2$'); axes[1].set_ylabel('throws'); axes[1].legend(fontsize=9)
fig.tight_layout(); fig.savefig(OUT / 'poisson_null.pdf', bbox_inches='tight')
print(pd.DataFrame(rows).round(3).to_string(index=False))
''')

md(r"""
## 3. How much does each channel constrain the systematics?

With the CV as fake data (no pulls expected), the width of each nuisance parameter's posterior measures how well the
data constrain it; a width of 1 is the prior. Below, the posterior widths of the nuisance parameters in the
nuisance-only fit for the three configurations.
""")

code(r'''
def widths(cfg, stem='ctrl', variant='fix'):
    f = fd.read_fit(cfg, variant, stem)
    return pd.Series({n: f['sigma'].get(n, np.nan) for n in fd.nuisance_names(f)})


W = pd.DataFrame({c: widths(c) for c in CFG if fd.available(c, 'fix', 'ctrl')})
print('nuisance parameters:', len(W))
for c in [c for c in ('b', '1p', 'Np') if c in W]:
    print(f'  {CFG.get(c, c):20s} constrained below 0.5 sigma: {(W[c] < 0.5).sum():3d}   below 0.8 sigma: {(W[c] < 0.8).sum():3d}   median width {W[c].median():.2f}')
REF = 'b' if 'b' in W else W.columns[0]
top = W.sort_values(REF).head(25)
with pd.option_context('display.width', 200, 'display.float_format', lambda v: f'{v:.2f}'):
    print('\nMost constrained parameters (posterior width, prior = 1):')
    print(top[[c for c in ('b', '1p', 'Np') if c in W]].rename(columns=CFG).to_string())
''')

code(r'''
fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
ax = axes[0]
ax.scatter(W['1p'], W['Np'], s=14, color='C2')
ax.plot([0, 1.2], [0, 1.2], 'k--', lw=0.8)
for n in W.assign(d=(W['1p'] - W['Np']).abs()).sort_values('d').tail(7).index:
    ax.annotate(n, (W.loc[n, '1p'], W.loc[n, 'Np']), fontsize=7, xytext=(3, 3), textcoords='offset points')
ax.set_xlabel('posterior width, single-proton only'); ax.set_ylabel('posterior width, multi-proton only')
ax.set_title('Nuisance constraints: one channel vs the other')
ax = axes[1]
order = W[REF].sort_values().index
x = np.arange(len(order))
for c in [c for c in ('b', '1p', 'Np') if c in W]:
    ax.plot(x, W.loc[order, c].values, '.', color=CFG_COLOR[c], label=CFG[c], ms=5)
ax.axhline(1, color='0.5', ls='--', lw=0.8)
ax.set_xlabel('nuisance parameter, sorted by two-channel constraint'); ax.set_ylabel('posterior width [prior sigma]')
ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(OUT / 'constraint_widths.pdf', bbox_inches='tight')
''')

md(r"""
## 4. Which fake data can each analysis see?

$\chi^2_{\rm null}$ (the fake data against the best nuisance-only fit) and $\Delta\chi^2$ (the extra the oscillation fit
buys) for every injection and configuration. Larger means more detectable.
""")

code(r'''
def matrix(kind):
    M = np.full((len(INJ), 3), np.nan)
    for i, s in enumerate(INJ):
        for j, c in enumerate(CFG):
            if kind == 'dchi2':
                M[i, j] = fd.delta_chi2(c, s)
            elif fd.available(c, 'fix', s):
                M[i, j] = fd.read_fit(c, 'fix', s)['chi2']
    return M


fig, axes = plt.subplots(1, 2, figsize=(14, max(6, 0.3 * len(INJ))), sharey=True)
for ax, kind, title in zip(axes, ('chi2null', 'dchi2'), (r'null-fit $\chi^2$', r'$\Delta\chi^2$ (null $-$ oscillation fit)')):
    M = matrix(kind)
    im = ax.imshow(np.log10(np.clip(M, 1e-3, None)), aspect='auto', cmap='viridis', vmin=-2, vmax=2.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f'{M[i, j]:.2g}', ha='center', va='center', fontsize=7, color='w' if np.log10(max(M[i, j], 1e-3)) < 1 else 'k')
    ax.set_xticks(range(3)); ax.set_xticklabels([CFG[c] for c in CFG], rotation=20, ha='right', fontsize=9)
    ax.set_title(title)
axes[0].set_yticks(range(len(INJ))); axes[0].set_yticklabels([label(s) for s in INJ], fontsize=7)
fig.tight_layout(); fig.savefig(OUT / 'detectability_heatmap.pdf', bbox_inches='tight')
''')

md(r"""
## 5. Recovering the injected oscillation point

For the native IceCube injections the oscillation fit should return the injected $(\Delta m^2_{41},\sin^2 2\theta_{\mu\mu})$
with $\chi^2\approx0$. The same three points injected through the reBruce dial with a fixed 110 m baseline, and with the
event's own baseline, show how much the baseline treatment matters.
""")

code(r'''
inj_points = {'phys_lo': (2.07, 0.116), 'phys_best': (3.5, 0.16), 'phys_hi': (12.1, 0.187),
              'osc_ic2024_pebl_dm2lo': (2.07, 0.116), 'osc_ic2024_pebl_bestfit': (3.5, 0.16), 'osc_ic2024_pebl_dm2hi': (12.1, 0.187),
              'osc_ic2024_dm2lo': (2.07, 0.116), 'osc_ic2024_bestfit': (3.5, 0.16), 'osc_ic2024_dm2hi': (12.1, 0.187)}
rows = []
for s, (dm2, s22) in inj_points.items():
    for c in CFG:
        if fd.available(c, 'free', s) and fd.available(c, 'fix', s):
            o = fd.read_fit(c, 'free', s)
            rows.append({'injection': label(s), 'config': CFG[c], 'injected dm2': dm2, 'fit dm2': o['dm2'], 'injected s22': s22,
                         'fit s22': o['s22'], 'chi2 osc': o['chi2'], 'dchi2': fd.delta_chi2(c, s)})
print(pd.DataFrame(rows).round(3).to_string(index=False))
''')

md(r"""
## 6. Pulls: what the fit does to the nuisance parameters

Post-fit nuisance values in units of the prior sigma (error bars: MCMC posterior width), for the nuisance-only fit,
overlaying the three configurations. Where the two single-channel analyses disagree on a pull, the fake data would bias that
systematic differently depending on which channel is used.
""")

code(r'''
def pull_overlay(stem, variant='fix', save=True, show=False, ntop=6):
    avail = [c for c in CFG if fd.available(c, variant, stem)]
    if not avail:
        return
    fits = {c: fd.read_fit(c, variant, stem) for c in avail}
    names = fd.nuisance_names(fits[avail[0]])
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(18, 5.2))
    ax.axhspan(-1, 1, color='0.9', zorder=0)
    ax.axhline(0, color='k', lw=0.8)
    for c, off in zip(avail, np.linspace(-0.22, 0.22, len(avail))):
        y = np.array([fits[c]['values'].get(n, np.nan) for n in names])
        e = np.array([fits[c]['sigma'].get(n, np.nan) for n in names])
        ax.errorbar(x + off, y, yerr=e, fmt='o', ms=3, lw=0.9, color=CFG_COLOR[c], label=f"{CFG[c]} ($\\chi^2$ = {fits[c]['chi2']:.1f})")
    ymax = np.nanmax([np.abs(fits[c]['values'].get(n, 0.0)) for c in avail for n in names])
    for c in avail:
        yy = np.array([fits[c]['values'].get(n, 0.0) for n in names])
        for i in np.argsort(-np.abs(yy))[:2]:
            ax.annotate(names[i], (x[i], yy[i]), xytext=(0, 7 if yy[i] >= 0 else -12), textcoords='offset points', ha='center', fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=90, fontsize=6.5); ax.set_xlim(-1, len(names))
    ax.set_ylabel(r'post-fit value [prior $\sigma$]')
    ax.set_title(f"{label(stem)}  ({'oscillation fixed' if variant == 'fix' else 'oscillation floating'})", fontsize=12)
    ax.legend(loc='upper right', fontsize=9, ncol=3)
    if save:
        (OUT / 'pulls').mkdir(exist_ok=True)
        fig.savefig(OUT / 'pulls' / f'pulls_{variant}_{stem}.pdf', bbox_inches='tight')
    if not show:
        plt.close(fig)
    return fits


_ = pull_overlay('phys_best', show=True)
''')

code(r'''
_ = pull_overlay('horn', show=True)
''')

code(r'''
for s in INJ:
    for v in ('fix', 'free'):
        pull_overlay(s, variant=v)
print('pull overlays written to', OUT / 'pulls')
''')

md(r"""
### Single-channel pull bias

For each injection, the post-fit pull of every nuisance parameter in the single-proton-only analysis against the
multi-proton-only analysis (nuisance-only fits). Points off the diagonal are systematics whose best-fit value depends on
which channel is used; the table lists the largest differences in units of the quadrature sum of the posterior widths.
""")

code(r'''
def pull_pair(stem, variant='fix'):
    if not (fd.available('1p', variant, stem) and fd.available('Np', variant, stem)):
        return None
    a, b = fd.read_fit('1p', variant, stem), fd.read_fit('Np', variant, stem)
    names = fd.nuisance_names(a)
    df = pd.DataFrame({'1p': [a['values'].get(n, np.nan) for n in names], 'Np': [b['values'].get(n, np.nan) for n in names],
                       's1p': [a['sigma'].get(n, np.nan) for n in names], 'sNp': [b['sigma'].get(n, np.nan) for n in names]}, index=names)
    df['diff'] = (df['1p'] - df['Np']) / np.hypot(df['s1p'], df['sNp'])
    return df


sel = [s for s in ('phys_best', 'phys_lo', 'phys_hi', 'horn', 'mec_bdt', 'qe_zexp_mva_to_lqcd', 'flux_hadprod_sw_piplus_u387', 'ub_cc1p0pi_dpt_dat') if pull_pair(s) is not None]
fig, axes = plt.subplots(2, 4, figsize=(18, 8.5))
for ax, s in zip(axes.ravel(), sel):
    df = pull_pair(s)
    ax.scatter(df['1p'], df['Np'], s=14, color='C4', alpha=0.8)
    lim = max(1.2, np.nanmax(np.abs(df[['1p', 'Np']].values)) * 1.1)
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=0.8)
    for n in df['diff'].abs().sort_values().tail(3).index:
        ax.annotate(n, (df.loc[n, '1p'], df.loc[n, 'Np']), fontsize=9, xytext=(4, 4), textcoords='offset points')
    ax.set_xlabel('pull, single-proton only'); ax.set_ylabel('pull, multi-proton only'); ax.set_title(label(s), fontsize=11)
fig.tight_layout(); fig.savefig(OUT / 'pull_bias_scatter.pdf', bbox_inches='tight')
rows = []
for s in INJ:
    df = pull_pair(s)
    if df is None:
        continue
    j = df['diff'].abs().idxmax()
    rows.append({'injection': s, 'max |1p-Np| [sigma]': abs(df.loc[j, 'diff']), 'parameter': j, 'pull 1p': df.loc[j, '1p'], 'pull Np': df.loc[j, 'Np'],
                 'max|pull| 1p': np.nanmax(np.abs(df['1p'])), 'max|pull| Np': np.nanmax(np.abs(df['Np']))})
tab = pd.DataFrame(rows)
with pd.option_context('display.width', 220, 'display.max_rows', 100, 'display.float_format', lambda v: f'{v:.3g}'):
    print(tab.sort_values('max |1p-Np| [sigma]', ascending=False).to_string(index=False))
''')

md(r"""
## 7. Summary

Detectability of each fake dataset (nuisance-only fit $\chi^2$ and $\Delta\chi^2$) and the largest post-fit pull, per analysis
configuration.
""")

code(r'''
rows = []
for s in INJ:
    r = {'injection': s}
    for c in CFG:
        if fd.available(c, 'fix', s) and fd.available(c, 'free', s):
            f = fd.read_fit(c, 'fix', s)
            r[f'chi2 null {c}'] = f['chi2']
            r[f'dchi2 {c}'] = fd.delta_chi2(c, s)
            r[f'max|pull| {c}'] = max(abs(f['values'][n]) for n in fd.nuisance_names(f))
    rows.append(r)
summary = pd.DataFrame(rows)
summary.to_csv(OUT / 'summary.csv', index=False)
with pd.option_context('display.width', 250, 'display.max_rows', 100, 'display.float_format', lambda v: f'{v:.3g}'):
    print(summary.to_string(index=False))
''')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-exec', action='store_true')
    args = ap.parse_args()
    nb = new_notebook()
    nb.cells = [new_markdown_cell(t) if k == 'md' else new_code_cell(t) for k, t in CELLS]
    if not args.no_exec:
        execute(nb)
    nbformat.write(nb, NB)
    print('wrote', NB)


if __name__ == '__main__':
    main()
