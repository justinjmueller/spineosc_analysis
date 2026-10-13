#!/usr/bin/env python3
"""Build (and execute) plotting/fakedata_2d_studies.ipynb from the 2D fake-data campaign (runs/fakedata_2d).

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/build_fakedata_2d_notebook.py [--no-exec]

Writes figures to plotting/fakedata_2d_figures/, tables to plotting/fakedata_2d_figures/tables/ and the numbers quoted in the
report to tables/headline.json.
"""

import argparse
import sys
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_fakedata_injection_nb import execute  # noqa: E402  (shared cell executor)

NB = HERE.parent / 'fakedata_2d_studies.ipynb'
CELLS = []


def md(t):
    CELLS.append(('md', t.strip('\n')))


def code(t):
    CELLS.append(('code', t.strip('\n')))


md(r"""
# Fake-data fits in the 2D configuration: absorption by the nuisance parameters, with and without the oscillation parameters,
# for the combined analysis, each channel alone, and a single merged histogram

Analyses (all `xml/contours/contours_base.xml`, full systematics suite, 2D `RecoE;pn` in the single-proton channel):
**b** = single-proton and multi-proton combined (both channels fit together), **1p** = single-proton channel only,
**Np** = multi-proton channel only, **m** = merged (one inclusive 2D `RecoE;pn` histogram, no split by proton multiplicity). Each fake dataset is fit with the oscillation parameters fixed (null) and floating
(oscillation), and $\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm osc}$. Fake data are MC only (no real beam data anywhere):
the CV MC reweighted by a reBruce dial (separate data XML, so the nuisance parameters have no knob for it), or PROfit's own
fake data (`--use-fake-data`) for the IceCube points, Poisson throws and the technote-style QE-interference and
spectral-function injections. See `docs/fakedata_review/STUDY_DESIGN.md`.
""")

code(r'''
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, '/Users/mueller/Projects/spineosc/plotting/scripts')
import fd2d as f2
from mutools.plotting import use_style

use_style('rootlike')
REPO = Path('/Users/mueller/Projects/spineosc')
OUT = REPO / 'plotting' / 'fakedata_2d_figures'
TAB = OUT / 'tables'
OUT.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
CFG, CC = f2.CONFIGS, f2.CFG_COLOR
HEAD = {}                               # numbers quoted in the report

jobs = [l.split() for l in (f2.RUN / 'jobs_2d.txt').read_text().splitlines()]
DS = [j[2] for j in jobs if j[0] == 'ds']
DS = list(dict.fromkeys(DS))
FD = list(dict.fromkeys(j[3] for j in jobs if j[0] == 'fd'))
POIS = [s for s in FD if s.startswith('pois')]
PHYS = ['phys_lo', 'phys_best', 'phys_hi']
TECH = [s for s in FD if s.startswith(('qeint', 'specf'))]
DIALS = [s for s in DS if s != 'ctrl']
ICE = {'phys_lo': (2.07, 0.116), 'phys_best': (3.5, 0.16), 'phys_hi': (12.1, 0.187)}
TECHNOTE = {'poisson': (120.8392, 119.5703, 1.2689), 'qeint_p1': (0.6064, 0.6059, 4.51e-4), 'qeint_m1': (0.7188, 0.6306, 0.0882),
            'specf_p1': (1.33e-4, 1.22e-4, 1.12e-5), 'specf_m1': (1.23e-4, 1.07e-4, 1.63e-5)}
TECHNOTE_ICE = {'phys_lo': 26.85, 'phys_best': 21.57, 'phys_hi': 33.14}
print('datasets via data XML:', len(DS), '| fd stems:', len(FD), '| expected job lines:', len(jobs))
''')

md(r"""
## 0. Integrity of the campaign

Every job's return codes (code 99 = an injected systematic was not confirmed in the log, i.e. the run silently fell back to
plain Asimov). A fit with $\Delta\chi^2<0$ means the oscillation fit found a worse minimum than the null fit (minimizer failure).
""")

code(r'''
st = f2.job_status()
bad = f2.failures()
print('completed job markers:', len(st), 'of', len(jobs), '| with a nonzero return code:', len(bad))
if len(bad):
    print(bad.to_string(index=False))
neg = []
for route, stems in (('fd', FD), ('ds', DS)):
    for cfg in CFG:
        for s in stems:
            d = f2.delta_chi2(route, cfg, s)
            if np.isfinite(d) and d < -1e-3:
                neg.append((route, cfg, s, d))
print('fits with Delta chi2 < 0 (oscillation fit worse than null):', len(neg))
for n in neg[:20]:
    print('  ', n)
HEAD['n_jobs_done'] = int(len(st)); HEAD['n_jobs_expected'] = int(len(jobs)); HEAD['n_failed'] = int(len(bad)); HEAD['n_negative_dchi2'] = len(neg)
''')

md(r"""
## 1. Closure and cross-validation of the two routes

* The no-injection control must give $\chi^2\approx0$ in both routes and all analyses.
* PROfit's native IceCube injection (`--use-fake-data --inject`) and the per-event-baseline dial `osc_ic2024_pebl_*` (separate data
  XML) describe the same physics and should agree. The fixed-110 m dial `osc_ic2024_*` is expected to differ (ICARUS).
""")

code(r'''
rows = []
for cfg in CFG:
    for route, stem in (('fd', 'ctrl'), ('ds', 'ctrl')):
        p = f2.pair(route, cfg, stem)
        if p:
            rows.append({'route': route, 'config': CFG[cfg], 'chi2 null': p[0]['chi2'], 'chi2 osc': p[1]['chi2']})
clos = pd.DataFrame(rows)
print(clos.to_string(index=False)); clos.to_csv(TAB / 'closure.csv', index=False)

pairs = [('phys_lo', 'osc_ic2024_pebl_dm2lo', 'osc_ic2024_dm2lo'), ('phys_best', 'osc_ic2024_pebl_bestfit', 'osc_ic2024_bestfit'),
         ('phys_hi', 'osc_ic2024_pebl_dm2hi', 'osc_ic2024_dm2hi')]
rows = []
for nat, pebl, fixed in pairs:
    for cfg in CFG:
        r = {'point': nat, 'config': CFG[cfg], 'dchi2 native (fd)': f2.delta_chi2('fd', cfg, nat),
             'dchi2 per-event-L dial (ds)': f2.delta_chi2('ds', cfg, pebl), 'dchi2 fixed-110m dial (ds)': f2.delta_chi2('ds', cfg, fixed)}
        pn, pp = f2.pair('fd', cfg, nat), f2.pair('ds', cfg, pebl)
        r['chi2 null native'] = pn[0]['chi2'] if pn else np.nan
        r['chi2 null per-event-L dial'] = pp[0]['chi2'] if pp else np.nan
        rows.append(r)
xr = pd.DataFrame(rows)
with pd.option_context('display.width', 250, 'display.float_format', lambda v: f'{v:.4g}'):
    print(xr.to_string(index=False))
xr.to_csv(TAB / 'route_crosscheck.csv', index=False)
HEAD['route_ratio_b'] = [float(r['dchi2 per-event-L dial (ds)'] / r['dchi2 native (fd)']) for _, r in xr[xr.config == CFG['b']].iterrows() if np.isfinite(r['dchi2 native (fd)']) and r['dchi2 native (fd)'] > 0]
''')

md(r"""
## 2. Comparison with the technote's fake-data studies (combined analysis)

The technote quotes one Poisson throw, QE interference $\pm1\sigma$ and spectral function $\pm1\sigma$ with the injected knob left
in the fit, and the three IceCube points. Below, the same studies in the current production (v1.1.0 selection; the technote
numbers predate it), in the combined analysis, plus the variants with the injected knob excluded from the fit.
Poisson: mean and standard deviation over ten throws (seeds 1-10) rather than one throw.
""")

code(r'''
def line(label, cfg, stem, route='fd'):
    p = f2.pair(route, cfg, stem)
    if p is None:
        return None
    return {'study': label, 'chi2 null': p[0]['chi2'], 'chi2 osc': p[1]['chi2'], 'dchi2': p[0]['chi2'] - p[1]['chi2']}


rows = []
for cfg in ('b',):
    pois = [f2.pair('fd', cfg, s) for s in POIS]
    pois = [p for p in pois if p]
    if pois:
        n0 = np.array([p[0]['chi2'] for p in pois]); n1 = np.array([p[1]['chi2'] for p in pois])
        rows.append({'study': f'Poisson throws (n={len(pois)}, mean)', 'chi2 null': n0.mean(), 'chi2 osc': n1.mean(), 'dchi2': (n0 - n1).mean(),
                     'technote chi2 null': TECHNOTE['poisson'][0], 'technote dchi2': TECHNOTE['poisson'][2]})
    for s, lab in (('qeint_p1', 'QE interference +1 sigma'), ('qeint_m1', 'QE interference -1 sigma'), ('specf_p1', 'Spectral function +1 sigma'), ('specf_m1', 'Spectral function -1 sigma')):
        for var, tag in (('infit', 'knob in fit (as technote)'), ('excl', 'knob excluded')):
            r = line(f'{lab}, {tag}', cfg, f'{s}_{var}')
            if r:
                if var == 'infit':
                    r['technote chi2 null'], r['technote dchi2'] = TECHNOTE[s][0], TECHNOTE[s][2]
                rows.append(r)
    for s, lab in (('phys_lo', 'IceCube dm2 = 2.07, s22 = 0.116'), ('phys_best', 'IceCube dm2 = 3.5, s22 = 0.16'), ('phys_hi', 'IceCube dm2 = 12.1, s22 = 0.187')):
        r = line(lab, cfg, s)
        if r:
            r['technote dchi2'] = TECHNOTE_ICE[s]
            rows.append(r)
tn = pd.DataFrame(rows)
with pd.option_context('display.width', 250, 'display.float_format', lambda v: f'{v:.4g}'):
    print(tn.to_string(index=False))
tn.to_csv(TAB / 'technote_comparison_b.csv', index=False)
for s in ICE:
    HEAD[f'dchi2_{s}'] = {c: float(f2.delta_chi2('fd', c, s)) for c in CFG}
''')

md(r"""
## 3. Poisson fluctuations: the null distribution of $\Delta\chi^2$

Ten throws of the CV (seeds 1-10) for each analysis. Under Wilks' theorem $\Delta\chi^2\sim\chi^2_{2}$ (CDF dotted).
""")

code(r'''
fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
rows = []
for cfg in CFG:
    ps = [f2.pair('fd', cfg, s) for s in POIS]
    ps = [p for p in ps if p]
    if not ps:
        continue
    d = np.array([p[0]['chi2'] - p[1]['chi2'] for p in ps]); n0 = np.array([p[0]['chi2'] for p in ps])
    xs = np.sort(d)
    axes[0].step(xs, np.arange(1, len(xs) + 1) / len(xs), where='post', color=CC[cfg], label=f'{CFG[cfg]} (n = {len(xs)})')
    axes[1].hist(n0, bins=8, histtype='step', color=CC[cfg], label=CFG[cfg])
    ks = stats.kstest(d, stats.chi2(2).cdf)
    rows.append({'config': CFG[cfg], 'throws': len(d), 'mean dchi2': d.mean(), 'std dchi2': d.std(ddof=1), 'max dchi2': d.max(),
                 'fraction with Wilks p<0.05': float(np.mean([f2.wilks_p(x) < 0.05 for x in d])), 'mean chi2 null': n0.mean(),
                 'bins': f2.NBINS[cfg], 'KS p vs chi2(2)': ks.pvalue})
xx = np.linspace(0, 12, 200)
axes[0].plot(xx, stats.chi2(2).cdf(xx), 'k:', label=r'$\chi^2$, 2 dof')
axes[0].set_xlabel(r'$\Delta\chi^2=\chi^2_{\rm null}-\chi^2_{\rm osc}$'); axes[0].set_ylabel('cumulative fraction'); axes[0].legend(fontsize=9)
axes[1].set_xlabel(r'null-fit $\chi^2$'); axes[1].set_ylabel('throws'); axes[1].legend(fontsize=9)
fig.tight_layout(); fig.savefig(OUT / 'poisson_null.pdf', bbox_inches='tight')
pn = pd.DataFrame(rows)
print(pn.round(3).to_string(index=False)); pn.to_csv(TAB / 'poisson_null.csv', index=False)
HEAD['poisson'] = pn.to_dict(orient='records')
''')

md(r"""
## 4. Can the nuisance parameters absorb each dial? (data XML route, dial excluded from the fit)

$\chi^2_{\rm null}$ is how far the fake data sit from the best nuisance-only description (0 = fully absorbed);
$\Delta\chi^2$ is what the oscillation parameters add. Both for each analysis.
""")

code(r'''
def table(route, stems):
    rows = []
    for s in stems:
        r = {'dataset': s}
        for cfg in CFG:
            p = f2.pair(route, cfg, s)
            r[f'chi2 null {cfg}'] = p[0]['chi2'] if p else np.nan
            r[f'dchi2 {cfg}'] = (p[0]['chi2'] - p[1]['chi2']) if p else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


dt = table('ds', DS)
with pd.option_context('display.width', 250, 'display.max_rows', 100, 'display.float_format', lambda v: f'{v:.4g}'):
    print(dt.to_string(index=False))
dt.to_csv(TAB / 'dials_absorption.csv', index=False)
HEAD['dials_top_b'] = dt.sort_values('chi2 null b', ascending=False).head(8)[['dataset', 'chi2 null b', 'dchi2 b', 'chi2 null 1p', 'chi2 null Np']].to_dict(orient='records')
HEAD['n_dials_absorbed_b'] = int((dt[dt.dataset != 'ctrl']['chi2 null b'] < 1.0).sum()); HEAD['n_dials'] = int(len(DIALS))
''')

code(r'''
def heat(ax, M, title, rows_, cols_):
    im = ax.imshow(np.log10(np.clip(M, 1e-3, None)), aspect='auto', cmap='viridis', vmin=-2, vmax=2.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f'{M[i, j]:.2g}', ha='center', va='center', fontsize=6, color='w' if np.log10(max(M[i, j], 1e-3)) < 1 else 'k')
    ax.set_xticks(range(len(cols_))); ax.set_xticklabels(cols_, rotation=20, ha='right', fontsize=8); ax.set_title(title, fontsize=11)


order = [s for s in DS if s != 'ctrl']
fig, axes = plt.subplots(1, 2, figsize=(11, max(7, 0.28 * len(order))), sharey=True)
for ax, kind, title in zip(axes, ('chi2 null', 'dchi2'), (r'null-fit $\chi^2$', r'$\Delta\chi^2$ (null $-$ oscillation fit)')):
    M = np.array([[dt.set_index('dataset').loc[s, f'{kind} {c}'] for c in CFG] for s in order])
    heat(ax, M, title, order, [CFG[c] for c in CFG])
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([s.replace('_', ' ') for s in order], fontsize=7)
fig.tight_layout(); fig.savefig(OUT / 'dials_heatmap.pdf', bbox_inches='tight')
''')

md(r"""
## 5. What fitting one channel changes

For every dataset: how large is $\chi^2_{\rm null}$ (and $\Delta\chi^2$) in the single-channel analyses relative to the combined one,
and is the combined analysis more than the sum of its parts (which would signal that combining channels breaks degeneracies
between nuisance parameters)?
""")

code(r'''
allrows = []
for route, stems in (('fd', [s for s in PHYS + TECH + ['ctrl'] if s != 'ctrl']), ('ds', DIALS)):
    for s in stems:
        pb, p1, pN = (f2.pair(route, c, s) for c in ('b', '1p', 'Np'))
        if not (pb and p1 and pN):
            continue
        cb, c1, cN = pb[0]['chi2'], p1[0]['chi2'], pN[0]['chi2']
        db, d1, dN = (p[0]['chi2'] - p[1]['chi2'] for p in (pb, p1, pN))
        pm = f2.pair(route, 'm', s)
        allrows.append({'chi2 null m': pm[0]['chi2'] if pm else np.nan, 'dchi2 m': (pm[0]['chi2'] - pm[1]['chi2']) if pm else np.nan,'route': route, 'dataset': s, 'chi2 null b': cb, 'chi2 null 1p': c1, 'chi2 null Np': cN,
                        '(1p+Np)/b': (c1 + cN) / cb if cb > 1e-6 else np.nan, 'm/b': (pm[0]['chi2'] / cb) if (pm and cb > 1e-6) else np.nan, 'frac 1p of 1p+Np': c1 / (c1 + cN) if c1 + cN > 1e-9 else np.nan,
                        'dchi2 b': db, 'dchi2 1p': d1, 'dchi2 Np': dN})
ch = pd.DataFrame(allrows)
with pd.option_context('display.width', 250, 'display.max_rows', 120, 'display.float_format', lambda v: f'{v:.3g}'):
    print(ch.to_string(index=False))
ch.to_csv(TAB / 'channel_comparison.csv', index=False)
sig = ch[ch['chi2 null b'] > 1.0]
HEAD['channel_ratio_median'] = float(sig['(1p+Np)/b'].median()) if len(sig) else None
HEAD['merged_over_combined_median'] = float(sig['m/b'].median()) if len(sig) and sig['m/b'].notna().any() else None
HEAD['channel_1p_dominates'] = [r['dataset'] for _, r in sig.iterrows() if r['frac 1p of 1p+Np'] > 0.7]
HEAD['channel_Np_dominates'] = [r['dataset'] for _, r in sig.iterrows() if r['frac 1p of 1p+Np'] < 0.3]

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.scatter(ch['chi2 null 1p'] + 1e-4, ch['chi2 null Np'] + 1e-4, c='C2', s=18)
lim = [1e-4, max(100, ch[['chi2 null 1p', 'chi2 null Np']].max().max() * 2)]
ax.plot(lim, lim, 'k--', lw=0.8); ax.set_xscale('log'); ax.set_yscale('log')
for _, r in ch.sort_values('chi2 null b', ascending=False).head(10).iterrows():
    ax.annotate(r['dataset'], (r['chi2 null 1p'] + 1e-4, r['chi2 null Np'] + 1e-4), fontsize=7, xytext=(3, 3), textcoords='offset points')
ax.set_xlabel(r'null-fit $\chi^2$, single-proton only'); ax.set_ylabel(r'null-fit $\chi^2$, multi-proton only')
ax = axes[1]
ax.scatter(ch['chi2 null b'] + 1e-4, ch['chi2 null 1p'] + ch['chi2 null Np'] + 1e-4, c='C4', s=18)
ax.plot(lim, lim, 'k--', lw=0.8); ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'null-fit $\chi^2$, combined'); ax.set_ylabel(r'$\chi^2$ single-proton + $\chi^2$ multi-proton')
fig.tight_layout(); fig.savefig(OUT / 'channel_split.pdf', bbox_inches='tight')
''')

md(r"""
## 6. Constraints on the nuisance parameters

With the CV as fake data, the posterior width of each nuisance parameter measures the constraint (prior = 1). **Caveat found in
independent review:** the MCMC chains are short (20000 steps) and highly autocorrelated (lag-1 autocorrelation about 0.995), so a single
run's width is noisy. We therefore use the four runs that have the same fake data (the control, in both routes, null and oscillation fit) per
analysis, quote their mean, and only call a difference between analyses real if it exceeds three standard errors of the
repeatability. Posterior *means* are not used; pulls are the optimizer's best-fit values.
""")

code(r'''
RUNS4 = [('fd', 'free'), ('fd', 'fix'), ('ds', 'free'), ('ds', 'fix')]
WR = {}
for cfg in CFG:
    cols = []
    for route, var in RUNS4:
        if f2.available(route, cfg, var, 'ctrl'):
            f = f2.read_fit(route, cfg, var, 'ctrl')
            cols.append(pd.Series({n: f['sigma'].get(n, np.nan) for n in f2.nuisance_names(f)}))
    WR[cfg] = pd.concat(cols, axis=1)
M = pd.DataFrame({c: WR[c].mean(axis=1) for c in CFG})
SD = pd.DataFrame({c: WR[c].std(axis=1) for c in CFG})
rep = {c: float(SD[c].mean()) for c in CFG}
print('repeatability: mean over parameters of the std across the 4 control runs:', {c: round(v, 3) for c, v in rep.items()})
for c in CFG:
    print(f'  {CFG[c]:34s} median width {M[c].median():.3f}   below 0.5 sigma: {(M[c] < 0.5).sum():3d}   below 0.8 sigma: {(M[c] < 0.8).sum():3d}')
se = np.sqrt((SD['1p'] ** 2 + SD['Np'] ** 2) / 4)
d = M['1p'] - M['Np']
sig = d[(d.abs() > 0.15) & (d.abs() > 3 * se)].sort_values()
print('single-proton vs multi-proton width difference > 0.15 and > 3 s.e.:'); print(sig.round(3).to_string())
comb = (M['b'] < M[['1p', 'Np']].min(axis=1) - 0.05)
print('combining tightens by > 0.05 relative to the better channel:', int(comb.sum()), '(repeatability ~%.2f, so treat as indicative)' % rep['b'])
M.to_csv(TAB / 'constraint_widths.csv'); SD.to_csv(TAB / 'constraint_width_repeatability.csv')
HEAD['width_repeatability'] = rep
HEAD['width_median'] = {c: float(M[c].median()) for c in CFG}; HEAD['width_lt05'] = {c: int((M[c] < 0.5).sum()) for c in CFG}
HEAD['width_sig_1p_tighter'] = [n for n, v in sig.items() if v < 0]; HEAD['width_sig_Np_tighter'] = [n for n, v in sig.items() if v > 0]
fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
ax = axes[0]
ax.errorbar(M['1p'], M['Np'], xerr=SD['1p'], yerr=SD['Np'], fmt='o', ms=3, lw=0.5, color='C2')
ax.plot([0, 1.2], [0, 1.2], 'k--', lw=0.8)
for n in sig.index:
    ax.annotate(n, (M.loc[n, '1p'], M.loc[n, 'Np']), fontsize=8, xytext=(3, 3), textcoords='offset points')
ax.set_xlabel('mean posterior width, single-proton only'); ax.set_ylabel('mean posterior width, multi-proton only')
ax = axes[1]
order = M['b'].sort_values().index
for c in CFG:
    ax.errorbar(np.arange(len(order)), M.loc[order, c].values, yerr=SD.loc[order, c].values, fmt='.', ms=4, lw=0.4, color=CC[c], label=CFG[c])
ax.axhline(1, color='0.5', ls='--', lw=0.8); ax.set_xlabel('nuisance parameter, sorted by combined constraint'); ax.set_ylabel('posterior width [prior sigma]'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / 'constraint_widths.pdf', bbox_inches='tight')
''')

md(r"""
## 7. Recovering the injected oscillation point

For the native IceCube injections the oscillation fit should return the injected $(\Delta m^2,\sin^2 2\theta_{\mu\mu})$ and $\chi^2\approx0$.
""")

code(r'''
rows = []
for s, (dm2, s22) in ICE.items():
    for cfg in CFG:
        p = f2.pair('fd', cfg, s)
        if p:
            rows.append({'injection': s, 'config': CFG[cfg], 'dm2 injected': dm2, 'dm2 fit': p[1]['dm2'], 's22 injected': s22, 's22 fit': p[1]['s22'],
                         'chi2 osc': p[1]['chi2'], 'chi2 null': p[0]['chi2'], 'dchi2': p[0]['chi2'] - p[1]['chi2']})
rec = pd.DataFrame(rows)
print(rec.round(4).to_string(index=False)); rec.to_csv(TAB / 'icecube_recovery.csv', index=False)
HEAD['icecube_recovery_max_rel_err'] = float(max(abs(r['dm2 fit'] / r['dm2 injected'] - 1) for _, r in rec.iterrows())) if len(rec) else None
''')

md(r"""
## 8. Pulls

Post-fit nuisance values (prior sigma) with the MCMC posterior width as error bar, nuisance-only fits, overlaying the
combined, single-proton and multi-proton analyses; and, for each dataset, the largest pull in each analysis.
""")

code(r'''
def pull_overlay(route, stem, label_, variant='fix', save=True):
    avail = [c for c in CFG if f2.available(route, c, variant, stem)]
    if not avail:
        return None
    fits = {c: f2.read_fit(route, c, variant, stem) for c in avail}
    names = f2.nuisance_names(fits[avail[0]])
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(18, 5.2))
    ax.axhspan(-1, 1, color='0.9', zorder=0); ax.axhline(0, color='k', lw=0.8)
    for c, off in zip(avail, np.linspace(-0.22, 0.22, len(avail))):
        y = np.array([fits[c]['values'].get(n, np.nan) for n in names]); e = np.array([fits[c]['sigma'].get(n, np.nan) for n in names])
        ax.errorbar(x + off, y, yerr=e, fmt='o', ms=3, lw=0.9, color=CC[c], label=f"{CFG[c]} ($\\chi^2$ = {fits[c]['chi2']:.1f})")
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=90, fontsize=6.5); ax.set_xlim(-1, len(names))
    ax.set_ylabel(r'post-fit value [prior $\sigma$]'); ax.set_title(f"{label_} ({'oscillation fixed' if variant == 'fix' else 'oscillation floating'})", fontsize=12)
    ax.legend(loc='upper right', fontsize=9, ncol=3)
    if save:
        (OUT / 'pulls').mkdir(exist_ok=True)
        fig.savefig(OUT / 'pulls' / f'pulls_{variant}_{route}_{stem}.pdf', bbox_inches='tight')
    plt.close(fig)
    return fits


for route, stem in [('fd', s) for s in PHYS + TECH] + [('ds', s) for s in DIALS]:
    for v in ('fix', 'free'):
        pull_overlay(route, stem, stem.replace('_', ' '), v)
rows = []
for route, stem in [('fd', s) for s in PHYS + TECH] + [('ds', s) for s in DIALS]:
    r = {'route': route, 'dataset': stem}
    for c in CFG:
        p = f2.pair(route, c, stem)
        if p:
            nm = f2.nuisance_names(p[0])
            j = max(nm, key=lambda n: abs(p[0]['values'][n]))
            r[f'max|pull| {c}'] = abs(p[0]['values'][j]); r[f'param {c}'] = j
    rows.append(r)
pl = pd.DataFrame(rows)
pl.to_csv(TAB / 'max_pulls.csv', index=False)
with pd.option_context('display.width', 250, 'display.max_rows', 120, 'display.float_format', lambda v: f'{v:.3g}'):
    print(pl.sort_values('max|pull| b', ascending=False).head(25).to_string(index=False))
''')

md(r"""
## 9. Pull of the injected knob (technote-style injections, knob in the fit)

With the knob left in the fit the null fit should pull it to the injected value ($\pm1$ sigma).
""")

code(r'''
KNOB = {'qeint': 'QEIntf_dial0', 'specf': 'SF_q0bin0'}
rows = []
for s in TECH:
    if not s.endswith('infit'):
        continue
    nm = s.split('_')[0]; sign = 1 if '_p1_' in s else -1
    for cfg in CFG:
        p = f2.pair('fd', cfg, s)
        if p:
            for var, fit in (('null', p[0]), ('osc', p[1])):
                v = fit['values'].get(KNOB[nm], np.nan)
                rows.append({'injection': s, 'config': CFG[cfg], 'fit': var, 'injected': sign, 'fitted pull': v, 'posterior width': fit['sigma'].get(KNOB[nm], np.nan)})
kp = pd.DataFrame(rows)
print(kp.round(3).to_string(index=False)); kp.to_csv(TAB / 'knob_pulls.csv', index=False)
''')

md(r"""
## 10. Size of the injected effect, and what "absorbed" means

$\chi^2_{\rm null}$ alone is not normalised to how large the injection is (the stat-only distortions span about three orders of magnitude). Below, the
stat-only $\chi^2_{\rm stat}=\sum(D-C)^2/D$ of the fake data $D$ against the CV prediction $C$ in the projected spectra that PROfit stores (the 16+14 reco-energy
bins per detector, i.e. not the 2D bins that are fit, so shape effects in $p_n$ are under-counted), against the null-fit $\chi^2$. The ratio is only indicative.
""")

code(r'''
def stat_chi2(route, cfg, stem):
    tag = f2._paths(route, cfg, 'fix', stem)[0]
    p = f2.RUN / f'{tag}_PROglobal_histograms.root'
    if not p.exists():
        return np.nan
    import uproot
    h = uproot.open(p)['hist1d'].arrays(library='np')
    pre = np.asarray(h['prefix']).astype(str)
    tot = 0.0
    for d_ in range(3):
        for c_ in range(2):
            def get(pf):
                m = (pre == pf) & (h['detector'] == d_) & (h['channel'] == c_)
                _, o = np.unique(h['bin_index'][m], return_index=True)
                # PROfit was run with --scale-by-width: bin_content is events per unit RecoE; multiply back by the bin width to get counts
                w = (h['bin_high_edge'][m] - h['bin_low_edge'][m])[o]
                return h['bin_content'][m][o] * w
            D, C = get('DATA'), get('CV')
            if len(D):
                tot += float(np.sum((D - C) ** 2 / np.maximum(D, 1e-9)))
    return tot


rows = []
for route, stems in (('fd', PHYS + [s for s in TECH]), ('ds', DIALS)):
    for s in stems:
        S = stat_chi2(route, 'b', s)
        p = f2.pair(route, 'b', s)
        if p and np.isfinite(S):
            rows.append({'route': route, 'dataset': s, 'stat chi2 (projected)': S, 'chi2 null b': p[0]['chi2'],
                         'null/stat': p[0]['chi2'] / S if S > 0 else np.nan})
sz = pd.DataFrame(rows)
with pd.option_context('display.width', 220, 'display.max_rows', 100, 'display.float_format', lambda v: f'{v:.4g}'):
    print(sz.sort_values('stat chi2 (projected)', ascending=False).to_string(index=False))
sz.to_csv(TAB / 'injection_size.csv', index=False)
HEAD['stat_chi2_range'] = [float(sz['stat chi2 (projected)'].min()), float(sz['stat chi2 (projected)'].max())] if len(sz) else None
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(sz['stat chi2 (projected)'], sz['chi2 null b'] + 1e-4, c=['C3' if r == 'fd' else 'C0' for r in sz['route']], s=20)
lim = [1e-4, 1e5]; ax.plot(lim, lim, 'k--', lw=0.8); ax.set_xscale('log'); ax.set_yscale('log')
for _, r in sz.sort_values('chi2 null b', ascending=False).head(9).iterrows():
    ax.annotate(r['dataset'], (r['stat chi2 (projected)'], r['chi2 null b'] + 1e-4), fontsize=7, xytext=(3, 3), textcoords='offset points')
ax.set_xlabel(r'stat-only $\chi^2$ of the injection (projected spectra)'); ax.set_ylabel(r'null-fit $\chi^2$ (combined analysis)')
fig.tight_layout(); fig.savefig(OUT / 'injection_size.pdf', bbox_inches='tight')
''')

md(r"""
## 11. Nuisance parameters at their bounds

Several detector parameters have bounded (flat) priors on $[-1,1]$. If the null fit pushes them to the bound, the remaining $\chi^2$ is limited by the prior envelope
and not by the flexibility of the model; list them for every dataset (combined analysis, null fit).
""")

code(r'''
rows = []
for route, stems in (('fd', PHYS + TECH + POIS), ('ds', DS)):
    for s in stems:
        p = f2.pair(route, 'b', s)
        if p:
            nm = f2.nuisance_names(p[0])
            peg = [n for n in nm if abs(abs(p[0]['values'][n]) - 1) < 0.01]
            rows.append({'dataset': s, 'n at bound': len(peg), 'parameters': '; '.join(peg), 'chi2 null': p[0]['chi2']})
pg = pd.DataFrame(rows)
print(pg[pg['n at bound'] > 0].to_string(index=False)); pg.to_csv(TAB / 'bounded_parameters.csv', index=False)
HEAD['n_datasets_with_bounded'] = int((pg['n at bound'] > 0).sum()); HEAD['n_datasets'] = int(len(pg))
''')

md(r"""
## 12. QE interference with the whole dial family excluded

In the technote-style `excl` runs only `QEIntf_dial0` is excluded from the fit; dials 1--5 of the same family stay in, so they can mimic it. Here all six are excluded
(the injection is still `QEIntf_dial0` at $\pm1\sigma$).
""")

code(r'''
rows = []
for sg in ('p1', 'm1'):
    for cfg in CFG:
        a_, b_ = (f2.pair('fd', cfg, f'qeint_{sg}_{v}') for v in ('excl', 'exclall'))
        if a_ and b_:
            rows.append({'injection': f'qeint_{sg}', 'config': CFG[cfg], 'chi2 null, dial0 excluded': a_[0]['chi2'], 'chi2 null, family excluded': b_[0]['chi2'],
                         'dchi2 family excluded': b_[0]['chi2'] - b_[1]['chi2']})
fam = pd.DataFrame(rows)
print(fam.round(4).to_string(index=False)); fam.to_csv(TAB / 'qeint_family_excluded.csv', index=False)
''')

md(r"""
## 13. Where the oscillation fit lands for non-oscillation fake data

A non-oscillation mismodelling that the oscillation parameters partly absorb would appear as a preference for a sterile-like signal. For every non-oscillation dataset with
$\Delta\chi^2>1$ in the combined analysis: the fitted point and the nominal (Wilks, 2 dof, **not validated here**) Gaussian-equivalent two-sided significance.
""")

code(r'''
osc_like = ('phys_', 'osc_ic2024_pebl')
rows = []
for route, stems in (('fd', PHYS + TECH), ('ds', DIALS)):
    for s in stems:
        if s.startswith(osc_like):
            continue
        for cfg in CFG:
            p = f2.pair(route, cfg, s)
            if p and (p[0]['chi2'] - p[1]['chi2']) > 1.0:
                d_ = p[0]['chi2'] - p[1]['chi2']
                rows.append({'dataset': s, 'config': CFG[cfg], 'dchi2': d_, 'fit dm2 [eV2]': p[1]['dm2'], 'fit sin2 2th': p[1]['s22'],
                             'nominal sigma (Wilks)': stats.norm.isf(f2.wilks_p(d_) / 2)})
fo = pd.DataFrame(rows)
with pd.option_context('display.width', 220, 'display.max_rows', 100, 'display.float_format', lambda v: f'{v:.3g}'):
    print(fo.sort_values('dchi2', ascending=False).to_string(index=False))
fo.to_csv(TAB / 'false_oscillation_signal.csv', index=False)
HEAD['false_signal_b'] = fo[fo.config == CFG['b']].sort_values('dchi2', ascending=False).to_dict(orient='records')
''')

md(r"""
## 14. Minimiser robustness: refits with different PROfit seeds

The same data-XML fake datasets (combined analysis) refit with two extra RNG seeds (`-s 7`, `-s 13`); the default run is the reference. Note that repeating a run with the
default seed reproduced the earlier $\chi^2$ exactly, so a re-run alone is not a robustness test.
""")

code(r'''
def seed_chi2(stem, seed, var):
    log = f2.RUN / 'logs' / f'global_d2b_{stem}_s{seed}_{var}.out'
    if not log.exists():
        return np.nan
    v = np.nan
    for line in log.read_text(errors='ignore').splitlines():
        if 'Global Best Fit chi^2:' in line:
            v = float(line.split('chi^2:')[1].split()[0])
    return v


rows = []
for stem in ('mec_bdt', 'minerva_3dqelike_bnb', 'ub_cc1p0pi_dpt_dat_midW', 'osc_ic2024_dm2hi'):
    p = f2.pair('ds', 'b', stem)
    if not p:
        continue
    r = {'dataset': stem, 'null (default seed)': p[0]['chi2'], 'osc (default seed)': p[1]['chi2']}
    for seed in (7, 13):
        r[f'null s{seed}'] = seed_chi2(stem, seed, 'fix'); r[f'osc s{seed}'] = seed_chi2(stem, seed, 'free')
    vals_n = [r['null (default seed)'], r['null s7'], r['null s13']]; vals_o = [r['osc (default seed)'], r['osc s7'], r['osc s13']]
    r['max |dnull|'] = float(np.nanmax(np.abs(np.array(vals_n) - vals_n[0]))); r['max |dosc|'] = float(np.nanmax(np.abs(np.array(vals_o) - vals_o[0])))
    rows.append(r)
sr = pd.DataFrame(rows)
with pd.option_context('display.width', 250, 'display.float_format', lambda v: f'{v:.5g}'):
    print(sr.to_string(index=False))
sr.to_csv(TAB / 'seed_robustness.csv', index=False)
HEAD['seed_robustness_max_diff'] = float(np.nanmax(sr[['max |dnull|', 'max |dosc|']].values)) if len(sr) else None
''')

code(r'''
Path(TAB / 'headline.json').write_text(json.dumps(HEAD, indent=1, default=float))
print(json.dumps(HEAD, indent=1, default=float)[:3000])
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
