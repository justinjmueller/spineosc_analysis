#!/usr/bin/env python3
"""Assemble docs/fakedata_review/report.tex (+ report.pdf) from the campaign tables and figures.

Tables and figures come from plotting/fakedata_2d_figures/ (written by build_fakedata_2d_notebook.py). The prose is kept in
docs/fakedata_review/prose/<section>.tex so it can be written after the numbers are known; a missing file leaves a visible
"[prose to be written]" marker. Reviewer reports, if present, are appended verbatim.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/build_fakedata_report.py
"""

import json
import re
import subprocess
from pathlib import Path

import pandas as pd

REPO = Path('/Users/mueller/Projects/spineosc')
FIG = REPO / 'plotting' / 'fakedata_2d_figures'
TAB = FIG / 'tables'
DOC = REPO / 'docs' / 'fakedata_review'
PROSE = DOC / 'prose'
PROSE.mkdir(exist_ok=True)


def esc(s):
    s = str(s)
    for a, b in (('\\', r'\textbackslash{}'), ('&', r'\&'), ('%', r'\%'), ('_', r'\_'), ('#', r'\#'), ('$', r'\$')):
        s = s.replace(a, b)
    return s


def fmt(v, digits=3):
    if isinstance(v, str):
        return esc(v)
    try:
        if pd.isna(v):
            return '--'
        a = abs(float(v))
    except Exception:
        return esc(v)
    if a != 0 and (a < 1e-3 or a >= 1e4):
        m, e = f'{float(v):.2e}'.split('e')
        return f'${m}\\times10^{{{int(e)}}}$'
    return f'{float(v):.{digits}g}'


def table(df, cols=None, rename=None, caption='', label='', digits=3, small=True, max_rows=None):
    if df is None or len(df) == 0:
        return '\\emph{[table not available]}\n'
    cols = cols or list(df.columns)
    d = df[cols].copy()
    if max_rows:
        d = d.head(max_rows)
    head = ' & '.join((rename[c] if rename and c in rename else esc(c)) for c in cols)   # renames are authored LaTeX
    body = '\n'.join(' & '.join(fmt(r[c], digits) for c in cols) + r' \\' for _, r in d.iterrows())
    # long free-text columns get a fixed-width paragraph column; the first column is left-aligned, the rest right-aligned
    spec = ''.join(('p{0.36\\linewidth}' if c in ('parameters', 'study') else ('l' if i == 0 else 'r')) for i, c in enumerate(cols))
    size = '\\footnotesize\n' if len(cols) > 5 else '\\small\n'
    return ('\\begin{table}[htbp]\n\\centering\n' + size + '\\adjustbox{max width=\\linewidth}{%\n' + f'\\begin{{tabular}}{{{spec}}}\n\\hline\n{head} \\\\\n\\hline\n{body}\n'
            f'\\hline\n\\end{{tabular}}}}\n\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{table}}\n')


def fig(path, caption, label, width='\\textwidth'):
    p = FIG / path
    if not p.exists():
        return f'\\emph{{[figure {esc(path)} not available]}}\n'
    return (f'\\begin{{figure}}[htbp]\n\\centering\n\\includegraphics[width={width}]{{{p}}}\n\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{figure}}\n')


def prose(name):
    p = PROSE / f'{name}.tex'
    if not p.exists():
        return '\\textbf{[prose to be written]}\n\n'
    t = p.read_text().replace('@FIG@', str(FIG)).replace('@GUMPLE_TABLE@', gumple_table())
    return t + '\n'


def gumple_table():
    """LaTeX longtable from docs/fakedata_review/gumple_vs_spine_tests.md (the markdown is the source of truth)."""
    f = DOC / 'gumple_vs_spine_tests.md'
    if not f.exists():
        return ''
    rows = [[c.strip() for c in ln.strip().strip('|').split('|')] for ln in f.read_text().splitlines() if ln.startswith('|') and not ln.startswith('|---')]
    head, body = rows[0], rows[1:]
    spec = '@{}r p{0.22\\linewidth} p{0.07\\linewidth} p{0.11\\linewidth} p{0.22\\linewidth} p{0.17\\linewidth}@{}'
    def cell(c):
        return esc(c).replace('<', '$<$').replace('>', '$>$').replace('Delta chi2', '$\\Delta\\chi^2$').replace('chi2', '$\\chi^2$')
    out = ['{\\scriptsize\\setlength{\\tabcolsep}{3pt}', f'\\begin{{longtable}}{{{spec}}}',
           '\\caption{Quantitative fake-data tests, GUMP/GUMPLE vs SPINE. Y = done, P = partial or planned, N = not done, ? = not shown in the sources read. Source: \\texttt{docs/fakedata\\_review/gumple\\_vs\\_spine\\_tests.md}.}\\label{tab:gumple_tests}\\\\',
           '\\toprule', ' & '.join(cell(h) for h in head) + ' \\\\ \\midrule\\endfirsthead', '\\toprule', ' & '.join(cell(h) for h in head) + ' \\\\ \\midrule\\endhead']
    out += [' & '.join(cell(c) for c in r) + ' \\\\' for r in body]
    out += ['\\bottomrule', '\\end{longtable}}']
    return '\n'.join(out)


SPECTRA_APPENDIX = [("cpi_subthreshold_minus100", "sub-threshold charged pion $-100\\%$"), ("cpi_subthreshold_plus100", "sub-threshold charged pion $+100\\%$"),
                    ("mec_bdt", "asymmetric MEC (BDT)"), ("jaesung_lowq2_pi_enhancement_prefsi", "low-$Q^2$ pion enhancement, pre-FSI"),
                    ("jaesung_lowq2_pi_enhancement_postfsi", "low-$Q^2$ pion enhancement, post-FSI"), ("pi_fsi_ha2025", "pion FSI hA2025"),
                    ("ub_cc1p0pi_dpt_dat_midW", "MicroBooNE CC1p0$\\pi$, mid-$W$"), ("ub_cc1p0pi_dpt_dat", "MicroBooNE CC1p0$\\pi$, all $W$"),
                    ("t2k_nc1pi_p_costh", "T2K NC1$\\pi$"), ("minerva_3dqelike_bnb", "MINERvA 3D QE-like"), ("qe_zexp_mva_to_lqcd", "LQCD")]


def spectra_appendix():
    """One figure per (dataset, channel): rows = SBND, ICARUS Run 2, Run 4; left = null fit, right = oscillation-floating fit (combined analysis).
    Pages of spectra/<dataset>_{null,osc}.pdf: 1-3 single-proton (SBND, Run 2, Run 4), 4-6 multi-proton."""
    out = ['\\clearpage\\newgeometry{margin=0.55in}\\section{Best-fit spectra of the fake datasets}\\label{app:spectra}\n',
           'Combined analysis, data-XML route, nominal exposure. Each figure shows one channel for the three detectors (rows: SBND, ICARUS Run 2, ICARUS Run 4); left: null fit (oscillation fixed), right: oscillation-floating fit. '
           'Stacked: the prefit prediction by true subchannel with its prefit uncertainty band; points: the injected fake data; red line: the best fit; lower panel: data / best fit with the best-fit uncertainty as the band. '
           'The right-hand watermark is the fake dataset name. The single-proton panels are RecoE projections of the 2D (RecoE $\\times p_n$) fit. The single-channel analyses are in the same files of \\path{plotting/fakedata_2d_figures/spectra/}.\n']
    for stem, label in SPECTRA_APPENDIX:
        n, o = FIG / 'spectra' / f'{stem}_null.pdf', FIG / 'spectra' / f'{stem}_osc.pdf'
        if not (n.exists() and o.exists()):
            continue
        for ch, chname, first in ((0, 'single-proton', 1), (1, 'multi-proton', 4)):
            rows = '\\\\\n'.join(f'\\includegraphics[page={first + d},width=0.495\\textwidth]{{{n}}}\\hfill\\includegraphics[page={first + d},width=0.495\\textwidth]{{{o}}}' for d in range(3))
            out.append('\\begin{figure}[p]\\centering\n' + rows + f'\n\\caption{{{label}, {chname} channel (dataset \\texttt{{{esc(stem)}}}): null fit (left) and oscillation-floating fit (right) for SBND, ICARUS Run 2 and ICARUS Run 4 (top to bottom).}}\n\\end{{figure}}\n')
    out.append('\\clearpage\\restoregeometry')
    return '\n'.join(out)


def pulls_appendix():
    out = ['\\section{Largest nuisance-parameter pulls}\\label{app:pulls}\n',
           'For each fake dataset in the combined analysis, the ten parameters with the largest post-fit value (in prior $\\sigma$) in the null fit: the CV closure fit (grey), the fake data with the oscillation fixed (filled) and floating (open); '
           'error bars are the MCMC posterior widths (noisy, see Sec.~\\ref{sec:limits}) and the shaded band is the prior $\\pm1\\sigma$. Values are the optimiser\'s best fit. Numbers in \\path{plotting/fakedata_2d_figures/tables/pull_top10.csv}; '
           'the full matrix is Fig.~\\ref{fig:pullmatrix}.\n']
    for i in (1, 2):
        f = FIG / f'pull_top10_{i}.pdf'
        if f.exists():
            out.append(f'\\begin{{figure}}[p]\\centering\n\\includegraphics[width=\\textwidth]{{{f}}}\n\\caption{{Top-ten pulls, page {i}.}}\n\\end{{figure}}\n')
    return '\n'.join(out)


def csv(name):
    p = TAB / name
    if not p.exists() or p.stat().st_size < 3:
        return None
    try:
        return pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return None


def verbatim_file(path):
    p = Path(path)
    if not p.exists():
        return '\\emph{[not available]}\n'
    import textwrap
    lines = []
    for ln in p.read_text().replace('\\end{verbatim}', '\\end {verbatim}').splitlines():
        lines += textwrap.wrap(ln, width=92, subsequent_indent='    ', break_long_words=True, break_on_hyphens=False) or ['']
    return '\\begin{footnotesize}\n\\begin{verbatim}\n' + '\n'.join(lines) + '\n\\end{verbatim}\n\\end{footnotesize}\n'


def main():
    head = json.loads((TAB / 'headline.json').read_text()) if (TAB / 'headline.json').exists() else {}
    tex = []
    tex.append(r'''\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}\usepackage{caption}
\usepackage{graphicx,amsmath,booktabs,float,hyperref,xcolor,adjustbox,enumitem,url,longtable,array}\setlist{itemsep=2pt,topsep=3pt}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}
\title{Fake-data fits in the 2D configuration:\\ absorption of injected effects by the nuisance parameters,\\ with and without oscillation, for combined and single-channel fits}
\author{Analysis report (draft with independent review)}
\date{\today}
\begin{document}
\maketitle
\tableofcontents
\clearpage''')
    tex.append('\\section{Key takeaways}\n' + prose('00_takeaways'))
    tex.append('\\section{Summary of results}\n' + prose('01_summary'))
    tex.append('\\section{Studies and methods}\n' + prose('02_methods'))
    tex.append('\\section{Validation of the set-up}\n' + prose('03_validation'))
    tex.append(table(csv('closure.csv'), caption='No-injection control (chi2 should be zero).', label='tab:closure'))
    xr = csv('route_crosscheck.csv')
    tex.append(table(xr, cols=['point', 'config', 'dchi2 native (fd)', 'dchi2 per-event-L dial (ds)', 'dchi2 fixed-110m dial (ds)'],
                     rename={'dchi2 native (fd)': r'$\Delta\chi^2$ native', 'dchi2 per-event-L dial (ds)': r'$\Delta\chi^2$ dial (event $L$)',
                             'dchi2 fixed-110m dial (ds)': r'$\Delta\chi^2$ dial (110 m)'},
                     caption='IceCube points: PROfit native injection vs the per-event-baseline dial and the fixed-110~m dial.', label='tab:route'))
    tex.append('\\section{Comparison with the technote}\n' + prose('04_technote'))
    tex.append(table(csv('technote_comparison_b.csv'), caption='Technote-style studies, combined analysis (both channels): this production vs the technote.', label='tab:technote'))
    tex.append(table(csv('icecube_recovery.csv'), cols=['injection', 'config', 'dm2 injected', 'dm2 fit', 's22 injected', 's22 fit', 'chi2 osc', 'dchi2'],
                     caption='Native IceCube injections: oscillation-fit result.', label='tab:recovery'))
    tex.append('\\section{Poisson fluctuations}\n' + prose('05_poisson'))
    tex.append(table(csv('poisson_null.csv'), caption='Poisson throws (seeds 1--10): distribution of $\\Delta\\chi^2$.', label='tab:poisson'))
    tex.append(fig('poisson_null.pdf', 'Left: cumulative distribution of $\\Delta\\chi^2$ for ten Poisson throws per analysis vs $\\chi^2_2$. Right: null-fit $\\chi^2$.', 'fig:poisson'))
    tex.append('\\section{Absorption of the injected effects (dials)}\\label{sec:dials}\n' + prose('06_dials'))
    dt = csv('dials_absorption.csv')
    if dt is not None:
        top = dt[dt.dataset != 'ctrl'].sort_values('chi2 null b', ascending=False).head(15)
        tex.append(table(top, caption='The 15 dials with the largest null-fit $\\chi^2$ in the combined analysis.', label='tab:dials'))
    tex.append(fig('dials_heatmap.pdf', 'Null-fit $\\chi^2$ and $\\Delta\\chi^2$ for every dial, in the combined, single-proton and multi-proton analyses.', 'fig:heat', '0.85\\textwidth'))
    tex.append(fig('injection_size.pdf', 'Null-fit $\\chi^2$ against the stat-only $\\chi^2$ of each injection (combined analysis).', 'fig:size', '0.6\\textwidth'))
    tex.append(table(csv('bounded_parameters.csv'), cols=['dataset', 'n at bound', 'parameters', 'chi2 null'], caption='Datasets in which the combined-analysis null fit puts nuisance parameters at their $\\pm1$ bound (rows with none omitted below).', label='tab:bounded', max_rows=None) if False else '')
    tex.append('\\section{Fitting a single channel}\n' + prose('07_channels'))
    tex.append(fig('channel_split.pdf', 'Left: null-fit $\\chi^2$ in the single-proton vs multi-proton analysis. Right: combined vs the sum of the two.', 'fig:split'))
    tex.append('\\section{Constraints on the nuisance parameters and pulls}\n' + prose('08_constraints'))
    tex.append(fig('constraint_widths.pdf', 'Nuisance-parameter posterior widths with the CV as fake data.', 'fig:widths'))
    tex.append(table(csv('knob_pulls.csv'), caption='Technote-style injections with the knob in the fit: the injected knob\'s fitted pull.', label='tab:knob'))
    bp = csv('bounded_parameters.csv')
    tex.append(table(bp[bp['n at bound'] > 0] if bp is not None else None, cols=['dataset', 'n at bound', 'parameters', 'chi2 null'], caption='Datasets in which the combined-analysis null fit puts nuisance parameters at their $\\pm1$ bound.', label='tab:bounded'))
    tex.append(table(csv('false_oscillation_signal.csv'), cols=['dataset', 'config', 'dchi2', 'fit dm2 [eV2]', 'fit sin2 2th', 'nominal sigma (Wilks)'], caption='Non-oscillation fake datasets for which the oscillation fit gains $\\Delta\\chi^2>1$: fitted point and nominal (Wilks, not validated) significance.', label='tab:falsesig'))
    tex.append('\\label{sec:seed}\n')
    tex.append(table(csv('seed_robustness.csv'), caption='Minimiser robustness: global-fit $\\chi^2$ of four dial datasets (combined analysis) with the default seed and two other PROfit RNG seeds.', label='tab:seed'))
    tex.append(table(csv('qeint_family_excluded.csv'), caption='QE-interference injection: null-fit $\\chi^2$ with only QEIntf\\_dial0 excluded from the fit and with the whole QEIntf family excluded.', label='tab:family'))
    tex.append('\\section{Comparison with the GUMP/GUMPLE studies and further tests}\\label{sec:gumple}\n' + prose('11_gumple'))
    tex.append('\\section{Limitations}\\label{sec:limits}\n' + prose('09_limitations'))
    tex.append('\\section{Independent review and response}\n' + prose('10_review_response'))
    tex.append('\\appendix\n\\section{Reproducibility}\n' + prose('A_reproducibility'))
    tex.append(spectra_appendix())
    tex.append(pulls_appendix())
    tex.append('\\section{Independent review, stage 1 (verbatim)}\\label{app:r1}\n' + verbatim_file(DOC / 'review_stage1.md'))
    tex.append('\\section{Independent review, stage 2 (verbatim)}\n' + verbatim_file(DOC / 'review_stage2.md'))
    tex.append('\\section{Independent review of the summary pages (verbatim)}\\label{app:rs}\n' + verbatim_file(DOC / 'review_summary.md'))
    tex.append('\\end{document}\n')
    out = DOC / 'report.tex'
    out.write_text('\n\n'.join(tex))
    for _ in range(2):
        r = subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'report.tex'], cwd=DOC, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        raise SystemExit('pdflatex failed')
    print('wrote', DOC / 'report.pdf')


if __name__ == '__main__':
    main()
