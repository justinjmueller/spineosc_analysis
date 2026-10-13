"""Readers and tables for the 2D fake-data campaign (runs/fakedata_2d, see docs/fakedata_review/STUDY_DESIGN.md).

Two routes:
  fd: PROfit --use-fake-data on the unmodified main XML      tag m2<cfg>,        out g<free|fix>_<stem>
  ds: separate data XML (-d), dial not in the main XML       tag d2<cfg>_<stem>, out free|fix
cfg: b (combined: single- and multi-proton channels fit together), 1p (single-proton, 2D RecoE;pn), Np (multi-proton, 1D RecoE),
     m (merged: one inclusive 2D RecoE;pn histogram, no proton-multiplicity split).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import uproot
from scipy import stats

RUN = Path('/Users/mueller/Projects/spineosc/runs/fakedata_2d')
CONFIGS = {'b': 'combined (1p+Np)', '1p': 'single-proton only', 'Np': 'multi-proton only', 'm': 'merged (one inclusive histogram)'}
CFG_COLOR = {'b': 'k', '1p': 'C0', 'Np': 'C3', 'm': 'C2'}
PHYSICS = ['dmsq', 'sinsq2thmm']      # fitted as log10 values

# number of bins in each analysis configuration, for reference in the text (3 detectors)
NBINS = {'b': 3 * (48 + 14), '1p': 3 * 48, 'Np': 3 * 14, 'm': 3 * 72}


def _paths(route, cfg, variant, stem):
    if route == 'fd':
        tag = f'm2{cfg}_g{variant}_{stem}'
        log = RUN / 'logs' / f'fd_{cfg}_g{variant}_{stem}.out'
        marker = RUN / 'done' / f'fd_{cfg}_{variant}_{stem}'
    else:
        tag = f'd2{cfg}_{stem}_{variant}'
        log = RUN / 'logs' / f'global_d2{cfg}_{stem}_{variant}.out'
        marker = RUN / 'done' / f'ds_{cfg}_{stem}'
    return tag, log, marker


def available(route, cfg, variant, stem):
    tag, log, marker = _paths(route, cfg, variant, stem)
    return marker.exists() and (RUN / f'{tag}_global_fit.txt').exists() and log.exists()


_cache = {}


def read_fit(route, cfg, variant, stem):
    """dict(values, sigma, chi2, dm2, s22, warnings) for one global fit; cached."""
    key = (route, cfg, variant, stem)
    if key in _cache:
        return _cache[key]
    tag, log, _ = _paths(route, cfg, variant, stem)
    values = {}
    for line in (RUN / f'{tag}_global_fit.txt').read_text().splitlines()[1:]:
        name, _, val = line.partition(' : ')
        if val.strip():
            values[name.strip()] = float(val)
    sigma = {}
    with uproot.open(RUN / f'{tag}_PROglobal.root') as f:
        for k in f.keys():
            if k.endswith('_posterior;1'):
                counts, edges = f[k].to_numpy()
                x = 0.5 * (edges[1:] + edges[:-1])
                mean = np.sum(counts * x) / counts.sum()
                sigma[k[:-len('_posterior;1')]] = float(np.sqrt(np.sum(counts * (x - mean) ** 2) / counts.sum()))
    chi2, warn = np.nan, 0
    for line in log.read_text(errors='ignore').splitlines():
        if 'Global Best Fit chi^2:' in line:
            chi2 = float(line.split('chi^2:')[1].split()[0])
        if line.startswith('1 ') or ' WARNING' in line:
            warn += 1
    out = {'values': values, 'sigma': sigma, 'chi2': chi2, 'warnings': warn,
           'dm2': 10 ** values.get('dmsq', np.nan), 's22': 10 ** values.get('sinsq2thmm', np.nan)}
    _cache[key] = out
    return out


def nuisance_names(fit):
    return [n for n in fit['values'] if n not in PHYSICS]


def pair(route, cfg, stem):
    """(null fit, oscillation fit) or None."""
    if available(route, cfg, 'fix', stem) and available(route, cfg, 'free', stem):
        return read_fit(route, cfg, 'fix', stem), read_fit(route, cfg, 'free', stem)
    return None


def delta_chi2(route, cfg, stem):
    p = pair(route, cfg, stem)
    return np.nan if p is None else p[0]['chi2'] - p[1]['chi2']


def job_status():
    """DataFrame of every completed job marker with its return codes."""
    rows = []
    for m in sorted((RUN / 'done').glob('*')):
        rows.append({'job': m.name, 'rc': m.read_text().strip()})
    return pd.DataFrame(rows)


def failures():
    st = job_status()
    if st.empty:
        return st
    bad = st[st.rc.map(lambda s: any(t.split('=')[-1] not in ('0',) for t in s.split()))]
    return bad


def wilks_p(dchi2, dof=2):
    return stats.chi2.sf(max(dchi2, 0.0), dof)
