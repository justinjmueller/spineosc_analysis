"""Readers for the global fits of the fake-data studies (runs/fakedata_full, produced by fitjob.sh).

File naming (config tag fd<cfg>, fit tag g<variant>_<stem>):
    fd<cfg>_g<variant>_<stem>_global_fit.txt     best-fit vector (physics log10 values, then nuisance parameters)
    fd<cfg>_g<variant>_<stem>_PROglobal.root     MCMC posterior of every nuisance parameter
    logs/global_<cfg>_g<variant>_<stem>.out      log (global best fit chi^2)
cfg:     b = both channels, 1p = single-proton only, Np = multi-proton only
variant: free = oscillation parameters floating, fix = nuisance parameters only (--syst-only)
"""

from pathlib import Path

import numpy as np
import uproot

RUN = Path('/Users/mueller/Projects/spineosc/runs/fakedata_full')
CONFIGS = {'b': 'both channels', '1p': 'single-proton only', 'Np': 'multi-proton only'}
PHYSICS = ['dmsq', 'sinsq2thmm']      # fitted as log10 values


def tag(cfg, variant, stem):
    return f'fd{cfg}_g{variant}_{stem}'


def available(cfg, variant, stem):
    return (RUN / 'done' / f'{cfg}_{variant}_{stem}').exists() and \
        (RUN / f'{tag(cfg, variant, stem)}_global_fit.txt').exists()


_cache = {}


def read_fit(cfg, variant, stem):
    """{'values', 'sigma', 'chi2', 'dm2', 's22'} of one global fit (cached)."""
    key = (cfg, variant, stem)
    if key in _cache:
        return _cache[key]
    t = tag(cfg, variant, stem)
    values = {}
    for line in (RUN / f'{t}_global_fit.txt').read_text().splitlines()[1:]:
        name, _, val = line.partition(' : ')
        if val.strip():
            values[name.strip()] = float(val)
    sigma = {}
    with uproot.open(RUN / f'{t}_PROglobal.root') as f:
        for k in f.keys():
            if k.endswith('_posterior;1'):
                counts, edges = f[k].to_numpy()
                x = 0.5 * (edges[1:] + edges[:-1])
                mean = np.sum(counts * x) / counts.sum()
                sigma[k[:-len('_posterior;1')]] = float(np.sqrt(np.sum(counts * (x - mean) ** 2) / counts.sum()))
    chi2 = np.nan
    log = RUN / 'logs' / f'global_{cfg}_g{variant}_{stem}.out'
    for line in log.read_text(errors='ignore').splitlines():
        if 'Global Best Fit chi^2:' in line:
            chi2 = float(line.split('chi^2:')[1].split()[0])
    out = {'values': values, 'sigma': sigma, 'chi2': chi2,
           'dm2': 10 ** values.get('dmsq', np.nan), 's22': 10 ** values.get('sinsq2thmm', np.nan)}
    _cache[key] = out
    return out


def nuisance_names(fit):
    return [n for n in fit['values'] if n not in PHYSICS]


def delta_chi2(cfg, stem):
    """chi2_null - chi2_osc for one injection and channel configuration (nan if a fit is missing)."""
    if not (available(cfg, 'free', stem) and available(cfg, 'fix', stem)):
        return np.nan
    return read_fit(cfg, 'fix', stem)['chi2'] - read_fit(cfg, 'free', stem)['chi2']
