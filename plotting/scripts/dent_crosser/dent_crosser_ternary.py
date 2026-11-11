import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
from spineosc_paths import DATA_DIR, REPO  # noqa: E402
import uproot
import numpy as np
import matplotlib.pyplot as plt

from mutools.plotting import ternary as T
from mutools.plotting import use_style
from mutools.plotting.helpers import mark_axis

use_style('rootlike')

DETVAR = str(DATA_DIR / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root")
OUTDIR = str(REPO / "figures/dent_crosser_ternary")

import os
os.makedirs(OUTDIR, exist_ok=True)

BRANCHES = [
    'reco_leading_muon_muon_softmax', 'reco_leading_muon_pion_softmax', 'reco_leading_muon_proton_softmax',
    'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz',
    'reco_vertex_x', 'reco_leading_muon_end_x',
    'true_cc', 'true_pdg',
]

f = uproot.open(DETVAR)


def load(key):
    t = f[f'events/{key}/selected']
    d = t.arrays(BRANCHES, library='np')
    base = (d['reco_fiducialize_cathode'] == 1) & (d['reco_veto_sbnd_highy_highz'] == 1)
    # "true muon present" proxy: true numu CC (no per-object truth-match flag
    # exists for the leading-muon slot in this ntuple) -- see conversation.
    true_muon = (d['true_cc'] == 1) & (np.abs(d['true_pdg']) == 14)
    crosser = (d['reco_vertex_x'] * d['reco_leading_muon_end_x']) < 0
    keep = base & true_muon
    return {
        'muon': d['reco_leading_muon_muon_softmax'][keep],
        'pion': d['reco_leading_muon_pion_softmax'][keep],
        'proton': d['reco_leading_muon_proton_softmax'][keep],
        'crosser': crosser[keep],
    }


cv = load('sbnd_detvar_newcv')
dent = load('sbnd_detvar_newdent')

print(f"CV: {len(cv['muon'])} true-muon rows ({cv['crosser'].sum()} crossers, {(~cv['crosser']).sum()} non-crossers)")
print(f"DENT: {len(dent['muon'])} true-muon rows ({dent['crosser'].sum()} crossers, {(~dent['crosser']).sum()} non-crossers)")

GRID = T.triangle_grid(140)
LABELS = (r'$\mu$', r'$\pi$', r'$p$')

#: Cap on scattered points per panel -- the full samples (60k-190k rows)
#: render fine, but a plain scatter that large is slow and oversaturates;
#: a random subsample keeps the picture legible while the per-point colour
#: still comes from the full-sample density (see `density_at`'s docstring).
MAX_POINTS = 8000
RNG = np.random.default_rng(0)

PANELS = [
    ('cv', 'CV', cv),
    ('dent', 'DENT', dent),
]
POPULATIONS = [
    ('stay', 'Non-crosser', False),
    ('xing', 'Cathode crosser', True),
]

for sample_key, sample_label, data in PANELS:
    for pop_key, pop_label, want_crosser in POPULATIONS:
        mask = data['crosser'] if want_crosser else ~data['crosser']
        x, y, _ = T.to_ternary(data['muon'][mask], data['pion'][mask], data['proton'][mask])
        dens = T.density(x, y, *GRID)
        colors = T.density_at(x, y, dens, GRID)

        n_total = len(x)
        if n_total > MAX_POINTS:
            sel = RNG.choice(n_total, MAX_POINTS, replace=False)
        else:
            sel = slice(None)

        fig, ax = plt.subplots(figsize=(6.5, 4.8))
        sc = ax.scatter(x[sel], y[sel], c=colors[sel], cmap='viridis', s=4,
                         linewidths=0, rasterized=True)
        T.draw_frame(ax, labels=LABELS)
        fig.colorbar(sc, ax=ax, label='Density (normalized)', fraction=0.046, pad=0.04)
        mark_axis(ax, rf'$\bf{{{sample_label}}}$ -- {pop_label} (n={n_total}, {len(np.atleast_1d(x[sel]))} shown)', vadj=0.08)

        out = f'{OUTDIR}/{sample_key}_{pop_key}.png'
        fig.savefig(out, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print('wrote', out)
