"""
Variant A: CV vs DENT using the genuine per-particle true-particle file
(v1.0.10 production, old DENT variation -- the only CV/DENT-paired
particle-level file on disk). No displacement arrows: this tree has no
per-particle index (only Run/Subrun/Evt, not unique per-particle), so
particles can't be safely paired across samples -- density/ratio only.
"""
import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

from mutools.plotting import ternary as T
from mutools.plotting import use_style

use_style('rootlike')

PARTICLES = '/Users/mueller/data/spineosc/v1.0.10/spineosc_numu_sbnd_release-v1.0.10_particles_nosyst.root'
OUTDIR = '/Users/mueller/Projects/spineosc/figures/dent_crosser_ternary'
os.makedirs(OUTDIR, exist_ok=True)

BRANCHES = [
    'true_particle_pid',
    'reco_particle_muon_softmax', 'reco_particle_pion_softmax', 'reco_particle_proton_softmax',
    'reco_particle_start_x', 'reco_particle_end_x',
]
MUON_PID = 2  # pidtools.io.SPECIES[2] == 'muon'

f = uproot.open(PARTICLES)


def load(key):
    t = f[f'events/{key}/sideband_particles_reco']
    d = t.arrays(BRANCHES, library='np')
    keep = d['true_particle_pid'] == MUON_PID
    crosser = (d['reco_particle_start_x'] * d['reco_particle_end_x']) < 0
    return {
        'muon': d['reco_particle_muon_softmax'][keep],
        'pion': d['reco_particle_pion_softmax'][keep],
        'proton': d['reco_particle_proton_softmax'][keep],
        'crosser': crosser[keep],
    }


print('loading true-particle CV/DENT samples (v1.0.10, old DENT variation)...')
cv = load('sbnd_detvar_cv')
dent = load('sbnd_detvar_dent')
print(f"  CV: {len(cv['muon'])} true muons ({cv['crosser'].sum()} crossers)")
print(f"  DENT: {len(dent['muon'])} true muons ({dent['crosser'].sum()} crossers)")

GRID = T.triangle_grid(140)
gx, gy, inside = GRID
LABELS = (r'$\mu$', r'$\pi$', r'$p$')
POPULATIONS = [('stay', 'Non-crosser', False), ('xing', 'Cathode crosser', True)]

for pop_key, pop_label, want_crosser in POPULATIONS:
    cv_mask = cv['crosser'] if want_crosser else ~cv['crosser']
    dent_mask = dent['crosser'] if want_crosser else ~dent['crosser']

    xcv, ycv, _ = T.to_ternary(cv['muon'][cv_mask], cv['pion'][cv_mask], cv['proton'][cv_mask])
    xdent, ydent, _ = T.to_ternary(dent['muon'][dent_mask], dent['pion'][dent_mask], dent['proton'][dent_mask])

    reference = T.density(xcv, ycv, *GRID)
    comparison = T.density(xdent, ydent, *GRID)
    ratio = T.log_ratio(comparison, reference)

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    im = ax.pcolormesh(gx, gy, ratio, cmap='RdBu_r', vmin=-1.5, vmax=1.5, shading='auto', zorder=1)
    ax.contour(gx, gy, reference / np.nanmax(reference), levels=[0.01, 0.1, 0.5],
               colors='0.30', linewidths=0.7, zorder=3)
    T.draw_frame(ax, labels=LABELS)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r'$\log_2$ (DENT density / CV density)')

    fig.text(0.5, 0.965, rf'True-particle file -- {pop_label}', ha='center', fontsize=14, color='#d67a11')
    fig.text(0.5, 0.015,
             f'v1.0.10 sideband_particles_reco, true_particle_pid==muon  |  '
             f'CV n={cv_mask.sum()}, DENT n={dent_mask.sum()} (unpaired -- no arrows, see script docstring)  |  '
             'contours: CV density at 1%, 10%, 50% of peak',
             ha='center', fontsize=8, color='0.35', wrap=True)

    out = f'{OUTDIR}/true_particles_ratio_{pop_key}.png'
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print('wrote', out)
