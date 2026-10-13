import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

from mutools.plotting import ternary as T
from mutools.plotting import use_style

use_style('rootlike')

DETVAR = '/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root'
OUTDIR = '/Users/mueller/Projects/spineosc/figures/dent_crosser_ternary'
os.makedirs(OUTDIR, exist_ok=True)

ID_BRANCHES = ['Run', 'Subrun', 'Evt', 'event_max_neutrino_energy']
SCORE_BRANCHES = ['reco_leading_muon_muon_softmax', 'reco_leading_muon_pion_softmax',
                   'reco_leading_muon_proton_softmax']
BRANCHES = ID_BRANCHES + SCORE_BRANCHES + [
    'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz',
    'reco_vertex_x', 'reco_leading_muon_end_x',
    'true_cc', 'true_pdg',
]

f = uproot.open(DETVAR)


def keytuples(d):
    return list(zip(d['Run'].tolist(), d['Subrun'].tolist(), d['Evt'].tolist(),
                     d['event_max_neutrino_energy'].tolist()))


def base_mask(d):
    fid = (d['reco_fiducialize_cathode'] == 1) & (d['reco_veto_sbnd_highy_highz'] == 1)
    true_muon = (d['true_cc'] == 1) & (np.abs(d['true_pdg']) == 14)
    return fid & true_muon


def crosser_mask(d):
    return (d['reco_vertex_x'] * d['reco_leading_muon_end_x']) < 0


def load_full(key):
    """Full (unpaired) baseline+true-muon selected sample -- same as the scatter plots."""
    t = f[f'events/{key}/signal']
    d = t.arrays(BRANCHES, library='np')
    keep = base_mask(d)
    return {
        'muon': d['reco_leading_muon_muon_softmax'][keep],
        'pion': d['reco_leading_muon_pion_softmax'][keep],
        'proton': d['reco_leading_muon_proton_softmax'][keep],
        'crosser': crosser_mask(d)[keep],
    }


def load_event_keys(key):
    t = f[f'events/{key}/events']
    d = t.arrays(ID_BRANCHES, library='np')
    return set(keytuples(d))


print('loading full (unpaired) samples for the ratio background...')
cv_full = load_full('sbnd_detvar_newcv')
dent_full = load_full('sbnd_detvar_newdent')
print(f"  CV: {len(cv_full['muon'])}   DENT: {len(dent_full['muon'])}")

print('loading common-processed events for the paired displacement field...')
cv_ev_keys = load_event_keys('sbnd_detvar_newcv')
dent_ev_keys = load_event_keys('sbnd_detvar_newdent')
common_processed = cv_ev_keys & dent_ev_keys
print(f'  common processed: {len(common_processed)}')


def load_selected_restricted(key):
    t = f[f'events/{key}/signal']
    d = t.arrays(BRANCHES, library='np')
    keys = keytuples(d)
    mask = np.array([k in common_processed for k in keys])
    return {b: v[mask] for b, v in d.items()}, [k for k, m in zip(keys, mask) if m]


cv_c, cv_keys = load_selected_restricted('sbnd_detvar_newcv')
dent_c, dent_keys = load_selected_restricted('sbnd_detvar_newdent')

cv_pass = base_mask(cv_c)
dent_pass = base_mask(dent_c)
cv_crosser = crosser_mask(cv_c)

cv_pass_keyset = set(k for k, p in zip(cv_keys, cv_pass) if p)
dent_pass_keyset = set(k for k, p in zip(dent_keys, dent_pass) if p)
both_keys = cv_pass_keyset & dent_pass_keyset
print(f'  paired (selected + baseline on both sides): {len(both_keys)}')

cv_lookup = {k: i for i, k in enumerate(cv_keys)}
dent_lookup = {k: i for i, k in enumerate(dent_keys)}
both_list = list(both_keys)
cv_idx = np.array([cv_lookup[k] for k in both_list])
dent_idx = np.array([dent_lookup[k] for k in both_list])

paired_cv_crosser = cv_crosser[cv_idx]
paired = {
    'cv_muon': cv_c['reco_leading_muon_muon_softmax'][cv_idx],
    'cv_pion': cv_c['reco_leading_muon_pion_softmax'][cv_idx],
    'cv_proton': cv_c['reco_leading_muon_proton_softmax'][cv_idx],
    'dent_muon': dent_c['reco_leading_muon_muon_softmax'][dent_idx],
    'dent_pion': dent_c['reco_leading_muon_pion_softmax'][dent_idx],
    'dent_proton': dent_c['reco_leading_muon_proton_softmax'][dent_idx],
    'crosser': paired_cv_crosser,
}

GRID = T.triangle_grid(140)
gx, gy, inside = GRID
LABELS = (r'$\mu$', r'$\pi$', r'$p$')

POPULATIONS = [('stay', 'Non-crosser', False), ('xing', 'Cathode crosser', True)]

for pop_key, pop_label, want_crosser in POPULATIONS:
    cv_mask = cv_full['crosser'] if want_crosser else ~cv_full['crosser']
    dent_mask = dent_full['crosser'] if want_crosser else ~dent_full['crosser']

    xcv, ycv, _ = T.to_ternary(cv_full['muon'][cv_mask], cv_full['pion'][cv_mask], cv_full['proton'][cv_mask])
    xdent, ydent, _ = T.to_ternary(dent_full['muon'][dent_mask], dent_full['pion'][dent_mask], dent_full['proton'][dent_mask])

    reference = T.density(xcv, ycv, *GRID)
    comparison = T.density(xdent, ydent, *GRID)
    ratio = T.log_ratio(comparison, reference)

    pmask = paired['crosser'] == want_crosser
    start = (paired['cv_muon'][pmask], paired['cv_pion'][pmask], paired['cv_proton'][pmask])
    end = (paired['dent_muon'][pmask], paired['dent_pion'][pmask], paired['dent_proton'][pmask])
    tx, ty, dx, dy, weights = T.displacement_field(start, end)

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    im = ax.pcolormesh(gx, gy, ratio, cmap='RdBu_r', vmin=-1.5, vmax=1.5, shading='auto', zorder=1)

    # Density contours at 1%, 10%, 50% of the CV peak, so a loud ratio in
    # the sparse tails (e.g. the mu-p edge) isn't mistaken for a large
    # effect on a well-populated region -- same convention as
    # spineosc_pid/ternary_migration.py.
    ax.contour(gx, gy, reference / np.nanmax(reference), levels=[0.01, 0.1, 0.5],
               colors='0.30', linewidths=0.7, zorder=3)

    T.draw_frame(ax, labels=LABELS)

    if len(tx):
        shade = T.arrow_shades(weights)
        ax.quiver(tx, ty, dx, dy, angles='xy', scale_units='xy', scale=1,
                  color=shade, width=0.005, headwidth=4, zorder=4)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r'$\log_2$ (DENT density / CV density)')

    fig.text(0.5, 0.965, rf'Signal tree -- {pop_label}', ha='center', fontsize=14, color='#d67a11')
    fig.text(0.5, 0.015,
             f'CV n={cv_mask.sum()}, DENT n={dent_mask.sum()}, paired n={pmask.sum()}  |  '
             '(signal tree, not selected) arrows: mean per-muon CV->DENT displacement (paired), faded by cell population  |  '
             'contours: CV density at 1%, 10%, 50% of peak',
             ha='center', fontsize=8, color='0.35', wrap=True)

    out = f'{OUTDIR}/signal_ratio_displacement_{pop_key}.png'
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print('wrote', out)
