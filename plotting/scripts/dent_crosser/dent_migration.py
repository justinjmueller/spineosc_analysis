import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
from spineosc_paths import DATA_DIR  # noqa: E402
import uproot, numpy as np

DETVAR = str(DATA_DIR / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root")
f = uproot.open(DETVAR)

id_branches = ['Run', 'Subrun', 'Evt', 'event_max_neutrino_energy']
branches = id_branches + ['reco_vertex_x', 'reco_leading_muon_end_x',
                           'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz',
                           'reco_leading_muon_muon_softmax']


def keytuples(d):
    return list(zip(d['Run'].tolist(), d['Subrun'].tolist(), d['Evt'].tolist(), d['event_max_neutrino_energy'].tolist()))


def load_event_keys(key):
    t = f[f'events/{key}/events']
    d = t.arrays(id_branches, library='np')
    return set(keytuples(d))


cv_ev_keys = load_event_keys('sbnd_detvar_newcv')
var_ev_keys = load_event_keys('sbnd_detvar_newdent')
common_processed = cv_ev_keys & var_ev_keys


def load_selected(key):
    t = f[f'events/{key}/selected']
    return t.arrays(branches, library='np')


cv_sel = load_selected('sbnd_detvar_newcv')
var_sel = load_selected('sbnd_detvar_newdent')

cv_keys = keytuples(cv_sel)
var_keys = keytuples(var_sel)
cv_common_mask = np.array([k in common_processed for k in cv_keys])
var_common_mask = np.array([k in common_processed for k in var_keys])
cv_c = {b: v[cv_common_mask] for b, v in cv_sel.items()}
var_c = {b: v[var_common_mask] for b, v in var_sel.items()}


def base_mask(d):
    return (d['reco_fiducialize_cathode'] == 1) & (d['reco_veto_sbnd_highy_highz'] == 1)


def crosser_mask(d):
    return (d['reco_vertex_x'] * d['reco_leading_muon_end_x']) < 0


cv_base = base_mask(cv_c)
var_base = base_mask(var_c)
cv_crosser = crosser_mask(cv_c)
var_crosser = crosser_mask(var_c)

cv_key_arr = keytuples(cv_c)
var_key_arr = keytuples(var_c)

cv_pass_keyset = set(k for k, p in zip(cv_key_arr, cv_base) if p)
var_pass_keyset = set(k for k, p in zip(var_key_arr, var_base) if p)
both_keys = cv_pass_keyset & var_pass_keyset
print(f'both (selected + baseline on both sides): {len(both_keys)}')

# per-event score lookup, restricted to "both"
cv_score_by_key = {k: s for k, s, p in zip(cv_key_arr, cv_c['reco_leading_muon_muon_softmax'], cv_base) if p and k in both_keys}
var_score_by_key = {k: s for k, s, p in zip(var_key_arr, var_c['reco_leading_muon_muon_softmax'], var_base) if p and k in both_keys}
cv_crosser_by_key = {k: c for k, c, p in zip(cv_key_arr, cv_crosser, cv_base) if p and k in both_keys}

common_keys = list(both_keys)
cv_scores = np.array([cv_score_by_key[k] for k in common_keys])
var_scores = np.array([var_score_by_key[k] for k in common_keys])
is_crosser = np.array([cv_crosser_by_key[k] for k in common_keys])
displacement = var_scores - cv_scores

print(f'\n=== ALL "both" events (n={len(common_keys)}) ===')
print(f'mean displacement = {displacement.mean():.5f}  std = {displacement.std():.5f}')
print(f'frac with |displacement| > 0.01 : {np.mean(np.abs(displacement) > 0.01):.4f}')
print(f'frac with |displacement| > 0.05 : {np.mean(np.abs(displacement) > 0.05):.4f}')
print(f'frac with |displacement| > 0.10 : {np.mean(np.abs(displacement) > 0.10):.4f}')

for label, mask in [('CROSSERS', is_crosser), ('NON-CROSSERS', ~is_crosser)]:
    d = displacement[mask]
    cvs = cv_scores[mask]
    print(f'\n=== {label} (n={mask.sum()}) ===')
    print(f'  mean displacement = {d.mean():.5f}  std = {d.std():.5f}')
    print(f'  frac |d|>0.01: {np.mean(np.abs(d) > 0.01):.4f}   frac |d|>0.05: {np.mean(np.abs(d) > 0.05):.4f}   frac |d|>0.10: {np.mean(np.abs(d) > 0.10):.4f}')
    # bucket by starting (CV) score to see if low-score events show correlated migration
    edges = [0.0, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 1.0]
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        m = (cvs >= lo) & (cvs < hi)
        if m.sum() == 0:
            continue
        dd = d[m]
        print(f'    CV score in [{lo:.2f},{hi:.2f}) n={m.sum():7d}  mean_disp={dd.mean():+.5f}  std={dd.std():.5f}  frac|d|>0.05={np.mean(np.abs(dd)>0.05):.3f}')

np.save('/tmp/dent_migration_arrays.npy', {'cv_scores': cv_scores, 'var_scores': var_scores, 'is_crosser': is_crosser}, allow_pickle=True)
