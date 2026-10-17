import pathlib as _pl, sys as _sys; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
from spineosc_paths import DATA_DIR  # noqa: E402
import uproot, numpy as np

DETVAR = str(DATA_DIR / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root")
WSYST = str(DATA_DIR / "spineosc_numu_sbnd_v1.1.0_wsyst.root")

id_branches = ['Run', 'Subrun', 'Evt', 'event_max_neutrino_energy']
branches = id_branches + ['reco_vertex_x', 'reco_leading_muon_end_x',
                           'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz']

f = uproot.open(DETVAR)


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
    # True if the muon track crosses the cathode (vertex and endpoint on opposite sides of x=0)
    return (d['reco_vertex_x'] * d['reco_leading_muon_end_x']) < 0


cv_base = base_mask(cv_c)
var_base = base_mask(var_c)
n_cv_base = cv_base.sum()
n_var_base = var_base.sum()
print(f"Baseline: CV={n_cv_base} VAR={n_var_base} ratio={n_var_base / n_cv_base:.4f}")

cv_crosser = crosser_mask(cv_c)
var_crosser = crosser_mask(var_c)
print(f"Cathode-crosser fraction among CV baseline-passing events: {cv_crosser[cv_base].mean():.4f}")
print(f"Cathode-crosser fraction among VAR baseline-passing events: {var_crosser[var_base].mean():.4f}")

# with the cathode-crosser cut stacked on top of baseline (keep only non-crossers)
cv_m = cv_base & (~cv_crosser)
var_m = var_base & (~var_crosser)
n_cv, n_var = cv_m.sum(), var_m.sum()
print(f"\nWith cathode-crosser cut: CV={n_cv} VAR={n_var} ratio={n_var / n_cv:.4f}  (baseline ratio {n_var_base / n_cv_base:.4f})")
print(f"DENT-sample stats kept: {100 * n_cv / n_cv_base:.2f}%")

# Now: is the crosser fraction elevated among the lost/gained populations specifically?
cv_key_arr = keytuples(cv_c)
var_key_arr = keytuples(var_c)
cv_pass_keyset = set(k for k, p in zip(cv_key_arr, cv_base) if p)
var_pass_keyset = set(k for k, p in zip(var_key_arr, var_base) if p)

lost_keys = cv_pass_keyset - var_pass_keyset
gained_keys = var_pass_keyset - cv_pass_keyset
both_keys = cv_pass_keyset & var_pass_keyset

lost_mask = np.array([k in lost_keys for k in cv_key_arr]) & cv_base
gained_mask = np.array([k in gained_keys for k in var_key_arr]) & var_base
both_cv_mask = np.array([k in both_keys for k in cv_key_arr]) & cv_base

print(f"\nCathode-crosser fraction among LOST (pass CV, fail VAR): {cv_crosser[lost_mask].mean():.4f}  (n={lost_mask.sum()})")
print(f"Cathode-crosser fraction among GAINED (pass VAR, fail CV): {var_crosser[gained_mask].mean():.4f}  (n={gained_mask.sum()})")
print(f"Cathode-crosser fraction among BOTH (pass both, CV side): {cv_crosser[both_cv_mask].mean():.4f}  (n={both_cv_mask.sum()})")

# Stats cost on the full-statistics CV production sample
fw = uproot.open(WSYST)
tw = fw['events/sbnd_large/selected']
wbr = tw.arrays(['reco_vertex_x', 'reco_leading_muon_end_x', 'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz'], library='np')
wbase = (wbr['reco_fiducialize_cathode'] == 1) & (wbr['reco_veto_sbnd_highy_highz'] == 1)
wcrosser = (wbr['reco_vertex_x'] * wbr['reco_leading_muon_end_x']) < 0
n_wbase = wbase.sum()
n_wkept = (wbase & (~wcrosser)).sum()
print(f"\nFull-stats CV (sbnd_large) baseline: {n_wbase}")
print(f"Cathode-crosser fraction (full-stats CV): {wcrosser[wbase].mean():.4f}")
print(f"After removing cathode crossers: kept={n_wkept}  reduction={100 * (1 - n_wkept / n_wbase):.2f}%")
