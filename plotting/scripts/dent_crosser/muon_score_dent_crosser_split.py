import uproot, numpy as np

DETVAR = '/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root'

id_branches = ['Run', 'Subrun', 'Evt', 'event_max_neutrino_energy']
branches = id_branches + ['reco_vertex_x', 'reco_leading_muon_end_x',
                           'reco_fiducialize_cathode', 'reco_veto_sbnd_highy_highz',
                           'reco_leading_muon_muon_softmax']

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
    return (d['reco_vertex_x'] * d['reco_leading_muon_end_x']) < 0


cv_base = base_mask(cv_c)
var_base = base_mask(var_c)
cv_crosser = crosser_mask(cv_c)
var_crosser = crosser_mask(var_c)

edges = np.arange(0.0, 1.0001, 0.05)


def report(name, cv_mask, var_mask):
    cv_vals = cv_c['reco_leading_muon_muon_softmax'][cv_mask]
    var_vals = var_c['reco_leading_muon_muon_softmax'][var_mask]
    cv_hist, _ = np.histogram(cv_vals, bins=edges)
    var_hist, _ = np.histogram(var_vals, bins=edges)
    print(f"\n=== {name} (CV n={cv_mask.sum()}, VAR n={var_mask.sum()}) ===")
    overall = var_hist.sum() / cv_hist.sum() if cv_hist.sum() else float('nan')
    print(f"  overall ratio = {overall:.4f}")
    for i in range(len(cv_hist)):
        r = var_hist[i] / cv_hist[i] if cv_hist[i] > 0 else float('nan')
        print(f"  [{edges[i]:>5.2f},{edges[i+1]:>5.2f})  cv_n={cv_hist[i]:>7d}  ratio={r:.4f}")


report("ALL baseline-selected", cv_base, var_base)
report("NON-crossers only", cv_base & (~cv_crosser), var_base & (~var_crosser))
report("CROSSERS only", cv_base & cv_crosser, var_base & var_crosser)
