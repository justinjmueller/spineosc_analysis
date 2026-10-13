import uproot, numpy as np

DETVAR = '/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root'
f = uproot.open(DETVAR)

id_branches = ['Run','Subrun','Evt','event_max_neutrino_energy']
rich_branches = id_branches + [
    'true_vertex_x','true_vertex_y','true_vertex_z','true_category','true_cc',
    'reco_vertex_x','reco_fiducialize_cathode','reco_veto_sbnd_highy_highz',
    'reco_leading_muon_end_x','reco_proton_multiplicity','event_largest_flash_pe',
]

def keytuples(d):
    return list(zip(d['Run'].tolist(), d['Subrun'].tolist(), d['Evt'].tolist(), d['event_max_neutrino_energy'].tolist()))

# 1) common *processed* events (from the 'events' tree, the superset)
def load_event_keys(key):
    t = f[f'events/{key}/events']
    d = t.arrays(id_branches, library='np')
    return set(keytuples(d))

cv_ev_keys = load_event_keys('sbnd_detvar_newcv')
var_ev_keys = load_event_keys('sbnd_detvar_newdent')
common_processed = cv_ev_keys & var_ev_keys
print('common processed events:', len(common_processed))

# 2) load the rich 'selected' branches for both sides
def load_selected(key):
    t = f[f'events/{key}/selected']
    return t.arrays(rich_branches, library='np')

cv_sel = load_selected('sbnd_detvar_newcv')
var_sel = load_selected('sbnd_detvar_newdent')

cv_keys = keytuples(cv_sel)
var_keys = keytuples(var_sel)

# 3) restrict each side's selected rows to those whose key is in common_processed
cv_common_mask = np.array([k in common_processed for k in cv_keys])
var_common_mask = np.array([k in common_processed for k in var_keys])
cv_c = {b: v[cv_common_mask] for b, v in cv_sel.items()}
var_c = {b: v[var_common_mask] for b, v in var_sel.items()}
print('cv common-restricted selected rows:', len(cv_c['Run']))
print('var common-restricted selected rows:', len(var_c['Run']))

# 4) apply the actual base selection (fiducialize_cathode + veto) used in the real analysis
def base_mask(d):
    return (d['reco_fiducialize_cathode'] == 1) & (d['reco_veto_sbnd_highy_highz'] == 1)

cv_pass = base_mask(cv_c)
var_pass = base_mask(var_c)
cv_pass_keys = set(k for k, p in zip(keytuples(cv_c), cv_pass) if p)
var_pass_keys = set(k for k, p in zip(keytuples(var_c), var_pass) if p)

n_cv_base = len(cv_pass_keys)
n_var_base = len(var_pass_keys)
print(f"baseline: CV={n_cv_base} VAR={n_var_base} ratio={n_var_base/n_cv_base:.4f}")

lost_keys = cv_pass_keys - var_pass_keys   # pass in CV, fail in variation
gained_keys = var_pass_keys - cv_pass_keys  # pass in variation, fail in CV
both_keys = cv_pass_keys & var_pass_keys
print(f"lost={len(lost_keys)}  gained={len(gained_keys)}  both={len(both_keys)}  net change={len(gained_keys)-len(lost_keys)}")
print(f"net change / baseline = {(len(gained_keys)-len(lost_keys))/n_cv_base:.4f}  (should match baseline ratio - 1 = {n_var_base/n_cv_base - 1:.4f})")

cv_key_arr = keytuples(cv_c)
var_key_arr = keytuples(var_c)

lost_mask = np.array([k in lost_keys for k in cv_key_arr])
gained_mask = np.array([k in gained_keys for k in var_key_arr])
cv_base_mask_all = np.array([k in cv_pass_keys for k in cv_key_arr])  # all CV events passing baseline (denominator)

lost = {b: v[lost_mask] for b, v in cv_c.items()}
gained = {b: v[gained_mask] for b, v in var_c.items()}
cv_baseline_full = {b: v[cv_base_mask_all] for b, v in cv_c.items()}

def stat(name, d):
    print(f"\n--- {name} (n={len(d['Run'])}) ---")
    ax = np.abs(d['true_vertex_x'])
    rx = np.abs(d['reco_vertex_x'])
    mex = np.abs(d['reco_leading_muon_end_x'])
    pe = d['event_largest_flash_pe']
    pe_ok = pe[np.isfinite(pe)]
    print(f"  |true_vertex_x|  mean={np.mean(ax):.2f} median={np.median(ax):.2f}  frac<10cm={np.mean(ax<10):.3f}  frac<20cm={np.mean(ax<20):.3f}")
    print(f"  |reco_vertex_x|  mean={np.mean(rx):.2f} median={np.median(rx):.2f}  frac<10cm={np.mean(rx<10):.3f}")
    print(f"  |muon_end_x|     mean={np.mean(mex):.2f} median={np.median(mex):.2f}  frac<10cm={np.mean(mex<10):.3f}")
    print(f"  flash_pe finite frac={len(pe_ok)/len(pe):.3f}  mean={np.mean(pe_ok) if len(pe_ok) else float('nan'):.1f}  median={np.median(pe_ok) if len(pe_ok) else float('nan'):.1f}")
    tc = d['true_category']
    tc = tc[np.isfinite(tc)]
    cats, counts = np.unique(tc, return_counts=True)
    order = np.argsort(-counts)
    total = counts.sum()
    print("  true_category frac:", {int(cats[i]): f"{counts[i]/total:.3f}" for i in order})
    print(f"  true_cc frac: {np.mean(d['true_cc']==1):.3f}")

stat("LOST under DENT (pass base selection in CV, fail in variation)", lost)
stat("GAINED under DENT (pass in variation, fail in CV)", gained)
stat("Full CV baseline (denominator)", cv_baseline_full)

print("\n=== Per-category net effect (ratio of var-pass / cv-pass counts, among common events) ===")
cats_all = np.unique(cv_c['true_category'][np.isfinite(cv_c['true_category'])])
for cat in cats_all:
    cv_n = np.sum((cv_c['true_category']==cat) & cv_pass)
    var_n = np.sum((var_c['true_category']==cat) & var_pass)
    if cv_n == 0: continue
    frac_of_total = cv_n/n_cv_base
    print(f"  category {int(cat):>2d}: CV={cv_n:7d} ({100*frac_of_total:4.1f}% of baseline)  VAR={var_n:7d}  ratio={var_n/cv_n:.4f}")
