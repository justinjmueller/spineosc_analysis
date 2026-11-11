"""Shared helpers for the per-systematic study (ps_*.py).  cvwgt is not used.  Inputs: runs/syst_study/pairs/<det>__<var>.pkl (see build_pairs.py)."""
import sys, json, pathlib
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from ana import *          # load, binidx, prep, counts, boot_ratio, flows, BINS, NBIN, CTR, FIG, OUT, REPO, prod_hist
_load_orig = load
def load(name, v2=False):
    """as ana.load, but uses the de-duplicated table (ps_fixdup.py) where one exists"""
    f = OUT / "pairs_fixdup" / f"{name}.pkl"
    return pd.read_pickle(f) if (f.exists() and not v2) else _load_orig(name, v2)
PSFIG = REPO / "docs" / "syst_study" / "persyst"; PSFIG.mkdir(parents=True, exist_ok=True)
PS_JSON = PSFIG / "summary.json"

# the production histogram name for each pair
def prod_name(pair):
    det, var = pair.split("__"); run = {"sbnd_run1": "sbnd_Run1", "icarus_run2": "icarus_Run2", "icarus_run4": "icarus_Run4"}[det]
    return f"{var}_{run}"
PAIRS = ["sbnd_run1__SCE0x", "sbnd_run1__SCE2x", "sbnd_run1__WireModxThetaXW", "sbnd_run1__WireModYZ", "sbnd_run1__MattModYZ", "sbnd_run1__DENT",
         "icarus_run2__SCE0x", "icarus_run2__SCE2x", "icarus_run2__WireModxThetaXW", "icarus_run2__WireModYZ",
         "icarus_run4__SCE0x", "icarus_run4__SCE2x", "icarus_run4__WireModxThetaXW", "icarus_run4__WireModYZ"]

def col(df, side, name):
    c = f"{side}_{name}"; return df[c].to_numpy(float) if c in df else np.full(len(df), np.nan)

def truth_class(cc, pdg, cat):
    """0 signal-like numu CC (category 0/2), 1 other numu CC, 2 NC, 3 nue CC, 4 other/cosmic/unknown(-1 = no truth)"""
    out = np.full(len(cc), -1, int); ok = np.isfinite(cc) & np.isfinite(cat)
    nu_mu = ok & (cc == 1) & (np.abs(pdg) == 14); out[nu_mu] = 1; out[nu_mu & np.isin(cat, [0, 2])] = 0
    out[ok & (cc == 0)] = 2; out[ok & (cc == 1) & (np.abs(pdg) == 12)] = 3
    out[ok & (cc == 1) & ~np.isin(np.abs(pdg), [12, 14])] = 4; out[ok & (cat == 9)] = 4
    return out
CLASSES = ["signal numuCC(0,2)", "other numuCC", "NC", "nueCC", "cosmic/other"]

def drift_dist(det, x):
    """distance to the cathode (cm). SBND cathode x=0; ICARUS cathodes at |x|=210.2"""
    return np.abs(x) if det.startswith("sbnd") else np.abs(np.abs(x) - 210.2)

def species_margins(df, side):
    """argmax-with-2x-muon-boost margins.  muon candidate: 2*s_mu - max(s_pi,s_p);  proton candidate: s_p - max(2*s_mu, s_pi)."""
    smu, spi, sp = col(df, side, "reco_leading_muon_muon_softmax"), col(df, side, "reco_leading_muon_pion_softmax"), col(df, side, "reco_leading_muon_proton_softmax")
    pmu, ppi, pp = col(df, side, "reco_leading_proton_muon_softmax"), col(df, side, "reco_leading_proton_pion_softmax"), col(df, side, "reco_leading_proton_proton_softmax")
    return 2 * smu - np.maximum(spi, sp), pp - np.maximum(2 * pmu, ppi)

def feature_table(df, side):
    """dict of named per-event arrays on one side"""
    f = {}
    for n, c in [("E_vis", "reco_visible_energy"), ("mu_ke", "reco_leading_muon_ke"), ("mu_len", "reco_leading_muon_length"), ("mu_score", "reco_leading_muon_muon_softmax"),
                 ("mu_pion_score", "reco_leading_muon_pion_softmax"), ("mu_proton_score", "reco_leading_muon_proton_softmax"), ("mu_primary", "reco_leading_muon_primary_softmax"),
                 ("mu_calo_csda", "reco_leading_muon_calo_csda_diff"), ("mu_mcs_csda", "reco_leading_muon_mcs_csda_diff"),
                 ("p_ke", "reco_leading_proton_ke"), ("p_len", "reco_leading_proton_length"), ("p_score", "reco_leading_proton_proton_softmax"), ("p_mu_score", "reco_leading_proton_muon_softmax"),
                 ("p_pion_score", "reco_leading_proton_pion_softmax"), ("p_primary", "reco_leading_proton_primary_softmax"), ("p_calo_csda", "reco_leading_proton_calo_csda_diff"),
                 ("p_mip", "reco_leading_proton_mip_softmax"), ("opening", "reco_bivar_muon_proton_opening_angle"), ("Eqel", "reco_energy_qel"),
                 ("vx", "reco_vertex_x"), ("vy", "reco_vertex_y"), ("vz", "reco_vertex_z"), ("mu_dirx", "reco_leading_muon_start_dir_x"), ("mu_diry", "reco_leading_muon_start_dir_y"),
                 ("mu_dirz", "reco_leading_muon_start_dir_z"), ("mu_endx", "reco_leading_muon_end_x"), ("mu_endy", "reco_leading_muon_end_y"), ("mu_endz", "reco_leading_muon_end_z"),
                 ("p_endx", "reco_leading_proton_end_x"), ("p_endy", "reco_leading_proton_end_y"), ("p_endz", "reco_leading_proton_end_z"),
                 ("n_p", "reco_proton_multiplicity"), ("flash_pe", "event_largest_flash_pe")]:
        f[n] = col(df, side, c)
    f["mu_margin"], f["p_margin"] = species_margins(df, side)
    f["theta_xz"] = np.degrees(np.arctan2(np.abs(f["mu_dirx"]), np.abs(f["mu_dirz"])))   # angle of the muon in the drift-beam plane: 0 = along z (perpendicular to drift), 90 = along drift
    return f

def pm_err(k, n):
    """binomial error of k/n"""
    n = np.maximum(n, 1); p = k / n; return np.sqrt(np.maximum(p * (1 - p), 0) / n)

def cond_frac(x, num_mask, den_mask, edges):
    """fraction num/den in bins of x (both masks restricted to events with finite x), with binomial errors"""
    ok = np.isfinite(x); b = np.digitize(x, edges) - 1; nb = len(edges) - 1; b[~ok | (b < 0) | (b >= nb)] = -1
    den = np.bincount(b[den_mask & (b >= 0)], minlength=nb).astype(float); num = np.bincount(b[num_mask & den_mask & (b >= 0)], minlength=nb).astype(float)
    return num / np.maximum(den, 1), pm_err(num, den), den

def smd(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 5: return np.nan
    s = np.sqrt(0.5 * (a.var() + b.var())); return float((a.mean() - b.mean()) / s) if s > 0 else 0.0

def boot_channels(df, nb=200, seed=3):
    """paired Poisson bootstrap of combined, 1p and Mp ratios per bin (channel defined on each side by its own reco_proton_multiplicity)."""
    bn, bv = prep(df); mn = col(df, "n", "reco_proton_multiplicity"); mv = col(df, "v", "reco_proton_multiplicity")
    cn = np.where(np.isfinite(mn), (mn > 1).astype(int), -1); cv_ = np.where(np.isfinite(mv), (mv > 1).astype(int), -1)
    rng = np.random.default_rng(seed); out = {}
    def rat(w):
        res = []
        for sel in (None, 0, 1):
            mN = (bn >= 0) if sel is None else (bn >= 0) & (cn == sel); mV = (bv >= 0) if sel is None else (bv >= 0) & (cv_ == sel)
            n = np.bincount(bn[mN], weights=w[mN], minlength=NBIN); v = np.bincount(bv[mV], weights=w[mV], minlength=NBIN); res.append((n, v))
        return res
    r0 = rat(np.ones(len(bn))); R = np.empty((nb, 3, NBIN))
    for i in range(nb):
        w = rng.poisson(1.0, len(bn)).astype(float); rr = rat(w)
        for j in range(3): R[i, j] = np.where(rr[j][0] > 0, rr[j][1] / np.maximum(rr[j][0], 1e-9), np.nan)
    nom = np.array([r[0] for r in r0]); var = np.array([r[1] for r in r0]); ratio = np.where(nom > 0, var / np.maximum(nom, 1e-9), np.nan)
    return dict(ratio=ratio, err=np.nanstd(R, axis=0), boots=R, nom=nom, var=var)

def savefig(fig, name):
    fig.savefig(PSFIG / name, dpi=110, bbox_inches="tight"); import matplotlib.pyplot as plt; plt.close(fig)

TRUTH_COLS = ["true_cc", "true_pdg", "true_category", "true_neutrino_energy", "true_leading_muon_ke", "true_leading_proton_ke"]
def truth_map(det):
    """Event-keyed truth lookup collected from every pair table of this detector (nominal side and varied side where filled).  Needed because for some variations the varied tree has no truth."""
    f = OUT / f"ps_truthmap_{det}.pkl"
    if f.exists(): return pd.read_pickle(f)
    parts = []
    for p in sorted((OUT / "pairs").glob(f"{det}__*.pkl")):
        d = pd.read_pickle(p)
        for s in ("n", "v"):
            c = f"{s}_true_category"
            if c not in d: continue
            m = d[c].notna().to_numpy()
            if m.sum() == 0: continue
            t = d.loc[m, KEYS + ["Ei"]].copy()
            for tc in TRUTH_COLS: t[tc] = d.loc[m, f"{s}_{tc}"].to_numpy() if f"{s}_{tc}" in d else np.nan
            parts.append(t)
        del d
    t = pd.concat(parts).drop_duplicates(KEYS + ["Ei"]); t.to_pickle(f); return t
def attach_truth(df, det):
    """returns arrays (cc, pdg, cat, Enu, mu_ke, p_ke) aligned with df from the truth map (all events, selected on either side or not)"""
    tm = truth_map(det); m = df[KEYS + ["Ei"]].merge(tm, on=KEYS + ["Ei"], how="left")
    return {c: m[c].to_numpy(float) for c in TRUTH_COLS}

def boot_int_channels(df, nb=300, seed=5):
    """paired Poisson bootstrap of the integrated (all 17 bins) ratios: all, 1p, Mp, and Mp-1p"""
    bn, bv = prep(df); mn = col(df, "n", "reco_proton_multiplicity"); mv = col(df, "v", "reco_proton_multiplicity"); rng = np.random.default_rng(seed)
    an, av = bn >= 0, bv >= 0; cn = (mn > 1); cv = (mv > 1); R = np.empty((nb, 3))
    for i in range(nb):
        w = rng.poisson(1.0, len(bn)).astype(float)
        R[i] = [(w * av).sum() / (w * an).sum() - 1, (w * (av & ~cv)).sum() / (w * (an & ~cn)).sum() - 1, (w * (av & cv)).sum() / (w * (an & cn)).sum() - 1]
    return dict(e_all=float(R[:, 0].std()), e_1p=float(R[:, 1].std()), e_Mp=float(R[:, 2].std()), e_diff=float((R[:, 2] - R[:, 1]).std()))
