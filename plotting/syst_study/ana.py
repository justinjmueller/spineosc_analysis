"""Analysis helpers for the detector-systematics study."""
import json, sys, numpy as np, pandas as pd
from common import *
NBIN = len(BINS) - 1; CTR = 0.5 * (BINS[1:] + BINS[:-1])
def load(name, v2=False):
    return pd.read_pickle((OUT / ("pairs2" if v2 else "pairs")) / f"{name}.pkl")
def binidx(E):
    b = np.digitize(np.asarray(E, float), BINS) - 1; b[(b < 0) | (b >= NBIN) | ~np.isfinite(E)] = -1; return b
def prep(df):
    """bin indices of the nominal and varied candidates (-1 = not selected / out of range)"""
    bn = np.where(df["in_n"].to_numpy(), binidx(df["n_reco_visible_energy"].to_numpy()), -1); bv = np.where(df["in_v"].to_numpy(), binidx(df["v_reco_visible_energy"].to_numpy()), -1)
    return bn, bv
def counts(bn, bv, w=None):
    w = np.ones(len(bn)) if w is None else w
    return np.bincount(bn[bn >= 0], weights=w[bn >= 0], minlength=NBIN), np.bincount(bv[bv >= 0], weights=w[bv >= 0], minlength=NBIN)
def boot_ratio(df, nb=400, seed=1):
    bn, bv = prep(df); rng = np.random.default_rng(seed); n0, v0 = counts(bn, bv); r0 = v0 / n0; R = np.empty((nb, NBIN))
    for i in range(nb):
        w = rng.poisson(1.0, len(bn)).astype(float); n, v = counts(bn, bv, w); R[i] = v / n
    return r0, R.std(0), R, n0, v0
def flows(df):
    """per bin (as fraction of nominal): lost, gained, migration (both selected, moved bin), by nominal/varied bin assignment"""
    bn, bv = prep(df); both = (bn >= -1) & df["in_n"].to_numpy() & df["in_v"].to_numpy()
    n0 = np.bincount(bn[bn >= 0], minlength=NBIN).astype(float)
    lost = np.bincount(bn[(bn >= 0) & ~df["in_v"].to_numpy()], minlength=NBIN) / n0
    gained = np.bincount(bv[(bv >= 0) & ~df["in_n"].to_numpy()], minlength=NBIN) / n0
    inb_n = np.bincount(bn[both & (bn >= 0)], minlength=NBIN); inb_v = np.bincount(bv[both & (bv >= 0)], minlength=NBIN)
    mig = (inb_v - inb_n) / n0
    return dict(lost=-lost, gained=gained, migration=mig, n_nom=n0)
def pairs_available(v2=False):
    d = OUT / ("pairs2" if v2 else "pairs"); return sorted(p.stem for p in d.glob("*.pkl"))
def prod_hist(hist):
    import uproot; return uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root")[hist].values()
