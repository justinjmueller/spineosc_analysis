"""Q2 library: statistics of the detector-variation ratios r_b = V_b/N_b (V: varied, N: nominal selected counts per bin of reco visible energy).
Events are the independent units. Per-event contribution to r_b: d_ib = 1[bv_i=b] - r_b 1[bn_i=b] (delta method), Var(r)_bb' = sum_i d_ib d_ib' /(N_b N_b').
Bootstrap weights are Poisson(1) generated from a hash of the event key and the replicate index, so the SAME event gets the same weight in every pair table (joint resampling
across variations that share events/nominal). cvwgt is not used."""
import sys, json, numpy as np, scipy.sparse as sp
sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from common import BINS, OUT, FIG
from q2_cache import SCR, keyhash
NB = len(BINS) - 1; CTR = 0.5 * (BINS[1:] + BINS[:-1]); LOGX = (np.log(CTR) - np.log(CTR).mean()); LOGX = LOGX / np.abs(LOGX).max()   # in [-1,1]-ish
Q2 = FIG / "q2"; Q2.mkdir(parents=True, exist_ok=True); RES = SCR / "res"; RES.mkdir(exist_ok=True)
SBND = ["ChargeScale", "ChargeScale15", "SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ", "MattModYZ", "GainHigh", "GainLow", "DENT"]
ICAR = ["ChargeScale", "SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ"]
NAMES = [f"sbnd_run1__{v}" for v in SBND] + [f"icarus_run2__{v}" for v in ICAR] + [f"icarus_run4__{v}" for v in ICAR]
EXTRA = ["icarus_run4__WireModxThetaXWclean"]; ALL = NAMES + EXTRA   # EXTRA: duplicate-event runs removed (see q2_clean.py)
def det_of(n): return n.split("__")[0]
def var_of(n): return n.split("__")[1]
def load(name):
    d = dict(np.load(SCR / f"{name}.npz")); d["bn"] = d["bn"].astype(np.int64); d["bv"] = d["bv"].astype(np.int64); return d
def counts(d, w=None, mask=None):
    bn, bv = d["bn"], d["bv"]
    if mask is not None: bn, bv = bn[mask], bv[mask]
    if w is not None:
        w = w[mask] if mask is not None else w; wn = w[bn >= 0]; wv = w[bv >= 0]
    else: wn = wv = None
    return np.bincount(bn[bn >= 0], weights=wn, minlength=NB).astype(float), np.bincount(bv[bv >= 0], weights=wv, minlength=NB).astype(float)
def ratio(d, w=None, mask=None):
    n, v = counts(d, w, mask); return v / n
M64 = 2 ** 64 - 1
def split64(x):
    x = x.copy(); x += np.uint64(0x9E3779B97F4A7C15); x ^= x >> np.uint64(30); x *= np.uint64(0xBF58476D1CE4E5B9); x ^= x >> np.uint64(27); x *= np.uint64(0x94D049BB133111EB); x ^= x >> np.uint64(31); return x
import math
_CDF = np.cumsum([math.exp(-1) / math.factorial(k) for k in range(12)])
def pois_w(key, k, salt=0):
    """Poisson(1) weights, deterministic in (event key, replicate k, salt)."""
    u = split64(key ^ split64(np.full(1, ((k + 1) * 0x2545F4914F6CDD1D + salt * 0x9E3779B1) & M64, dtype=np.uint64))[0]) >> np.uint64(11)
    return np.searchsorted(_CDF, u.astype(np.float64) / 2.0 ** 53).astype(np.float64)
def uni(key, salt):
    return (split64(key ^ split64(np.full(1, (salt * 0x2545F4914F6CDD1D + 12345) & M64, dtype=np.uint64))[0]) >> np.uint64(11)).astype(np.float64) / 2.0 ** 53
def boot_counts(d, nb=1000, salt=0, mask=None):
    N = np.empty((nb, NB)); V = np.empty((nb, NB))
    for k in range(nb): N[k], V[k] = counts(d, pois_w(d["key"], k, salt), mask)
    return N, V
def design(d, r=None, w=None):
    """sparse D (rows x NB) with d_ib; also N. w: optional event weights (e.g. a half-split indicator)."""
    bn, bv = d["bn"], d["bv"]; n = len(bn); w = np.ones(n) if w is None else w
    N, V = counts(d, w); r = V / N if r is None else r; rows = np.arange(n)
    mv = bv >= 0; mn = bn >= 0
    Dv = sp.csr_matrix((w[mv], (rows[mv], bv[mv])), shape=(n, NB)); Dn = sp.csr_matrix((w[mn] * r[bn[mn]], (rows[mn], bn[mn])), shape=(n, NB))
    return Dv - Dn, N, V, r
def cov_analytic(d, w=None):
    D, N, V, r = design(d, w=w); C = (D.T @ D).toarray() / np.outer(N, N); return r, C, N, V
def cross_cov(da, db):
    """analytic covariance between ratio vectors of two pair tables (exact joint treatment through the shared event keys)."""
    Da, Na, Va, ra = design(da); Db, Nb, Vb, rb = design(db)
    common, ia, ib = np.intersect1d(da["key"], db["key"], return_indices=True)
    return (Da[ia].T @ Db[ib]).toarray() / np.outer(Na, Nb)
def naive_sigma(N, V, r): return r * np.sqrt(1.0 / N + 1.0 / V)
