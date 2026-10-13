"""Q2: smooth models for the binned ratio (GLS Legendre polynomials in ln E; Whittaker/P-spline-like penalised smoother with 2nd-divided-difference penalty)."""
import numpy as np
from scipy import stats
from numpy.polynomial import legendre as L
from q2_lib import NB, LOGX
def basis(deg, x=LOGX): return np.stack([L.legval(x, np.eye(deg + 1)[k]) for k in range(deg + 1)], 1) if deg >= 0 else np.zeros((len(x), 0))
def whiten(S):
    ev, U = np.linalg.eigh(S); return U @ np.diag(ev ** -0.5) @ U.T   # symmetric Sigma^{-1/2}
def gls(y, S, deg, x=LOGX):
    W = whiten(S); X = basis(deg, x); yw = W @ y
    if deg < 0: return np.zeros_like(y), dict(chi2=float(yw @ yw), k=0, loo=float(np.mean(yw ** 2)), H=np.zeros((len(y), len(y))), beta=np.zeros(0))
    Xw = W @ X; Hw = Xw @ np.linalg.solve(Xw.T @ Xw, Xw.T); fw = Hw @ yw; e = yw - fw
    f = X @ np.linalg.solve(Xw.T @ Xw, Xw.T @ yw); loo = float(np.mean((e / (1 - np.diag(Hw))) ** 2))
    H = X @ np.linalg.solve(Xw.T @ Xw, Xw.T @ W @ W)   # f = H y
    return f, dict(chi2=float(e @ e), k=deg + 1, loo=loo, H=H, beta=np.linalg.solve(Xw.T @ Xw, Xw.T @ yw))
def penalty(x=LOGX):
    n = len(x); D = np.zeros((n - 2, n)); h = np.mean(np.diff(x))
    for i in range(1, n - 1):
        a, b = x[i] - x[i - 1], x[i + 1] - x[i]; D[i - 1, i - 1] = 2 / (a * (a + b)); D[i - 1, i] = -2 / (a * b); D[i - 1, i + 1] = 2 / (b * (a + b))
    D *= h ** 2; return D.T @ D
K0 = penalty()
LAMS = 10.0 ** np.arange(-3, 6.01, 0.25)
def whit(y, S, lam, K=K0):
    Si = np.linalg.inv(S); H = np.linalg.solve(Si + lam * K, Si); f = H @ y; r = y - f; W = whiten(S)
    Hw = W @ H @ np.linalg.inv(W); e = W @ r; loo = float(np.mean((e / (1 - np.diag(Hw))) ** 2)); return f, dict(chi2=float(e @ e), edf=float(np.trace(H)), loo=loo, H=H, lam=lam)
def whit_cv(y, S, K=K0, crit="loo"):
    best = None
    for lam in LAMS:
        f, i = whit(y, S, lam, K); c = i["loo"] if crit == "loo" else i["chi2"] + 2 * i["edf"]
        if best is None or c < best[0]: best = (c, f, i)
    return best[1], best[2]
def pval(chi2, ndf): return float(stats.chi2.sf(chi2, ndf)) if ndf > 0 else np.nan
