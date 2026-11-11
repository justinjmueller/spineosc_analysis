"""Q2: vectorised smoothers operating on batches of ratio vectors Y (K x n) with a common covariance S.
 poly(BIC): GLS Legendre polynomials in ln E, degree chosen by BIC among {zero (optional), 0..dmax}.
 whittaker(LOO-CV): GLS penalised smoother with 2nd-divided-difference penalty, lambda from leave-one-out CV on the whitened residuals.
 edge-aware variants: fit bins 0..n-2 and treat the last bin separately (raw, or precision-weighted with the extrapolation + empirical edge variance tau^2)."""
import numpy as np
from numpy.polynomial import legendre as L
from q2_fit import whiten, penalty, LAMS
def _basis(deg, x): return np.stack([L.legval(x, np.eye(deg + 1)[k]) for k in range(deg + 1)], 1)
class Smoother:
    def __init__(self, S, x, dmax=4, allow_zero=True, lams=LAMS):
        n = len(x); self.n, self.x, self.S = n, x, S; W = whiten(S); self.W = W; self.degs = ([-1] if allow_zero else []) + list(range(0, dmax + 1)); self.A = {}; self.H = {}; self.k = {}; self.P = {}
        for dg in self.degs:
            if dg < 0: self.H[dg] = np.zeros((n, n)); self.A[dg] = W; self.k[dg] = 0; self.P[dg] = np.zeros((n, 0)); continue
            X = _basis(dg, x); Xw = W @ X; M = np.linalg.inv(Xw.T @ Xw); self.H[dg] = X @ M @ Xw.T @ W; self.A[dg] = W @ (np.eye(n) - self.H[dg]); self.k[dg] = dg + 1
            self.P[dg] = X @ M   # beta = M Xw^T W y  -> not needed
        K = penalty(x); Si = np.linalg.inv(S); self.lams = lams; self.WA = []; self.WH = []; self.WD = []; self.edf = []
        for lam in lams:
            H = np.linalg.solve(Si + lam * K, Si); Hw = W @ H @ np.linalg.inv(W); self.WH.append(H); self.WA.append(W @ (np.eye(n) - H)); self.WD.append(1 / (1 - np.diag(Hw))); self.edf.append(np.trace(H))
    def poly_bic(self, Y, dmin=-1, dmax=4):
        degs = [d for d in self.degs if dmin <= d <= dmax]; bic = np.stack([((Y @ self.A[d].T) ** 2).sum(1) + self.k[d] * np.log(self.n) for d in degs], 1); sel = np.argmin(bic, 1)
        F = np.zeros_like(Y)
        for j, d in enumerate(degs): m = sel == j; F[m] = Y[m] @ self.H[d].T
        return F, np.array(degs)[sel]
    def poly_fixed(self, Y, d): return Y @ self.H[d].T
    def whit_cv(self, Y):
        loo = np.stack([(((Y @ A.T) * D) ** 2).mean(1) for A, D in zip(self.WA, self.WD)], 1); sel = np.argmin(loo, 1); F = np.zeros_like(Y)
        for j in np.unique(sel): m = sel == j; F[m] = Y[m] @ self.WH[j].T
        return F, sel
def edge_extrap(Y, S, x, deg=2):
    """fit bins 0..n-2 with GLS poly(deg), return prediction and its variance at the last bin"""
    n = len(x); keep = np.arange(n - 1); Sk = S[np.ix_(keep, keep)]; W = whiten(Sk); X = _basis(deg, x[keep]); Xw = W @ X; M = np.linalg.inv(Xw.T @ Xw)
    Xh = _basis(deg, x[[n - 1]]); coef = Xh @ M @ Xw.T @ W; pred = Y[:, keep] @ coef.T; return pred[:, 0], float((Xh @ M @ Xh.T)[0, 0]), coef[0]
class Estimator:
    """Complete recommended-class estimator. mode: 'poly', 'whit'; last: 'smooth' (last bin treated like the others), 'raw', 'shrink' (precision weighted with extrapolation + tau2)."""
    def __init__(self, S, x, mode="poly", last="shrink", tau2=0.0, dmax=3, allow_zero=True):
        self.S, self.x, self.mode, self.last, self.tau2, self.dmax, self.az = S, x, mode, last, tau2, dmax, allow_zero; n = len(x); self.n = n
        if last == "smooth": self.sm = Smoother(S, x, dmax, allow_zero)
        elif last == "inflate":
            S2 = S.copy(); S2[-1, -1] += tau2; self.sm = Smoother(S2, x, dmax, allow_zero)
        else:
            keep = np.arange(n - 1); self.sm = Smoother(S[np.ix_(keep, keep)], x[keep], dmax, allow_zero)
            self.deg_ext = max(1, min(dmax, 2))
    def __call__(self, Y):
        n = self.n
        if self.last in ("smooth", "inflate"):
            F, sel = self.sm.poly_bic(Y, -1 if self.az else 0, self.dmax) if self.mode == "poly" else self.sm.whit_cv(Y); return F, sel
        F0, sel = self.sm.poly_bic(Y[:, :n - 1], -1 if self.az else 0, self.dmax) if self.mode == "poly" else self.sm.whit_cv(Y[:, :n - 1])
        out = np.zeros_like(Y); out[:, :n - 1] = F0
        if self.last == "raw": out[:, n - 1] = Y[:, n - 1]
        elif self.last == "pretest":
            pred, vp, _ = edge_extrap(Y, self.S, self.x, self.deg_ext); v0 = self.S[n - 1, n - 1]; vx = vp + self.tau2; shr = (Y[:, n - 1] / v0 + pred / vx) / (1 / v0 + 1 / vx)
            pull = (Y[:, n - 1] - pred) / np.sqrt(v0 + vp); out[:, n - 1] = np.where(np.abs(pull) < 2.0, shr, Y[:, n - 1])
        else:
            pred, vp, _ = edge_extrap(Y, self.S, self.x, self.deg_ext)
            if self.mode == "poly" and self.az:   # if the bulk model chose "zero", the extrapolation is zero as well
                pass
            v0 = self.S[n - 1, n - 1]; vx = vp + self.tau2; out[:, n - 1] = (Y[:, n - 1] / v0 + pred / vx) / (1 / v0 + 1 / vx)
        return out, sel
