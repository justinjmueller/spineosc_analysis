"""Replica of PROfit's shape-only chi2 for single-channel studies, validated against PROfit (see scripts/dent_study/export_covariance.py).
    chi2 = (n - p)^T [ diag(n) + F o (p p^T) ]^-1 (n - p),   p = prediction area-normalised to the data integral,
    F = PROfit's collapsed total fractional covariance (already shape-projected under --shapeonly), n = data counts.
Its use here: (1) hold the covariance fixed while the DENT strength changes (PROfit's --inject-systs-cv re-derives it around the shifted
prediction), and (2) run many pseudo-experiments (Poisson + Gaussian systematic throws through F) cheaply, so p-values are not limited
to the 50 PROfit throws. Nothing is fitted: DENT strength is only ever 0, 0.5 or 1.0."""
import pathlib, numpy as np
RDIR = pathlib.Path(__file__).resolve().parents[2] / "runs" / "dent_study" / "replica"
THETAS = (0.0, 0.5, 1.0)

def load(obs):
    z = np.load(RDIR / f"{obs}.npz"); return {k: z[k] for k in z.files}

def chi2(n, p, F, p_cov=None):
    """p_cov: the prediction used for the covariance outer product (default p itself)."""
    pc = p if p_cov is None else p_cov
    C = np.diag(n) + F * np.outer(pc, pc)
    r = n - p
    return float(r @ np.linalg.solve(C, r))

def validate(d):
    """replica minus PROfit chi2 at each strength (PROfit's own covariance at that strength)."""
    return {t: chi2(d["n"], d[f"p_{t:g}"], d[f"F_{t:g}"]) - float(d[f"chi_profit_{t:g}"]) for t in THETAS}

def fixed_cov_dchi2(d, n=None, scale=None):
    """chi2(theta=0) - chi2(theta=1) with the covariance HELD at its theta=0 form: F0 o (p0' p0'^T) for both predictions,
    where p' are rescaled to the (possibly fake) data integral. Positive: the DENT-shifted prediction describes the data better."""
    n = d["n"] if n is None else n; s = n.sum() / d["n"].sum()
    p0, p1 = d["p_0"] * s, d["p_1"] * s
    return chi2(n, p0, d["F_0"], p_cov=p0) - chi2(n, p1, d["F_0"], p_cov=p0)

def pseudo_experiments(d, truth, nthrow=20000, seed=7):
    """Fake data: expectation = prediction at strength `truth`, plus a Gaussian systematic shift through F_0 o (mu mu^T), plus Poisson."""
    rng = np.random.default_rng(seed); mu = d[f"p_{truth:g}"]; S = d["F_0"] * np.outer(mu, mu)
    w, V = np.linalg.eigh((S + S.T) / 2); L = V * np.sqrt(np.clip(w, 0, None))
    out = np.empty(nthrow)
    for i in range(nthrow):
        shifted = np.clip(mu + L @ rng.standard_normal(len(mu)), 0, None)
        out[i] = fixed_cov_dchi2(d, n=rng.poisson(shifted).astype(float))
    return out

def summarize(d, nthrow=20000):
    obs = fixed_cov_dchi2(d); res = dict(dchi2_fixed_data=obs)
    for t in THETAS:
        x = pseudo_experiments(d, t, nthrow); lo = (np.sum(x <= obs) + 1) / (len(x) + 1); hi = (np.sum(x >= obs) + 1) / (len(x) + 1)
        res[f"mean_{t:g}"] = float(x.mean()); res[f"sd_{t:g}"] = float(x.std(ddof=1)); res[f"p2_{t:g}"] = float(min(1.0, 2 * min(lo, hi)))
    res["dprime"] = (res["mean_1"] - res["mean_0"]) / np.sqrt(0.5 * (res["sd_0"] ** 2 + res["sd_1"] ** 2))
    return res
