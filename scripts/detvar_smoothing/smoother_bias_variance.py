#!/usr/bin/env python3
"""Bias and variance of smoothers for detector-systematic ratios, against a known truth.

A smoother can only be judged against a truth, and a real ratio has none. This study therefore
builds scenarios from a known truth ratio and noise at the real per-bin level: for each of the four
systematics the per-bin statistical sigma is taken from the paired bootstrap
(docs/detvar_stat/bootstrap.csv), scaled to a one-third sample (the size of a part in
split_closure.py), and Gaussian noise of that size is added to each truth shape. Each noisy ratio is
smoothed with every method; bias and variance per bin follow exactly from the repeated draws.

Truth shapes (analysis reco-energy binning, 17 bins): ``null`` (ratio 1), ``linear``, ``quadratic``,
``rising`` (a smooth rise like WireModxThetaXW), ``edge`` (a linear shape with a -10% step in the
last bin, like ChargeScale at 1.5-2 GeV: real structure a smoother must not erase) and ``wave``.

Smoothers (weighted by 1 / sigma^2): constant, polynomials of order 1-3, local-linear (tricube) and
Gaussian-kernel fits with fixed bandwidths, the smoothing spline with the penalty chosen by
generalized cross-validation, and cross-validated variants (leave-one-out) that choose the
polynomial order or the bandwidth for each noisy ratio, so nothing is tuned by eye.

    python scripts/detvar_smoothing/smoother_bias_variance.py [--repeats 400]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import make_smoothing_spline

REPO = Path(__file__).resolve().parents[2]
BOOTSTRAP = REPO / "docs" / "detvar_stat" / "bootstrap.csv"
SYSTEMATICS = ["DENT", "ChargeScale15", "WireModxThetaXW", "MattModYZ"]
PART_FACTOR = np.sqrt(3.0)          # a one-third sample has sqrt(3) times the sigma of the full one
BANDWIDTHS = (1.0, 1.5, 2.0, 3.0, 5.0)


def truths(n: int) -> dict[str, np.ndarray]:
    """The truth ratio shapes over n bins."""
    x = np.linspace(0.0, 1.0, n)
    edge = 1 + 0.04 * (x - 0.5)
    edge[-1] -= 0.10
    return {"null": np.ones(n), "linear": 1 + 0.04 * (x - 0.5),
            "quadratic": 1 + 0.06 * (x - 0.5) ** 2 * 4 - 0.01,
            "rising": 1 + 0.12 * x ** 2, "edge": edge, "wave": 1 + 0.03 * np.sin(2 * np.pi * x)}


def poly_fit(x: np.ndarray, y: np.ndarray, w: np.ndarray, order: int, at: np.ndarray) -> np.ndarray:
    """Weighted polynomial fit (order 0 is a weighted mean) evaluated at ``at``."""
    scale = np.ptp(x)
    center = x.mean()
    series = np.polynomial.polynomial
    coefficients = series.polyfit((x - center) / scale, y, order, w=np.sqrt(w))
    return series.polyval((at - center) / scale, coefficients)


def kernel_fit(index: np.ndarray, y: np.ndarray, w: np.ndarray, bandwidth: float, local: bool,
               skip: int = -1) -> np.ndarray:
    """Gaussian-kernel average or local-linear (tricube) fit at every bin; ``skip`` leaves one out.
    """
    out = np.empty(len(y))
    for i in range(len(y)):
        distance = np.abs(index - i) / bandwidth
        kernel = np.where(distance < 1, (1 - distance ** 3) ** 3, 0.0) if local \
            else np.exp(-0.5 * distance ** 2)
        weights = kernel * w
        if skip >= 0:
            weights = np.where(index == skip, 0.0, weights)
        if weights.sum() <= 0:
            out[i] = np.nan if skip >= 0 else y[i]       # nothing left to predict a left-out bin
        elif local and (weights > 0).sum() >= 3:
            design = np.vstack([np.ones(len(y)), index - i]).T * np.sqrt(weights)[:, None]
            out[i] = np.linalg.lstsq(design, y * np.sqrt(weights), rcond=None)[0][0]
        else:
            out[i] = (weights * y).sum() / weights.sum()
    return out


def rolling_fit(y: np.ndarray, w: np.ndarray, half: int, median: bool = False) -> np.ndarray:
    """Rolling window of 2 * half + 1 bins: the weighted mean (or median), truncated at the edges.

    A bin at the edge of the range averages only the bins that exist, so the window is
    asymmetric there and the average is pulled towards the neighbours.
    """
    out = np.empty(len(y))
    for i in range(len(y)):
        low, high = max(0, i - half), min(len(y), i + half + 1)
        out[i] = np.median(y[low:high]) if median else np.average(y[low:high], weights=w[low:high])
    return out


def loo_poly(x: np.ndarray, y: np.ndarray, w: np.ndarray, order: int) -> float:
    """Weighted leave-one-out error of a polynomial of the given order."""
    total = 0.0
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        prediction = poly_fit(x[keep], y[keep], w[keep], order, x[i:i + 1])[0]
        total += w[i] * (y[i] - prediction) ** 2
    return total


def loo_kernel(x: np.ndarray, y: np.ndarray, w: np.ndarray, bandwidth: float, local: bool) -> float:
    """Weighted leave-one-out error of a kernel smoother (infinite if a bin cannot be predicted)."""
    total = 0.0
    for i in range(len(y)):
        prediction = kernel_fit(x, y, w, bandwidth, local, skip=i)[i]
        total += np.inf if np.isnan(prediction) else w[i] * (y[i] - prediction) ** 2
    return total


def smooth_all(y: np.ndarray, sigma: np.ndarray) -> dict[str, np.ndarray]:
    """Every smoothed version of one noisy ratio."""
    x = np.arange(len(y), dtype=float)
    w = 1.0 / sigma ** 2
    out = {"none": y, "constant": poly_fit(x, y, w, 0, x)}
    for order in (1, 2, 3):
        out[f"poly{order}"] = poly_fit(x, y, w, order, x)
    for bandwidth in (1.5, 3.0):
        out[f"local{bandwidth:g}"] = kernel_fit(x, y, w, bandwidth, True)
        out[f"gauss{bandwidth:g}"] = kernel_fit(x, y, w, bandwidth, False)
    out["roll3"] = rolling_fit(y, w, 1)
    out["roll5"] = rolling_fit(y, w, 2)
    out["median3"] = rolling_fit(y, w, 1, median=True)
    out["spline_gcv"] = make_smoothing_spline(x, y, w=np.sqrt(w))(x)
    errors = {order: loo_poly(x, y, w, order) for order in (0, 1, 2, 3)}
    out["poly_cv"] = poly_fit(x, y, w, min(errors, key=errors.get), x)
    for name, local in (("local_cv", True), ("gauss_cv", False)):
        scores = {h: loo_kernel(x, y, w, h, local) for h in BANDWIDTHS}
        out[name] = kernel_fit(x, y, w, min(scores, key=scores.get), local)
    return out


def sigma_profile(syst: str) -> np.ndarray:
    """Per-bin statistical sigma of a one-third-sample ratio, from the paired bootstrap."""
    table = pd.read_csv(BOOTSTRAP)
    return table[table.systematic == syst].sort_values("bin").sigma_stat.to_numpy() * PART_FACTOR


def run(repeats: int, seed: int = 20261015) -> pd.DataFrame:
    """Bias, variance and relative error of every smoother in every scenario."""
    rng = np.random.default_rng(seed)
    rows = []
    for syst in SYSTEMATICS:
        sigma = sigma_profile(syst)
        for shape, truth in truths(len(sigma)).items():
            draws = {}
            for _ in range(repeats):
                noisy = truth + rng.normal(0.0, sigma)
                for name, values in smooth_all(noisy, sigma).items():
                    draws.setdefault(name, []).append(values)
            baseline = (sigma ** 2).sum()
            for name, values in draws.items():
                values = np.array(values)
                bias = values.mean(axis=0) - truth
                variance = values.var(axis=0, ddof=1)
                rows.append({"syst": syst, "truth": shape, "method": name,
                             "bias2": (bias ** 2).sum(), "variance": variance.sum(),
                             "relative_mse": ((bias ** 2).sum() + variance.sum()) / baseline,
                             "bias_last_bin": bias[-1], "sigma_last_bin": sigma[-1]})
        print(f"done {syst}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    """Run the scenarios and print the summaries."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--repeats", type=int, default=400)
    args = parser.parse_args()
    table = run(args.repeats)
    table.to_csv(REPO / "docs" / "detvar_smoothing" / "bias_variance.csv", index=False)
    pd.set_option("display.width", 200)
    print("\nrelative MSE (smoothed / unsmoothed; below 1 is a gain), mean over truth shapes:")
    print(table.groupby(["method", "syst"]).relative_mse.mean().unstack().round(2).to_string())
    print("\nworst case over truth shapes:")
    print(table.groupby(["method", "syst"]).relative_mse.max().unstack().round(2).to_string())
    edge = table[table.truth == "edge"].groupby(["method", "syst"]).bias_last_bin.mean().unstack()
    print("\nbias in the last bin for the edge truth (-10% step), in %:")
    print((100 * edge).round(1).to_string())


if __name__ == "__main__":
    sys.exit(main())
