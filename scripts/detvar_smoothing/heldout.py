#!/usr/bin/env python3
"""Held-out predictive error of smoothed detector-systematic ratios, on the real samples.

The three-way event split gives three independent ratios of each systematic (parts 0, 1, 2). A
smoother is applied to the ratio of part i and the result is compared with the *raw* ratio of a
different part j, per bin, in units of that raw ratio's statistical sigma:

    chi2(i, j) = sum_b (smoothed_i,b - raw_j,b)^2 / sigma_j,b^2

If the smoothed ratio were exactly the truth the expectation is the number of bins n (the noise of
part j alone). Noise left in part i adds to it (an unsmoothed ratio gives about 2 n) and so does any
real structure the smoother removes. No fit and no truth are needed, so this works for any
observable, and it is much more sensitive than a global fit with every systematic floating.
All six ordered pairs of parts are used; the spread of chi2 across them is the uncertainty on a
mean, and paired differences relative to the unsmoothed ratio are reported.

    python scripts/detvar_smoothing/heldout.py muon_score_inclusive_s2 reco_energy_inclusive_s2
"""
import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
import smoother_bias_variance as bv  # noqa: E402
import split_closure as sc  # noqa: E402

REPO = sc.REPO
METHODS = ["none", "constant", "poly1", "poly2", "poly3", "poly_cv", "roll3", "roll5", "median3",
           "local1.5", "local3", "gauss1.5", "gauss3", "spline_gcv"]


def run(study: str) -> pd.DataFrame:
    """Held-out chi-square of every method, systematic and ordered pair of parts."""
    variable, source = sc.STUDIES[study]
    original = uproot.open(REPO / source)
    window = sc.fit_range(study)
    rows = []
    for syst in sc.PAIRS:
        edges = original[f"{syst}_sbnd_Run1"].to_numpy()[1]
        centers = 0.5 * (edges[1:] + edges[:-1])
        inside = (centers > window[0]) & (centers < window[1])
        parts = {part: sc.ratio_from_part(syst, variable, edges, part) for part in range(3)}
        for i, j in itertools.permutations(range(3), 2):
            ratio_i, sigma_i = (v[inside] for v in parts[i])
            ratio_j, sigma_j = (v[inside] for v in parts[j])
            smoothed = bv.smooth_all(ratio_i, np.maximum(sigma_i, 1e-4))
            for method in METHODS:
                chi2 = (((smoothed[method] - ratio_j) / sigma_j) ** 2).sum()
                rows.append({"study": study, "syst": syst, "method": method, "i": i, "j": j,
                             "bins": int(inside.sum()), "chi2": chi2})
        print(f"done {study} {syst}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    """Run each study and print the held-out chi-square per bin and the paired change."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("studies", nargs="+")
    args = parser.parse_args()
    pd.set_option("display.width", 220)
    for study in args.studies:
        table = run(study)
        table.to_csv(REPO / "docs" / "detvar_smoothing" / f"heldout_{study}.csv", index=False)
        table["chi2_per_bin"] = table.chi2 / table.bins
        mean = table.groupby(["method", "syst"]).chi2_per_bin.mean().unstack()
        mean["all"] = mean.mean(axis=1)
        print(f"\n{study}: held-out chi2 per bin (1 = a perfect smoother, ~2 = unsmoothed)")
        print(mean.loc[METHODS].round(2).to_string())
        base = table[table.method == "none"].set_index(["syst", "i", "j"]).chi2_per_bin
        gain = []
        for method in METHODS[1:]:
            current = table[table.method == method].set_index(["syst", "i", "j"]).chi2_per_bin
            gain.append((current - base).rename(method))
        gain = pd.concat(gain, axis=1)
        summary = gain.groupby(level="syst").mean().T
        summary["all"] = summary.mean(axis=1)
        summary["se_all"] = gain.groupby(level=["i", "j"]).mean().std().to_numpy() / np.sqrt(6)
        print("\npaired change in chi2 per bin relative to unsmoothed (negative = better):")
        print(summary.round(2).to_string())


if __name__ == "__main__":
    sys.exit(main())
