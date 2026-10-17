#!/usr/bin/env python3
"""Tier A power check in drift-distance regions: can data resolve the DENT response shift locally?

The stage-1 map shows the DENT response shift depends on drift distance and differs between the two
halves of the detector, but not on (y, z). This repeats the simulation-only fake-data power check of
tier_a_power.py in regions defined by the half (sign of the reco vertex x) and the drift distance
``|x|``, for the DENT pair and the selection S2. No data are read.

    python scripts/dent_improvement/tier_a_power_regions.py [--pot 4.5e18 1e20]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tier_a_power as ta  # noqa: E402

DRIFT_EDGES = [5, 50, 100, 150, 205]


def regions(frame: pd.DataFrame) -> pd.Series:
    """Region label (half, drift bin) of each event from the reco vertex x."""
    drift = pd.cut(frame.reco_vertex_x.abs(), DRIFT_EDGES, labels=False)
    half = np.where(frame.reco_vertex_x < 0, "x<0", "x>0")
    labels = [f"{h}, |x| {DRIFT_EDGES[int(d)]}-{DRIFT_EDGES[int(d) + 1]}" if d == d else "out"
              for h, d in zip(half, drift, strict=True)]
    return pd.Series(labels, index=frame.index)


def main() -> None:
    """Compute d' per observable and region."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--pot", type=float, nargs="+", default=[4.5e18, 1.0e20])
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--throws", type=int, default=200)
    args = parser.parse_args()
    rng = np.random.default_rng(20261015)
    rows = []
    for study in ta.STUDIES:
        variable, edges = ta.read_binning(study)
        branches = [variable, "reco_vertex_x"]
        cv, pot_cv = ta.load("sbnd_detvar_newcv", branches)
        dent, _ = ta.load("sbnd_detvar_newdent", branches)
        cv, dent = cv[ta.selection(cv, "s2")], dent[ta.selection(dent, "s2")]
        cv_region, dent_region = regions(cv), regions(dent)
        for name in sorted(set(cv_region) - {"out"}):
            samples = {"cv": cv.loc[cv_region == name, variable].to_numpy(),
                       "dent": dent.loc[dent_region == name, variable].to_numpy()}
            for pot in args.pot:
                n_expected = len(samples["cv"]) * pot / pot_cv
                deltas = []
                for _ in range(args.repeats):
                    counts = {}
                    for sample, values in samples.items():
                        in_a = rng.random(len(values)) < 0.5
                        for half, chosen in (("A", values[in_a]), ("B", values[~in_a])):
                            hist = np.histogram(chosen, edges)[0]
                            counts[f"{sample}_{half}"] = hist.astype(float)
                    deltas.append(ta.separation(counts, n_expected, rng, args.throws))
                mean_cv, mean_dent, sd = np.mean(deltas, axis=0)
                rows.append({"study": study, "region": name, "pot": pot,
                             "expected_events": round(n_expected),
                             "d_prime": (mean_dent - mean_cv) / sd})
        print(f"done {study}", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(REPO / "docs" / "dent_improvement" / "tier_a_power_regions.csv", index=False)
    small = table[np.isclose(table.pot, min(args.pot))]
    print(f"\nd' at {min(args.pot):.1e} POT:")
    print(small.pivot_table(index="study", columns="region", values="d_prime").round(2).to_string())


if __name__ == "__main__":
    main()
