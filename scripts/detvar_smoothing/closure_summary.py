#!/usr/bin/env python3
"""Summarize the multi-systematic sample closure across ratio variants (paired by split).

For each variant (unsmoothed or a smoother) and systematic: the fitted pull of the parameter under
test in the fit to the variation sample, the same pull in the null fit (pseudo-data from the
reference sample), the excess chi-square (varied minus null), and their changes relative to the
unsmoothed ratios split by split. A smoother helps if it keeps the pull of a real effect, lowers the
spurious pull and the excess chi-square of noise-dominated cases, and shrinks the spread across
splits.

    python scripts/detvar_smoothing/closure_summary.py [closure2m_muon_score_inclusive_s2.csv]
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SYSTS = ["ChargeScale15", "MattModYZ", "DENT", "WireModxThetaXW"]


def per_split(table: pd.DataFrame) -> pd.DataFrame:
    """One row per variant, systematic and split: pulls, chi-squares and the excess."""
    rows = []
    for (variant, perm), group in table.groupby(["variant", "perm"]):
        samples = group.set_index("data")
        for syst in SYSTS:
            reference = samples.loc["null_dent" if syst == "DENT" else "null_nominal"]
            varied = samples.loc[f"varied_{syst}"]
            rows.append({"variant": variant, "perm": perm, "syst": syst,
                         "pull_varied": varied[f"fit_{syst}"],
                         "pull_null": reference[f"fit_{syst}"],
                         "excess_chi2": varied.chi2_total - reference.chi2_total,
                         "chi2_varied": varied.chi2_total, "chi2_null": reference.chi2_total})
    return pd.DataFrame(rows)


def main() -> None:
    """Print the per-variant table and the paired changes relative to unsmoothed."""
    name = sys.argv[1] if len(sys.argv) > 1 else "closure2m_muon_score_inclusive_s2.csv"
    table = pd.read_csv(REPO / "docs" / "detvar_smoothing" / name, keep_default_na=False)
    splits = per_split(table)
    splits.to_csv(REPO / "docs" / "detvar_smoothing" / name.replace("closure2m_", "by_split_"),
                  index=False)
    summary = splits.groupby(["syst", "variant"]).agg(
        pull_varied=("pull_varied", "mean"), pull_sd=("pull_varied", "std"),
        pull_null=("pull_null", "mean"), excess_mean=("excess_chi2", "mean"),
        excess_sd=("excess_chi2", "std"), chi2_null=("chi2_null", "mean"), n=("perm", "count"))
    pd.set_option("display.width", 200)
    print(summary.round(2).to_string())
    base = splits[splits.variant == "unsmoothed"].set_index(["syst", "perm"])
    changes = []
    for variant in sorted(set(splits.variant) - {"unsmoothed"}):
        current = splits[splits.variant == variant].set_index(["syst", "perm"])
        delta = pd.DataFrame({"d_pull_varied": current.pull_varied - base.pull_varied,
                              "d_excess_chi2": current.excess_chi2 - base.excess_chi2,
                              "d_chi2_null": current.chi2_null - base.chi2_null}).dropna()
        delta["variant"] = variant
        changes.append(delta.reset_index())
    if changes:
        paired = pd.concat(changes).groupby(["syst", "variant"]).agg(["mean", "std"])
        print("\npaired change relative to unsmoothed (same split):")
        print(paired.round(2).to_string())


if __name__ == "__main__":
    main()
