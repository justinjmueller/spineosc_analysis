#!/usr/bin/env python3
"""Exports compact result tables from the untracked run directory to docs/dent_improvement/results/.

* ``tier_b_power.csv``: the fake-data power check of the 12 Tier B studies. Only the
  pseudo-experiment (``tierC_pe``) field of each result file is read; no data chi-square is read
  or written.
* ``regional_chi2.csv``: the data chi-square of the regional Tier A comparison (muon calorimetric
  difference, S2), per region and sample. These are the values behind
  ``regional_contrast_data.csv``.

    python scripts/dent_improvement/export_results.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import regional_contrast as rc  # noqa: E402

RESULTS = REPO / "runs" / "dent_study" / "results"
OUT = REPO / "docs" / "dent_improvement" / "results"
TIER_B = ["muon_score_inclusive", "muon_primary_score_inclusive", "proton_score_inclusive",
          "proton_primary_score_inclusive", "muon_calo_csda", "proton_calo_csda"]


def tier_b_power() -> pd.DataFrame:
    """Mean dchi2 for fake data from CV, half DENT and full DENT, and the separation d'."""
    rows = []
    for study in TIER_B:
        for selection in ("s0", "s2"):
            result = json.loads((RESULTS / f"dimp_{study}_{selection}.json").read_text())
            pe = result["tierC_pe"]          # the only field read
            chi2 = {(e["truth"], e["seed"], e["test"]): float(np.sum(e["chi2"])) for e in pe}
            delta = {truth: np.array([chi2[truth, s, 0.0] - chi2[truth, s, 1.0]
                                      for s in range(1, 51)]) for truth in ("cv", "half", "dent")}
            sd = np.sqrt(0.5 * (delta["cv"].var(ddof=1) + delta["dent"].var(ddof=1)))
            rows.append({"study": study, "selection": selection, "seeds": len(delta["cv"]),
                         "mean_dchi2_truth_cv": delta["cv"].mean(),
                         "mean_dchi2_truth_half": delta["half"].mean(),
                         "mean_dchi2_truth_dent": delta["dent"].mean(), "sd": sd,
                         "d_prime": (delta["dent"].mean() - delta["cv"].mean()) / sd})
    return pd.DataFrame(rows)


def regional_chi2() -> pd.DataFrame:
    """Data chi-square of every regional Tier A sample (muon calorimetric difference, S2)."""
    rows = []
    for region in rc.region_names():
        prefix = f"dimpAR_{rc.STUDY}_s2_{region}_"
        for path in sorted(RESULTS.glob(f"{prefix}*.json")):
            sample = path.stem.removeprefix(prefix)
            rows.append({"region": region, "sample": sample, "chi2": rc.chi2(region, sample),
                         "ndf": 19})
    return pd.DataFrame(rows)


def main() -> None:
    """Write both tables."""
    OUT.mkdir(parents=True, exist_ok=True)
    tier_b_power().to_csv(OUT / "tier_b_power.csv", index=False)
    table = regional_chi2()
    table.to_csv(OUT / "regional_chi2.csv", index=False)
    print(f"wrote {OUT / 'tier_b_power.csv'} (12 rows) and {OUT / 'regional_chi2.csv'} "
          f"({len(table)} rows)")


if __name__ == "__main__":
    main()
