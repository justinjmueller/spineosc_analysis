#!/usr/bin/env python3
"""Tier A summary: does DENT improve the data agreement more than the other variations?

For each study and selection that passes the pre-registered power floor, reads the shape-only
chi-square of the detector-access data against each Tier A sample and forms, for every variation,
the observed ``dchi2 = chi2(reference) - chi2(varied)`` (positive: the variation agrees better
with the data). The reference is the paired CV sample for DENT (and for the older DENT pair) and
the production nominal sample for the other variations (a nominal without the trigger block for
the WireMod pair, whose samples have no truth or flash information). The simulated expectation of
``dchi2`` if the reference or the variation were the truth, from tier_a_power.py, is shown next
to it.

Only the studies passing the floor are read; the others are not examined.

    python scripts/dent_improvement/tier_a_summary.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "runs" / "dent_study" / "results"
POT = 4.5e18          # the conservative assumed sample size used for the power floor
FLOOR = 0.5
REFERENCE = {"dent": "cv", "olddent": "oldcv", "wiremodxtheta": "nominal_noflash",
             "wiremodyz": "nominal_noflash"}      # all others use "nominal"


def chi2(study: str, selection: str, sample: str) -> tuple[float, int]:
    """Data chi-square and ndf of one Tier A sample from its result file."""
    result = json.loads((RESULTS / f"dimpA_{study}_{selection}_{sample}.json").read_text())
    rows = result["chi2"]["v1"]
    if not result.get("shapeonly_live") or len(rows) != 1:
        raise ValueError(f"unexpected result for {study} {selection} {sample}")
    return rows[0]["chi2"], rows[0]["ndf"]


def main() -> None:
    """Print the table for every study and selection passing the power floor."""
    power = pd.read_csv(REPO / "docs" / "dent_improvement" / "tier_a_power_all.csv")
    power = power[np.isclose(power.pot, POT)]
    dent = power[power.variation == "dent"].set_index(["study", "selection"]).d_prime
    passed = [key for key, value in dent.items() if value >= FLOOR]
    print(f"studies passing d' >= {FLOOR} for DENT at {POT:.1e} POT: "
          f"{len(passed)} of {len(dent)}\n")
    rows = []
    for study, selection in passed:
        for variation in power.variation.unique():
            reference = REFERENCE.get(variation, "nominal")
            sample = {"dent": "dent", "olddent": "olddent"}.get(variation, variation)
            c_ref, ndf = chi2(study, selection, reference)
            c_var, _ = chi2(study, selection, sample)
            match = ((power.study == study) & (power.selection == selection)
                     & (power.variation == variation))
            p = power[match].iloc[0]
            observed = c_ref - c_var
            rows.append({"study": study, "sel": selection, "variation": variation, "ndf": ndf,
                         "chi2_ref": c_ref, "chi2_var": c_var, "dchi2_obs": observed,
                         "exp_if_ref": p.mean_dchi2_ref, "exp_if_var": p.mean_dchi2_var,
                         "sd": p.sd, "d_prime": p.d_prime,
                         "obs_over_expected_shift":
                             observed / (p.mean_dchi2_var - p.mean_dchi2_ref)})
    table = pd.DataFrame(rows)
    table.to_csv(REPO / "docs" / "dent_improvement" / "tier_a_summary.csv", index=False)
    show = table.round({"chi2_ref": 0, "chi2_var": 0, "dchi2_obs": 0, "exp_if_ref": 1,
                        "exp_if_var": 1, "sd": 1, "d_prime": 1, "obs_over_expected_shift": 1})
    for (study, selection), group in show.groupby(["study", "sel"], sort=False):
        print(f"== {study}  {selection}  (ndf {group.ndf.iloc[0]})")
        columns = ["variation", "chi2_ref", "chi2_var", "dchi2_obs", "exp_if_ref", "exp_if_var",
                   "d_prime", "obs_over_expected_shift"]
        print(group[columns].to_string(index=False))
        better = (group.dchi2_obs > 0).sum()
        dent_value = group[group.variation == "dent"].dchi2_obs.iloc[0]
        rank = int((group.dchi2_obs > dent_value).sum()) + 1
        print(f"   variations that agree better with the data than their reference: {better} of "
              f"{len(group)}; DENT ranks {rank} of {len(group)} by dchi2\n")


if __name__ == "__main__":
    main()
