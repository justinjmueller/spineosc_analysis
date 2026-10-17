#!/usr/bin/env python3
"""Stage 3: the half-asymmetry contrast of the regional Tier A comparison.

The contrast is ``C = sum over drift bins 50-100, 100-150, 150-205 of [dchi2(x>0) - dchi2(x<0)]``
with ``dchi2 = chi2(reference) - chi2(variation)`` per region, for the muon calorimetric difference
under S2.
DENT is predicted to give C > 0 (its response shift persists through the drift in the x>0 half and
decays in the x<0 half); half-symmetric variations should give C near 0.

    python scripts/dent_improvement/regional_contrast.py simulate     # expectations, no data
    python scripts/dent_improvement/regional_contrast.py data         # observed C next to them

``simulate`` thrown fake data only (the half-split procedure of tier_a_power.py, region by region;
regions are disjoint, so their draws are independent). ``data`` reads the regional Tier A results.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tier_a_power as ta  # noqa: E402

STUDY = "muon_calo_csda"
DRIFT_EDGES = [5, 50, 100, 150, 205]
HALVES = ("xneg", "xpos")
CONTRAST_BINS = [(50, 100), (100, 150), (150, 205)]
EXPECTATION_FILE = REPO / "docs" / "dent_improvement" / "regional_contrast_expectation.csv"
REFERENCE_SAMPLE = {"dent": "cv", "olddent": "oldcv", "wiremodxtheta": "nominal_noflash",
                    "wiremodyz": "nominal_noflash"}      # every other variation uses "nominal"
RESULTS = REPO / "runs" / "dent_study" / "results"


def region_names() -> list[str]:
    """Names of the eight regions, matching build_xmls_regions.py."""
    return [f"{h}_{lo}_{hi}" for h in HALVES for lo, hi in zip(DRIFT_EDGES[:-1], DRIFT_EDGES[1:],
                                                               strict=True)]


def assign(frame: pd.DataFrame) -> pd.Series:
    """Region name of each event from the reco vertex x."""
    drift = pd.cut(frame.reco_vertex_x.abs(), DRIFT_EDGES, labels=False)
    half = np.where(frame.reco_vertex_x < 0, "xneg", "xpos")
    names = [f"{h}_{DRIFT_EDGES[int(d)]}_{DRIFT_EDGES[int(d) + 1]}" if d == d else "out"
             for h, d in zip(half, drift, strict=True)]
    return pd.Series(names, index=frame.index)


def contrast_weights() -> dict[str, int]:
    """+1 for the x>0 regions of the contrast, -1 for the x<0 regions."""
    weights = {}
    for lo, hi in CONTRAST_BINS:
        weights[f"xpos_{lo}_{hi}"] = 1
        weights[f"xneg_{lo}_{hi}"] = -1
    return weights


def simulate(pots: list[float], repeats: int, throws: int) -> None:
    """Expected regional dchi2 and contrast for every variation, from fake data only."""
    rng = np.random.default_rng(20261016)
    variable, edges = ta.read_binning(STUDY)
    cache: dict = {}
    rows = []
    for label, ((ref_tree, ref_path), var_tree) in ta.PAIRS.items():
        flash = label not in ta.NO_FLASH
        ref, pot_ref = ta.cached_load(cache, ref_tree, ref_path, variable, ["reco_vertex_x"])
        var, _ = ta.cached_load(cache, var_tree, ta.cs.DET, variable, ["reco_vertex_x"])
        ref, var = ref[ta.selection(ref, "s2", flash)], var[ta.selection(var, "s2", flash)]
        ref_region, var_region = assign(ref), assign(var)
        for pot in pots:
            for region in region_names():
                samples = {"cv": ref.loc[ref_region == region, variable].to_numpy(),
                           "dent": var.loc[var_region == region, variable].to_numpy()}
                n_expected = len(samples["cv"]) * pot / pot_ref
                deltas = []
                for _ in range(repeats):
                    counts = {}
                    for sample, values in samples.items():
                        in_a = rng.random(len(values)) < 0.5
                        for half, chosen in (("A", values[in_a]), ("B", values[~in_a])):
                            hist = np.histogram(chosen, edges)[0]
                            counts[f"{sample}_{half}"] = hist.astype(float)
                    deltas.append(ta.separation(counts, n_expected, rng, throws))
                mean_ref, mean_var, sd = np.mean(deltas, axis=0)
                rows.append({"variation": label, "pot": pot, "region": region,
                             "expected_events": round(n_expected), "exp_if_ref": mean_ref,
                             "exp_if_var": mean_var, "sd": sd})
        print(f"done {label}", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(EXPECTATION_FILE, index=False)
    print(contrast_table(table).round(1).to_string())


def contrast_table(table: pd.DataFrame) -> pd.DataFrame:
    """Expected contrast if the reference / the variation is true, and its sd, per variation."""
    weights = contrast_weights()
    rows = []
    for (variation, pot), group in table.groupby(["variation", "pot"], sort=False):
        group = group[group.region.isin(weights)]
        sign = group.region.map(weights)
        rows.append({"variation": variation, "pot": pot,
                     "C_if_ref": float((sign * group.exp_if_ref).sum()),
                     "C_if_var": float((sign * group.exp_if_var).sum()),
                     "sd_C": float(np.sqrt((group.sd ** 2).sum()))})
    return pd.DataFrame(rows)


def chi2(region: str, sample: str) -> float:
    """Data chi-square of one regional Tier A result."""
    result = json.loads((RESULTS / f"dimpAR_{STUDY}_s2_{region}_{sample}.json").read_text())
    rows = result["chi2"]["v1"]
    if not result["steps"]["v1"]["ok"] or not result.get("shapeonly_live") or len(rows) != 1:
        raise ValueError(f"unexpected result for {region} {sample}")
    return rows[0]["chi2"]


def data(pot: float) -> None:
    """Observed regional dchi2 and the contrast next to their expectations."""
    expectation = pd.read_csv(EXPECTATION_FILE)
    expectation = expectation[np.isclose(expectation.pot, pot)]
    weights = contrast_weights()
    expected = contrast_table(expectation).set_index("variation")
    rows, regional = [], []
    for variation in expected.index:
        reference = REFERENCE_SAMPLE.get(variation, "nominal")
        observed = {region: chi2(region, reference) - chi2(region, variation)
                    for region in region_names()}
        regional.append({"variation": variation, **{r: round(v) for r, v in observed.items()}})
        c_obs = sum(weights[r] * observed[r] for r in weights)
        e = expected.loc[variation]
        rows.append({"variation": variation, "C_obs": c_obs, "C_if_ref": e.C_if_ref,
                     "C_if_var": e.C_if_var, "sd_C": e.sd_C,
                     "z_vs_ref": (c_obs - e.C_if_ref) / e.sd_C,
                     "z_vs_var": (c_obs - e.C_if_var) / e.sd_C})
    pd.set_option("display.width", 220)
    print("observed regional dchi2 (reference - variation):")
    print(pd.DataFrame(regional).to_string(index=False))
    table = pd.DataFrame(rows)
    table.to_csv(REPO / "docs" / "dent_improvement" / "regional_contrast_data.csv", index=False)
    print(f"\ncontrast C at the {pot:.1e} POT assumption:")
    print(table.round(1).to_string(index=False))


def main() -> None:
    """Run the requested part."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("part", choices=["simulate", "data"])
    parser.add_argument("--pot", type=float, nargs="+", default=[4.5e18, 1.0e20])
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--throws", type=int, default=200)
    args = parser.parse_args()
    if args.part == "simulate":
        simulate(args.pot, args.repeats, args.throws)
    else:
        data(args.pot[0])


if __name__ == "__main__":
    main()
