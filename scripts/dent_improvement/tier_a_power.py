#!/usr/bin/env python3
"""Power of the Tier A comparison: can data tell the detvar DENT sample from the detvar CV sample?

Tier A uses the detvar CV and DENT samples directly as the MC, with MC-statistics covariance only,
so the DENT shift is not absorbed by a systematic covariance. This script estimates, from simulation
alone and before any data are compared, how well that test separates DENT from CV.

For one observable and selection, each detvar sample is split at random into two halves. Fake data
are Poisson throws of the shape of half A; the predictions are the shapes of half B of the CV and of
the DENT sample (independent of the truth, so the MC-statistics noise of a prediction is realistic).
Each prediction is area-normalised to the fake data (shape only) and compared with the Pearson
chi-square ``sum((d - m)^2 / (m + var_m))`` where ``var_m`` is the MC-statistics variance of the
prediction.
``dchi2 = chi2(CV prediction) - chi2(DENT prediction)``; positive favours DENT. Its distribution
for truth CV and for truth DENT gives the separation d'.

The expected number of selected data events is the selected MC count scaled to a stated POT. No data
are read.

    python scripts/dent_improvement/tier_a_power.py [--pot 1e20 4.5e18]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plotting" / "dent_study"))
import cathode_scan as cs  # noqa: E402

XML = REPO / "xml" / "dent_improvement" / "sbnd_run1"
STUDIES = ["muon_score_inclusive", "muon_primary_score_inclusive", "proton_score_inclusive",
           "proton_primary_score_inclusive", "muon_calo_csda", "proton_calo_csda"]
FLAGS = ["reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
         "reco_cathode_crosser", "reco_leading_muon_end_x", "event_largest_flash_pe"]


def read_binning(study: str) -> tuple[str, np.ndarray]:
    """Fit variable (a branch name) and bin edges of a study, read from its S0 MC XML."""
    text = (XML / f"{study}_s0.xml").read_text()
    variable = re.search(r"<variable>\s*([A-Za-z0-9_]+)\s*</variable>", text).group(1)
    line = re.search(r'<bins unit="[^"]*"[^>]*/>', text).group(0)
    if "edges=" in line:
        edges = np.array(re.search(r'edges="([^"]*)"', line).group(1).split(), dtype=float)
    else:
        lo, hi, n = (float(re.search(f'{key}="([^"]*)"', line).group(1))
                     for key in ("min", "max", "nbins"))
        edges = np.linspace(lo, hi, int(n) + 1)
    return variable, edges


# variation label -> (reference tree and file, varied tree); the paired DENT samples and the other
# detector variations, which the detsys calculator compares with the production nominal sample.
NOMINAL = ("sbnd_large", cs.NOS)
PAIRS = {
    "dent": (("sbnd_detvar_newcv", cs.DET), "sbnd_detvar_newdent"),
    "olddent": (("sbnd_detvar_cv", cs.DET), "sbnd_detvar_dent"),
    "chargescale": (NOMINAL, "sbnd_detvar_charge_scale"),
    "smear15": (NOMINAL, "sbnd_detvar_smear15"),
    "sce0x": (NOMINAL, "sbnd_detvar_sce0x"),
    "sce2x": (NOMINAL, "sbnd_detvar_sce2x"),
    "wiremodxtheta": (NOMINAL, "sbnd_detvar_wiremod_xthetaxw"),
    "wiremodyz": (NOMINAL, "sbnd_detvar_wiremod_yz"),
    "mattmodyz": (NOMINAL, "sbnd_detvar_mattmod_yz"),
    "gainhigh": (NOMINAL, "sbnd_detvar_p5_gain"),
    "gainlow": (NOMINAL, "sbnd_detvar_m5_gain"),
}


def load(tree: str, branches: list[str], path: Path = cs.DET) -> tuple[pd.DataFrame, float]:
    """Selected-event columns of one sample and its POT."""
    f = uproot.open(path)
    frame = f[f"events/{tree}/selected"].arrays(branches + FLAGS, library="pd")
    return frame, float(f[f"events/{tree}/POT"].values()[0])


def cached_load(cache: dict, tree: str, path: Path, variable: str,
                extra: list[str] | None = None) -> tuple[pd.DataFrame, float]:
    """``load`` with a per-observable cache, so each sample is read once."""
    if (tree, str(path)) not in cache:
        cache[tree, str(path)] = load(tree, [variable, *(extra or [])], path)
    return cache[tree, str(path)]


NO_FLASH = {"wiremodxtheta", "wiremodyz"}     # samples with no flash information


def selection(frame: pd.DataFrame, selection_name: str, flash: bool = True) -> pd.Series:
    """The study's MC selection: S0, or S0 plus the crosser veto and the 4 cm end-point cut.

    The flash-threshold (trigger emulation) term is omitted for the WireMod pairs.
    """
    keep = ((frame.reco_proton_multiplicity >= 1) & (frame.reco_fiducialize_cathode == 1)
            & (frame.reco_veto_sbnd_highy_highz == 1))
    if flash:
        keep &= 0.647 * frame.event_largest_flash_pe > 2000
    if selection_name == "s2":
        keep &= (frame.reco_cathode_crosser == 0) & (frame.reco_leading_muon_end_x.abs() >= 4)
    return keep


def separation(counts: dict[str, np.ndarray], n_expected: float, rng: np.random.Generator,
               n_throws: int) -> tuple[float, float, float]:
    """Mean dchi2 for truth CV, for truth DENT, and the pooled sd, from half-split counts."""
    def chi2(data: np.ndarray, pred: np.ndarray) -> np.ndarray:
        total = data.sum(axis=1, keepdims=True)
        scale = total / pred.sum()
        model = scale * pred
        variance = model + scale ** 2 * pred
        return ((data - model) ** 2 / np.maximum(variance, 1e-9)).sum(axis=1)

    out = {}
    for truth in ("cv", "dent"):
        shape = counts[f"{truth}_A"] / counts[f"{truth}_A"].sum()
        data = rng.poisson(n_expected * shape, size=(n_throws, len(shape))).astype(float)
        out[truth] = chi2(data, counts["cv_B"]) - chi2(data, counts["dent_B"])
    pooled = np.sqrt(0.5 * (out["cv"].var(ddof=1) + out["dent"].var(ddof=1)))
    return out["cv"].mean(), out["dent"].mean(), pooled


def main() -> None:
    """Compute and print d' for every study, selection and assumed POT."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--pot", type=float, nargs="+", default=[1.0e20, 4.5e18])
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--throws", type=int, default=200)
    args = parser.parse_args()
    rng = np.random.default_rng(20261014)
    rows = []
    for study in STUDIES:
        variable, edges = read_binning(study)
        cache: dict[tuple[str, str], tuple[pd.DataFrame, float]] = {}
        for label, ((ref_tree, ref_path), var_tree) in PAIRS.items():
            ref, pot_ref = cached_load(cache, ref_tree, ref_path, variable)
            var, _ = cached_load(cache, var_tree, cs.DET, variable)
            for sel in ("s0", "s2"):
                flash = label not in NO_FLASH
                samples = {name: frame[selection(frame, sel, flash)][variable].to_numpy()
                           for name, frame in (("cv", ref), ("dent", var))}
                for pot in args.pot:
                    n_expected = len(samples["cv"]) * pot / pot_ref
                    deltas = []
                    for _ in range(args.repeats):
                        counts = {}
                        for name, values in samples.items():
                            in_a = rng.random(len(values)) < 0.5
                            for half, chosen in (("A", values[in_a]), ("B", values[~in_a])):
                                hist = np.histogram(chosen, edges)[0]
                                counts[f"{name}_{half}"] = hist.astype(float)
                        deltas.append(separation(counts, n_expected, rng, args.throws))
                    mean_cv, mean_dent, sd = np.mean(deltas, axis=0)
                    rows.append({"study": study, "variation": label, "selection": sel, "pot": pot,
                                 "expected_events": round(n_expected),
                                 "mean_dchi2_ref": mean_cv, "mean_dchi2_var": mean_dent, "sd": sd,
                                 "d_prime": (mean_dent - mean_cv) / sd})
        print(f"done {study}", flush=True)
    table = pd.DataFrame(rows)
    out = REPO / "docs" / "dent_improvement" / "tier_a_power_all.csv"
    table.to_csv(out, index=False)
    print(table.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
