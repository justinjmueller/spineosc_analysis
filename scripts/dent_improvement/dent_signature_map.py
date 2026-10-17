#!/usr/bin/env python3
"""Where does the DENT variation change the simulation? A position-resolved, simulation-only map.

The detvar CV and DENT samples contain the same events. Events are matched on the key (Run, Subrun,
Evt, event_max_neutrino_energy) and assigned to a region by the CV vertex: the half of the detector
(sign of the reco vertex x), a drift-distance bin ``|x|``, and a coarse (y, z) grid. For each region
the script reports the mean DENT-minus-CV difference, with its paired statistical error, of:

* the true and the reco muon end point (x, y, z), since DENT moves the true positions as well;
* the muon and proton calorimetric difference and the muon and proton PID and primary scores;

and the DENT/CV ratio of selected events and the fraction of events selected in only one sample.
No data are read.

    python scripts/dent_improvement/dent_signature_map.py [--selection s2]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "plotting" / "dent_study"))
import cathode_scan as cs  # noqa: E402

DRIFT_EDGES = [5, 50, 100, 150, 205]
Y_EDGES = [-200, -100, 0, 100, 200]
Z_EDGES = [0, 125, 250, 375, 500]
DIFFERENCES = [
    "true_leading_muon_end_x", "true_leading_muon_end_y", "true_leading_muon_end_z",
    "reco_leading_muon_end_x", "reco_leading_muon_end_y", "reco_leading_muon_end_z",
    "reco_leading_muon_calo_csda_diff", "reco_leading_proton_calo_csda_diff",
    "reco_leading_muon_muon_softmax", "reco_leading_proton_proton_softmax",
    "reco_leading_muon_primary_softmax", "reco_leading_proton_primary_softmax",
]
COLUMNS = [*cs.K, "reco_vertex_x", "reco_vertex_y", "reco_vertex_z", "reco_proton_multiplicity",
           "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser",
           "reco_leading_muon_end_x", "event_largest_flash_pe", *DIFFERENCES]


def common_events() -> pd.MultiIndex:
    """Keys of the events present in both samples, so that selection rates are comparable.

    The DENT sample has less POT than the CV sample; events absent from it are not selection churn.
    """
    return cs.keys(cs.DET, "sbnd_detvar_newcv").intersection(cs.keys(cs.DET, "sbnd_detvar_newdent"))


def load(tree: str, selection: str, common: pd.MultiIndex) -> pd.DataFrame:
    """Selected events of one sample under S0 or S2, indexed by event key."""
    frame = uproot.open(cs.DET)[f"events/{tree}/selected"].arrays(
        sorted(set(COLUMNS)), library="pd")
    frame = frame.drop_duplicates(cs.K, keep=False).set_index(cs.K)
    frame = frame[frame.index.isin(common)]
    keep = ((frame.reco_proton_multiplicity >= 1) & (frame.reco_fiducialize_cathode == 1)
            & (frame.reco_veto_sbnd_highy_highz == 1)
            & (0.647 * frame.event_largest_flash_pe > 2000))
    if selection == "s2":
        keep &= (frame.reco_cathode_crosser == 0) & (frame.reco_leading_muon_end_x.abs() >= 4)
    return frame[keep]


def region(frame: pd.DataFrame) -> pd.DataFrame:
    """Region labels from the vertex: half, drift-distance bin, y bin and z bin."""
    out = pd.DataFrame(index=frame.index)
    out["half"] = np.where(frame.reco_vertex_x < 0, "x<0", "x>0")
    out["drift"] = pd.cut(frame.reco_vertex_x.abs(), DRIFT_EDGES, labels=False)
    out["ybin"] = pd.cut(frame.reco_vertex_y, Y_EDGES, labels=False)
    out["zbin"] = pd.cut(frame.reco_vertex_z, Z_EDGES, labels=False)
    return out


def summarize(cv: pd.DataFrame, dent: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Per-region mean differences (DENT - CV) with paired errors, and selection counts."""
    both = cv.index.intersection(dent.index)
    labels = region(cv)
    rows = []
    for key, idx in labels.groupby(by).groups.items():
        key = key if isinstance(key, tuple) else (key,)
        in_both = idx.intersection(both)
        row = dict(zip(by, key, strict=True))
        row["n_cv"] = len(idx)
        row["n_dent_same_region"] = int(dent.index.isin(idx).sum())
        row["only_cv_frac"] = 1 - len(in_both) / max(len(idx), 1)
        for column in DIFFERENCES:
            delta = (dent.loc[in_both, column] - cv.loc[in_both, column]).dropna()
            row[f"d_{column}"] = delta.mean() if len(delta) else np.nan
            row[f"e_{column}"] = delta.std() / np.sqrt(len(delta)) if len(delta) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def uniformity(table: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """Chi-square of region-to-region constancy of each mean difference, per half.

    ``chi2 / ndf`` well above 1 means the effect varies with position by more than statistics.
    """
    rows = []
    for half, group in table.groupby("half"):
        for column in DIFFERENCES:
            value, error = group[f"d_{column}"], group[f"e_{column}"]
            ok = value.notna() & (error > 0)
            weight = 1 / error[ok] ** 2
            mean = (value[ok] * weight).sum() / weight.sum()
            chi2 = (((value[ok] - mean) / error[ok]) ** 2).sum()
            rows.append({"half": half, "binning": "/".join(by[1:]), "quantity": column,
                         "mean": mean, "chi2": chi2, "ndf": int(ok.sum()) - 1,
                         "chi2_per_ndf": chi2 / max(ok.sum() - 1, 1)})
    return pd.DataFrame(rows)


def main() -> None:
    """Build the maps and print the summaries."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--selection", default="s2", choices=["s0", "s2"])
    selection = parser.parse_args().selection
    common = common_events()
    cv = load("sbnd_detvar_newcv", selection, common)
    dent = load("sbnd_detvar_newdent", selection, common)
    print(f"selected ({selection}): CV {len(cv)}, DENT {len(dent)}, "
          f"matched in both {len(cv.index.intersection(dent.index))}")
    out = REPO / "docs" / "dent_improvement"
    by_drift = summarize(cv, dent, ["half", "drift"])
    by_yz = summarize(cv, dent, ["half", "ybin", "zbin"])
    by_drift.to_csv(out / f"signature_drift_{selection}.csv", index=False)
    by_yz.to_csv(out / f"signature_yz_{selection}.csv", index=False)
    show = ["half", "drift", "n_cv", "only_cv_frac", "d_true_leading_muon_end_x",
            "d_reco_leading_muon_calo_csda_diff", "d_reco_leading_muon_muon_softmax",
            "d_reco_leading_proton_calo_csda_diff"]
    pd.set_option("display.width", 200)
    print("\nBy drift distance and half (mean DENT - CV):")
    print(by_drift[show].round(4).to_string(index=False))
    columns = ["half", "ybin", "zbin", "n_cv", "d_true_leading_muon_end_x",
               "d_reco_leading_muon_calo_csda_diff", "d_reco_leading_muon_muon_softmax"]
    tests = pd.concat([uniformity(by_drift, ["half", "drift"]),
                       uniformity(by_yz, ["half", "ybin", "zbin"])])
    tests.to_csv(out / f"signature_uniformity_{selection}.csv", index=False)
    print("\nChi-square of constancy across regions (per half); chi2/ndf >> 1 = real variation:")
    print(tests.pivot_table(index="quantity", columns=["binning", "half"], values="chi2_per_ndf")
          .round(1).to_string())
    print(columns[:0] and "")


if __name__ == "__main__":
    main()
