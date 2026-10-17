#!/usr/bin/env python3
"""
Which selection criteria cause events to enter or leave the selection in a detector variation?

For a detector variation (default: ``ChargeScale``) paired event by event with the nominal sample,
this script finds the events that are selected in only one of the two samples and attributes each to
the reconstruction criteria that fail in the sample where it is *not* selected. The per-criterion
flags come from the ``signal`` tree, so the attribution covers true-signal events. Selected events
that are not true signal have no per-criterion flags; they are broken down by true category instead,
which says what kind of background they are but not which cut removed them.

Criteria (all must pass for an event to be selected under the analysis selection S2): fiducial
volume, cathode fiducial, high-y/z veto, containment, single muon, at least one proton, no charged
pions, no electrons, no photons, not a cathode crosser, and a muon end point at least 4 cm from the
cathode.

Usage
-----
    python plotting/studies/detvar_churn.py [--variation sbnd_detvar_charge_scale]
                                            [--min-energy 1.5]

Outputs go to ``docs/detvar_churn/``: a CSV of counts per criterion and energy bin, and a figure.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, REPO, configure_style, save_figure, set_provenance  # noqa: E402

DATA = DATA_DIR
NOMINAL_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_nosyst.root"
DETVAR_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"
KEY = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
END_CUT = 4.0
EDGES = [0.3, 0.7, 1.0, 1.3, 1.5, 2.0, 3.0]

# criterion name -> (branch, test) ; the test returns True where the criterion is PASSED
CRITERIA = {
    "fiducial volume": ("reco_fiducial_cut_osc", lambda x: x == 1),
    "cathode fiducial": ("reco_fiducialize_cathode", lambda x: x == 1),
    "high-y/z veto": ("reco_veto_sbnd_highy_highz", lambda x: x == 1),
    "containment": ("reco_user_containment_cut", lambda x: x == 1),
    "single muon": ("reco_single_muon", lambda x: x == 1),
    "at least one proton": ("reco_proton_multiplicity", lambda x: x >= 1),
    "no charged pions": ("reco_no_charged_pions", lambda x: x == 1),
    "no electrons": ("reco_no_electrons", lambda x: x == 1),
    "no photons": ("reco_no_photons", lambda x: x == 1),
    "not a crosser": ("reco_cathode_crosser", lambda x: x == 0),
    "end point >= 4 cm": ("reco_leading_muon_end_x", lambda x: np.abs(x) >= END_CUT),
}
# true event categories of the ntuples (truth-level definitions from the category configuration)
CATEGORIES = {
    0: "contained, 1 muon, 1 proton (signal)",
    1: "non-contained, 1 muon, 1 proton",
    2: "contained, 1 muon, N>1 protons (signal)",
    3: "non-contained, 1 muon, N>1 protons",
    4: "contained, 1 muon, any protons, >0 charged pions",
    5: "non-contained, 1 muon, any protons, >0 charged pions",
    6: "fiducial CC, non-signal",
    7: "non-fiducial CC",
    8: "NC",
    9: "cosmic",
}
SELECTED_COLUMNS = [*KEY, "reco_visible_energy", "true_category", "reco_fiducialize_cathode",
                    "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x"]


def load_selected(path: Path, tree: str) -> pd.DataFrame:
    """Selected events under S2, indexed by event key, duplicate keys dropped."""
    frame = uproot.open(path)[f"events/{tree}/selected"].arrays(SELECTED_COLUMNS, library="pd")
    frame = frame.drop_duplicates(KEY, keep=False).set_index(KEY)
    keep = ((frame.reco_fiducialize_cathode == 1) & (frame.reco_veto_sbnd_highy_highz == 1)
            & (frame.reco_cathode_crosser == 0) & (frame.reco_leading_muon_end_x.abs() >= END_CUT))
    return frame[keep]


def load_flags(path: Path, tree: str) -> pd.DataFrame:
    """Per-criterion inputs for every true-signal event, indexed by event key."""
    branches = sorted({branch for branch, _ in CRITERIA.values()})
    frame = uproot.open(path)[f"events/{tree}/signal"].arrays([*KEY, *branches], library="pd")
    return frame.drop_duplicates(KEY, keep=False).set_index(KEY)


def attribute(events: pd.Index, flags: pd.DataFrame) -> pd.DataFrame:
    """Boolean table (events x criteria): True where the event FAILS the criterion.

    Events absent from the signal tree are marked in the ``not signal`` column instead.
    """
    present = events.isin(flags.index)
    failed = pd.DataFrame(False, index=events, columns=[*CRITERIA, "not signal"])
    failed["not signal"] = ~present
    sub = flags.reindex(events[present])
    for name, (branch, passes) in CRITERIA.items():
        failed.loc[sub.index, name] = ~passes(sub[branch].to_numpy())
    return failed


def check_flags_reproduce_selection(selected: pd.DataFrame, flags: pd.DataFrame) -> float:
    """Fraction of selected true-signal events for which every criterion passes."""
    both = selected.index[selected.index.isin(flags.index)]
    return float((~attribute(both, flags).drop(columns="not signal").any(axis=1)).mean())


def breakdown(lost: pd.Index, flags: pd.DataFrame) -> pd.DataFrame:
    """Counts per criterion: failing at all, and failing as the only criterion."""
    failed = attribute(lost, flags)
    criteria = failed.drop(columns="not signal")
    n_fail = criteria.sum(axis=1)
    rows = [{"criterion": c, "fails": int(criteria[c].sum()),
             "only_this": int((criteria[c] & (n_fail == 1)).sum())} for c in CRITERIA]
    rows.append({"criterion": "not signal (unattributed)", "fails": int(failed["not signal"].sum()),
                 "only_this": int(failed["not signal"].sum())})
    rows.append({"criterion": "signal, no failing criterion", "fails": int(((n_fail == 0)
                 & ~failed["not signal"]).sum()), "only_this": int(((n_fail == 0)
                 & ~failed["not signal"]).sum())})
    return pd.DataFrame(rows)


def category_table(nominal: pd.DataFrame, varied: pd.DataFrame, flags_nom: pd.DataFrame,
                   flags_var: pd.DataFrame, min_energy: float) -> pd.DataFrame:
    """True-category counts of the unattributable lost and gained events above ``min_energy``."""
    high = nominal.index[nominal.reco_visible_energy >= min_energy]
    lost = high.intersection(nominal.index.difference(varied.index))
    lost = lost[~lost.isin(flags_var.index)]
    high_var = varied.index[varied.reco_visible_energy >= min_energy]
    gained = high_var.intersection(varied.index.difference(nominal.index))
    gained = gained[~gained.isin(flags_nom.index)]
    table = pd.DataFrame({
        "selected_nominal": nominal.loc[high].true_category.value_counts(),
        "lost": nominal.loc[lost].true_category.value_counts(),
        "gained": varied.loc[gained].true_category.value_counts(),
    }).fillna(0).astype(int).sort_index()
    table["lost_fraction"] = (table.lost / table.selected_nominal).round(3)
    table.insert(0, "meaning", [CATEGORIES.get(int(c), "?") for c in table.index])
    return table


def run(variation: str, min_energy: float, output: Path) -> None:
    """Compute the churn breakdown for one variation and write the CSV and figure."""
    nominal = load_selected(NOMINAL_FILE, "sbnd_large")
    varied = load_selected(DETVAR_FILE, variation)
    flags_nom = load_flags(NOMINAL_FILE, "sbnd_large")
    flags_var = load_flags(DETVAR_FILE, variation)
    print(f"selected (S2): nominal {len(nominal)}, {variation} {len(varied)}")
    print("flags reproduce the selection for selected signal events: "
          f"nominal {check_flags_reproduce_selection(nominal, flags_nom):.4f}, "
          f"variation {check_flags_reproduce_selection(varied, flags_var):.4f}")

    energy_nom = nominal.reco_visible_energy
    lost_all = nominal.index.difference(varied.index)       # selected in nominal only
    gained_all = varied.index.difference(nominal.index)     # selected in the variation only
    tables = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:], strict=True):
        in_bin_nom = nominal.index[(energy_nom >= lo) & (energy_nom < hi)]
        lost = lost_all.intersection(in_bin_nom)
        energy_var = varied.reco_visible_energy
        in_bin_var = varied.index[(energy_var >= lo) & (energy_var < hi)]
        gained = gained_all.intersection(in_bin_var)
        for direction, events, flags in (("lost", lost, flags_var), ("gained", gained, flags_nom)):
            table = breakdown(events, flags)
            table.insert(0, "direction", direction)
            table.insert(1, "bin", f"{lo}-{hi}")
            table["n_events"] = len(events)
            tables.append(table)
        print(f"E in [{lo}, {hi}): nominal {len(in_bin_nom):6d}, lost {len(lost):5d}, "
              f"gained {len(gained):5d}")
    result = pd.concat(tables, ignore_index=True)
    output.mkdir(parents=True, exist_ok=True)
    result.to_csv(output / f"churn_{variation}.csv", index=False)

    categories = category_table(nominal, varied, flags_nom, flags_var, min_energy)
    categories.to_csv(output / f"churn_categories_{variation}.csv")
    print(f"\nE >= {min_energy} GeV, events with no per-criterion flags (not true signal), "
          "by true category and the fraction of that category's selected events involved:")
    print(categories.to_string())
    window = result[result.bin.map(lambda b: float(b.split("-")[0]) >= min_energy)]
    window = window.groupby(["direction", "criterion"], sort=False)[["fails", "only_this"]].sum()
    print(f"\nE >= {min_energy} GeV: events not selected in one sample, by failing criterion")
    print(window.to_string())
    plot(result, variation, output)


def plot(result: pd.DataFrame, variation: str, output: Path) -> None:
    """Failing-criterion fractions of lost events, per energy bin."""
    import matplotlib.pyplot as plt

    configure_style()
    lost = result[(result.direction == "lost") & ~result.criterion.str.startswith("signal,")]
    bins = list(dict.fromkeys(lost.bin))
    criteria = [c for c in lost.criterion.unique() if lost[lost.criterion == c].fails.sum() > 0]
    fig, ax = plt.subplots(figsize=(9, 5))
    bottom = np.zeros(len(bins))
    for name in criteria:
        values = np.array([lost[(lost.bin == b) & (lost.criterion == name)].only_this.sum()
                           for b in bins])
        ax.bar(bins, values, bottom=bottom, label=name)
        bottom += values
    ax.set_xlabel("Nominal reconstructed visible energy [GeV]")
    ax.set_ylabel("Events selected in nominal only (sole failing criterion)")
    ax.legend(fontsize=8, ncol=2)
    save_figure(fig, output / f"churn_{variation}.pdf")


def main() -> None:
    """Parse arguments and run."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--variation", default="sbnd_detvar_charge_scale")
    parser.add_argument("--min-energy", type=float, default=1.5)
    parser.add_argument("--output", type=Path, default=REPO / "docs" / "detvar_churn")
    args = parser.parse_args()
    set_provenance(__file__, [NOMINAL_FILE, DETVAR_FILE],
                   {"variation": args.variation, "min_energy": args.min_energy})
    run(args.variation, args.min_energy, args.output)


if __name__ == "__main__":
    main()
