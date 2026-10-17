"""
Per-event study of track breaking at the SBND cathode in the DENT variation.

The DENT detector variation shifts the leading-muon PID scores of cathode-crossing tracks. This
script asks how much of that shift comes from two different causes, using event-by-event matching
of the nominal (CV) and varied simulations:

* **break**: a true crosser that is reconstructed as a crosser (one stitched track) in the CV but
  as a non-crosser in the variation, i.e. the track is broken at the cathode and not stitched;
* **gap**: a true crosser that is a crosser in both samples (stitched) but whose start and end
  points are unchanged while the path length shrinks, i.e. a gap in an otherwise joined track.

Every other stitched track is "intact", and true non-crossers are a cathode-free reference. The
same decomposition is repeated for variations with no cathode-specific physics (the controls), so
that the effect of simply re-simulating the events can be separated from the effect of DENT.

The plan, class definitions and the rule used to call a cause "present" were written down before
the first run: ``docs/cathode_break/PLAN.md``.

Inputs
------
``spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root`` (signal trees of the detvar samples) and
``spineosc_numu_sbnd_v1.1.0_nosyst.root`` (``sbnd_large``). Unweighted MC; no data is used.

Outputs
-------
``<output>/transitions.csv``, ``classes.csv``, ``decomposition.csv``, ``tolerance_scan.csv``,
``summary.json`` and the figures ``crossing_transitions``, ``pid_decomposition``,
``gap_signature`` and ``tolerance_scan``.

Usage::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/dent_cathode_break.py \
        --output docs/cathode_break
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

from common import DATA_DIR, REPO, base_parser, save_figure, set_provenance

#: Directory holding the v1.1.0 ntuples.
DATA = DATA_DIR

#: Event key used to match the same generated event across samples.
KEY = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]

#: Nominal/varied sample pairs: label -> (file, nominal tree, varied tree, role).
#: Controls carry no cathode-specific physics; the pairing follows detsys/calculator.
PAIRS = {
    "DENT": ("detvar", "sbnd_detvar_newcv", "sbnd_detvar_newdent", "primary"),
    "DENT (replication)": ("detvar", "sbnd_detvar_cv", "sbnd_detvar_dent", "replication"),
    "ChargeScale": ("nosyst", "sbnd_large", "sbnd_detvar_charge_scale", "control"),
    "ChargeScale15": ("nosyst", "sbnd_large", "sbnd_detvar_smear15", "control"),
}

#: Which file the trees of each pair side live in (the controls mix the two files).
FILES = {
    "sbnd_large": "spineosc_numu_sbnd_v1.1.0_nosyst.root",
}
DETVAR_FILE = "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"

#: Columns read from the ``signal`` tree.
COLUMNS = KEY + [
    "true_cathode_crosser",
    "reco_cathode_crosser",
    "reco_fiducialize_cathode",
    "reco_veto_sbnd_highy_highz",
    "reco_leading_muon_length",
    "reco_leading_muon_muon_softmax",
    "reco_leading_muon_primary_softmax",
    "reco_leading_muon_mip_softmax",
    "reco_leading_muon_calo_csda_diff",
    "reco_leading_muon_start_x",
    "reco_leading_muon_start_y",
    "reco_leading_muon_start_z",
    "reco_leading_muon_end_x",
    "reco_leading_muon_end_y",
    "reco_leading_muon_end_z",
    "reco_leading_muon_start_dir_x",
    "reco_vertex_x",
]

#: Analysis end-point cut on the muon end position, in cm (the production SBND selection).
END_CUT = 4.0

#: Bin edges for the break-rate dependence studies.
VERTEX_EDGES = [0, 5, 10, 20, 40, 70, 120, 200]
DIRECTION_EDGES = [0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]

#: Quantities whose variation-minus-nominal shift is decomposed, with display labels.
QUANTITIES = {
    "reco_leading_muon_muon_softmax": "Muon softmax",
    "reco_leading_muon_primary_softmax": "Primary softmax",
    "reco_leading_muon_mip_softmax": "MIP softmax",
    "reco_leading_muon_length": "Length [cm]",
    "reco_leading_muon_calo_csda_diff": "Calo-CSDA diff",
}

#: Classes of true-crosser events, in plotting order, plus the reference class.
CLASSES = ["break", "gap", "intact", "gain", "non-crosser"]

#: Default tolerances (cm). Same end points to within TOL_POINT with a path-length loss larger than
#: TOL_LENGTH is called a gap. Both are scanned in :func:`tolerance_scan`.
TOL_POINT = 5.0
TOL_LENGTH = 3.0

#: Bootstrap resamples (Poisson weights) for every uncertainty.
N_BOOT = 500

#: Seed for the bootstrap.
SEED = 20261014

#: Tolerance grids for the stability scan (cm).
SCAN_POINT = [2.0, 3.0, 5.0, 8.0, 12.0]
SCAN_LENGTH = [1.0, 2.0, 3.0, 5.0, 8.0]


def load(tree: str, data: Path = DATA) -> pd.DataFrame:
    """
    Read the ``signal`` tree of one sample, indexed by event key.

    Events whose key is duplicated are dropped, as in the detsys calculator. Only events passing
    the cathode-fiducial and high-y/z vetoes are kept, matching the selection used elsewhere.

    Parameters
    ----------
    tree : str
        Sample name, e.g. ``"sbnd_detvar_newdent"``.
    data : Path, optional
        Directory containing the ntuples.

    Returns
    -------
    pd.DataFrame
        One row per event, indexed by :data:`KEY`.
    """
    path = data / FILES.get(tree, DETVAR_FILE)
    frame = uproot.open(path)[f"events/{tree}/signal"].arrays(COLUMNS, library="pd")
    frame = frame.drop_duplicates(KEY, keep=False).set_index(KEY)
    keep = (frame.reco_fiducialize_cathode == 1) & (frame.reco_veto_sbnd_highy_highz == 1)
    return frame[keep]


def selected_index(tree: str, data: Path = DATA) -> pd.MultiIndex:
    """
    Event keys present in the ``selected`` tree of one sample (duplicated keys dropped).

    Parameters
    ----------
    tree : str
        Sample name.
    data : Path, optional
        Directory containing the ntuples.

    Returns
    -------
    pd.MultiIndex
        Keys of selected events.
    """
    path = data / FILES.get(tree, DETVAR_FILE)
    frame = uproot.open(path)[f"events/{tree}/selected"].arrays(KEY, library="pd")
    return pd.MultiIndex.from_frame(frame.drop_duplicates(KEY, keep=False))


def pair(
    nominal: pd.DataFrame, varied: pd.DataFrame, require_reco: bool = True
) -> pd.DataFrame:
    """
    Join the two samples event by event.

    Parameters
    ----------
    nominal, varied : pd.DataFrame
        Output of :func:`load`.
    require_reco : bool, optional
        Keep only events with a reconstructed leading muon in both samples. The selection-fate
        study sets this to False so that events lost by one sample's reconstruction are visible.

    Returns
    -------
    pd.DataFrame
        Columns carry the suffixes ``_a`` (nominal) and ``_b`` (varied).
    """
    joined = nominal.join(varied, how="inner", lsuffix="_a", rsuffix="_b")
    if not require_reco:
        return joined
    return joined.dropna(subset=["reco_leading_muon_length_a", "reco_leading_muon_length_b"])


def analysis_selected(joined: pd.DataFrame, side: str) -> pd.Series:
    """
    Whether each event passes the analysis crosser veto and the end-point cut in one sample.

    Parameters
    ----------
    joined : pd.DataFrame
        Output of :func:`pair`.
    side : str
        ``"a"`` (nominal) or ``"b"`` (varied).

    Returns
    -------
    pd.Series
        Boolean per event; events without a reconstructed muon fail.
    """
    crosser = joined[f"reco_cathode_crosser_{side}"] == 1
    end = joined[f"reco_leading_muon_end_x_{side}"].abs()
    return (~crosser) & (end >= END_CUT)


def fate_table(
    joined: pd.DataFrame, sel_a: pd.MultiIndex, sel_b: pd.MultiIndex
) -> pd.DataFrame:
    """
    Selection fate of every class, for the ``selected`` tree and for the analysis cuts.

    For each class the events selected in the union of the two samples are split into selected in
    both, nominal only and variation only. The analysis version additionally requires
    ``reco_cathode_crosser == 0`` and ``|reco_leading_muon_end_x| >= END_CUT`` in each sample.
    Classes are defined on events reconstructed in both samples; the rest are ``"unreco"``.

    Parameters
    ----------
    joined : pd.DataFrame
        Output of :func:`pair` with ``require_reco=False``.
    sel_a, sel_b : pd.MultiIndex
        Keys of the selected events of the nominal and varied samples.

    Returns
    -------
    pd.DataFrame
        One row per (definition, class) with counts and shares of the class's union.
    """
    reco = joined[["reco_leading_muon_length_a", "reco_leading_muon_length_b"]].notna().all(axis=1)
    cls = pd.Series("unreco", index=joined.index, dtype=object)
    cls[reco] = classify(joined[reco])
    rows = []
    for definition in ("selected tree", "analysis cuts"):
        in_a = pd.Series(joined.index.isin(sel_a), index=joined.index)
        in_b = pd.Series(joined.index.isin(sel_b), index=joined.index)
        if definition == "analysis cuts":
            in_a &= analysis_selected(joined, "a")
            in_b &= analysis_selected(joined, "b")
        for name in [*CLASSES, "unreco"]:
            mask = cls == name
            both = int((in_a & in_b & mask).sum())
            only_a = int((in_a & ~in_b & mask).sum())
            only_b = int((~in_a & in_b & mask).sum())
            union = max(both + only_a + only_b, 1)
            rows.append({"definition": definition, "class": name, "n_union": both + only_a + only_b,
                         "both": both, "nominal_only": only_a, "variation_only": only_b,
                         "share_both": both / union, "share_nominal_only": only_a / union,
                         "share_variation_only": only_b / union})
    return pd.DataFrame(rows)


def break_rate(
    joined: pd.DataFrame, values: pd.Series, edges: list[float], weights: np.ndarray
) -> pd.DataFrame:
    """
    Break rate of nominal reco crossers that are true crossers, in bins of a variable.

    Parameters
    ----------
    joined : pd.DataFrame
        Output of :func:`pair`.
    values : pd.Series
        Variable to bin in (same index as ``joined``), taken from the nominal sample.
    edges : list[float]
        Bin edges.
    weights : np.ndarray
        Bootstrap weights, row 0 nominal.

    Returns
    -------
    pd.DataFrame
        Bin edges, break rate and bootstrap error, and the number of events per bin.
    """
    pop = ((joined.true_cathode_crosser_a == 1) & (joined.reco_cathode_crosser_a == 1)).to_numpy()
    broke = (joined.reco_cathode_crosser_b == 0).to_numpy()
    bins = np.digitize(values.to_numpy(), edges) - 1
    rows = []
    for k in range(len(edges) - 1):
        mask = pop & (bins == k)
        denom = weights[:, mask].sum(axis=1)
        num = weights[:, mask & broke].sum(axis=1)
        rate = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)
        rows.append({"lo": edges[k], "hi": edges[k + 1], "rate": float(rate[0]),
                     "error": float(rate[1:].std()), "n": int(mask.sum())})
    return pd.DataFrame(rows)


def classify(
    joined: pd.DataFrame, tol_point: float = TOL_POINT, tol_length: float = TOL_LENGTH
) -> pd.Series:
    """
    Assign every paired event to a class.

    Parameters
    ----------
    joined : pd.DataFrame
        Output of :func:`pair`.
    tol_point : float, optional
        Largest start/end-point movement (cm) still counted as "same end points".
    tol_length : float, optional
        Path-length loss (cm) above which a stitched track with unchanged end points is a gap.

    Returns
    -------
    pd.Series
        One of :data:`CLASSES` per event.
    """
    truth = joined.true_cathode_crosser_a == 1
    reco_a = joined.reco_cathode_crosser_a == 1
    reco_b = joined.reco_cathode_crosser_b == 1

    def move(end: str) -> np.ndarray:
        stem = f"reco_leading_muon_{end}"
        delta = [joined[f"{stem}_{axis}_b"] - joined[f"{stem}_{axis}_a"] for axis in "xyz"]
        return np.sqrt(sum(d**2 for d in delta))

    same_points = (move("start") < tol_point) & (move("end") < tol_point)
    lost = (joined.reco_leading_muon_length_a - joined.reco_leading_muon_length_b) > tol_length

    cls = pd.Series("non-crosser", index=joined.index, dtype=object)
    cls[truth & reco_a & ~reco_b] = "break"
    cls[truth & reco_a & reco_b] = "intact"
    cls[truth & reco_a & reco_b & same_points & lost] = "gap"
    cls[truth & ~reco_a & reco_b] = "gain"
    # true crossers reconstructed as non-crossers in both samples carry no cathode information
    cls[truth & ~reco_a & ~reco_b] = "non-crosser"
    return cls


def poisson_weights(n: int, n_boot: int = N_BOOT, seed: int = SEED) -> np.ndarray:
    """Poisson(1) bootstrap weights of shape ``(n_boot + 1, n)``; row 0 is the nominal sample."""
    rng = np.random.default_rng(seed)
    weights = rng.poisson(1.0, size=(n_boot + 1, n)).astype(np.float32)
    weights[0] = 1.0
    return weights


def contributions(
    delta: np.ndarray, cls: np.ndarray, weights: np.ndarray
) -> dict[str, tuple[float, float, float, float]]:
    """
    Class contributions to the total mean shift of one quantity.

    The mean shift of the whole population is split as ``sum_c f_c * m_c`` with ``f_c`` the class
    fraction and ``m_c`` the mean shift inside the class, so the contributions sum to the total.

    Parameters
    ----------
    delta : np.ndarray
        Per-event variation-minus-nominal shift.
    cls : np.ndarray
        Class label per event.
    weights : np.ndarray
        Bootstrap weights, row 0 nominal.

    Returns
    -------
    dict
        class (and ``"total"``) -> (fraction, mean shift in class, contribution, contribution
        bootstrap error).
    """
    total = weights.sum(axis=1)
    out = {}
    for name in [*CLASSES, "total"]:
        mask = np.ones(len(delta), dtype=bool) if name == "total" else cls == name
        w = weights[:, mask]
        count = w.sum(axis=1)
        summed = w @ np.nan_to_num(delta[mask])
        contrib = summed / total
        mean = np.divide(summed, count, out=np.zeros_like(summed), where=count > 0)
        out[name] = (
            float(count[0] / total[0]),
            float(mean[0]),
            float(contrib[0]),
            float(contrib[1:].std()),
        )
    return out


def decompose(joined: pd.DataFrame, cls: pd.Series, weights: np.ndarray) -> pd.DataFrame:
    """Class contributions to the shift of every quantity in :data:`QUANTITIES`."""
    rows = []
    for column, label in QUANTITIES.items():
        delta = (joined[f"{column}_b"] - joined[f"{column}_a"]).to_numpy()
        for name, (frac, mean, contrib, err) in contributions(
            delta, cls.to_numpy(), weights
        ).items():
            rows.append(
                {"quantity": label, "class": name, "fraction": frac, "mean_shift": mean,
                 "contribution": contrib, "contribution_err": err}
            )
    return pd.DataFrame(rows)


def transitions(joined: pd.DataFrame, weights: np.ndarray) -> pd.DataFrame:
    """Reco crossing-flag transition rates among true crossers, with bootstrap errors."""
    sel = (joined.true_cathode_crosser_a == 1).to_numpy()
    a = (joined.reco_cathode_crosser_a == 1).to_numpy()[sel]
    b = (joined.reco_cathode_crosser_b == 1).to_numpy()[sel]
    w = weights[:, sel]
    rows = []
    for name, mask in (("X -> X", a & b), ("X -> NX (break)", a & ~b),
                       ("NX -> X", ~a & b), ("NX -> NX", ~a & ~b)):
        frac = (w @ mask) / w.sum(axis=1)
        rows.append({"transition": name, "fraction": float(frac[0]),
                     "error": float(frac[1:].std()), "n": int(mask.sum())})
    return pd.DataFrame(rows)


def tolerance_scan(joined: pd.DataFrame) -> pd.DataFrame:
    """Gap fraction of crossing tracks for every (point, length) tolerance pair."""
    rows = []
    for tp in SCAN_POINT:
        for tl in SCAN_LENGTH:
            cls = classify(joined, tp, tl)
            stitched = cls.isin(["intact", "gap"])
            rows.append({"tol_point": tp, "tol_length": tl,
                         "gap_fraction": float((cls == "gap").sum() / max(stitched.sum(), 1))})
    return pd.DataFrame(rows)


def plot_decomposition(table: pd.DataFrame, output: Path, fmt: str) -> None:
    """Stacked class contributions to the PID-score shifts, one panel per sample pair."""
    pairs = list(table.pair.unique())
    quantities = ["Muon softmax", "Primary softmax", "MIP softmax"]
    colors = {"break": "#D55E00", "gap": "#E69F00", "intact": "#56B4E9",
              "gain": "#009E73", "non-crosser": "0.7"}
    fig, axes = plt.subplots(1, len(quantities), figsize=(13, 4.6), sharey=False)
    for ax, quantity in zip(axes, quantities, strict=True):
        x = np.arange(len(pairs))
        pos = np.zeros(len(pairs))
        neg = np.zeros(len(pairs))
        for name in CLASSES:
            vals = np.array([
                table[(table.pair == p) & (table.quantity == quantity) & (table["class"] == name)]
                .contribution.iloc[0] for p in pairs
            ])
            bottom = np.where(vals >= 0, pos, neg)
            ax.bar(x, vals, bottom=bottom, color=colors[name], label=name, edgecolor="k", lw=0.5)
            pos += np.where(vals >= 0, vals, 0)
            neg += np.where(vals < 0, vals, 0)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(pairs, rotation=20, ha="right", fontsize=9)
        ax.set_title(quantity)
        ax.set_ylabel("Contribution to mean shift (var - nom)")
    axes[0].legend(fontsize=9, frameon=False)
    fig.tight_layout()
    save_figure(fig, output / f"pid_decomposition.{fmt}")


def plot_transitions(table: pd.DataFrame, output: Path, fmt: str) -> None:
    """Break and gain rates of true crossers for every sample pair."""
    pairs = list(table.pair.unique())
    fig, ax = plt.subplots(figsize=(7, 4.4))
    x = np.arange(len(pairs))
    for k, (name, color) in enumerate((("X -> NX (break)", "#D55E00"), ("NX -> X", "#009E73"))):
        sub = [table[(table.pair == p) & (table.transition == name)].iloc[0] for p in pairs]
        ax.bar(x + (k - 0.5) * 0.35, [100 * s.fraction for s in sub], 0.35,
               yerr=[100 * s.error for s in sub], color=color, label=name, edgecolor="k", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("True crossers [%]")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, output / f"crossing_transitions.{fmt}")


def plot_gap_signature(joins: dict[str, pd.DataFrame], output: Path, fmt: str) -> None:
    """Path-length change versus end-point movement for stitched tracks, DENT against a control."""
    show = [name for name in ("DENT", "ChargeScale") if name in joins]
    fig, axes = plt.subplots(1, len(show), figsize=(6 * len(show), 4.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, name in zip(axes, show, strict=True):
        j = joins[name]
        stitched = (
            (j.true_cathode_crosser_a == 1)
            & (j.reco_cathode_crosser_a == 1)
            & (j.reco_cathode_crosser_b == 1)
        )
        j = j[stitched]
        dend = np.sqrt(sum(
            (j[f"reco_leading_muon_end_{ax_}_b"] - j[f"reco_leading_muon_end_{ax_}_a"]) ** 2
            for ax_ in "xyz"
        ))
        dlen = j.reco_leading_muon_length_b - j.reco_leading_muon_length_a
        hist = ax.hist2d(np.clip(dend, 0, 30), np.clip(dlen, -40, 40), bins=[60, 80],
                         range=[[0, 30], [-40, 40]], cmap="viridis", norm="log")
        ax.axhline(-TOL_LENGTH, color="w", ls="--", lw=0.8)
        ax.axvline(TOL_POINT, color="w", ls="--", lw=0.8)
        ax.set_xlabel("End-point movement [cm]")
        ax.set_title(f"{name}: stitched true crossers")
        fig.colorbar(hist[3], ax=ax)
    axes[0].set_ylabel("Length change, var - nom [cm]")
    fig.tight_layout()
    save_figure(fig, output / f"gap_signature.{fmt}")


def plot_scan(scans: dict[str, pd.DataFrame], output: Path, fmt: str) -> None:
    """Gap fraction against the length tolerance, one curve per point tolerance."""
    fig, axes = plt.subplots(1, len(scans), figsize=(5.5 * len(scans), 4.2), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (name, scan) in zip(axes, scans.items(), strict=True):
        for tp, sub in scan.groupby("tol_point"):
            ax.plot(sub.tol_length, 100 * sub.gap_fraction, marker="o", label=f"{tp:g} cm")
        ax.set_xlabel("Length-loss tolerance [cm]")
        ax.set_title(name)
        ax.legend(title="End-point tolerance", fontsize=8, frameon=False)
    axes[0].set_ylabel("Gap fraction of stitched tracks [%]")
    fig.tight_layout()
    save_figure(fig, output / f"tolerance_scan.{fmt}")


def plot_fate(fates: pd.DataFrame, output: Path, fmt: str) -> None:
    """Selection fate of each class for DENT, under the selected tree and the analysis cuts."""
    sub = fates[fates.pair == "DENT"]
    colors = {"share_both": "0.7", "share_nominal_only": "#D55E00",
              "share_variation_only": "#0072B2"}
    labels = {"share_both": "selected in both", "share_nominal_only": "nominal only",
              "share_variation_only": "DENT only"}
    classes = [*CLASSES, "unreco"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    for ax, definition in zip(axes, ("selected tree", "analysis cuts"), strict=True):
        d = sub[sub.definition == definition].set_index("class").loc[classes]
        bottom = np.zeros(len(classes))
        for key, color in colors.items():
            ax.bar(classes, 100 * d[key].to_numpy(), bottom=bottom, color=color,
                   label=labels[key], edgecolor="k", lw=0.5)
            bottom += 100 * d[key].to_numpy()
        ax.set_title(definition)
        ax.tick_params(axis="x", rotation=20)
    axes[0].set_ylabel("Share of events selected in either sample [%]")
    axes[0].legend(fontsize=9, frameon=False, loc="lower left")
    fig.tight_layout()
    save_figure(fig, output / f"selection_fate.{fmt}")


def plot_dependence(rates: pd.DataFrame, output: Path, fmt: str) -> None:
    """Break rate against vertex distance to the cathode and muon direction."""
    titles = {"vertex_x": "|reco vertex x| [cm] (cathode at 0)",
              "direction_x": "|start direction x| (1 = along the drift)"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, (variable, title) in zip(axes, titles.items(), strict=True):
        for label in rates.pair.unique():
            d = rates[(rates.pair == label) & (rates.variable == variable)]
            centers = 0.5 * (d.lo + d.hi)
            ax.errorbar(centers, 100 * d.rate, yerr=100 * d.error, marker="o", label=label,
                        capsize=2)
        ax.set_xlabel(title)
        ax.set_ylabel("Break rate of true crossers [%]")
        ax.set_yscale("symlog", linthresh=0.1)
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    save_figure(fig, output / f"break_rate_dependence.{fmt}")


def run(pairs: list[str], data: Path, output: Path, fmt: str, tol_point: float,
        tol_length: float) -> dict:
    """
    Run the study for the chosen sample pairs and write tables, figures and a JSON summary.

    Returns
    -------
    dict
        The summary written to ``summary.json``.
    """
    output.mkdir(parents=True, exist_ok=True)
    trans, deco, joins, scans, fates, rates = [], [], {}, {}, [], []
    summary: dict = {"tol_point": tol_point, "tol_length": tol_length, "pairs": {}}
    for label in pairs:
        _, nom, var, role = PAIRS[label]
        print(f"{label}: {nom} -> {var}", flush=True)
        frame_a, frame_b = load(nom, data), load(var, data)
        joined = pair(frame_a, frame_b)
        weights = poisson_weights(len(joined))
        fate = fate_table(pair(frame_a, frame_b, require_reco=False),
                          selected_index(nom, data), selected_index(var, data))
        fates.append(fate.assign(pair=label))
        for name, values, edges in (
            ("vertex_x", joined.reco_vertex_x_a.abs(), VERTEX_EDGES),
            ("direction_x", joined.reco_leading_muon_start_dir_x_a.abs(), DIRECTION_EDGES),
        ):
            rates.append(break_rate(joined, values, edges, weights).assign(pair=label,
                                                                           variable=name))
        cls = classify(joined, tol_point, tol_length)
        t = transitions(joined, weights)
        d = decompose(joined, cls, weights)
        t.insert(0, "pair", label)
        d.insert(0, "pair", label)
        trans.append(t)
        deco.append(d)
        joins[label] = joined
        scans[label] = tolerance_scan(joined)
        summary["pairs"][label] = {
            "role": role,
            "n_matched": int(len(joined)),
            "class_counts": cls.value_counts().to_dict(),
        }
    trans = pd.concat(trans)
    deco = pd.concat(deco)
    trans.to_csv(output / "transitions.csv", index=False)
    deco.to_csv(output / "decomposition.csv", index=False)
    pd.concat([s.assign(pair=k) for k, s in scans.items()]).to_csv(
        output / "tolerance_scan.csv", index=False)
    fates = pd.concat(fates)
    rates = pd.concat(rates)
    fates.to_csv(output / "selection_fate.csv", index=False)
    rates.to_csv(output / "break_rate_dependence.csv", index=False)
    plot_decomposition(deco, output, fmt)
    plot_fate(fates, output, fmt)
    plot_dependence(rates, output, fmt)
    plot_transitions(trans, output, fmt)
    plot_gap_signature(joins, output, fmt)
    plot_scan({k: v for k, v in scans.items() if k in ("DENT", "ChargeScale")}, output, fmt)
    (output / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(trans.to_string(index=False))
    print(deco[deco.quantity == "Muon softmax"].to_string(index=False))
    print(fates[fates.pair == "DENT"].to_string(index=False))
    print(rates[rates.pair == "DENT"].to_string(index=False))
    return summary


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    parser = base_parser(__doc__.split("\n\n")[0])
    parser.set_defaults(output=REPO / "docs" / "cathode_break")
    parser.add_argument("--pair", action="append", choices=list(PAIRS),
                        help="sample pair to run (repeatable; default all)")
    parser.add_argument("--data", type=Path, default=DATA, help="directory with the ntuples")
    parser.add_argument("--tol-point", type=float, default=TOL_POINT,
                        help="end-point tolerance [cm]")
    parser.add_argument("--tol-length", type=float, default=TOL_LENGTH,
                        help="path-length loss defining a gap [cm]")
    args: argparse.Namespace = parser.parse_args(argv)
    set_provenance(__file__, [args.data / DETVAR_FILE], {"pairs": args.pair or list(PAIRS)})
    run(args.pair or list(PAIRS), args.data, args.output, args.format, args.tol_point,
        args.tol_length)


if __name__ == "__main__":
    main()
