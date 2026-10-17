"""
The SBND APA-gap deficit in the detector data access: where it is, how large, and what it correlates with.

Compares the SBND release detector-access data with the production simulation (``cvrw`` with ``cvwgt``), both with the
detector-access selection of ``xml/detector_data_access/sbnd_run1/*_inclusive.xml`` (at least one proton, cathode
fiducialisation, high-y/z veto, no cathode crosser, muon end point at least 4 cm from the cathode, and the flash-PE trigger
emulation on the simulation). **Area-normalised comparisons only, with no sense of the size of the data**: the simulation is
scaled to the data of the selection (the only normalisation), every histogram is shown and stored as a fraction of the
selected events, and no event count, total, rate or file size of the data is printed, written to a figure or stored in a
sidecar (the data file is recorded by hash only). Ratios, errors and pulls are what carry the comparison.

Only the variables the detector data access allows are read from the data file: the vertex, the muon and proton start and
end points, their calorimetric (E_calo - E_range) / E_range, their PID and primary softmax scores and the selection flags
(plus Run, which labels the data-taking period). No energy, no length, no polar angle and no direction components; the
allowed-branch list is enforced when the data are loaded. Two-dimensional views combine two allowed coordinates, never a
difference of them.

Outputs (``--output``, default ``docs/apa_gap``)::

    gap_vertex_z.<fmt>        data/MC in 0.5 cm bins around the gap, with the 4 cm PROfit systematic band for scale
    gap_coordinates.<fmt>     data/MC ratio of the vertex, muon and proton start and end z in 2 cm bins
    gap_proton_end_vs_vertex.<fmt>   data/MC ratio in (vertex z, proton end z)
    gap_splits.<fmt>          ratio inside the gap window against x, y, run and the proton score
    gap_summary.json          the ratios, pulls and fractions quoted in docs/apa_gap/REPORT.md (no counts)

Usage::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/apagap.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Release detector-access data and production simulation (``<DATA_DIR>/...``).
DATA_FILE = "spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root"
MC_FILE = "spineosc_numu_sbnd_v1.1.0_wsyst_cvrw.root"
DATA_TREE = "events/sbnd_bnblight_full/selected"
MC_TREE = "events/sbnd_large/selected"

#: PROfit detector-access output for the 4 cm vertex-z comparison (statistics + systematic band).
PROFIT_FILE = common.RUNS / "da_v2" / "da2_sbnd_run1_vertexz_apagap_inclusive_v1_PROplot_histograms.root"

#: What the detector data access allows (technical note, "Detector data access").
_ALLOWED_EXACT = {
    "reco_vertex_x", "reco_vertex_y", "reco_vertex_z", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
    "reco_cathode_crosser", "reco_proton_multiplicity", "Run", "Subrun", "Evt",
}
_ALLOWED_PATTERNS = [
    r"reco_leading_(muon|proton)_(start|end)_[xyz]",
    r"reco_leading_(muon|proton)_calo_csda_diff",
    r"reco_leading_(muon|proton)_(muon|pion|proton|mip|electron|photon|secondary|hadron)_softmax",
    r"reco_leading_(muon|proton)_primary_softmax",
]

#: Data branches read (all allowed).
POSITIONS = [f"reco_leading_{p}_{k}_{c}" for p in ("muon", "proton") for k in ("start", "end") for c in "xyz"]
VERTEX = [f"reco_vertex_{c}" for c in "xyz"]
SCORES = [
    "reco_leading_muon_muon_softmax", "reco_leading_muon_primary_softmax", "reco_leading_proton_proton_softmax",
    "reco_leading_proton_primary_softmax",
]
FLAGS = ["reco_proton_multiplicity", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser"]
DATA_BRANCHES = VERTEX + POSITIONS + SCORES + FLAGS + ["Run"]

#: Simulation-only branches (truth categories, trigger emulation, weight).
MC_EXTRA = ["true_cc", "true_pdg", "true_neutrino_id", "true_category", "true_detector", "event_largest_flash_pe", "cvwgt"]

#: Gap window in vertex z [cm] used for the split studies and the headline numbers.
WINDOW = (251.5, 254.5)

#: Smallest share of the selected events a (vertex z, proton end z) cell must hold to be shown.
MIN_CELL_SHARE = 2.0e-4

#: Colours of the data/MC curves.
DATA_COLOR, MC_COLOR = "k", "C0"


def check_allowed(branches: list[str]) -> None:
    """Raise if a branch is outside the detector-data-access whitelist."""
    for name in branches:
        if name not in _ALLOWED_EXACT and not any(re.fullmatch(p, name) for p in _ALLOWED_PATTERNS):
            raise ValueError(f"{name} is not allowed in the detector data access sample")


def selection(df: pd.DataFrame) -> pd.Series:
    """The detector-access (inclusive) selection."""
    return (
        (df.reco_proton_multiplicity >= 1) & (df.reco_fiducialize_cathode == 1)
        & (df.reco_veto_sbnd_highy_highz == 1) & (df.reco_cathode_crosser == 0)
        & (df.reco_leading_muon_end_x.abs() >= 4)
    )


def load_data(data_dir: Path = common.DATA_DIR) -> pd.DataFrame:
    """Selected data events, allowed branches only."""
    import uproot

    check_allowed(DATA_BRANCHES)
    with uproot.open(data_dir / DATA_FILE) as handle:
        frame = handle[DATA_TREE].arrays(DATA_BRANCHES, library="pd")
    return frame[selection(frame)].reset_index(drop=True)


def load_mc(data_dir: Path = common.DATA_DIR) -> pd.DataFrame:
    """Selected simulation events (neutrino categories, flash threshold) with weight ``w`` = cvwgt."""
    import uproot

    branches = [b for b in DATA_BRANCHES if b != "Run"] + MC_EXTRA
    with uproot.open(data_dir / MC_FILE) as handle:
        frame = handle[MC_TREE].arrays(branches, library="pd")
    keep = selection(frame) & (0.647 * frame.event_largest_flash_pe > 2000) & (frame.true_neutrino_id >= 0)
    frame = frame[keep].reset_index(drop=True)
    frame["w"] = frame.cvwgt
    return frame


def normalise(data: pd.DataFrame, mc: pd.DataFrame) -> pd.DataFrame:
    """Scale the simulation to the number of data events (the only normalisation used)."""
    mc = mc.copy()
    mc["w"] = mc.cvwgt * len(data) / mc.cvwgt.sum()
    return mc


def binned(data_values, mc_values, mc_weights, edges):
    """Data counts, MC prediction and the MC statistical variance in the given bins."""
    nd, _ = np.histogram(data_values, edges)
    nm, _ = np.histogram(mc_values, edges, weights=mc_weights)
    var, _ = np.histogram(mc_values, edges, weights=mc_weights**2)
    return nd, nm, var


def pull(nd, nm, var):
    """(data - MC) / sqrt(MC + MC statistical variance), the Poisson-plus-MC-statistics pull."""
    return (nd - nm) / np.sqrt(np.maximum(nm, 1e-9) + var)


def in_window(frame: pd.DataFrame, window=WINDOW) -> pd.Series:
    """Events whose vertex z lies in the gap window."""
    return (frame.reco_vertex_z >= window[0]) & (frame.reco_vertex_z < window[1])


def window_stats(data: pd.DataFrame, mc: pd.DataFrame, mask_data=None, mask_mc=None, window=WINDOW) -> dict:
    """
    Data/MC ratio, its error and the pull in the gap window, optionally within a subset of the events.

    The counts are used here and not returned: the result holds only the ratio, its statistical error and the pull. (A
    share of the selected events would be a count in disguise: with the ratio error it fixes the size of the sample.)
    """
    sd = in_window(data, window) if mask_data is None else in_window(data, window) & mask_data
    sm = in_window(mc, window) if mask_mc is None else in_window(mc, window) & mask_mc
    nd, nm = int(sd.sum()), float(mc.w[sm].sum())
    err = float(np.sqrt(nm + (mc.w[sm] ** 2).sum()))
    return {
        "ratio": nd / nm if nm else float("nan"), "error": float(np.sqrt(nd) / nm) if nm else float("nan"),
        "pull": (nd - nm) / err if err else float("nan"),
    }


def scan_windows(data: pd.DataFrame, mc: pd.DataFrame, width=3.0, lo=226.0, hi=280.0, step=0.5) -> pd.DataFrame:
    """Pull of every window of the given width over a vertex-z range, to put the gap in context."""
    rows = []
    for start in np.arange(lo, hi - width + 1e-9, step):
        rows.append({"z_lo": start, **window_stats(data, mc, window=(start, start + width))})
    return pd.DataFrame(rows)


def profit_four_cm_bins() -> pd.DataFrame | None:
    """Data/CV and the relative systematic band of the existing 4 cm PROfit comparison, if its output exists."""
    if not PROFIT_FILE.exists():
        return None
    import uproot

    total = 4294967295
    with uproot.open(PROFIT_FILE) as handle:
        hist = {k: handle["hist1d"][k].array(library="np") for k in
                ("variable", "subchannel", "prefix", "bin_low_edge", "bin_high_edge", "bin_content")}
        band = {k: handle["errorband"][k].array(library="np") for k in
                ("variable", "subchannel", "x_value", "y_value", "error_y_low", "error_y_high")}
    cv = (hist["variable"] == 0) & (hist["subchannel"] == total) & (hist["prefix"] == "CV")
    da = (hist["variable"] == 0) & (hist["subchannel"] == total) & (hist["prefix"] == "DATA")
    index = np.where((band["variable"] == 0) & (band["subchannel"] == total))[0]
    rows = []
    # PROfit writes the CV rows of the total subchannel twice (identical); keep the first of each bin
    first = np.sort(np.unique(hist["bin_low_edge"][cv], return_index=True)[1])
    for lo, hi, mc, obs in zip(hist["bin_low_edge"][cv][first], hist["bin_high_edge"][cv][first],
                               hist["bin_content"][cv][first], hist["bin_content"][da], strict=True):
        j = index[np.argmin(np.abs(band["x_value"][index] - (lo + hi) / 2))]
        rows.append({"lo": lo, "hi": hi, "ratio": obs / mc, "band_low": band["error_y_low"][j] / band["y_value"][j],
                     "band_high": band["error_y_high"][j] / band["y_value"][j]})
    return pd.DataFrame(rows)


def plot_vertex(data, mc, path: Path, profit: pd.DataFrame | None):
    """Data and MC in 0.5 cm bins around the gap, ratio below, 4 cm systematic band for scale."""
    import matplotlib.pyplot as plt

    edges = np.arange(244.0, 262.01, 0.5)
    centers = (edges[:-1] + edges[1:]) / 2
    nd, nm, var = binned(data.reco_vertex_z, mc.reco_vertex_z, mc.w, edges)
    total = len(data)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    top.step(edges[:-1], nm / total, where="post", color=MC_COLOR, label="Simulation")
    top.errorbar(centers, nd / total, yerr=np.sqrt(nd) / total, fmt="o", color=DATA_COLOR, ms=3.5, label="Data")
    top.axvspan(*WINDOW, color="C3", alpha=0.12, label=f"Gap window {WINDOW[0]}-{WINDOW[1]} cm")
    top.set_ylabel("Area normalized / 0.5 cm")
    top.set_yticklabels([])  # shape only: the fractions plus the error bars would fix the size of the sample
    top.legend(loc="upper left", fontsize=9)
    ratio = nd / nm
    bottom.errorbar(centers, ratio, yerr=np.sqrt(nd) / nm, fmt="o", color=DATA_COLOR, ms=3.5)
    bottom.axhline(1.0, color="0.5", lw=0.8)
    bottom.axvspan(*WINDOW, color="C3", alpha=0.12)
    if profit is not None:
        shown = profit[(profit.hi > 244) & (profit.lo < 262)]
        for k, (_, row) in enumerate(shown.iterrows()):
            bottom.fill_between([max(row.lo, 244), min(row.hi, 262)], 1 - row.band_low, 1 + row.band_high, color="C1",
                                alpha=0.35, lw=0, label="PROfit systematic band (4 cm bins)" if k == 0 else None)
        bottom.legend(loc="lower left", fontsize=8)
    bottom.set_ylim(0.2, 1.9)
    bottom.set_ylabel("Data / MC")
    bottom.set_xlabel("Reconstructed vertex z [cm]")
    common.save_figure(fig, path, arguments={"kind": "vertex_z"})


def plot_coordinates(data, mc, path: Path):
    """Data/MC ratio of the vertex and the muon and proton start and end z."""
    import matplotlib.pyplot as plt

    edges = np.arange(238.0, 268.01, 2.0)
    centers = (edges[:-1] + edges[1:]) / 2
    names = {
        "reco_vertex_z": "vertex", "reco_leading_muon_start_z": "muon start", "reco_leading_proton_start_z": "proton start",
        "reco_leading_muon_end_z": "muon end", "reco_leading_proton_end_z": "proton end",
    }
    fig, axes = plt.subplots(len(names), 1, figsize=(8, 9), sharex=True)
    for ax, (column, label) in zip(axes, names.items(), strict=True):
        nd, nm, var = binned(data[column], mc[column], mc.w, edges)
        ax.errorbar(centers, nd / nm, yerr=np.sqrt(nd) / nm, fmt="o", color=DATA_COLOR, ms=3.5)
        ax.axhline(1.0, color="0.5", lw=0.8)
        ax.axvspan(*WINDOW, color="C3", alpha=0.12)
        ax.set_ylim(0.3, 1.7)
        ax.set_ylabel(f"{label}\nData / MC")
    axes[-1].set_xlabel("z [cm]")
    common.save_figure(fig, path, arguments={"kind": "coordinates"})


def plot_end_vs_vertex(data, mc, path: Path):
    """Data/MC ratio in (vertex z, proton end z), for cells holding at least ``MIN_CELL_SHARE`` of the selected events."""
    import matplotlib.pyplot as plt

    edges = np.array([240, 246, 250, 252, 254, 256, 258, 264, 270], dtype=float)
    nd, _, _ = np.histogram2d(data.reco_vertex_z, data.reco_leading_proton_end_z, [edges, edges])
    nm, _, _ = np.histogram2d(mc.reco_vertex_z, mc.reco_leading_proton_end_z, [edges, edges], weights=mc.w)
    var, _, _ = np.histogram2d(mc.reco_vertex_z, mc.reco_leading_proton_end_z, [edges, edges], weights=mc.w**2)
    populated = nm >= MIN_CELL_SHARE * len(data)
    ratio = np.where(populated, nd / np.where(nm > 0, nm, 1), np.nan)
    significance = np.where(populated, (nd - nm) / np.sqrt(nm + var + 1e-9), 0.0)
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(ratio, origin="lower", cmap="RdBu_r", vmin=0.4, vmax=1.6)
    for i in range(len(edges) - 1):
        for j in range(len(edges) - 1):
            if populated[i, j]:
                mark = "*" if abs(significance[i, j]) > 2.5 else ""
                ax.text(j, i, f"{ratio[i, j]:.2f}{mark}", ha="center", va="center", fontsize=8)
    labels = [f"{a:.0f}-{b:.0f}" for a, b in zip(edges[:-1], edges[1:], strict=True)]
    ax.set_xticks(range(len(labels)), labels, rotation=45)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("Leading-proton end z [cm]")
    ax.set_ylabel("Vertex z [cm]")
    fig.colorbar(image, label="Data / MC (* = beyond 2.5 sigma)")
    common.save_figure(fig, path, arguments={"kind": "proton_end_vs_vertex"})


def split_table(data, mc) -> dict:
    """The gap-window ratio against x, y, run and the proton score (the outside ratio gives the reference)."""
    outside_d = (data.reco_vertex_z.sub(sum(WINDOW) / 2).abs() < 20) & ~in_window(data)
    outside_m = (mc.reco_vertex_z.sub(sum(WINDOW) / 2).abs() < 20) & ~in_window(mc)
    out = {}

    def slices(name, column, edges):
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:], strict=True):
            sd, sm = (data[column] >= lo) & (data[column] < hi), (mc[column] >= lo) & (mc[column] < hi)
            inside = window_stats(data, mc, sd, sm)
            ref_d, ref_m = int((outside_d & sd).sum()), float(mc.w[outside_m & sm].sum())
            rows.append({"lo": float(lo), "hi": float(hi), "ratio": inside["ratio"], "error": inside["error"],
                         "pull": inside["pull"], "outside_ratio": ref_d / ref_m if ref_m else float("nan")})
        out[name] = rows

    slices("vertex_x", "reco_vertex_x", [-200, -120, -60, -4, 4, 60, 120, 200])
    slices("vertex_y", "reco_vertex_y", [-200, -100, -50, 0, 50, 100, 170])
    slices("proton_score", "reco_leading_proton_proton_softmax", [0, 0.5, 0.9, 0.99, 1.01])
    slices("proton_primary", "reco_leading_proton_primary_softmax", [0, 0.3, 0.7, 0.9, 1.01])
    slices("muon_score", "reco_leading_muon_muon_softmax", [0, 0.5, 0.9, 1.01])
    runs = np.quantile(data.Run, [0, 0.25, 0.5, 0.75, 1.0])
    rows = []
    for lo, hi in zip(runs[:-1], runs[1:], strict=True):
        sd = (data.Run >= lo) & (data.Run <= hi)
        keep = (data.reco_vertex_z >= 226) & (data.reco_vertex_z < 280) & ~in_window(data) & sd
        k = keep.sum() / mc.w[(mc.reco_vertex_z >= 226) & (mc.reco_vertex_z < 280) & ~in_window(mc)].sum()
        nd = int((in_window(data) & sd).sum())
        nm = float(mc.w[in_window(mc)].sum() * k)
        rows.append({"run_lo": int(lo), "run_hi": int(hi), "ratio": nd / nm, "error": np.sqrt(nd) / nm})
    out["run"] = rows
    return out


def plot_splits(table: dict, path: Path):
    """Ratio in the gap window (red) and in the surrounding positions (grey) against each split variable."""
    import matplotlib.pyplot as plt

    names = ["vertex_x", "vertex_y", "proton_score", "muon_score", "proton_primary", "run"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
    for ax, name in zip(axes.flat, names, strict=True):
        rows = table[name]
        x = np.arange(len(rows))
        ratio = [r["ratio"] for r in rows]
        error = [r["error"] for r in rows]
        ax.errorbar(x, ratio, yerr=error, fmt="o", color="C3", label="gap window")
        if "outside_ratio" in rows[0]:
            ax.plot(x, [r["outside_ratio"] for r in rows], "s", color="0.5", label="outside, same slice")
        ax.axhline(1.0, color="0.5", lw=0.8)
        labels = ([f"{r['run_lo']}-{r['run_hi']}" for r in rows] if name == "run"
                  else [f"{r['lo']:g}-{r['hi']:g}" for r in rows])
        ax.set_xticks(x, labels, rotation=40, fontsize=7)
        ax.set_title(name.replace("_", " "), fontsize=10)
        ax.set_ylim(0, 1.8)
        ax.set_ylabel("Data / MC")
    axes[0, 0].legend(fontsize=8)
    common.save_figure(fig, path, arguments={"kind": "splits"})


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    parser = common.base_parser(__doc__)
    parser.set_defaults(output=common.DOCS / "apa_gap")
    args = parser.parse_args(argv)
    common.configure_style(fmt=args.format)
    inputs = [common.DATA_DIR / DATA_FILE, common.DATA_DIR / MC_FILE, PROFIT_FILE]
    common.set_provenance(__file__, inputs, {"window": list(WINDOW)}, restricted=[common.DATA_DIR / DATA_FILE])

    data, mc = load_data(), load_mc()
    mc = normalise(data, mc)
    profit = profit_four_cm_bins()
    out = args.output
    fmt = args.format

    plot_vertex(data, mc, out / f"gap_vertex_z.{fmt}", profit)
    plot_coordinates(data, mc, out / f"gap_coordinates.{fmt}")
    plot_end_vs_vertex(data, mc, out / f"gap_proton_end_vs_vertex.{fmt}")
    table = split_table(data, mc)
    plot_splits(table, out / f"gap_splits.{fmt}")

    scan = scan_windows(data, mc)
    worst = scan.loc[scan.pull.idxmin()]
    n_independent = (280 - 226) / 3.0
    from math import erfc, sqrt

    local_p = 0.5 * erfc(-worst.pull / sqrt(2))
    summary = {
        "window": list(WINDOW),
        "gap_window": window_stats(data, mc),
        "worst_3cm_window": {"z_lo": float(worst.z_lo), "pull": float(worst.pull), "ratio": float(worst.ratio),
                             "local_p": local_p, "global_p": 1 - (1 - local_p) ** n_independent},
        "profit_4cm": None if profit is None else profit[(profit.lo >= 244) & (profit.hi <= 262)].round(4).to_dict("records"),
        "splits": table,
    }
    for lo, hi, label in ((WINDOW[1], WINDOW[1] + 2, "just_above"), (WINDOW[0] - 2, WINDOW[0], "just_below")):
        summary[label] = window_stats(data, mc, window=(lo, hi))
    (out / "gap_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps({k: summary[k] for k in ("gap_window", "worst_3cm_window", "just_above", "just_below")},
                     indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1:])
