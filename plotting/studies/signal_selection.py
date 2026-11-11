"""
Signal-definition and selection figures for the v1.1.0 technote.

Replaces ``plotting/signal.ipynb`` and the efficiency cells of
``plotting/notebook.ipynb``. Two families of figures are made, both for
SBND, ICARUS Run 2 and ICARUS Run 4, and staged in
``docs/signal_selection_v1.1.0/signal_and_selection`` for review before
they are copied into the technote.

``histograms``
    Stacked true-neutrino-energy spectra drawn from the PROfit
    ``*_PROplot_histograms.root`` files in ``runs/legacy_root`` with the
    mutools TOML dispatcher (statistical errors only, systematics
    disabled, events per GeV):

    * signal split by interaction mode (CC QE / MEC / RES / DIS), for the
      inclusive channel and the exclusive single- and multi-proton
      channels (``signal_{inclusive,exclusive}_v1``);
    * selected events split by final state (target and other nu_mu CC,
      nu_e CC, NC, out-of-time and in-time cosmics), same channels
      (``selected_{inclusive,exclusive}_v1``). The legend carries the
      fraction of the total in each category (``counter_index = 0``, the
      selected-event count); for the signal plots it is the fraction of
      the CC QE/MEC/RES/DIS split (``counter_index = 3``).

``efficiency``
    Cut-by-cut signal efficiency read directly from the ``signal`` tree
    of the ``nosyst`` ntuples (``<data>/v1.1.0``). Cuts are applied
    cumulatively: TPC preselection, shower veto, pion veto, proton
    requirements, muon requirements. The efficiency is binomial, so the
    error bars are Clopper-Pearson intervals (the default of ROOT's
    TEfficiency, and well behaved close to 1 where the normal
    approximation gives zero error). The legend quotes the efficiency
    integrated over the *whole* signal sample, not just the plotted
    range, so the same cut reports the same number on every plot.
    Efficiencies are drawn versus true neutrino energy, muon azimuthal
    angle, leading-proton kinetic energy and muon kinetic energy, for the
    inclusive (Np), single-proton (1p, CCQE-like) and multi-proton
    (Mp, n-proton) selections.

The signal definition includes the flash-PE trigger emulation and the
truth-level geometry: the SBND-only cathode-crosser cut (ICARUS has two
cryostats and sets the flag for about 40% of its events), the SBND muon
end point at least 4 cm from the cathode, and, for ICARUS Run 2, the
truth-level EE filter containment. The reco-level TPC preselection
mirrors the ``selected/*.xml`` cuts.

Inputs: ``runs/legacy_root/{signal,selected}_{inclusive,exclusive}_v1_PROplot_histograms.root``
and ``<data>/v1.1.0/spineosc_numu_<detector>_v1.1.0_nosyst.root``.

Figures: ``signal_and_selection/{inclusive,exclusive_ccqe,exclusive_multiproton}_
{signal_bymode,selected_byfs}_<detector>`` and
``signal_and_selection/efficiency_<detector>_<variable>[_singlep|_multip]``.

Usage::

    python plotting/studies/signal_selection.py --part histograms
    python plotting/studies/signal_selection.py --part efficiency --detector sbnd --variable muon_ke
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import provenance  # noqa: E402

import common  # noqa: E402

#: Default location of the v1.1.0 ntuples.
DATA = common.DATA_ROOT

#: Sub-directory (of the staging directory and of ``--output``) that holds the figures.
SUBDIR = "signal_and_selection"

#: Staging directory for reviewed figures.
STAGE = common.DOCS / "signal_selection_v1.1.0"

#: Version of the medulla selection stamped under the efficiency legend.
MEDULLA_VERSION = "v1.1.0"

#: Watermark and legend title colour of the efficiency plots.
LEGEND_TITLE_COLOR = "#d67a11"

#: Interaction-mode split of the signal histograms, and the final-state
#: split (with matching colours) of the selected-event histograms. mutools
#: pairs the first label with the last subchannel stored in the file.
MODE_SUBCHANNELS = ["CC QE", "CC MEC", "CC RES", "CC DIS"]
FS_SUBCHANNELS = [
    r"Target $\nu_\mu$ CC",
    r"Other $\nu_\mu$ CC",
    r"$\nu_e$ CC",
    r"$\nu$ NC",
    "Out-of-time Cosmics",
    "In-time Cosmics",
]
FS_COLORS = ["C0", "C1", "C2", "C4", "C3", "C5"]
EXCLUSIVE_CHANNELS = ["Single-Proton", "Multi-Proton"]


@dataclass(frozen=True)
class Detector:
    """
    One detector / run configuration.

    Parameters
    ----------
    name : str
        Short name used in file names.
    title : str
        Detector label used as histogram and legend title.
    pot : str
        Target POT stamped on the histograms (from the XML ``<detector>`` blocks).
    ntuple : str
        File name of the ``nosyst`` ntuple, relative to ``<data>/v1.1.0``.
    tree : str
        Name of the detector directory inside ``events/`` in the ntuple.
    scale, sigma, threshold : float
        Trigger-emulation calibration: PE-to-threshold scale, its
        1-sigma uncertainty (documentation only) and the PE-equivalent
        threshold. Supplied directly, not derivable from the ntuples;
        same values as the XMLs.
    """

    name: str
    title: str
    pot: str
    ntuple: str
    tree: str
    scale: float
    sigma: float
    threshold: float


#: Detector table; the position is the detector index in the PROfit files.
DETECTORS = [
    Detector("sbnd", "SBND Simulation", "1.0e20 POT",
             "spineosc_numu_sbnd_v1.1.0_nosyst.root", "sbnd_large", 0.647, 0.004, 2000),
    Detector("icarus_run2", "ICARUS Run 2 Simulation", "2.0e20 POT",
             "spineosc_numu_icarus_run2_v1.1.0_nosyst.root", "icarus_run2", 0.615, 0.018, 6000),
    Detector("icarus_run4", "ICARUS Run 4 Simulation", "3.0e20 POT",
             "spineosc_numu_icarus_run4_v1.1.0_nosyst.root", "icarus_run4", 0.361, 0.013, 3000),
]


@dataclass(frozen=True)
class HistogramSet:
    """
    One stacked-histogram figure per detector from one PROfit file.

    Parameters
    ----------
    prefix : str
        Output name before the detector suffix.
    source : str
        Stem of the ``*_PROplot_histograms.root`` file.
    subchannels : list of str
        Legend labels of the stacked categories.
    channel : int
        Channel index of the variable (0: inclusive / single-proton, 1: multi-proton).
    xlim : tuple of float
        Energy range, GeV.
    counter_index : int
        Variable whose counts give the legend fractions.
    colors : list of str, optional
        Stack colours; mutools defaults when empty.
    channels : list of str, optional
        Channel labels, for the exclusive files.
    """

    prefix: str
    source: str
    subchannels: list[str]
    channel: int
    xlim: tuple[float, float]
    counter_index: int
    colors: list[str] = field(default_factory=list)
    channels: list[str] = field(default_factory=list)


#: The six histogram figures (times three detectors). The xlim of the
#: exclusive channels is narrower for single-proton than multi-proton.
HISTOGRAM_SETS = [
    HistogramSet("inclusive_signal_bymode", "signal_inclusive_v1", MODE_SUBCHANNELS,
                 0, (0.3, 2.0), 3),
    HistogramSet("exclusive_ccqe_signal_bymode", "signal_exclusive_v1", MODE_SUBCHANNELS,
                 0, (0.3, 1.5), 3, channels=EXCLUSIVE_CHANNELS),
    HistogramSet("exclusive_multiproton_signal_bymode", "signal_exclusive_v1", MODE_SUBCHANNELS,
                 1, (0.4, 2.0), 3, channels=EXCLUSIVE_CHANNELS),
    HistogramSet("inclusive_selected_byfs", "selected_inclusive_v1", FS_SUBCHANNELS,
                 0, (0.3, 2.0), 0, colors=FS_COLORS),
    HistogramSet("exclusive_ccqe_selected_byfs", "selected_exclusive_v1", FS_SUBCHANNELS,
                 0, (0.3, 1.5), 0, colors=FS_COLORS, channels=EXCLUSIVE_CHANNELS),
    HistogramSet("exclusive_multiproton_selected_byfs", "selected_exclusive_v1", FS_SUBCHANNELS,
                 1, (0.4, 2.0), 0, colors=FS_COLORS, channels=EXCLUSIVE_CHANNELS),
]

#: Edges of the true-energy efficiency (the project's energy edges, every
#: other one kept for per-bin statistics).
ENERGY_EDGES = np.array(
    [0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.25, 1.5, 2.0]
)


@dataclass(frozen=True)
class EfficiencyVariable:
    """
    One variable the efficiency is plotted against.

    Parameters
    ----------
    name : str
        Short name used in the figure file name.
    column : str
        Column of the ``signal`` tree.
    xlabel : str
        Axis label.
    edges : numpy.ndarray
        Fixed bin edges.
    """

    name: str
    column: str
    xlabel: str
    edges: np.ndarray


#: Variables; file names are ``efficiency_<detector>_<name>[suffix].pdf``.
EFFICIENCY_VARIABLES = [
    EfficiencyVariable("true_neutrino_energy", "true_neutrino_energy",
                       "True Neutrino Energy [GeV]", ENERGY_EDGES[::2]),
    EfficiencyVariable("muon_azimuthal_angle", "true_leading_muon_azimuthal_angle",
                       "Muon Azimuthal Angle [rad]", np.linspace(-np.pi, np.pi, 15)),
    EfficiencyVariable("proton_ke", "true_leading_proton_ke",
                       "Leading Proton Kinetic Energy [MeV]", np.linspace(50, 800, 15)),
    EfficiencyVariable("muon_ke", "true_leading_muon_ke",
                       "Muon Kinetic Energy [MeV]", np.linspace(140, 1200, 15)),
]


# ----------------------------------------------------------------------
# Signal and selection definitions
# ----------------------------------------------------------------------
def detector_of(df) -> str | None:
    """
    Detector name recorded on an ntuple table.

    Parameters
    ----------
    df : pandas.DataFrame
        Table loaded by `load_signal`.

    Returns
    -------
    str or None
        The detector name stored in ``df.attrs``.
    """
    return df.attrs.get("detector")


def apply_trigger_emulation(df):
    """
    Flash-PE trigger emulation with the per-detector calibration.

    Parameters
    ----------
    df : pandas.DataFrame
        Table loaded by `load_signal`.

    Returns
    -------
    pandas.Series
        Boolean mask of events above the PE-equivalent threshold.
    """
    det = next(d for d in DETECTORS if d.name == detector_of(df))
    return det.scale * df["event_largest_flash_pe"] > det.threshold


def apply_signal_tpc(df):
    """
    Truth-level signal geometry, including the trigger emulation.

    The cathode crosser cut is SBND only. SBND also requires the true
    muon end point to be at least 4 cm from the cathode; Run 2
    additionally requires the truth-level EE filter containment. The
    trigger emulation is part of the signal definition, so it is applied
    here and not as a cut in the efficiency plot.

    Parameters
    ----------
    df : pandas.DataFrame
        Table loaded by `load_signal`.

    Returns
    -------
    pandas.Series
        Boolean mask.
    """
    mask = (
        apply_trigger_emulation(df)
        & (df["true_fiducialize_cathode"] == 1)
        & (df["true_veto_sbnd_highy_highz"] == 1)
    )
    if detector_of(df) == "sbnd":
        mask &= (df["true_cathode_crosser"] == 0) & (df["true_leading_muon_end_x"].abs() >= 4)
    if detector_of(df) == "icarus_run2":
        mask &= df["true_user_containment_cut_filterEE"] == 1
    return mask


def apply_signal_ccqe(df):
    """Signal with true category 0 (CCQE-like, single proton) inside the signal geometry."""
    return (df["true_category"] == 0) & apply_signal_tpc(df)


def apply_signal_nproton(df):
    """Signal with true category 2 (n-proton) inside the signal geometry."""
    return (df["true_category"] == 2) & apply_signal_tpc(df)


def apply_signal(df):
    """
    Full signal definition: true category 0 or 2 inside the signal geometry.

    The ``signal`` tree is already restricted to these categories, so the
    category part is a no-op there but kept explicit for broader tables.
    """
    return apply_signal_ccqe(df) | apply_signal_nproton(df)


def apply_tpc_cut(df):
    """
    TPC preselection, the project's base selection.

    Cathode fiducialization, oscillation fiducial and containment cuts,
    the high-y/high-z veto, the SBND cathode-crosser removal, the SBND
    muon end point at least 4 cm from the cathode and the ICARUS Run 2
    EE filter containment. Mirrors the ``selected/*.xml`` cuts.

    Parameters
    ----------
    df : pandas.DataFrame
        Table loaded by `load_signal`.

    Returns
    -------
    pandas.Series
        Boolean mask.
    """
    mask = (
        (df["reco_fiducialize_cathode"] == 1)
        & (df["reco_fiducial_cut_osc"])
        & (df["reco_user_containment_cut"] == 1)
        & (df["reco_veto_sbnd_highy_highz"] == 1)
    )
    if detector_of(df) == "sbnd":
        mask &= (df["reco_cathode_crosser"] == 0) & (df["reco_leading_muon_end_x"].abs() >= 4)
    if detector_of(df) == "icarus_run2":
        mask &= df["reco_user_containment_cut_filterEE"] == 1
    return mask


def apply_shower_multiplicity(df):
    """Shower veto: no photon and no electron."""
    return (df["reco_photon_multiplicity"] == 0) & (df["reco_electron_multiplicity"] == 0)


def apply_pion_multiplicity(df):
    """Pion veto: no pion."""
    return df["reco_pion_multiplicity"] == 0


def apply_proton_multiplicity_1p(df):
    """Exactly one proton."""
    return df["reco_proton_multiplicity"] == 1


def apply_proton_multiplicity_multip(df):
    """More than one proton."""
    return df["reco_proton_multiplicity"] > 1


def apply_proton_multiplicity(df):
    """At least one proton."""
    return apply_proton_multiplicity_1p(df) | apply_proton_multiplicity_multip(df)


def apply_muon_cut(df):
    """Exactly one muon."""
    return df["reco_muon_multiplicity"] == 1


@dataclass(frozen=True)
class EfficiencySelection:
    """
    One signal / proton-requirement pairing.

    Parameters
    ----------
    suffix : str
        File-name suffix (empty for the inclusive selection).
    signal : Callable
        Signal mask (the efficiency denominator).
    proton_cut : Callable
        Reco proton-multiplicity requirement.
    label : str
        Final-state label drawn at the top right.
    """

    suffix: str
    signal: Callable
    proton_cut: Callable
    label: str


EFFICIENCY_SELECTIONS = [
    EfficiencySelection("", apply_signal, apply_proton_multiplicity, r"1$\mu$Np0$\pi$"),
    EfficiencySelection("_singlep", apply_signal_ccqe, apply_proton_multiplicity_1p,
                        r"1$\mu$1p0$\pi$"),
    EfficiencySelection("_multip", apply_signal_nproton, apply_proton_multiplicity_multip,
                        r"1$\mu$(M>1)p0$\pi$"),
]


# ----------------------------------------------------------------------
# Efficiency computation
# ----------------------------------------------------------------------
def load_signal(detector: Detector, data: Path = DATA):
    """
    Read the ``signal`` tree of one detector's ``nosyst`` ntuple.

    Parameters
    ----------
    detector : Detector
        Detector to load.
    data : Path, optional
        Directory holding ``v1.1.0/``.

    Returns
    -------
    pandas.DataFrame
        The signal events, with ``attrs["detector"]`` set because several
        cuts depend on the detector.
    """
    import uproot

    path = data / "v1.1.0" / detector.ntuple
    with uproot.open(path) as rf:
        df = rf[f"events/{detector.tree}/signal"].arrays(library="pd")
    df.attrs["detector"] = detector.name
    print(f"{detector.name}: signal events = {len(df)}")
    return df


def equal_population_edges(values, lo: float, hi: float, nbins: int) -> np.ndarray:
    """
    Bin edges that place an equal number of parent-population entries in each bin.

    The parent population is the denominator of the efficiency (the
    signal sample), so each point carries comparable statistics.

    Parameters
    ----------
    values : array_like
        Parent-population values.
    lo, hi : float
        Range to bin over; entries outside it are ignored.
    nbins : int
        Number of bins requested.

    Returns
    -------
    numpy.ndarray
        Edges (``nbins + 1``, or fewer if the distribution cannot support
        that many distinct quantiles).

    Raises
    ------
    ValueError
        If no entries fall inside ``[lo, hi]``.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (v >= lo) & (v <= hi)]
    if v.size == 0:
        raise ValueError(f"no parent-population entries in [{lo}, {hi}]")
    edges = np.quantile(v, np.linspace(0.0, 1.0, nbins + 1))
    edges[0], edges[-1] = lo, hi
    # A spiked distribution can repeat a quantile; collapse duplicates.
    unique = np.unique(edges)
    if unique.size < edges.size:
        print(
            f"  note: {edges.size - unique.size} duplicate quantile edge(s) "
            f"collapsed; using {unique.size - 1} bins instead of {nbins}"
        )
    return unique


def make_counts(energy, mask0, mask1, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Bin a variable under two masks.

    Parameters
    ----------
    energy : pandas.Series
        The variable to bin.
    mask0 : pandas.Series
        Denominator mask (the signal sample).
    mask1 : pandas.Series
        Numerator mask (signal passing the cuts so far).
    edges : numpy.ndarray
        Bin edges.

    Returns
    -------
    passed, total : numpy.ndarray
        Per-bin numerator and denominator counts.
    """
    total, _ = np.histogram(energy[mask0], bins=edges)
    passed, _ = np.histogram(energy[mask1], bins=edges)
    return passed, total


def efficiency_errors(
    passed: np.ndarray,
    total: np.ndarray,
    method: str = "clopper-pearson",
    level: float = 0.6827,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Per-bin efficiency and its statistical uncertainty.

    The numerator is a subset of the denominator, so the content is
    binomial. ``"clopper-pearson"`` is the exact interval (asymmetric,
    inside [0, 1]); ``"normal"`` is sqrt(eff (1 - eff) / total).

    Parameters
    ----------
    passed, total : numpy.ndarray
        Per-bin numerator and denominator counts.
    method : str, optional
        ``"clopper-pearson"`` or ``"normal"``.
    level : float, optional
        Coverage of the interval; the default is one sigma.

    Returns
    -------
    eff : numpy.ndarray
        Efficiency per bin (NaN where the denominator is empty).
    yerr : numpy.ndarray
        ``(2, nbins)`` lower/upper errors for ``ax.errorbar``.

    Raises
    ------
    ValueError
        If ``method`` is not supported.
    """
    from scipy.stats import beta

    passed = np.asarray(passed, dtype=float)
    total = np.asarray(total, dtype=float)
    filled = total > 0
    eff = np.full(total.shape, np.nan)
    np.divide(passed, total, out=eff, where=filled)

    if method == "normal":
        err = np.full(total.shape, np.nan)
        np.divide(eff * (1.0 - eff), total, out=err, where=filled)
        err = np.sqrt(err)
        return eff, np.vstack([err, err])
    if method != "clopper-pearson":
        raise ValueError(f"unsupported method: {method!r}")

    alpha = 1.0 - level
    # beta.ppf is NaN at the boundaries, where the interval is one-sided.
    low = beta.ppf(alpha / 2.0, passed, total - passed + 1.0)
    high = beta.ppf(1.0 - alpha / 2.0, passed + 1.0, total - passed)
    low = np.where(np.isnan(low), 0.0, low)
    high = np.where(np.isnan(high), 1.0, high)
    return eff, np.vstack([eff - low, high - eff])


def add_efficiency(
    ax,
    passed: np.ndarray,
    total: np.ndarray,
    edges: np.ndarray,
    label: str | None = None,
    color: str | None = None,
    method: str = "clopper-pearson",
    integrated: float | None = None,
) -> None:
    """
    Draw one efficiency curve as points with error bars, x error spanning each bin.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    passed, total : numpy.ndarray
        Per-bin counts.
    edges : numpy.ndarray
        Bin edges.
    label : str, optional
        Legend label; the integrated efficiency is appended.
    color : str, optional
        Marker colour.
    method : str, optional
        Passed to `efficiency_errors`.
    integrated : float, optional
        Efficiency over the whole sample to quote in the legend. Computing
        it from the binned counts would cover only the plotted range, so
        the binned counts are used only when this is not supplied.
    """
    eff, yerr = efficiency_errors(passed, total, method=method)
    if integrated is None:
        integrated = passed.sum() / total.sum() if total.sum() else np.nan
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.errorbar(
        centers, eff, xerr=0.5 * np.diff(edges), yerr=yerr, fmt="o",
        markersize=4, linewidth=1.0, capsize=0, label=f"{label} ({integrated:.1%})", color=color,
    )


def build_efficiency_plot(
    df,
    apply_sig: Callable,
    apply_pcut: Callable,
    edges: np.ndarray,
    var_key: str = "true_neutrino_energy",
    xlabel: str = "True Neutrino Energy [GeV]",
    detector_label: str | None = None,
    channel_label: str | None = None,
    mark_fs: str | None = None,
    binning: str = "fixed",
    nbins: int | None = None,
    error_method: str = "clopper-pearson",
    savefig: Path | None = None,
):
    """
    Cut-by-cut efficiency versus ``var_key`` with binomial error bars.

    Parameters
    ----------
    df : pandas.DataFrame
        Table from `load_signal`.
    apply_sig : Callable
        Signal mask function (the denominator).
    apply_pcut : Callable
        Proton requirement applied as the fourth cut.
    edges : numpy.ndarray
        Bin edges.
    var_key : str, optional
        Column plotted on the x axis.
    xlabel : str, optional
        Axis label.
    detector_label, channel_label : str, optional
        Legend title lines.
    mark_fs : str, optional
        Final-state label at the top right.
    binning : str, optional
        ``"fixed"`` uses ``edges``; ``"equal-population"`` keeps only the
        range ``[edges[0], edges[-1]]`` and rebins it so each bin holds
        the same number of signal events.
    nbins : int, optional
        Number of bins for ``"equal-population"``; defaults to that of ``edges``.
    error_method : str, optional
        ``"clopper-pearson"`` or ``"normal"``.
    savefig : Path, optional
        Where to write the figure.

    Returns
    -------
    matplotlib.figure.Figure
        The (still open) figure.

    Raises
    ------
    ValueError
        If ``binning`` is not supported.
    """
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis

    fig, ax = plt.subplots(figsize=(8, 6))
    energy = df[var_key]
    mask_signal = apply_sig(df)

    edges = np.asarray(edges, dtype=float)
    if binning == "equal-population":
        edges = equal_population_edges(
            energy[mask_signal], edges[0], edges[-1], len(edges) - 1 if nbins is None else nbins
        )
    elif binning != "fixed":
        raise ValueError(f"unsupported binning: {binning!r}")

    cuts = {
        "TPC Preselection": apply_tpc_cut,
        "Shower Veto": apply_shower_multiplicity,
        "Pion Veto": apply_pion_multiplicity,
        "Proton Requirements": apply_pcut,
        "Muon Requirements": apply_muon_cut,
    }

    # Cuts are cumulative; one curve after each.
    selection = np.repeat(True, len(df))
    n_signal = int(np.count_nonzero(mask_signal))
    for cut_name, cut_func in cuts.items():
        selection &= cut_func(df)
        passed, total = make_counts(energy, mask_signal, mask_signal & selection, edges)
        # Integrated over the whole signal sample, not the plotted range.
        integrated = np.count_nonzero(mask_signal & selection) / n_signal if n_signal else np.nan
        add_efficiency(ax, passed, total, edges, label=cut_name, method=error_method,
                       integrated=integrated)

    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Efficiency")
    mark_axis(ax, common.WATERMARK)
    if mark_fs:
        mark_axis(ax, mark_fs, alignment="right")

    # Legend lower left, lifted so the version stamp fits above the x axis.
    legend = ax.legend(loc="lower left", bbox_to_anchor=(0.02, 0.045))
    title = detector_label if channel_label is None else f"{detector_label}\n{channel_label}"
    legend.set_title(title)
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_fontsize(14)
    legend.get_title().set_color(LEGEND_TITLE_COLOR)

    # Version stamp below the legend, small and faded.
    fig.canvas.draw()
    box = legend.get_window_extent().transformed(ax.transAxes.inverted())
    ax.text(0.5 * (box.x0 + box.x1), box.y0 - 0.004, f"medulla {MEDULLA_VERSION}",
            transform=ax.transAxes, fontsize=8, color="#555555", alpha=0.7, ha="center", va="top")

    if savefig is not None:
        Path(savefig).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savefig)
    return fig


def make_efficiency_figures(
    detectors: list[Detector],
    variables: list[EfficiencyVariable],
    stage: Path,
    data: Path = DATA,
    fmt: str = "pdf",
) -> list[Path]:
    """
    Draw every efficiency figure for the chosen detectors and variables.

    Parameters
    ----------
    detectors : list of Detector
        Detectors to process.
    variables : list of EfficiencyVariable
        Variables to plot against.
    stage : Path
        Directory the ``signal_and_selection`` folder is created in.
    data : Path, optional
        Directory holding ``v1.1.0/``.
    fmt : str, optional
        Figure file format.

    Returns
    -------
    list of Path
        Files written.
    """
    import matplotlib.pyplot as plt

    written = []
    for det in detectors:
        df = load_signal(det, data)
        for var in variables:
            for sel in EFFICIENCY_SELECTIONS:
                path = stage / SUBDIR / f"efficiency_{det.name}_{var.name}{sel.suffix}.{fmt}"
                fig = build_efficiency_plot(
                    df, sel.signal, sel.proton_cut, var.edges, var_key=var.column,
                    xlabel=var.xlabel, detector_label=det.title, mark_fs=sel.label, savefig=path,
                )
                plt.close(fig)
                provenance.record(path, inputs=[data / "v1.1.0" / det.ntuple], script=__file__,
                                  arguments={"detector": det.name, "variable": var.name,
                                             "selection": sel.label, "kind": "efficiency"},
                                  versions={"selection": common.SELECTION_VERSION})
                written.append(path)
                print(f"wrote {path}")
    return written


# ----------------------------------------------------------------------
# Stacked histograms from the PROfit files
# ----------------------------------------------------------------------
def histogram_config(hset: HistogramSet, detectors: list[Detector], runs: Path, output: Path):
    """
    Build the mutools dispatcher configuration for one histogram set.

    Parameters
    ----------
    hset : HistogramSet
        The set to draw.
    detectors : list of Detector
        Detectors to draw (one ``[[plot]]`` each).
    runs : Path
        Directory holding the PROfit histogram files.
    output : Path
        Directory the dispatcher writes to.

    Returns
    -------
    dict
        Configuration accepted by ``mutools.plotting.run``.
    """
    general = {
        "input": str(runs / f"{hset.source}_PROplot_histograms.root"),
        "output": str(output),
        "code_version": common.CODE_VERSION,
        "selection_version": common.SELECTION_VERSION,
        "subchannels": hset.subchannels,
        "detectors": [d.title for d in DETECTORS],
        "counter_index": hset.counter_index,
        "savefig": True,
    }
    if hset.channels:
        general["channels"] = hset.channels
    plots = []
    for det in detectors:
        plot = {
            "type": "histogram",
            "variable": 0,
            "channel": hset.channel,
            "detectors": [DETECTORS.index(det)],
            "xlabel": "True Neutrino Energy [GeV]",
            "ylabel": "Events / GeV",
            "xlim": list(hset.xlim),
            "watermark": common.WATERMARK,
            "disable_systematics": True,
            "counter_fmt": ".1%",
            "scale_by_width": "enabled",
            "right_watermark": det.pot,
        }
        if hset.colors:
            plot["colors"] = hset.colors
        plots.append(plot)
    return {"general": general, "plot": plots}


def make_histogram_figures(
    detectors: list[Detector],
    sets: list[HistogramSet],
    output: Path,
    stage: Path | None,
    runs: Path = common.LEGACY_ROOT,
    fmt: str = "pdf",
) -> list[Path]:
    """
    Draw the stacked histograms and move them to their staged names.

    The dispatcher names files ``hist_<detector>_<variable>_<version>``;
    these are renamed to ``<prefix>_<detector>`` under ``stage``. Missing
    outputs are reported rather than raised.

    Parameters
    ----------
    detectors : list of Detector
        Detectors to draw.
    sets : list of HistogramSet
        Figure sets to draw.
    output : Path
        Directory the dispatcher writes to.
    stage : Path or None
        Staging directory; ``None`` leaves the figures under their
        dispatcher names in ``output``.
    runs : Path, optional
        Directory with the PROfit histogram files.
    fmt : str, optional
        Figure file format.

    Returns
    -------
    list of Path
        Final paths of the figures.
    """
    from mutools.plotting import run

    written = []
    for hset in sets:
        run(histogram_config(hset, detectors, runs, output))
        for det in detectors:
            idx = DETECTORS.index(det)
            old = output / f"hist_{idx}_0_{common.CODE_VERSION}.{fmt}"
            if stage is None:
                written.append(old)
                continue
            new = stage / SUBDIR / f"{hset.prefix}_{det.name}.{fmt}"
            new.parent.mkdir(parents=True, exist_ok=True)
            if old.exists():
                old.rename(new)
                print(f"Renamed {old} to {new}")
                written.append(new)
            else:
                print(f"Warning: {old} does not exist and cannot be renamed to {new}")
    return written


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.add_argument("--part", choices=["histograms", "efficiency", "all"], default="all",
                        help="which family of figures to make")
    parser.add_argument("--detector", action="append", choices=[d.name for d in DETECTORS],
                        help="restrict to this detector (repeatable)")
    parser.add_argument("--variable", action="append",
                        choices=[v.name for v in EFFICIENCY_VARIABLES],
                        help="restrict the efficiency plots to this variable (repeatable)")
    parser.add_argument("--histogram", action="append",
                        choices=[h.prefix for h in HISTOGRAM_SETS],
                        help="restrict the stacked histograms to this set (repeatable)")
    parser.add_argument("--stage", type=Path, default=STAGE,
                        help="staging directory for the finished figures")
    parser.add_argument("--no-stage", action="store_true",
                        help="leave histograms under the dispatcher names in --output")
    parser.add_argument("--data", type=Path, default=DATA, help="directory holding v1.1.0 ntuples")
    parser.add_argument("--runs", type=Path, default=common.LEGACY_ROOT,
                        help="directory with the PROfit histogram files")
    args = parser.parse_args(argv)

    common.configure_style(fmt=args.format)
    detectors = [d for d in DETECTORS if not args.detector or d.name in args.detector]

    if args.part in ("histograms", "all"):
        sets = [h for h in HISTOGRAM_SETS if not args.histogram or h.prefix in args.histogram]
        make_histogram_figures(detectors, sets, args.output, None if args.no_stage else args.stage,
                               args.runs, args.format)
    if args.part in ("efficiency", "all"):
        variables = [
            v for v in EFFICIENCY_VARIABLES if not args.variable or v.name in args.variable
        ]
        make_efficiency_figures(detectors, variables, args.stage, args.data, args.format)


if __name__ == "__main__":
    main()
