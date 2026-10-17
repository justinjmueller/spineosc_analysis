"""
SBND data access request: data/MC chi-square tables and panels, phases 1-3
(plus phases 4 and 5, which live in ``data_access_phases.py``).

Replaces ``plotting/data_access_results.ipynb``. It reads finished PROfit
plot outputs from ``runs/legacy_root`` (``t<tier>_<variable>_<arm>_PROplot_histograms.root``)
and never reruns PROfit. Figures feed ``figures/sbnd_data_access`` of the
technote (``sbnd_data_access.tex``) as ``sbndaccess_<variable>_sbnd_<channel>_<arm>``
and ``pvalues_<arm>``.

Phases (``--phase``)
--------------------
``1``  Area-normalised, shape-only data/MC, by tier.
    The chi-square PROfit draws on its channel plots is also written to
    the ``chi2`` tree of each run's ``_PROplot_histograms.root``, so the
    numbers are read back, not recomputed. The runs use both ``area-norm``
    and ``shapeonly`` (set in every ``xml/sbnd_data_access/phase1-tier*/
    *_cvrw_plot.ini``): area-norm rescales the prediction onto the data
    integral and costs one degree of freedom (``ndf = nbins - 1``), and
    ``shapeonly`` also projects each systematic covariance onto
    per-channel shape. The ndf assertion proves area normalisation only,
    not that ``shapeonly`` was set (an earlier set of runs passed it
    without being shape-only); that is read off the ini.
    Several studies have more than one run on disk; the table names the
    file each number came from so a stale or rebinned run is visible.
    Caveats: the numbers carry whatever detector-systematics set and
    selection the runs used; the dirt samples are empty in the inputs
    (``events/sbnd_dirt/selected`` has no entries); the chi-square is
    whichever ``--chi2`` the run used (PROfit default: Neyman); the
    prediction is rescaled to the data integral so the data POT
    mismatch in the ``_data.xml`` files does not enter; and chi-square carries about
    +-0.6 of run-to-run scatter from the unpinned seed in the
    spline-to-covariance throws, and with 16 comparisons a single
    p-value near 0.01 is unremarkable.
``2``  Shape + normalisation: the same runs with ``cvrw_cvrw_phase2``.
    No area normalisation, so ``ndf = nbins`` (the phase-1 accessor would
    fail its assertion) and the MC/data POT matching now matters.
    The runs read the same ``t<N>_<var>_cvrw`` caches as phase 1.
``3``  Visual inspection: the phase-2 files with data overlaid on the
    absolutely normalised prediction, for what a p-value misses (a slope
    the band covers, one discrepant bin, structure at a threshold).
``4``, ``5``  See ``data_access_phases.py`` (conditional-constraint
    overlays; global null fit pulls). They are dispatched from here.

Tiers are gated. ``Phase.enabled_tiers`` controls what is tabulated and
plotted, and ``Phase.overlay_tiers`` separately which tiers may show data
points. Where the request withholds data overlays (phase 1 tiers 2 and 3,
phase 2 all tiers) the panels are MC plus uncertainty band with the
chi-square and p-value quoted but no points drawn, so a tier can be
assessed before its data is displayed. Each phase has its own gate, so
widening one does not drag the others along.

Input path note
---------------
The notebook defines ``INPUT = .../dp_test09142026_v1_PROplot.root`` but
never uses it (it was deleted with class F; see ``docs/RESIDUAL_TASKS.md``).
It is not ported. The nearest live file for that study is
``t1_dp_cvrw_cvrw_PROplot.root`` (``BASELINE_FILE``), which is what
phase 1 reads for delta p anyway.

Run::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/data_access.py \\
        --phase 1 --what table
"""

from __future__ import annotations

import shutil
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Directory of PROfit plot outputs read by every phase.
RUNS_DIR = common.LEGACY_ROOT

#: Live stand-in for the notebook's unused, deleted ``dp_test09142026_v1_PROplot.root``.
BASELINE_FILE = RUNS_DIR / "t1_dp_cvrw_cvrw_PROplot.root"

#: Study -> (tier, display label, channel indices). ``energy_qel`` is
#: single-proton only, by construction of its XML.
VARIABLES: dict[str, tuple[int, str, tuple[int, ...]]] = {
    "dp": (1, r"$\delta p$", (0, 1)),
    "dpt": (1, r"$\delta p_{T}$", (0, 1)),
    "dalphat": (1, r"$\delta\alpha_{T}$", (0, 1)),
    "dphit": (1, r"$\delta\phi_{T}$", (0, 1)),
    "muon_azimuth": (1, r"$\phi_{\mu}$", (0, 1)),
    "proton_ke": (1, r"$T_{p}$", (0, 1)),
    "muon_score": (1, "Muon PID score", (0, 1)),
    "proton_score": (1, "Proton PID score", (0, 1)),
    "muon_p": (2, r"$p_{\mu}$", (0, 1)),
    "muon_polar_angle": (2, r"$\theta_{\mu}$", (0, 1)),
    "opening_angle": (2, r"$\theta_{\mu p}$", (0,)),
    "vis_energy": (3, r"$E^{reco}_{vis}$", (0, 1)),
    "energy_qel": (3, r"$E^{reco}_{ccqe}$", (0,)),
}

#: Detector index (XML declaration order) of SBND, and the detectors drawn.
SBND = 0
DETECTORS = {0: "SBND"}

#: Channel indices -> XML channel names (the selection itself).
CHANNELS = {0: "ccqelike", 1: "nproton"}

#: Channel presentation, in plot legends and in file names
#: (``datamc_<detector>_<variable>_<singlep|multip>`` in the technote).
CHANNEL_LABELS = {0: "Single-Proton", 1: "Multi-Proton"}
CHANNEL_TOKENS = {0: "singlep", 1: "multip"}

#: Legend labels are the XML subchannel declaration order REVERSED (mutools
#: pairs the first label with the last subchannel); the colours are in the
#: XML order as written (intime, cos, nc, nuecc, othercc, target). Checked
#: against a rendered panel: Target is #99CCFF, not #FBBC04.
SUBCHANNELS = [
    "Target",
    r"Other $\nu$ CC",
    r"$\nu_e$ CC",
    r"$\nu$ NC",
    "Cosmics",
    "In-time Cosmics",
]
SUBCHANNEL_COLORS = ["C0", "C1", "C2", "C4", "C3", "C5"]

#: Study -> (x-axis label, xlim). Ranges match the var0 binning in the XMLs.
#: xlim is a pair, or ``{channel: pair}`` where the channels are binned
#: over different ranges (``vis_energy``).
PLOT_AXES: dict[str, tuple[str, list[float] | dict[int, list[float]]]] = {
    "dp": (r"Reconstructed $\delta p$ [MeV/c]", [0.0, 1000.0]),
    "dpt": (r"Reconstructed $\delta p_{T}$ [MeV/c]", [0.0, 1000.0]),
    "dalphat": (r"Reconstructed $\delta\alpha_{T}$ [rad]", [0.0, 3.14]),
    "dphit": (r"Reconstructed $\delta\phi_{T}$ [rad]", [0.0, 3.14]),
    "muon_score": ("Muon PID Score", [0.30, 1.00]),
    "proton_score": ("Proton PID Score", [0.40, 1.00]),
    "muon_azimuth": ("Muon Azimuthal Angle [rad]", [-3.14, 3.14]),
    "proton_ke": ("Leading Proton Kinetic Energy [MeV]", [40.0, 600.0]),
    "muon_p": ("Muon Momentum [MeV/c]", [225.0, 1225.0]),
    "muon_polar_angle": ("Muon Polar Angle [rad]", [0.0, 3.14]),
    "opening_angle": ("Muon-Proton Opening Angle [rad]", [0.0, 3.14]),
    "vis_energy": ("Reconstructed Visible Energy [GeV]", {0: [0.3, 1.5], 1: [0.4, 2.0]}),
    "energy_qel": (r"Reconstructed $E_{ccqe}$ [GeV]", [0.2, 1.5]),
}

#: Phase-1 y-extent per variable, shared by the channels (area
#: normalisation puts them on the same footing). Each is ~55% above the
#: tallest of CV, band top and (tier 1) data + error. Tier 2 and 3 entries
#: were set from CV and band only, so setting a y-scale never required
#: looking at their data.
PHASE1_YLIM: dict[str, list[float]] = {
    "dp": [0.0, 0.006],
    "dpt": [0.0, 0.006],
    "dalphat": [0.0, 1.5],
    "dphit": [0.0, 3.0],
    "muon_score": [0.0, 15.0],
    "proton_score": [0.0, 50.0],
    "muon_azimuth": [0.0, 0.50],
    "proton_ke": [0.0, 0.01],
    "muon_p": [0.0, 0.005],
    "muon_polar_angle": [0.0, 1.11],
    "opening_angle": [0.0, 1.44],
    "vis_energy": [0.0, 3.18],
    "energy_qel": [0.0, 2.44],
}

#: Phase-1 y-limit overrides keyed ``(variable,)``, ``(variable, detector)``
#: or ``(variable, detector, channel)``; the most specific wins.
PHASE1_YLIM_OVERRIDES: dict[tuple, list[float]] = {}

#: Phase-2 y-limits per variable and channel. Absolute normalisation makes
#: the channels differ in scale (single-proton peak 2.8 to 6.2 times the
#: multi-proton one, median 3.4), so one shared limit would flatten every
#: multi-proton panel. Measured from CV + band only (overlays are off):
#: peak plus ~55% headroom, rounded up to three significant figures.
PHASE2_YLIM: dict[str, dict[int, list[float]]] = {
    "dp": {0: [0.0, 20.0], 1: [0.0, 4.00]},
    "dpt": {0: [0.0, 22.0], 1: [0.0, 5.00]},
    "dalphat": {0: [0.0, 3800], 1: [0.0, 1200]},
    "dphit": {0: [0.0, 10000], 1: [0.0, 1600]},
    "muon_score": {0: [0.0, 45000], 1: [0.0, 13000]},
    "proton_score": {0: [0.0, 160000], 1: [0.0, 53000]},
    "muon_azimuth": {0: [0.0, 1400], 1: [0.0, 500]},
    "proton_ke": {0: [0.0, 35.0], 1: [0.0, 9.75]},
    "proton_p": {0: [0.0, 15.0], 1: [0.0, 4.97]},
    "muon_p": {0: [0.0, 17.0], 1: [0.0, 5.81]},
    "muon_polar_angle": {0: [0.0, 4500], 1: [0.0, 1600]},
    "opening_angle": {0: [0.0, 5500], 1: [0.0, 1400]},
    "muon_ke": {0: [0.0, 20.0], 1: [0.0, 7.00]},
    "vis_energy": {0: [0.0, 12500], 1: [0.0, 4000]},
    "energy_qel": {0: [0.0, 9500]},
}

#: Phase-3 y-limits, independent of phase 2. Phase 3 draws a ratio panel
#: (mutools uses ``height_ratios=[4, 1]``, leaving the main axis 4/5 as
#: tall) and the legend is sized in points, so it covers 1/0.8 as much of
#: the shorter axis. Each value started as the phase-2 limit x 1.159, which
#: holds the absolute clearance above the peak fixed, and was rounded up.
#: Nothing recomputes them; tune freely.
PHASE3_YLIM: dict[str, dict[int, list[float]]] = {
    "dp": {0: [0.0, 25], 1: [0.0, 7]},
    "dpt": {0: [0.0, 30], 1: [0.0, 6]},
    "dalphat": {0: [0.0, 4500], 1: [0.0, 1500]},
    "dphit": {0: [0.0, 12000], 1: [0.0, 2000]},
    "muon_score": {0: [0.0, 60000], 1: [0.0, 18000]},
    "proton_score": {0: [0.0, 200000], 1: [0.0, 70000]},
    "muon_azimuth": {0: [0.0, 1800], 1: [0.0, 600]},
    "proton_ke": {0: [0.0, 45], 1: [0.0, 12]},
    "proton_p": {0: [0.0, 18], 1: [0.0, 6]},
    "muon_p": {0: [0.0, 20], 1: [0.0, 7]},
    "muon_polar_angle": {0: [0.0, 6000], 1: [0.0, 2000]},
    "opening_angle": {0: [0.0, 8000], 1: [0.0, 2000]},
    "muon_ke": {0: [0.0, 25], 1: [0.0, 9]},
    "vis_energy": {0: [0.0, 18000], 1: [0.0, 5000]},
    "energy_qel": {0: [0.0, 14000]},
}

#: Per-variable legend placement ``(loc, ncol)`` and its default.
PLOT_LEGEND: dict[str, tuple[str, int]] = {}
DEFAULT_LEGEND = ("upper right", 2)

#: Telescope stamp on the p-value summary figures.
POT_LABEL = r"$4.51 \times 10^{18}$ POT"

#: Colours of the two channels and of the threshold line in the p-value figure.
PVALUE_COLORS = {"ccqelike": "C1", "nproton": "C2"}
THRESHOLD_COLOR = "#eb6834"
PVALUE_THRESHOLD = 0.05


def _ylim_phase1(study: str, detector: int, channel: int) -> list[float] | None:
    """
    Phase-1 y-limit, honouring `PHASE1_YLIM_OVERRIDES` (most specific wins).

    Parameters
    ----------
    study : str
        Variable name.
    detector, channel : int
        Panel indices.

    Returns
    -------
    list of float or None
        The limits, or None to autoscale.
    """
    for key in ((study, detector, channel), (study, detector), (study,)):
        if key in PHASE1_YLIM_OVERRIDES:
            return PHASE1_YLIM_OVERRIDES[key]
    return PHASE1_YLIM.get(study)


def _ylim_per_channel(table: dict[str, dict[int, list[float]]]) -> Callable:
    """
    Build a y-limit lookup from a per-channel table.

    Falls back to None (autoscale) rather than to another channel's limit:
    borrowing the single-proton scale for a multi-proton panel is the
    mistake the table exists to prevent, and an autoscaled panel is
    obvious where a silently squashed one is not.

    Parameters
    ----------
    table : dict
        Variable -> {channel: [low, high]}.

    Returns
    -------
    callable
        ``f(study, detector, channel)``.
    """

    def lookup(study: str, detector: int, channel: int) -> list[float] | None:
        return table.get(study, {}).get(channel)

    return lookup


@dataclass(frozen=True)
class Phase:
    """
    Configuration of one analysis phase, independent of the others.

    Attributes
    ----------
    key : str
        ``"1"``, ``"2"`` or ``"3"``.
    arm : str
        File-name arm of the PROfit runs read.
    chi2_label : str
        Legend label of the chi-square.
    ylabel : str
        y-axis label of the panels.
    pvalue_ylabel : str
        y-axis label of the p-value summary figure.
    area_norm : bool
        True: ``ndf == nbins - 1`` is asserted; False: ``ndf == nbins``.
    enabled_tiers : tuple of int
        Tiers tabulated and plotted.
    overlay_tiers : tuple of int
        Tiers whose panels carry data points (None/empty: none).
    ylim : callable
        ``f(study, detector, channel)`` giving the y-limits.
    figure_arm : str
        Suffix in the written panel file names.
    """

    key: str
    arm: str
    chi2_label: str
    ylabel: str
    pvalue_ylabel: str
    area_norm: bool
    enabled_tiers: tuple[int, ...]
    overlay_tiers: tuple[int, ...]
    ylim: Callable = field(compare=False)
    figure_arm: str

    def studies(self) -> list[tuple[str, int, str]]:
        """
        Variables in the enabled tiers.

        Returns
        -------
        list of (study, tier, label)
        """
        return [(s, v[0], v[1]) for s, v in VARIABLES.items() if v[0] in self.enabled_tiers]

    def file(self, study: str) -> Path:
        """
        PROfit output path of a study in this phase.

        Parameters
        ----------
        study : str
            Variable name.

        Returns
        -------
        Path
            ``t<tier>_<study>_<arm>_PROplot_histograms.root``.
        """
        return RUNS_DIR / f"t{VARIABLES[study][0]}_{study}_{self.arm}_PROplot_histograms.root"


#: Phase 1 (area-normalised, shape-only), 2 (absolute rates) and 3 (absolute
#: rates with data overlaid). All three tiers are enabled in every phase (the
#: notebook's ``ENABLED_TIERS = (1, 2, 3)``); which tiers may show data points
#: is set separately by ``overlay_tiers``.
PHASES: dict[str, Phase] = {
    "1": Phase(
        key="1",
        arm="cvrw_cvrw",
        chi2_label=r"$\chi^2_{area}$",
        ylabel="Area Normalized",
        pvalue_ylabel="Shape-only p-value",
        area_norm=True,
        enabled_tiers=(1, 2, 3),
        overlay_tiers=(1,),
        ylim=_ylim_phase1,
        figure_arm="cvrw_cvrw",
    ),
    "2": Phase(
        key="2",
        arm="cvrw_cvrw_phase2",
        chi2_label=r"$\chi^2_{abs}$",
        ylabel="Events / bin width",
        pvalue_ylabel="Absolutely-normalized p-value",
        area_norm=False,
        enabled_tiers=(1, 2, 3),
        overlay_tiers=(),
        ylim=_ylim_per_channel(PHASE2_YLIM),
        figure_arm="cvrw_cvrw_phase2",
    ),
    "3": Phase(
        key="3",
        arm="cvrw_cvrw_phase2",
        chi2_label=r"$\chi^2_{abs}$",
        ylabel="Events / bin width",
        pvalue_ylabel="Absolutely-normalized p-value",
        area_norm=False,
        enabled_tiers=(1, 2, 3),
        overlay_tiers=(1, 2, 3),
        ylim=_ylim_per_channel(PHASE3_YLIM),
        figure_arm="phase3",
    ),
}


def channels_of(study: str) -> tuple[int, ...]:
    """
    Channel indices a variable defines (``energy_qel`` has one).

    Parameters
    ----------
    study : str
        Variable name.

    Returns
    -------
    tuple of int
    """
    return VARIABLES[study][2]


def plot_axes(study: str, channel: int) -> tuple[str, list[float] | None]:
    """
    x label and range for one panel, resolving a per-channel xlim.

    Parameters
    ----------
    study : str
        Variable name.
    channel : int
        Channel index.

    Returns
    -------
    xlabel : str
    xlim : list of float or None
    """
    xlabel, xlim = PLOT_AXES[study]
    if isinstance(xlim, dict):
        xlim = xlim.get(channel)
    return xlabel, xlim


def plot_ready(path: Path) -> bool:
    """
    True only if the run has finished writing its trees.

    A job still in flight leaves a file that exists but has no keys, so
    ``Path.exists()`` is not enough.

    Parameters
    ----------
    path : Path
        A ``*_PROplot_histograms.root`` file.

    Returns
    -------
    bool
    """
    import uproot

    if not path.exists():
        return False
    try:
        with uproot.open(path) as handle:
            keys = {k.split(";")[0] for k in handle.keys(recursive=False)}
        return {"hist1d", "chi2"} <= keys
    except Exception:
        return False


def sbnd_chi2(path: Path, channel: int, area_norm: bool) -> tuple[float, int, int, float]:
    """
    chi-square, nbins, ndf and p-value for one SBND channel of one run.

    The assertion checks the normalisation convention: area normalisation
    costs one degree of freedom (``ndf == nbins - 1``); without it
    ``ndf == nbins``. A file failing the assertion is a mislabelled run.
    It says nothing about whether ``shapeonly`` was set; that is read off
    the ini.

    Parameters
    ----------
    path : Path
        PROfit plot histograms file.
    channel : int
        Channel index.
    area_norm : bool
        Which convention to assert.

    Returns
    -------
    value, nbins, ndf, pvalue : float, int, int, float

    Raises
    ------
    AssertionError
        If ndf does not match the expected convention.
    """
    from mutools.plotting.profit import ProfitPlotData
    from scipy.stats import chi2 as chi2_dist

    value, nbins, ndf = ProfitPlotData(path).get_chi2(variable=0, detector=SBND, channel=channel)
    expected = nbins - 1 if area_norm else nbins
    assert ndf == expected, (
        f"{path.name}: ndf {ndf} != {expected}; "
        f"{'not area-normalised?' if area_norm else 'is this really not area-normalised?'}"
    )
    return value, nbins, ndf, float(chi2_dist.sf(value, ndf))


def chi2_table(phase: Phase) -> tuple[pd.DataFrame, list[str]]:
    """
    chi-square, ndf and p-value per variable and channel.

    Parameters
    ----------
    phase : Phase
        Which phase's runs and tiers to read.

    Returns
    -------
    table : pandas.DataFrame
        Columns tier, variable, channel, file, chi2, nbins, ndf, chi2/ndf, p-value.
    missing : list of str
        Studies whose run is absent or still in progress.
    """
    rows, missing = [], []
    for study, tier, label in phase.studies():
        path = phase.file(study)
        if not plot_ready(path):
            missing.append(f"{study} (tier {tier})")
            continue
        for channel in channels_of(study):
            value, nbins, ndf, p = sbnd_chi2(path, channel, phase.area_norm)
            rows.append(
                {
                    "tier": tier,
                    "variable": label,
                    "channel": CHANNELS[channel],
                    "file": path.name,
                    "chi2": value,
                    "nbins": nbins,
                    "ndf": ndf,
                    "chi2/ndf": value / ndf,
                    "p-value": p,
                }
            )
    return pd.DataFrame(rows), missing


def print_table(phase: Phase) -> pd.DataFrame:
    """
    Print the chi-square table and the comparisons below p = 0.05.

    For phase 2/3 ``ndf = nbins``, so the numbers are not comparable
    bin-for-bin with the phase-1 table.

    Parameters
    ----------
    phase : Phase
        Phase to tabulate.

    Returns
    -------
    pandas.DataFrame
        The table.
    """
    table, missing = chi2_table(phase)
    print(f"phase {phase.key}: ARM={phase.arm!r}  tiers={phase.enabled_tiers}  "
          f"data overlaid on={phase.overlay_tiers or 'nothing'}")
    if table.empty:
        print("no outputs on disk yet")
    else:

        def fmt(v: float) -> str:
            return f"{v:.2e}" if 0 < abs(v) < 1e-3 else f"{v:.3f}"

        shown = table if not phase.area_norm else table.drop(columns="nbins")
        with pd.option_context("display.width", 200, "display.max_colwidth", 52):
            print(shown.to_string(index=False, float_format=fmt))
        below = table[table["p-value"] < PVALUE_THRESHOLD]
        print(f"\n{len(below)} of {len(table)} comparisons below the {PVALUE_THRESHOLD} threshold")
        for _, r in below.iterrows():
            print(f"  {r['variable']} / {r['channel']}: p = {r['p-value']:.4f}")
    if missing:
        print("\nnot on disk yet: " + ", ".join(missing))
        if not phase.area_norm:
            print(
                "  run:  ../profit-asymmetric/build/bin/PROfit "
                "--config xml/sbnd_data_access/phase1-tier<N>/<var>_cvrw_phase2_plot.ini plot"
            )
    return table


def plot_pvalues(table: pd.DataFrame, phase: Phase):
    """
    Per-variable p-values, one marker per channel, with the 0.05 threshold.

    The document's off-ramp is any p-value below 0.05, so the threshold is
    drawn rather than left to the reader. Phase 2 uses the same threshold
    and log scale, but carries the normalisation, so a variable can sit
    lower than in phase 1 without its shape having changed.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of `chi2_table`.
    phase : Phase
        Phase (sets the y label).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis

    fig, ax = plt.subplots(figsize=(10, 5))
    shown = table.dropna(subset=["p-value"])
    labels = list(dict.fromkeys(shown["variable"]))
    x = np.arange(len(labels))
    for name, key in (("Single-proton", "ccqelike"), ("Multi-proton", "nproton")):
        subset = shown[shown["channel"] == key].set_index("variable").reindex(labels)
        ax.scatter(x, subset["p-value"], s=45, color=PVALUE_COLORS[key], label=name, zorder=3)
    ax.axhline(
        PVALUE_THRESHOLD, color=THRESHOLD_COLOR, lw=1.2, ls="--", zorder=2, label="p = 0.05"
    )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel(phase.pvalue_ylabel)
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(min(shown["p-value"].min() / 5, 1e-2), 2.0)
    ax.legend(loc="lower right")
    mark_axis(ax, r"$\bf{SBND}$ Data Access", vadj=0.05)
    mark_axis(ax, POT_LABEL, vadj=0.05, alignment="right")
    return fig


def _data_top(trace: np.ndarray) -> np.ndarray:
    """
    Data bin content plus its error (Poisson if the trace has no error column).

    Parameters
    ----------
    trace : np.ndarray
        A ``HIST_CONTENTS`` trace; column 3 is content, 4 the error.

    Returns
    -------
    np.ndarray
    """
    err = trace[:, 4] if trace.shape[1] > 4 else np.sqrt(np.maximum(trace[:, 3], 0.0))
    return trace[:, 3] + err


def report_peaks(phase: Phase) -> None:
    """
    Print each panel's peak height next to the y-limit in use.

    Lets the y-limit tables be set to real numbers instead of autoscale.
    For tiers outside the overlay set the data trace is deliberately NOT
    read, so a y-scale can be chosen without looking at withheld data;
    those rows are marked ``CV+band`` and may clip if the data sits higher.

    Parameters
    ----------
    phase : Phase
        Phase whose runs and y-limits are reported.
    """
    from mutools.plotting.profit import ProfitPlotData, TraceType

    print(
        f"{'study':24s}{'tier':>5}{'det':>4}{'ch':>4}{'peak':>10}{'basis':>14}{'ylim in use':>18}"
    )
    for study, tier, _label in phase.studies():
        path = phase.file(study)
        if not plot_ready(path):
            state = "IN PROGRESS" if path.exists() else "MISSING"
            print(f"{study:24s}{tier:5d}{state:>36}")
            continue
        d = ProfitPlotData(path, True)
        use_data = tier in phase.overlay_tiers
        for det in DETECTORS:
            for ch in channels_of(study):
                key = f"0:0:{det}:{ch}:total"
                try:
                    cv = d.get_trace(f"{key}:CV", TraceType.HIST_CONTENTS)
                    band = d.get_trace(f"{key}:CV", TraceType.HIST_ERROR_BAND)
                    stack = [cv[:, 3], band[:, 1] + band[:, 3]]
                    if use_data:
                        stack.append(_data_top(d.get_trace(f"{key}:DATA", TraceType.HIST_CONTENTS)))
                except Exception as exc:
                    print(f"{study:24s}{tier:5d}{det:4d}{ch:4d}   trace error: {exc}")
                    continue
                peak = float(np.nanmax(np.concatenate(stack)))
                basis = "CV+band+data" if use_data else "CV+band"
                print(
                    f"{study:24s}{tier:5d}{det:4d}{ch:4d}{peak:10.4f}{basis:>14}"
                    f"{phase.ylim(study, det, ch)!s:>18}"
                )


def check_clipping(phase: Phase) -> list[tuple[str, int, float, float]]:
    """
    Find data points above their configured y-limit.

    Run before rendering: a clipped point is easy to miss in a finished
    panel and reads as a missing marker rather than a disagreement.

    Parameters
    ----------
    phase : Phase
        Phase to check (only tiers in its overlay set have data drawn).

    Returns
    -------
    list of (study, channel, data_top, ylim_top)
        The offending panels.
    """
    from mutools.plotting.profit import ProfitPlotData, TraceType

    clipped = []
    for study, tier, _label in phase.studies():
        path = phase.file(study)
        if not plot_ready(path):
            print(f"{study}: MISSING")
            continue
        if tier not in phase.overlay_tiers:
            continue
        d = ProfitPlotData(path, True)
        for ch in channels_of(study):
            trace = d.get_trace(f"0:0:0:{ch}:total:DATA", TraceType.HIST_CONTENTS)
            top = float(np.nanmax(_data_top(trace)))
            lim = phase.ylim(study, 0, ch)
            if lim is not None and top > lim[1]:
                clipped.append((study, ch, top, lim[1]))
    if clipped:
        print("Data exceeds the limit -- raise these y-limit entries:")
        for study, ch, top, lim in clipped:
            print(
                f"  {study:22s} ch{ch}  data top {top:.4g} > ylim {lim:.4g}  "
                f"(suggest {1.1 * top:.4g})"
            )
    else:
        print(
            f"no clipping: every drawn data point fits inside its y-limit "
            f"({len(phase.studies())} variables, tiers {phase.enabled_tiers})"
        )
    return clipped


def panel_config(
    phase: Phase, path: Path, study: str, det: int, ch: int, overlay: bool, output: Path
) -> dict:
    """
    mutools ``run`` configuration for one data/MC histogram panel.

    ``ratio='data'`` is the only thing that draws data points in mutools,
    so ``ratio=None`` gives MC plus band with no points and no ratio
    panel; the chi-square and p-value still come from the file's chi2
    tree. The full detector and channel lists are passed because
    dispatch indexes them with the absolute detector/channel index.

    Parameters
    ----------
    phase : Phase
        Phase (labels, y-limits).
    path : Path
        PROfit histograms file.
    study : str
        Variable name.
    det, ch : int
        Detector and channel indices.
    overlay : bool
        Draw data points and the ratio panel.
    output : Path
        Directory mutools writes its staged file to.

    Returns
    -------
    dict
        Configuration for ``mutools.plotting.run``.
    """
    legend_loc, legend_ncol = PLOT_LEGEND.get(study, DEFAULT_LEGEND)
    xlabel, xlim = plot_axes(study, ch)
    return {
        "general": {
            "input": str(path),
            "output": str(output),
            "savefig": True,
            "code_version": common.CODE_VERSION,
            "selection_version": common.SELECTION_VERSION,
            "subchannels": SUBCHANNELS,
            "detectors": list(DETECTORS.values()),
            "channels": [CHANNEL_LABELS[i] for i in sorted(CHANNEL_LABELS)],
            "scale_by_width": "enabled",
            "counter_index": 4,  # Count is var4 in the coverage convention
        },
        "plot": [
            {
                "type": "histogram",
                "variable": 0,
                "channel": ch,
                "detectors": [det],
                "colors": SUBCHANNEL_COLORS,
                "xlabel": xlabel,
                "ylabel": phase.ylabel,
                "xlim": xlim,
                "ylim": phase.ylim(study, det, ch),
                "rlim": [0.0, 2.0],
                "ratio": "data" if overlay else None,
                "scale_by_width": "enabled",
                "disable_systematics": False,
                "band_label": "Total Uncertainty",
                "counter_fmt": ".1%",  # percentages, never yields
                "watermark": common.WATERMARK,
                "show_chi2": True,
                "chi2_label": phase.chi2_label,
                "show_pvalue": True,
                "legend_ncol": legend_ncol,
                "legend_loc": legend_loc,
            }
        ],
    }


def render_panels(
    phase: Phase, output: Path, fmt: str, only: Sequence[str] | None = None
) -> int:
    """
    Render the data/MC panels of one phase, skipping missing runs.

    mutools picks the staged file name itself (``hist_<i>_<j>_<version>.<fmt>``;
    the indices do not map cleanly onto detector/channel), so the newest
    staged file is moved to
    ``sbndaccess_<variable>_<detector>_<channel>_<figure_arm>.<fmt>``.

    Parameters
    ----------
    phase : Phase
        Phase to render.
    output : Path
        Figure directory.
    fmt : str
        Figure format; must match the saver configuration.
    only : sequence of str, optional
        Restrict to these variable names.

    Returns
    -------
    int
        Number of panels written.

    Raises
    ------
    FileNotFoundError
        If mutools wrote no staged file.
    """
    from mutools.plotting import run

    output.mkdir(parents=True, exist_ok=True)
    made, skipped, suppressed = 0, [], []
    for study, tier, _label in phase.studies():
        if only and study not in only:
            continue
        path = phase.file(study)
        if not plot_ready(path):
            skipped.append(f"{study}/tier{tier}" + (" (in progress)" if path.exists() else ""))
            continue
        overlay = tier in phase.overlay_tiers
        if not overlay:
            suppressed.append(f"{study} (tier {tier})")
        for det, det_label in DETECTORS.items():
            for ch in channels_of(study):
                run(panel_config(phase, path, study, det, ch, overlay, output))
                staged = sorted(
                    output.glob(f"hist_*_{common.CODE_VERSION}.{fmt}"),
                    key=lambda q: q.stat().st_mtime,
                )
                if not staged:
                    raise FileNotFoundError(
                        f"{study} det{det} ch{ch}: no staged .{fmt} appeared"
                    )
                name = (
                    f"sbndaccess_{study}_{det_label.replace(' ', '').lower()}"
                    f"_{CHANNEL_TOKENS[ch]}_{phase.figure_arm}"
                )
                shutil.move(staged[-1], output / f"{name}.{fmt}")
                made += 1
    print(f"{made} phase-{phase.key} panels rendered for tiers {phase.enabled_tiers}")
    if suppressed:
        print("data points withheld: " + ", ".join(suppressed))
    if skipped:
        print("not produced yet: " + ", ".join(skipped))
    return made


def configure(fmt: str) -> None:
    """
    Configure style and the mutools saver for these panels.

    Parameters
    ----------
    fmt : str
        Figure format. The saver must use the same extension that
        `render_panels` globs for (they once disagreed, which silently
        produced PDFs named ``.png``).
    """
    import mutools.plotting as mp

    common.configure_style(fmt)
    mp.saver.configure(fmt=fmt, dpi=300, rasterized=True, ytick_precision=2, mark_hadj=0.055)


def main(argv: Sequence[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : sequence of str, optional
        Arguments; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.add_argument(
        "--phase",
        choices=("1", "2", "3", "4", "5", "all"),
        default="1",
        help="analysis phase (see the module docstring); 4 and 5 are in data_access_phases",
    )
    parser.add_argument(
        "--what",
        nargs="+",
        choices=("table", "pvalues", "peaks", "clipping", "panels", "all"),
        default=["all"],
        help="parts of phases 1-3 to run (default: all)",
    )
    parser.add_argument("--only", nargs="+", help="restrict panels to these variables")
    parser.add_argument(
        "--tiers", type=int, nargs="+", help="override the phase's enabled tiers"
    )
    parser.add_argument(
        "--overlay-tiers",
        type=int,
        nargs="*",
        help="override the tiers that draw data points (empty list: none)",
    )
    parser.add_argument("--flip", action="store_true", help="phase 4: partner constrains E_nu")
    parser.add_argument("--alpha", type=float, default=1.0, help="phase 4: exposure scale")
    args = parser.parse_args(argv)

    every = {"table", "pvalues", "peaks", "clipping", "panels"}
    what = every if "all" in args.what else set(args.what)
    phases = ("1", "2", "3", "4", "5") if args.phase == "all" else (args.phase,)
    configure(args.format)
    for key in phases:
        if key in ("4", "5"):
            import data_access_phases

            if key == "4":
                data_access_phases.run_phase4(args.output, args.format, args.flip, args.alpha)
            else:
                data_access_phases.run_phase5(args.output, args.format)
            continue
        phase = PHASES[key]
        changes = {}
        if args.tiers:
            changes["enabled_tiers"] = tuple(args.tiers)
        if args.overlay_tiers is not None:
            changes["overlay_tiers"] = tuple(args.overlay_tiers)
        if changes:
            from dataclasses import replace

            phase = replace(phase, **changes)
        if "table" in what or "pvalues" in what:
            table = print_table(phase) if "table" in what else chi2_table(phase)[0]
            if "pvalues" in what and not table.empty:
                fig = plot_pvalues(table, phase)
                print(common.save_figure(fig, args.output / f"pvalues_{phase.arm}.{args.format}"))
        if "peaks" in what:
            report_peaks(phase)
        if "clipping" in what and phase.key == "3":
            check_clipping(phase)
        if "panels" in what:
            render_panels(phase, args.output, args.format, args.only)


if __name__ == "__main__":
    main()
