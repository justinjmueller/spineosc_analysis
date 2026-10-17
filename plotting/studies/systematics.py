"""
Systematic-uncertainty breakdowns per detector, channel and energy variable.

Makes the fractional-uncertainty line plots of the systematics chapter
of the technote, in reconstructed and true neutrino energy, and the
near/far ratio versions of the same plots. Replaces
``plotting/systematics.ipynb``.

Inputs are PROfit ``plot`` outputs of ``jobs/systematics/systs_reco.ini``
(``--fit-variable 0``, Reco E) and ``systs_true.ini`` (``--fit-variable 2``,
True E), both reading ``xml/systematics/systematics_suite.xml``. That
XML is generated from ``xml/contours/contours_recoe.xml`` by
``scripts/systematics_suite_from_recoe.py``, so it carries the current
analysis (cuts, CV-reweighted inputs, systematic suite). PROfit only
writes the fractional-systematic tree for the fitting variable. The
histogram files are read from ``runs/<run>/<tag>_{recoE,trueE}_PROplot_histograms.root``.

Broad breakdown (``summary``, one line per type plus a quadrature total):

* Flux: Flux + G4 (beamline Geant4 reinteractions)
* Cross Section: every ``Cross-Section-*`` category plus the nuclear-model
  Binding Energy systematic (tag ``Nuclear-Model``)
* Uncorrelated Systematics: Detector Systematics + POT Normalization + CORSIKA
* Prediction Stat.: MC statistics

Narrow breakdowns (``flux``, ``xsec``, ``geant``, ``uncorrelated``) show
each type's components with their quadrature total. Prediction stat. is
a single component, so it has no narrow breakdown.

Totals are drawn from the plotted components (``show_total``), NOT PROfit's
TOTAL entry: the mcstat allowlist carries two tags ("stat,MC Stats") and
PROfit's TOTAL counts that one covariance twice.

The ratio figures (``--ratio``) show the uncertainty on the ratio of the
SBND (near) to the ICARUS Run 2 or Run 4 (far) prediction in each bin,
from the ``use_ratio`` entries (``detector``/``detector2``) of the PROfit
fractional-systematics tree. Only categories whose tags are common to both
detectors are drawn (the uncorrelated per-detector components have
different names and are summarized by the Uncorrelated Systematics line).
Cancellation between the detectors reduces the correlated flux and
interaction-model uncertainties; the uncorrelated ones do not cancel.

Y-axis ranges and legend placement are hand-set module constants
(``YLIM``, ``YLIM_OVERRIDES``, ``LEGEND``); run ``--report`` to print the
largest plotted value per category and choose ``YLIM`` from it.

Example::

    python plotting/studies/systematics.py --run systematics_v2 --ratio
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import provenance

import common
from common import CODE_VERSION, DOCS, RUNS, SELECTION_VERSION

#: Default run directory (under ``runs/``) holding the PROfit outputs of
#: ``jobs/systematics/systs_{reco,true}.ini``. Use another run to plot an
#: alternative scenario, e.g. the end-point-cut variant.
DEFAULT_RUN = "systematics_v2"

#: Default file-name prefix of the PROfit histogram files.
DEFAULT_TAG = "systs"

#: Default staging directory name (under ``docs/``). New figures are
#: written here first; copy to the technote once the set has been checked.
DEFAULT_STAGING = "systematics_v1.1.0v2"

#: Energy variables: token -> (file suffix, x label, [single-proton xlim,
#: multi-proton xlim]).
VARIABLES = {
    "recoe": (
        "recoE",
        "Reconstructed Neutrino Energy [GeV]",
        [[0.3, 1.5], [0.4, 2.0]],
    ),
    "truee": (
        "trueE",
        "True Neutrino Energy [GeV]",
        [[0.3, 3.0], [0.3, 3.0]],
    ),
}

#: Detectors, in PROfit index order: (legend label, file token, suffix of
#: this detector's systematic names).
DETECTORS = [
    ("SBND Simulation", "sbnd_run1", "sbnd_Run1"),
    ("ICARUS Run 2 Simulation", "icarus_run2", "icarus_Run2"),
    ("ICARUS Run 4 Simulation", "icarus_run4", "icarus_Run4"),
]

#: Analysis channels, in PROfit index order: (legend label, file token).
CHANNELS = [("Single-Proton", "singlep"), ("Multi-Proton", "multip")]

#: Subchannel labels passed to the mutools plot configuration.
SUBCHANNELS = [
    "Target $\\nu_\\mu$ CC",
    "Other $\\nu_\\mu$ CC",
    "$\\nu_e$ CC",
    "$\\nu$ NC",
    "Out-of-time Cosmics",
    "In-time Cosmics",
]

#: Cross-section tags -> legend labels. ``Nuclear-Model`` is the binding-energy
#: shift of the struck nucleons (hist1d, shared across detectors); it is grouped
#: with the interaction model. Move it into ``uncorrelated_components`` instead
#: if it should count as an uncorrelated systematic.
XSEC = {
    "Cross-Section-QE": "QE",
    "Cross-Section-MEC": "MEC",
    "Cross-Section-RES": "RES",
    "Cross-Section-DIS": "DIS",
    "Cross-Section-FSI": "FSI",
    "Cross-Section-COH": "COH",
    "Cross-Section-NCEL": "NC Elastic",
    "Cross-Section-NonRES-nu": r"Non-RES $\nu$",
    "Cross-Section-NonRES-nubar": r"Non-RES $\bar{\nu}$",
    "Nuclear-Model": "Binding Energy Shift",
}

#: Flux tags -> legend labels.
FLUX = {
    "Flux:piplus_Flux": r"$\pi^{+}$ Production",
    "Flux:piminus_Flux": r"$\pi^{-}$ Production",
    "Flux:kplus_Flux": r"$K^{+}$ Production",
    "Flux:kminus_Flux": r"$K^{-}$ Production",
    "Flux:kzero_Flux": r"$K^{0}$ Production",
    "Flux:nucleoninexsec_Flux": "Nucleon Inelastic XS",
    "Flux:nucleonqexsec_Flux": "Nucleon QE XS",
    "Flux:nucleontotxsec_Flux": "Nucleon Total XS",
    "Flux:pioninexsec_Flux": "Pion Inelastic XS",
    "Flux:pionqexsec_Flux": "Pion QE XS",
    "Flux:piontotxsec_Flux": "Pion Total XS",
    "Flux:horncurrent_Flux": "Horn Current",
    "Flux:expskin_Flux": "Skin Depth",
}

#: Beamline Geant4 reinteraction tags -> legend labels.
GEANT = {
    "G4:reinteractions_proton_Geant4": "G4 Proton Reinteractions",
    "G4:reinteractions_piplus_Geant4": r"G4 $\pi^{+}$ Reinteractions",
    "G4:reinteractions_piminus_Geant4": r"G4 $\pi^{-}$ Reinteractions",
}

#: Y-axis range ``[ymin, ymax]`` per breakdown category, applied to EVERY
#: detector and to BOTH energy variables so that plots which appear side by
#: side in the technote share an axis. The trailing comments give the largest
#: value found in the files when the ranges were set (see ``--report``).
YLIM = {
    "summary": [0.0, 0.80],  # 0.55
    "flux": [0.0, 0.50],  # 0.29
    "geant": [0.0, 0.25],  # 0.17
    "xsec": [0.0, 0.65],  # 0.44 (includes Binding Energy)
    "uncorrelated": [0.0, 0.24],  # 0.148 (SBND multi-proton Reco E; includes DENT)
}

#: Optional exceptions, for when one case genuinely should not share the
#: range in ``YLIM``. Keys are ``(category,)``, ``(category, variable)`` or
#: ``(category, variable, channel)``; the most specific matching key wins.
#: For example ``("uncorrelated", "truee"): [0.0, 0.06]`` or
#: ``("summary", "recoe", "singlep"): [0.0, 0.5]``.
YLIM_OVERRIDES: dict[tuple[str, ...], list[float]] = {}

#: Legend corner and column count per breakdown category.
LEGEND = {
    "summary": ("upper left", 1),
    "flux": ("upper left", 2),
    "geant": ("upper left", 1),
    "xsec": ("upper left", 2),
    "uncorrelated": ("upper left", 2),
}

#: Legend placement used by any category not listed in ``LEGEND``.
DEFAULT_LEGEND = ("upper left", 2)

#: Near/far pairs for the ratio figures: (numerator detector index,
#: denominator detector index, file token, legend label).
RATIO_PAIRS = [
    (0, 1, "sbnd_run1_over_icarus_run2", "SBND / ICARUS Run 2"),
    (0, 2, "sbnd_run1_over_icarus_run4", "SBND / ICARUS Run 4"),
]

#: Y-axis range per category for the ratio figures (only categories with
#: tags common to both detectors).
RATIO_YLIM = {
    "summary": [0.0, 0.5],
    "flux": [0.0, 0.3],
    "xsec": [0.0, 0.4],
    "geant": [0.0, 0.2],
}


def uncorrelated_components(suffix: str) -> dict[str, str]:
    """
    Detector-model, POT and (SBND only) CORSIKA components of one detector.

    Parameters
    ----------
    suffix : str
        The detector's systematic-name suffix, e.g. ``"sbnd_Run1"``.

    Returns
    -------
    dict[str, str]
        Systematic tag -> legend label.
    """
    det = "Detector Systematics"
    sbnd = suffix.startswith("sbnd")
    components = {}
    if sbnd:
        components[f"{det}:ChargeScale15_{suffix}"] = "Charge Scale / Resolution"
    else:
        components[f"{det}:ChargeScale_{suffix}"] = "Charge Scale / Resolution"
        components[f"{det}:SCE_{suffix}"] = "SCE"
    components[f"{det}:TriggerEmulation_{suffix}"] = "Trigger Emulation"
    components[f"{det}:WireModxThetaXW_{suffix}"] = r"TPC Signal Shape (WireMod X-$\theta_{xw}$)"
    if sbnd:
        # SBND uses the MattMod non-uniformity variation in place of WireMod
        # Y-Z (the latter is commented out of the XML so that one physical
        # effect is not counted twice); ICARUS keeps WireMod Y-Z.
        components[f"{det}:MattModYZ_{suffix}"] = "Non-uniformity"
        components[f"{det}:DENT_{suffix}"] = "DENT"
    else:
        components[f"{det}:WireModYZ_{suffix}"] = "Non-uniformity (WireMod Y-Z)"
    components[f"POT Normalization:ToroidCalibration_{suffix}"] = "Toroid Calibration"
    components[f"POT Normalization:ProtonMiss_{suffix}"] = "Protons Missing Target"
    if sbnd:
        components["Cosmic Background Normalization"] = "CORSIKA Normalization"
    return components


def breakdowns(suffix: str) -> dict[str, tuple[list[str], dict | None, dict[str, str]]]:
    """
    All breakdown categories for one detector.

    Parameters
    ----------
    suffix : str
        The detector's systematic-name suffix, e.g. ``"icarus_Run2"``.

    Returns
    -------
    dict
        Category name -> ``(tags, groups, label_map)`` as consumed by the
        mutools ``uncertainty`` plot; ``groups`` is None except for the
        ``summary`` category.
    """
    broad = {
        "Flux": ["Flux", "G4"],
        "Cross Section": list(XSEC),
        "Uncorrelated Systematics": [
            "Detector Systematics",
            "POT Normalization",
            "Cosmic Background Normalization",
        ],
    }
    unc = uncorrelated_components(suffix)
    return {
        "summary": (list(broad) + ["stat"], broad, {"stat": "Prediction Stat."}),
        "flux": (list(FLUX), None, dict(FLUX)),
        "xsec": (list(XSEC), None, dict(XSEC)),
        "geant": (list(GEANT), None, dict(GEANT)),
        "uncorrelated": (list(unc), None, dict(unc)),
    }


def ylim_for(category: str, variable: str, channel: str) -> list[float] | None:
    """
    Y-axis range for one plot: the most specific override, else ``YLIM``.

    Parameters
    ----------
    category : str
        Breakdown category.
    variable : str
        Energy-variable token (``"recoe"`` or ``"truee"``).
    channel : str
        Channel token (``"singlep"`` or ``"multip"``).

    Returns
    -------
    list[float] or None
        ``[ymin, ymax]``; None (matplotlib autoscale) when the category
        has no entry, in which case a warning is printed because its
        panels will not agree with each other.
    """
    for key in ((category, variable, channel), (category, variable), (category,)):
        if key in YLIM_OVERRIDES:
            return YLIM_OVERRIDES[key]
    if category not in YLIM:
        print(f"WARNING: no YLIM entry for {category!r}; autoscaling (ranges will not match)")
        return None
    return YLIM[category]


def input_file(run_dir: Path, tag: str, variable: str) -> Path:
    """
    Path of the PROfit histogram file for one energy variable.

    Parameters
    ----------
    run_dir : Path
        Run directory holding the PROfit outputs.
    tag : str
        File-name prefix of the PROfit output.
    variable : str
        Energy-variable token, a key of ``VARIABLES``.

    Returns
    -------
    Path
        ``<run_dir>/<tag>_<recoE|trueE>_PROplot_histograms.root``.
    """
    return run_dir / f"{tag}_{VARIABLES[variable][0]}_PROplot_histograms.root"


def category_maxima(run_dir: Path, tag: str) -> None:
    """
    Print the largest value plotted in each category.

    Covers components and their total over the x-range actually shown.
    Use it to choose the ``YLIM`` entries: pick a ymax comfortably above
    the "all" value so the legend has room.

    Parameters
    ----------
    run_dir : Path
        Run directory holding the PROfit outputs.
    tag : str
        File-name prefix of the PROfit output.
    """
    from mutools.plotting.profit import ProfitPlotData, TraceType

    rows: dict[str, dict[tuple[str, str], list[float]]] = {}
    for var_token, (_, _xlabel, xlims) in VARIABLES.items():
        data = ProfitPlotData(input_file(run_dir, tag, var_token))
        for detector, (_label, _token, suffix) in enumerate(DETECTORS):
            for channel, (_ch_label, ch_token) in enumerate(CHANNELS):
                lo, hi = xlims[channel]
                for name, (tags, groups, _lbl) in breakdowns(suffix).items():
                    values = []
                    for tag_name in tags:
                        members = groups[tag_name] if groups and tag_name in groups else [tag_name]
                        sq = None
                        for member in members:
                            t, _, sysname = member.partition(":")
                            trace = f"0:{detector}:{channel}:{detector}:{t}:{sysname or 'SUM'}"
                            a = data.get_trace(trace, TraceType.FRAC_SYST)
                            inx = (a[:, 2] > lo) & (a[:, 1] < hi)
                            v = a[inx, 3]
                            sq = v**2 if sq is None else sq + v**2
                        values.append(np.sqrt(sq))
                    total = np.sqrt(sum(v**2 for v in values))
                    peak = max(total.max(), max(v.max() for v in values))
                    rows.setdefault(name, {}).setdefault((var_token, ch_token), []).append(peak)
    for name, per_case in rows.items():
        overall = max(max(v) for v in per_case.values())
        print(f"{name:14s} all: {overall:.3f}   ylim in use: {YLIM.get(name, 'not set')}")
        for (var_token, ch_token), peaks in per_case.items():
            rounded = [round(float(p), 3) for p in peaks]
            print(f"    {var_token:6s} {ch_token:8s} per detector: {rounded}")


def make_config(
    infile: Path,
    output: Path,
    detector: int,
    channel: int,
    xlabel: str,
    xlim: list[float],
    ylim: list[float] | None,
    tags: list[str],
    groups: dict | None,
    label_map: dict[str, str],
    legend_loc: str,
    legend_ncol: int,
) -> dict:
    """
    Configuration dictionary for ONE uncertainty plot.

    Parameters
    ----------
    infile : Path
        PROfit histogram file.
    output : Path
        Directory the mutools plotter writes into.
    detector, channel : int
        Detector and channel index (see ``DETECTORS``, ``CHANNELS``).
    xlabel : str
        X-axis label.
    xlim, ylim : list[float]
        Axis ranges; ``ylim`` may be None for autoscale.
    tags : list[str]
        Systematic tags drawn as lines.
    groups : dict or None
        Tag -> member tags added in quadrature.
    label_map : dict[str, str]
        Tag -> legend label.
    legend_loc : str
        Legend corner.
    legend_ncol : int
        Legend column count.

    Returns
    -------
    dict
        Configuration accepted by ``mutools.plotting.run``.
    """
    from matplotlib import colormaps

    # Explicit qualitative palettes: the style's default cycle includes faint
    # greys that are hard to tell apart on a many-line plot.
    cmap = colormaps["tab10" if len(tags) <= 10 else "tab20"]
    plot = {
        "type": "uncertainty",
        "detectors": [detector],
        "channel": channel,
        "tags": tags,
        "show_total": True,
        "xlabel": xlabel,
        "ylabel": "Fractional Uncertainty",
        "xlim": xlim,
        "ylim": ylim,
        "watermark": r"$\bf{SBN}$ Internal",
        "legend_title": f"{DETECTORS[detector][0]}\n{CHANNELS[channel][0]}",
        "label_map": {**label_map, "Total": "Total Uncertainty"},
        "legend_loc": legend_loc,
        "legend_ncol": legend_ncol,
        "colors": [cmap(i) for i in range(len(tags))],
    }
    if groups:
        plot["groups"] = groups
    return {
        "general": {
            "input": str(infile),
            "output": str(output),
            "savefig": True,
            "code_version": CODE_VERSION,
            "selection_version": SELECTION_VERSION,
            "subchannels": SUBCHANNELS,
            "detectors": [d[0] for d in DETECTORS],
            "channels": [c[0] for c in CHANNELS],
        },
        "plot": [plot],
    }


def make_breakdowns(run_dir: Path, tag: str, out_dir: Path) -> int:
    """
    Render every breakdown: detectors x channels x variables x categories.

    Parameters
    ----------
    run_dir : Path
        Run directory holding the PROfit outputs.
    tag : str
        File-name prefix of the PROfit output.
    out_dir : Path
        Destination directory.

    Returns
    -------
    int
        Number of figures written.
    """
    from mutools.plotting import run

    out_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for var_token, (_, xlabel, xlims) in VARIABLES.items():
            infile = input_file(run_dir, tag, var_token)
            for detector, (_label, det_token, suffix) in enumerate(DETECTORS):
                for channel, (_ch_label, ch_token) in enumerate(CHANNELS):
                    for name, (tags, groups, label_map) in breakdowns(suffix).items():
                        loc, ncol = LEGEND.get(name, DEFAULT_LEGEND)
                        run(
                            make_config(
                                infile, tmp, detector, channel, xlabel, xlims[channel],
                                ylim_for(name, var_token, ch_token),
                                tags, groups, label_map, loc, ncol,
                            )
                        )
                        plt.close("all")
                        stem = f"systematics_{det_token}_{name}_{var_token}_{ch_token}"
                        target = out_dir / f"{stem}.pdf"
                        shutil.move(tmp / f"uncertainty_{detector}_{CODE_VERSION}.pdf", target)
                        provenance.record(target, [infile], __file__, versions={"selection": common.SELECTION_VERSION},
                                          arguments={"kind": "breakdown", "detector": det_token, "category": name,
                                                     "variable": var_token, "channel": ch_token})
                        count += 1
    return count


def make_ratio_breakdowns(run_dir: Path, tag: str, out_dir: Path) -> int:
    """
    Render the near/far ratio uncertainty figures.

    Parameters
    ----------
    run_dir : Path
        Run directory holding the PROfit outputs.
    tag : str
        File-name prefix of the PROfit output.
    out_dir : Path
        Destination directory (the ``ratio`` sub-directory of the staging area).

    Returns
    -------
    int
        Number of figures written.

    Raises
    ------
    RuntimeError
        If the plotter wrote no figure for a configuration.
    """
    from mutools.plotting import run

    out_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for var_token, (_, xlabel, xlims) in VARIABLES.items():
            infile = input_file(run_dir, tag, var_token)
            for num, den, token, pair_label in RATIO_PAIRS:
                for channel, (ch_label, ch_token) in enumerate(CHANNELS):
                    for name, ylim in RATIO_YLIM.items():
                        tags, groups, label_map = breakdowns(DETECTORS[num][2])[name]
                        loc, ncol = LEGEND.get(name, DEFAULT_LEGEND)
                        cfg = make_config(
                            infile, tmp, num, channel, xlabel, xlims[channel],
                            ylim, tags, groups, label_map, loc, ncol,
                        )
                        # PROfit's `detector2` (denominator) entries; the
                        # numerator is `detector`.
                        cfg["plot"][0]["use_ratio"] = [den]
                        cfg["plot"][0]["legend_title"] = f"{pair_label}\n{ch_label}"
                        cfg["plot"][0]["ylabel"] = "Fractional Uncertainty on the Ratio"
                        run(cfg)
                        plt.close("all")
                        produced = sorted(tmp.glob("uncertainty_*.pdf"))
                        if not produced:
                            raise RuntimeError(f"no figure for {token} {name} {var_token}")
                        stem = f"systematics_ratio_{token}_{name}_{var_token}_{ch_token}"
                        target = out_dir / f"{stem}.pdf"
                        shutil.move(produced[-1], target)
                        provenance.record(target, [infile], __file__, versions={"selection": common.SELECTION_VERSION},
                                          arguments={"kind": "ratio", "pair": token, "category": name,
                                                     "variable": var_token, "channel": ch_token})
                        count += 1
    return count


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str], optional
        Argument list; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.set_defaults(output=DOCS / DEFAULT_STAGING)
    parser.add_argument("--run", default=DEFAULT_RUN, help="run directory name under runs/")
    parser.add_argument("--tag", default=DEFAULT_TAG, help="PROfit output file-name prefix")
    parser.add_argument("--report", action="store_true",
                        help="print the largest plotted value per category and exit")
    parser.add_argument("--ratio", action="store_true",
                        help="also make the near/far ratio figures (in <output>/ratio)")
    parser.add_argument("--ratio-only", action="store_true",
                        help="make only the near/far ratio figures")
    args = parser.parse_args(argv)

    run_dir = RUNS / args.run
    if args.report:
        category_maxima(run_dir, args.tag)
        return

    common.configure_style(fmt=args.format)
    import mutools.plotting as mp

    mp.saver.configure(mark_hadj=0.005)
    if not args.ratio_only:
        n = make_breakdowns(run_dir, args.tag, args.output)
        print(f"wrote {n} figures to {args.output}")
    if args.ratio or args.ratio_only:
        n = make_ratio_breakdowns(run_dir, args.tag, args.output / "ratio")
        print(f"wrote {n} ratio figures to {args.output / 'ratio'}")


if __name__ == "__main__":
    main(sys.argv[1:])
