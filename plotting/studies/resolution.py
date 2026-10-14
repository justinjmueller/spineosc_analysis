"""
Energy and kinetic-energy resolution (relative bias) figures.

Replaces ``plotting/resolution.ipynb``. For SBND, ICARUS Run 2 and
ICARUS Run 4 it draws the distribution of the relative bias
``(reco - true) / true`` of the reconstructed quantities, from the PROfit
``resolution_{inclusive,exclusive}_v1_PROplot_histograms.root`` files in
``runs/legacy_root`` (MC only, statistical errors only).

Two kinds of figure:

``stacked``
    One figure per bias variable, detector and channel: the bias
    distribution stacked by final-state category (target and other nu_mu
    CC, nu_e CC, NC, out-of-time and in-time cosmics), drawn with the
    mutools TOML dispatcher. The y axis is sized to 1.5 times the tallest
    bin so the legend clears the stack. The variables are:

    * 4, 5: visible energy against true visible / true neutrino energy;
    * 6, 7: CCQE-estimator energy against true visible / true neutrino energy;
    * 8, 9: muon and leading-proton kinetic energy;
    * 10: CCQE-estimator against reconstructed visible energy. This one
      involves no truth, so it is computable in data too.

    The CCQE estimator assumes a single-proton final state, so variables
    6, 7 and 10 are dropped from the multi-proton channel.

``overlay``
    Several bias variables overlaid on one axis (translucent step
    histograms) so their shapes compare directly. Variables are grouped by
    scale: the energy biases (4, 5, 6, 7, 10) are broad, while the muon and
    proton kinetic-energy biases (8, 9) are far narrower and would sit as
    thin spikes beside them. ``overlay`` is called directly rather than
    through the dispatcher because the latter builds legend labels from the
    detector list instead of naming the variables.

Inputs: ``runs/legacy_root/resolution_<selection>_v1_PROplot_histograms.root``.

Figures: ``<output>/resolution/<selection>/ch<N>/hist_*`` (stacked) and
``<stage>/signal_and_selection/resolution_<group>_<channel>_<detector>``
(overlay), the latter staged for the technote.

Usage::

    python plotting/studies/resolution.py --selection exclusive --kind overlay
    python plotting/studies/resolution.py --selection inclusive --kind stacked --variable 4
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Selections: PROfit file stem and channel labels (``None`` for no channel split).
SELECTIONS: dict[str, tuple[str, list[str | None]]] = {
    "inclusive": ("resolution_inclusive_v1_PROplot_histograms.root", [None]),
    "exclusive": ("resolution_exclusive_v1_PROplot_histograms.root",
                  ["Single-Proton", "Multi-Proton"]),
}

#: Relative-bias variables, keyed by their variable index in the PROfit files.
BIAS_VARIABLES = {
    4: r"$(E^{reco}_{vis} - E^{true}_{vis}) / E^{true}_{vis}$",
    5: r"$(E^{reco}_{vis} - E^{true}_{\nu}) / E^{true}_{\nu}$",
    6: r"$(E^{reco}_{CCQE} - E^{true}_{vis}) / E^{true}_{vis}$",
    7: r"$(E^{reco}_{CCQE} - E^{true}_{\nu}) / E^{true}_{\nu}$",
    8: r"$(T^{reco}_{\mu} - T^{true}_{\mu}) / T^{true}_{\mu}$",
    9: r"$(T^{reco}_{p} - T^{true}_{p}) / T^{true}_{p}$",
    # Reco vs reco: no truth involved, so this one is computable in data too.
    10: r"$(E^{reco}_{CCQE} - E^{reco}_{vis}) / E^{reco}_{vis}$",
}

#: Stack labels. mutools pairs the FIRST label with the LAST subchannel in
#: the file, so this list runs opposite to the XML's declaration order of
#: intime, cos, nc, nuecc, othercc, target. Colours stay in XML order.
SUBCHANNELS = [r"Target $\nu_\mu$ CC", r"Other $\nu_\mu$ CC", r"$\nu_e$ CC", r"$\nu$ NC",
               "Out-of-time Cosmics", "In-time Cosmics"]
SUBCHANNEL_COLORS = ["C0", "C1", "C2", "C4", "C3", "C5"]

#: Biases built on the CCQE estimator (single-proton final state only).
CCQE_VARIABLES = {6, 7, 10}

#: Detector labels, file-name tags and target POT (from the XML ``<detector>``
#: blocks; these are MC-only configurations). Position is the detector index.
DETECTORS = ["SBND", "ICARUS Run 2", "ICARUS Run 4"]
DETECTOR_TAGS = ["sbnd", "icarus_run2", "icarus_run4"]
POT = ["1.00e20 POT", "2.00e20 POT", "3.00e20 POT"]

#: Channel tags used in overlay file names.
CHANNEL_TAGS = ["singleproton", "multiproton"]

#: Vertical headroom over the tallest bin, so the legend clears the stack.
HEADROOM = 1.5

#: Variables drawn together on one overlay axis, and their colours.
OVERLAY_GROUPS = {
    "energy": [4, 5, 6, 7, 10],
    "particle_ke": [8, 9],
}
OVERLAY_COLORS = {
    "energy": ["C0", "C1", "C2", "C3", "C4", "C5"],
    "particle_ke": ["C0", "C2"],
}

#: ``False`` draws outlines only, which keeps overlapping distributions
#: readable; ``OVERLAY_FILL_ALPHA`` sets the fill opacity otherwise.
OVERLAY_FILL = True
OVERLAY_FILL_ALPHA = 0.2

#: Watermark of the resolution plots (simulation only).
WATERMARK = r"$\bf{SBN}$ Simulation"

#: Staging directory for the overlay figures.
STAGE = common.DOCS / "signal_selection_v1.1.0" / "signal_and_selection"


def variables_for(channel_label: str | None, variables: list[int] | None = None) -> list[int]:
    """
    Bias variables to plot for one channel.

    The CCQE estimator assumes a single-proton final state, so its biases
    are meaningless in the multi-proton channel.

    Parameters
    ----------
    channel_label : str or None
        Channel label; ``"Multi-Proton"`` drops the CCQE biases.
    variables : list of int, optional
        Restrict to these variables; default is all of `BIAS_VARIABLES`.

    Returns
    -------
    list of int
        Variable indices.
    """
    chosen = list(BIAS_VARIABLES) if variables is None else list(variables)
    if channel_label == "Multi-Proton":
        return [v for v in chosen if v not in CCQE_VARIABLES]
    return chosen


def peak_content(contents, variable: int, detector: int, channel: int) -> float:
    """
    Tallest bin of the total CV histogram of one variable.

    Parameters
    ----------
    contents : mutools.plotting.profit.ProfitPlotData
        The opened histogram file.
    variable, detector, channel : int
        Indices into the file.

    Returns
    -------
    float
        Maximum bin content (column 3 of the stored trace).
    """
    from mutools.plotting.profit import TraceType

    trace = contents.get_trace(f"{variable}:0:{detector}:{channel}:total:CV",
                               TraceType.HIST_CONTENTS)
    return float(trace[:, 3].max())


def make_stacked(
    contents,
    input_file: Path,
    selection: str,
    channel_labels: list[str | None],
    detectors: list[int],
    variables: list[int] | None,
    output: Path,
) -> None:
    """
    Stacked bias distributions through the mutools dispatcher.

    Parameters
    ----------
    contents : mutools.plotting.profit.ProfitPlotData
        The opened histogram file, used to size each y axis.
    input_file : Path
        The PROfit histogram file.
    selection : str
        Selection name, used in the output sub-directory.
    channel_labels : list
        Channel labels of the file.
    detectors : list of int
        Detector indices to draw.
    variables : list of int, optional
        Restrict to these variables.
    output : Path
        Base directory; files go to ``<output>/<selection>/ch<N>``.
    """
    from mutools.plotting import run

    for channel, channel_label in enumerate(channel_labels):
        for variable in variables_for(channel_label, variables):
            for detector in detectors:
                run({
                    "general": {
                        "input": str(input_file),
                        "output": str(output / selection / f"ch{channel}"),
                        "savefig": True,
                        "code_version": common.CODE_VERSION,
                        "selection_version": common.SELECTION_VERSION,
                        "subchannels": SUBCHANNELS,
                        "detectors": DETECTORS,
                        "channels": channel_labels,
                        "counter_index": 3,  # the Count variable
                    },
                    "plot": [{
                        "type": "histogram",
                        "variable": variable,
                        "channel": channel,
                        "detectors": [detector],
                        "colors": SUBCHANNEL_COLORS,
                        "xlabel": BIAS_VARIABLES[variable],
                        "ylabel": "Events",
                        "xlim": [-1.0, 1.0],
                        "ylim": [0.0, HEADROOM * peak_content(contents, variable, detector,
                                                              channel)],
                        "disable_systematics": True,
                        "counter_fmt": ".1%",
                        "watermark": WATERMARK,
                        "right_watermark": POT[detector],
                        "legend_ncol": 2,
                        "legend_loc": "upper right",
                    }],
                })


def make_overlays(
    contents,
    channel_labels: list[str | None],
    detectors: list[int],
    variables: list[int] | None,
    stage: Path,
) -> list[Path]:
    """
    Overlay plots of several bias variables on one axis.

    ``overlay()`` names its file ``overlay_<channel>_<variables>_<version>``,
    which carries no detector, so every detector would write the same
    path. Each file is therefore renamed straight away, before the next
    call overwrites it.

    Parameters
    ----------
    contents : mutools.plotting.profit.ProfitPlotData
        The opened histogram file.
    channel_labels : list
        Channel labels of the file.
    detectors : list of int
        Detector indices to draw.
    variables : list of int, optional
        Restrict to these variables.
    stage : Path
        Directory for the renamed figures.

    Returns
    -------
    list of Path
        Files written.
    """
    import mutools.plotting as mp
    from mutools.plotting.profit import overlay

    written = []
    for group, group_variables in OVERLAY_GROUPS.items():
        for channel, channel_label in enumerate(channel_labels):
            allowed = set(variables_for(channel_label, variables))
            overlay_variables = [v for v in group_variables if v in allowed]
            if not overlay_variables:
                continue
            for detector in detectors:
                peak = max(peak_content(contents, v, detector, channel)
                           for v in overlay_variables)
                label = DETECTORS[detector]
                overlay(
                    contents,
                    variable=overlay_variables,
                    detectors=[detector] * len(overlay_variables),
                    channel=channel,
                    detector_labels=[BIAS_VARIABLES[v] for v in overlay_variables],
                    channel_label=label if channel_label is None else f"{label}\n{channel_label}",
                    xlabel="Relative Bias",
                    ylabel="Events",
                    xlim=[-1.0, 1.0],
                    ylim=[0.0, HEADROOM * peak],
                    colors=OVERLAY_COLORS[group][: len(overlay_variables)],
                    fill=OVERLAY_FILL,
                    fill_alpha=OVERLAY_FILL_ALPHA,
                    scale_by_width=False,
                    code_version=common.CODE_VERSION,
                    selection_version=common.SELECTION_VERSION,
                    watermark=WATERMARK,
                    right_watermark=POT[detector],
                    output=stage,
                )
                chname = CHANNEL_TAGS[channel]
                names = "-".join(str(v) for v in overlay_variables)
                auto = stage / f"overlay_{channel}_{names}_{common.CODE_VERSION}.{mp.saver.fmt}"
                tag = DETECTOR_TAGS[detector]
                new = stage / f"resolution_{group}_{chname}_{tag}.{mp.saver.fmt}"
                auto.rename(new)
                written.append(new)
                print(f"wrote {new}")
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
    parser.add_argument("--selection", choices=list(SELECTIONS), default="exclusive",
                        help="inclusive or exclusive (single/multi-proton) channels")
    parser.add_argument("--kind", choices=["stacked", "overlay", "all"], default="all",
                        help="which figures to make")
    parser.add_argument("--variable", type=int, action="append", choices=list(BIAS_VARIABLES),
                        help="restrict to this bias variable (repeatable)")
    parser.add_argument("--detector", type=int, action="append", choices=range(len(DETECTORS)),
                        help="restrict to this detector index (0 SBND, 1 Run 2, 2 Run 4)")
    parser.add_argument("--stage", type=Path, default=STAGE, help="directory for overlay figures")
    parser.add_argument("--runs", type=Path, default=common.LEGACY_ROOT,
                        help="directory with the PROfit histogram files")
    args = parser.parse_args(argv)

    import mutools.plotting as mp
    from mutools.plotting.profit import ProfitPlotData

    common.configure_style(fmt=args.format)
    mp.saver.configure(mark_hadj=0.055)

    stem, channel_labels = SELECTIONS[args.selection]
    input_file = args.runs / stem
    # The stored histograms are read once to size each y axis.
    contents = ProfitPlotData(input_file)
    detectors = args.detector if args.detector else list(range(len(DETECTORS)))

    if args.kind in ("stacked", "all"):
        make_stacked(contents, input_file, args.selection, channel_labels, detectors,
                     args.variable, args.output / "resolution")
    if args.kind in ("overlay", "all"):
        args.stage.mkdir(parents=True, exist_ok=True)
        make_overlays(contents, channel_labels, detectors, args.variable, args.stage)


if __name__ == "__main__":
    main()
