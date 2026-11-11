"""
SBND -> ICARUS conditional constraint on the systematic covariance.

How much of the ICARUS systematic uncertainty survives the SBND
measurement, in which directions, and why -- computed from PROfit's
PROplot output through ``mutools.profit``. Replaces the notebook
``plotting/conditional_constraint.ipynb``.

Conventions
-----------
Everything is *fractional*: a shift :math:`\\delta` means
:math:`N_i \\to N_i(1+\\delta_i)`, and the fractional statistical variance
of a Poisson bin with expectation :math:`N` is :math:`1/N`.

The model
---------
Nuisance knobs :math:`\\theta\\sim\\mathcal N(0,I)`; SBND fractional
residual :math:`r = A\\theta + \\varepsilon` with
:math:`\\varepsilon\\sim\\mathcal N(0,S)`, :math:`S=\\mathrm{diag}(1/N)`;
ICARUS shift :math:`B\\theta`. The surviving ICARUS covariance is the
Gaussian conditional

.. math::  R = C_{yy} - C_{yx}\\bigl(C_{xx} + S/\\alpha\\bigr)^{-1}C_{xy},

equivalently :math:`R = B\\Lambda B^{\\mathsf T}` with
:math:`\\Lambda = (I + \\alpha A^{\\mathsf T}S^{-1}A)^{-1}` the post-SBND
knob covariance -- so the conditional *is* the profiled linear fit.

Both ICARUS runs are conditioned together. Run 2 and Run 4 are separate
blocks of the covariance, but one near detector constrains them jointly;
conditioning them in a single step keeps the Run 2 / Run 4 correlations
instead of discarding them.

Input
-----
``runs/legacy_root/contours_recoe_covar_PROplot.root`` (override with
``--input``). It needs the ``Covariance`` directory, which PROfit writes
only when asked::

    PROfit -x xml/contours/contours_recoe.xml -t contours_recoe -o covar -n 8 \\
           plot --with-covar

The loader infers the covariance bin ordering by rebuilding the file's
own ``collapsed_total_frac_cov``, so a successful load is already a
consistency check. Channels may be binned differently --
``contours_recoe`` bins ``ccqelike`` in 16 reconstructed-energy bins and
``nproton`` in 14 -- which is why bin counts are per channel. The fit
variable must be 1D: for a 2D channel PROfit writes only the x
projection per subchannel, which cannot fill the covariance space
(loading a ``contours_base`` file fails with exactly that diagnostic).

Figures (one per ``--only`` name, written to ``--output``)
----------------------------------------------------------
``constraint``
    ICARUS fractional uncertainty per bin before and after SBND.
``family_bars``
    Each systematic family on its own, before and after the constraint.
    A family whose two bars are nearly equal is one the near detector
    cannot see; one whose second bar collapses is constrained away.
``family_removed``
    The full budget with one family removed; the gap to the total is
    what that family contributes *after* the constraint. Removing
    CORSIKA can lower the curve: it is near-detector noise, so dropping
    it sharpens the constraint.
``detector``
    Breakdown of the detector family, each systematic alone.
``exposure``
    Scan of :math:`\\alpha`, which scales the SBND statistics,
    :math:`S\\to S/\\alpha`. As :math:`\\alpha\\to\\infty` the posterior
    tends to the projector onto knob directions SBND cannot see at all,
    so the limiting curve is a degeneracy floor, not a statistical one:
    bins already on it want a new observable, not more POT.
``modes``
    The residual is low rank. Each :math:`\\sqrt{\\lambda_a}\\,v_a` is a
    one-sigma fractional shift of the ICARUS prediction and the
    participation ratio :math:`n_\\mathrm{eff}=(\\sum\\lambda)^2/\\sum\\lambda^2`
    counts the modes that matter.
``knobs``
    Per-knob budget. For one knob, Sherman-Morrison gives the ceiling
    :math:`\\varepsilon^\\mathrm{max}_k = g_k/(1+g_k)` with visibility
    :math:`g_k = a_k^{\\mathsf T}S^{-1}a_k` (the :math:`\\chi^2` with which
    SBND alone distinguishes a :math:`1\\sigma` excursion). The achieved
    efficiency in the joint fit is :math:`\\varepsilon_k = 1-\\Lambda_{kk}`;
    the gap is constraint lost to *degeneracy* with other knobs rather
    than to statistics. Building the responses eigen-decomposes every
    selected component, so this is the slow figure.

The systematic families
-----------------------
Fractional component covariances add linearly, so any ``select`` /
``exclude`` subset is an exact partition of the *covariance*. The
*constraint* is not linear in that partition -- conditioning inverts
:math:`C_{xx}+S/\\alpha` -- so a family in isolation and a family removed
from the total answer two different questions, and neither set of
numbers sums to the total; the summary table reports both. The families
partition every component in the file and :func:`check_partition`
asserts it, so a systematic added to the XML later cannot silently fall
outside the breakdown.

On the v1.1.0 file the constraint works on the families *shared*
between the detectors and does nothing to the rest:

* Cross section, flux and G4 are correlated near-to-far and SBND
  removes most of them (isolation ratios 0.15, 0.14, 0.20).
* Detector, POT normalisation and MC stat come out with ratio exactly
  1.000. Each of those components is declared per detector and run
  (``ChargeScale_icarus_Run2``, ``ToroidCalibration_sbnd_Run1``, ...),
  so there is no near/far correlation for the conditional to exploit.
* CORSIKA is a cosmic-rate normalisation on SBND alone, so its ICARUS
  prior is exactly zero. It enters only as near-detector noise, and
  removing it improves the total constraint by 4e-5 (negligible).

Since the detector family survives untouched, its breakdown is about
*size*. Two asymmetries come from the XML, not the physics: ``SCE``
carries only the two ICARUS entries and ``DENT`` only the SBND one (the
SBND SCE and Gain entries are commented out), and a SBND-only
systematic has zero ICARUS prior by construction -- DENT's curve sits on
zero for that reason, not because it is small.

Usage::

    python plotting/studies/conditional_constraint.py --only constraint family_bars
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Default PROplot file with the ``Covariance`` directory.
INPUT = common.LEGACY_ROOT / "contours_recoe_covar_PROplot.root"

#: Near detector and the far detectors that are conditioned together.
NEAR = "SBND"
FAR = ("ICARUSRun2", "ICARUSRun4")

#: Fixed series order, checked for colourblind separation. Assigned in
#: order and never cycled; the yellow slot is skipped because it is too
#: low-contrast for a thin line on white.
SERIES_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#e87ba4", "#008300")

#: Systematic families, as regexes over the component names. They must
#: tile the component list exactly (see :func:`check_partition`).
FAMILIES = {
    "Cross section": r"^(GENIEReWeight|ZExpPCA|CCQETemplate|MECq0q3|QEInterference)",
    "Flux": r"_Flux$",
    "Detector": r"^(ChargeScale|SCE|TriggerEmulation|WireMod|DENT|Gain)",
    "POT norm": r"^(ToroidCalibration|ProtonMiss)",
    "G4": r"^reinteractions",
    "CORSIKA": r"cos:",
    "MC stat": r"mcstat",
}

#: Families drawn in the "family removed" figure.
REMOVED_SHOWN = ("Cross section", "Flux", "Detector", "POT norm", "CORSIKA")

#: Individual detector systematics, as regexes over the component names.
DETECTOR_SOURCES = {
    "Charge scale": r"^ChargeScale",
    "Space charge": r"^SCE",
    "Trigger emulation": r"^TriggerEmulation",
    "WireMod x-thetaXW": r"^WireModxThetaXW",
    "WireMod YZ": r"^WireModYZ",
    "DENT": r"^DENT",
}

#: Statistics-scale factors of the exposure scan (the infinite-exposure
#: limit is added separately with ``statistics=False``).
ALPHAS = (0.25, 1.0, 4.0)

#: Number of leading eigenmodes drawn, and knobs shown in the budget.
N_MODES = 3
N_KNOBS = 20

#: Figure names, in the order they are made.
FIGURES = (
    "constraint",
    "family_bars",
    "family_removed",
    "detector",
    "exposure",
    "modes",
    "knobs",
)

BIN_LABEL = "ICARUS bin (Run 2 then Run 4)"
LEGEND_BOX = {"frameon": True, "framealpha": 0.85, "edgecolor": "none"}


def block_boundaries(proplot, detectors) -> tuple[list[int], list[tuple[float, str]]]:
    """
    Channel boundaries and labels within a concatenated detector block.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    detectors : sequence of str
        Detectors in the order they are concatenated.

    Returns
    -------
    starts : list of int
        Bin index at which each channel starts (for divider lines).
    labels : list of tuple
        A ``(centre, label)`` pair per channel, for annotation.
    """
    starts, labels = [], []
    position = 0
    for detector in detectors:
        for mode, channel in proplot.channels:
            width = proplot.channel_bins[(mode, channel)]
            starts.append(position)
            labels.append((position + 0.5 * width, f"{detector}\n{channel}"))
            position += width
    return starts, labels


def draw_blocks(ax, proplot, detectors, label_y: float = 0.92) -> None:
    """
    Mark the channel blocks of a concatenated detector block on an axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis whose x axis is the concatenated bin index.
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    detectors : sequence of str
        Detectors in the order they are concatenated.
    label_y : float, optional
        Axes-fraction height of the channel labels.
    """
    starts, labels = block_boundaries(proplot, detectors)
    for start in starts[1:]:
        ax.axvline(start, color="0.6", lw=0.8, ls=":", zorder=1)
    for centre, text in labels:
        ax.text(
            centre,
            label_y,
            text,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=7,
            color="0.35",
        )


def check_partition(proplot, families: dict[str, str]) -> dict[str, str]:
    """
    Assert that the families tile the component list exactly.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    families : dict
        Family label -> regex over component names.

    Returns
    -------
    dict
        Component name -> family label.

    Raises
    ------
    AssertionError
        If a component matches two families or none.
    """
    from mutools import profit

    assignment: dict[str, str] = {}
    for label, pattern in families.items():
        for name in profit.match_components(proplot.component_names, pattern, label):
            assert name not in assignment, f"{name} matches {label} and {assignment[name]}"
            assignment[name] = label
    missing = set(proplot.component_names) - set(assignment)
    assert not missing, f"components in no family: {sorted(missing)}"
    return assignment


def describe(proplot) -> None:
    """
    Print what the PROplot file holds.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    """
    print("detectors:  ", proplot.detectors)
    print("channels:   ", proplot.channels)
    print("bins/channel:", proplot.channel_bins)
    print("components: ", len(proplot.components))
    print(
        f"full space: {proplot.cv.shape[0]} bins -> collapsed {proplot.cv_collapsed.shape[0]}"
        f" ({proplot.detector_width()} per detector)"
    )


def family_table(proplot, result, assignment: dict[str, str]) -> list[tuple]:
    """
    Constrain each family alone and the total with it removed; print a table.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The constraint with every component included.
    assignment : dict
        Component name -> family label, from :func:`check_partition`.

    Returns
    -------
    list of tuple
        ``(label, n_components, alone_prior, alone_constrained, total_without)``
        per family, in :data:`FAMILIES` order.
    """
    from mutools import profit

    rows = []
    for label, pattern in FAMILIES.items():
        alone = profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR, select=pattern)
        removed = profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR, exclude=pattern)
        rows.append(
            (
                label,
                sum(v == label for v in assignment.values()),
                alone.prior_error.mean(),
                alone.constrained_error.mean(),
                removed.constrained_error.mean(),
            )
        )

    print(
        f"\n{'family':15s} {'n':>3s} {'alone prior':>12s} {'alone constr':>13s}"
        f" {'total w/o it':>13s}"
    )
    print(
        f"{'TOTAL':15s} {len(assignment):3d} {result.prior_error.mean():12.4f}"
        f" {result.constrained_error.mean():13.4f} {'-':>13s}"
    )
    for label, count, prior, constrained, without in sorted(rows, key=lambda r: -r[3]):
        print(f"{label:15s} {count:3d} {prior:12.4f} {constrained:13.4f} {without:13.4f}")
    return rows


def plot_constraint(proplot, result):
    """
    Per-bin ICARUS fractional uncertainty before and after the constraint.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The constraint with every component included.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools.plotting.helpers import mark_axis

    fig, ax = plt.subplots(figsize=(9, 4.5))
    bins = np.arange(result.residual.shape[0])
    ax.step(bins, result.prior_error, where="mid", color="0.35", lw=1.6, label="Before constraint")
    ax.step(
        bins,
        result.constrained_error,
        where="mid",
        color=SERIES_COLORS[0],
        lw=1.8,
        label="After SBND constraint",
    )
    draw_blocks(ax, proplot, FAR)
    ax.set_xlabel(BIN_LABEL)
    ax.set_ylabel("Fractional uncertainty")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, loc="lower left", **LEGEND_BOX)
    mark_axis(ax, "SBN Simulation")
    return fig


def plot_family_bars(rows: list[tuple]):
    """
    Each family on its own, before and after the SBND constraint.

    Parameters
    ----------
    rows : list of tuple
        Output of :func:`family_table`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools.plotting.helpers import mark_axis

    labels = [row[0] for row in rows]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.barh(
        y - 0.2,
        [row[2] for row in rows],
        height=0.38,
        color="0.75",
        label="Alone, before constraint",
    )
    ax.barh(
        y + 0.2,
        [row[3] for row in rows],
        height=0.38,
        color=SERIES_COLORS[0],
        label="Alone, after SBND constraint",
    )
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Mean fractional uncertainty on ICARUS")
    ax.legend(fontsize=8, loc="lower right", **LEGEND_BOX)
    mark_axis(ax, "SBN Simulation")
    return fig


def plot_family_removed(proplot, result):
    """
    Constrained uncertainty with one family removed at a time.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The constraint with every component included.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools import profit

    bins = np.arange(result.residual.shape[0])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.step(bins, result.constrained_error, where="mid", color="0.35", lw=1.8, label="Total")
    for color, label in zip(SERIES_COLORS, REMOVED_SHOWN, strict=False):
        res = profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR, exclude=FAMILIES[label])
        ax.step(
            bins, res.constrained_error, where="mid", color=color, lw=1.4, label=f"without {label}"
        )
    draw_blocks(ax, proplot, FAR)
    ax.set_xlabel(BIN_LABEL)
    ax.set_ylabel("Fractional uncertainty\nafter SBND constraint")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=7, loc="lower left", ncol=2, **LEGEND_BOX)
    return fig


def plot_detector(proplot, result):
    """
    Detector systematics, each alone, after the constraint.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The constraint with every component included (sets the x range).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools import profit

    bins = np.arange(result.residual.shape[0])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    detector_only = profit.constraint(
        proplot, nd_detector=NEAR, fd_detector=FAR, select=FAMILIES["Detector"]
    )
    ax.step(
        bins,
        detector_only.constrained_error,
        where="mid",
        color="0.35",
        lw=1.8,
        label="All detector systematics",
    )
    for color, (label, pattern) in zip(SERIES_COLORS, DETECTOR_SOURCES.items(), strict=False):
        names = profit.match_components(proplot.component_names, pattern, label)
        res = profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR, select=pattern)
        ax.step(
            bins,
            res.constrained_error,
            where="mid",
            color=color,
            lw=1.4,
            label=f"{label} ({len(names)})",
        )
        print(
            f"{label:20s} {len(names)} component(s), mean constrained "
            f"{res.constrained_error.mean():.4f}"
        )
    draw_blocks(ax, proplot, FAR)
    ax.set_xlabel(BIN_LABEL)
    ax.set_ylabel("Fractional uncertainty\nafter SBND constraint")
    # Headroom for the block labels, and the legend below the axes: seven
    # series do not fit inside the frame without covering the curves.
    ax.set_ylim(0, detector_only.constrained_error.max() * 1.35)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4)
    return fig


def plot_exposure(proplot, result):
    """
    Constrained uncertainty against SBND exposure, with the degeneracy floor.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The nominal constraint (sets the x range).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools import profit

    bins = np.arange(result.residual.shape[0])
    scan = {
        rf"$\alpha$ = {a:g}": profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR, alpha=a)
        for a in ALPHAS
    }
    scan[r"$\alpha \to \infty$"] = profit.constraint(
        proplot, nd_detector=NEAR, fd_detector=FAR, statistics=False
    )
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for color, (label, res) in zip(SERIES_COLORS, scan.items(), strict=False):
        ax.step(bins, res.constrained_error, where="mid", color=color, lw=1.6, label=label)
    draw_blocks(ax, proplot, FAR)
    ax.set_xlabel(BIN_LABEL)
    ax.set_ylabel("Fractional uncertainty\nafter SBND constraint")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, loc="lower left", **LEGEND_BOX)
    return fig


def plot_modes(proplot, result):
    """
    Leading eigenmodes of the surviving ICARUS covariance.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.
    result : mutools.profit.ConstraintResult
        The nominal constraint.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools import profit

    values, vectors, participation = profit.eigenmodes(result.residual)
    print(f"participation ratio: {participation:.2f} effective modes")
    print("leading eigenvalues: ", values[:5].round(6))

    bins = np.arange(result.residual.shape[0])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for index, color in enumerate(SERIES_COLORS[:N_MODES]):
        shape = np.sqrt(values[index]) * vectors[:, index]
        ax.step(
            bins,
            shape,
            where="mid",
            color=color,
            lw=1.6,
            label=rf"mode {index}, $\sqrt{{\lambda}}$ = {np.sqrt(values[index]):.3f}",
        )
    ax.axhline(0.0, color="0.6", lw=0.8)
    draw_blocks(ax, proplot, FAR)
    ax.set_xlabel(BIN_LABEL)
    ax.set_ylabel("One-sigma fractional shift")
    ax.legend(fontsize=8, loc="lower left", **LEGEND_BOX)
    return fig


def plot_knobs(proplot):
    """
    Per-knob constraint efficiency: statistical ceiling against achieved.

    Parameters
    ----------
    proplot : mutools.profit.ProPlot
        The loaded PROplot file.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools import profit

    model = profit.build_responses(proplot, nd_detector=NEAR, fd_detector=FAR)
    post = profit.posterior(model, proplot.fractional_stats(NEAR))
    print(f"{model.n_knobs} effective knobs from {len(model.captured)} components")

    order = np.argsort(post.efficiency)[::-1][:N_KNOBS]
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(y, post.efficiency_ceiling[order], color="0.85", label="Statistical ceiling")
    ax.barh(
        y,
        post.efficiency[order],
        height=0.55,
        color=SERIES_COLORS[0],
        label="Achieved in joint fit",
    )
    ax.set_yticks(y)
    ax.set_yticklabels([model.labels[i] for i in order], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Constraint efficiency")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    return fig


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.add_argument("--input", type=Path, default=INPUT, help="PROplot file with Covariance")
    parser.add_argument(
        "--only",
        nargs="+",
        choices=FIGURES,
        default=list(FIGURES),
        help="figures to make (default: all; 'knobs' is the slow one)",
    )
    args = parser.parse_args(argv)

    common.configure_style(fmt=args.format)
    common.set_provenance(__file__, [args.input])
    from mutools import profit

    proplot = profit.ProPlot(args.input)
    describe(proplot)
    result = profit.constraint(proplot, nd_detector=NEAR, fd_detector=FAR)
    print(f"far block: {result.residual.shape[0]} bins")
    print(f"prior       mean frac error: {result.prior_error.mean():.4f}")
    print(f"constrained mean frac error: {result.constrained_error.mean():.4f}")
    print(f"mean reduction:              {result.reduction.mean():.3f}")

    assignment = check_partition(proplot, FAMILIES)
    print(f"{len(assignment)} components in {len(FAMILIES)} families")

    rows = None
    if "family_bars" in args.only:
        rows = family_table(proplot, result, assignment)

    makers = {
        "constraint": lambda: plot_constraint(proplot, result),
        "family_bars": lambda: plot_family_bars(rows),
        "family_removed": lambda: plot_family_removed(proplot, result),
        "detector": lambda: plot_detector(proplot, result),
        "exposure": lambda: plot_exposure(proplot, result),
        "modes": lambda: plot_modes(proplot, result),
        "knobs": lambda: plot_knobs(proplot),
    }
    for name in FIGURES:
        if name in args.only:
            dest = args.output / f"conditional_{name}.{args.format}"
            path = common.save_figure(makers[name](), dest)
            print("wrote", path)


if __name__ == "__main__":
    main()
