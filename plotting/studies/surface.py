"""
Sensitivity contours (90% CL) from PROfit chi-square surfaces.

Overlays delta-chi2 contours of the sin^2(2 theta_mumu) - dm^2_41 plane
for sets of PROfit ``*_surf.root`` files. Replaces
``plotting/surface.ipynb``.

Inputs are the ``runs/legacy_root/contours_*_surf.root`` files written by
the contour jobs. Each holds a ``surf`` TH2D with the chi2 grid and its
axis binning, and a ``tree`` with (xbin, ybin, chi2). The tree stores
integer bin INDICES, not coordinates, so it cannot place the contour on
its own; the two agree exactly (max |surf - tree grid| == 0 on the
current files), so the TH2D is read directly. The axes are already in
linear units (sin^2(2theta) 0.01..1, dm^2 0.1..100, log-spaced).

Studies (select with ``--study``; default is all):

``baseline``
    Baseline sensitivity. Feeds the technote sensitivity-studies chapter
    (``sensitivity_baseline.pdf``).
``blinding``
    Fully blinded vs SBND absolute normalization vs full sensitivity.
    Feeds the SBND unblinding request (``surface.pdf``).
``sbnd_detvar``, ``icarus_detvar_trigger``, ``icarus_detvar_model``
    Leave-one-out detector-systematic comparisons against the full
    detector-variation set, for SBND (WireMod Y-Z, WireMod X-theta_xw,
    DENT, charge scale), ICARUS (SCE, trigger) and ICARUS (charge scale,
    WireMod). They show which detector systematic drives the loss of
    sensitivity. The notebook displayed these without saving.

Example::

    python plotting/studies/surface.py --study baseline --output figures/sensitivity_studies
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import provenance
from matplotlib.lines import Line2D

import common
from common import LEGACY_ROOT

#: Directory holding the PROfit ``*_surf.root`` files.
SURFACE_DIR = LEGACY_ROOT

#: Default contour level: delta chi2 = 1.6 is 90% CL in 1D.
LEVEL = 1.6

#: Number of grid points per axis for the interpolated contour.
GRID_POINTS = 125

#: Matplotlib colour-cycle indices used for successive curves.
COLOR_INDICES = [0, 2, 3, 4, 5]

#: Axis ranges of the sin^2(2 theta) - dm^2 plane.
XLIM = (1e-3, 1.0)
YLIM = (1e-2, 1e2)

#: Studies: name -> (output file name, [(surface file, legend label), ...]).
STUDIES: dict[str, tuple[str, list[tuple[str, str]]]] = {
    "baseline": (
        "sensitivity_baseline.pdf",
        [("contours_full_v1_surf.root", r"Baseline Sensitivity (90% CL, $\Delta \chi^2 = 1.6$)")],
    ),
    "blinding": (
        "surface.pdf",
        [
            ("contours_blind_v1_surf.root", r"Fully Blinded"),
            ("contours_sbndunblind_v1_surf.root", r"SBND Absolute Normalization"),
            ("contours_full_v1_surf.root", r"Full Sensitivity"),
        ],
    ),
    "sbnd_detvar": (
        "surface_sbnd_detvar.pdf",
        [
            ("contours_full_v1_surf.root", r"Full Detvar"),
            ("contours_full_nowmyz_surf.root", r"No SBND WireMod $Y-Z$"),
            ("contours_full_nowmxtxw_surf.root", r"No SBND WireMod $X-\theta_{xw}$"),
            ("contours_full_nodent_surf.root", r"No SBND DENT"),
            ("contours_full_nochargescale_surf.root", r"No SBND Charge Scale / Resolution"),
        ],
    ),
    "icarus_detvar_trigger": (
        "surface_icarus_detvar_trigger.pdf",
        [
            ("contours_full_v1_surf.root", r"Full Detvar"),
            ("contours_full_noicsce_surf.root", r"No ICARUS SCE"),
            ("contours_full_noictrigemu_surf.root", r"No ICARUS Trigger"),
        ],
    ),
    "icarus_detvar_model": (
        "surface_icarus_detvar_model.pdf",
        [
            ("contours_full_v1_surf.root", r"Full Detvar"),
            ("contours_full_noicchargescale_surf.root", r"No ICARUS Charge Scale / Resolution"),
            ("contours_full_noicwmxtxw_surf.root", r"No ICARUS WireMod $X-\theta_{xw}$"),
            ("contours_full_noicwmyz_surf.root", r"No ICARUS WireMod $Y-Z$"),
        ],
    ),
}


def add_surface(ax: plt.Axes, surface: Path, ci: int, level: float = LEVEL) -> Line2D:
    """
    Draw one chi2-surface contour on an axis.

    The log-spaced grid is interpolated (pchip) in log space, where it is
    evenly spaced, and bin centres are geometric: arithmetic ones would
    bias every contour toward the upper edge of each bin, increasingly so
    at low x.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to draw on.
    surface : Path
        PROfit ``*_surf.root`` file.
    ci : int
        Colour-cycle index of the contour line.
    level : float, optional
        Delta-chi2 level of the contour.

    Returns
    -------
    matplotlib.lines.Line2D
        A proxy line for the legend.
    """
    import uproot
    from scipy.interpolate import RegularGridInterpolator

    with uproot.open(surface) as rf:
        th2 = rf["surf"]
        chi2 = th2.values()  # (nx, ny), no under/overflow
        xedges = th2.axes[0].edges()
        yedges = th2.axes[1].edges()

    xc = np.sqrt(xedges[:-1] * xedges[1:])
    yc = np.sqrt(yedges[:-1] * yedges[1:])
    interp = RegularGridInterpolator(
        (np.log10(xc), np.log10(yc)), chi2 - chi2.min(), method="pchip"
    )
    xg = np.logspace(np.log10(xc.min()), np.log10(xc.max()), GRID_POINTS)
    yg = np.logspace(np.log10(yc.min()), np.log10(yc.max()), GRID_POINTS)
    xm, ym = np.meshgrid(xg, yg, indexing="ij")
    zg = interp((np.log10(xm), np.log10(ym)))

    ax.contour(xg, yg, zg.T, levels=[level], linewidths=2, colors=f"C{ci}")
    return Line2D([0], [0], color=f"C{ci}")


def plot_sensitivity(
    paths: list[Path],
    labels: list[str],
    level: float = LEVEL,
    savefig: Path | None = None,
) -> plt.Figure:
    """
    Plot sensitivity contours for several surfaces on one axis.

    Parameters
    ----------
    paths : list[Path]
        PROfit ``*_surf.root`` files.
    labels : list[str]
        Legend label of each curve.
    level : float, optional
        Delta-chi2 level of the contours.
    savefig : Path, optional
        Where to save the figure; not saved when None.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    from mutools.plotting.helpers import mark_axis

    figure = plt.figure(figsize=(8, 7))
    ax = figure.add_subplot()
    handles = [
        add_surface(ax, path, COLOR_INDICES[i], level=level) for i, path in enumerate(paths)
    ]

    ax.set_xlabel(r"$\sin^2 2 \theta_{\mu \mu}$")
    ax.set_ylabel(r"$\Delta m_{41}^2$ [eV$^2$]")
    ax.set(xscale="log", yscale="log")
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.legend(handles, labels, loc="lower left")
    mark_axis(ax, r"SBN $\nu_\mu$ Disappearance", hadj=-0.025, vadj=0.12)
    mark_axis(ax, r"PROfit v3.0.0", hadj=-0.015, vadj=0.12, alignment="right")
    if savefig is not None:
        savefig = Path(savefig)
        savefig.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(savefig)
        provenance.record(savefig, inputs=paths, script=__file__,
                          arguments={"labels": labels, "level": level},
                          versions={"selection": common.SELECTION_VERSION})
    return figure


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str], optional
        Argument list; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.add_argument("--study", action="append", choices=sorted(STUDIES),
                        help="study to draw (repeatable); default all")
    parser.add_argument("--surface-dir", type=Path, default=SURFACE_DIR,
                        help="directory holding the *_surf.root files")
    parser.add_argument("--level", type=float, default=LEVEL, help="delta-chi2 contour level")
    args = parser.parse_args(argv)

    common.configure_style(fmt=args.format)
    for name in args.study or list(STUDIES):
        filename, meta = STUDIES[name]
        fig = plot_sensitivity(
            [args.surface_dir / f for f, _ in meta],
            [label for _, label in meta],
            level=args.level,
            savefig=args.output / filename,
        )
        plt.close(fig)
        print(f"wrote {args.output / filename}")


if __name__ == "__main__":
    main(sys.argv[1:])
