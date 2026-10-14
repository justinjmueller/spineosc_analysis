"""
Injected oscillation parameter scan: spectra, rate deficit and SBND bins
outside the CV systematic envelope.

Replaces ``plotting/sbnd_data_access.ipynb``. The scan is a grid of
(dmsq, sinsq) points, one PROfit ``*_PROplot.root`` file per point, in
``runs/legacy_root/parameter_scan`` (900 files). Each file holds the
oscillated prediction (``Osc_hists/``, split by subchannel) and the
unoscillated CV (``ErrorBand/Var0/..._cv1d``, identical in every file).

Figures made (``--figure``):

``spectra``
    SBND and ICARUS spectra at maximal mixing for four dmsq values,
    against the unoscillated CV (two panels per detector pair).
``ratio``
    The same points as a ratio to the unoscillated CV, where the shape
    of the disappearance is easiest to read.
``deficit``
    Total rate lost relative to the unoscillated CV over the whole
    grid. The band of strongest disappearance sits at roughly 5-6x
    higher dmsq for SBND than ICARUS, the ratio of the baselines.
``outside``
    The comparison itself: for each injected point, how many of the 30
    SBND bins (16 ccqelike + 14 nproton) land outside CV - sigma_low to
    CV + sigma_high. The band is asymmetric and the two sides are read
    separately, not symmetrised. Written twice, auto-ranged and zoomed
    out (sinsq 1e-3..1, dmsq 1e-2..1e2). These two feed the SBND
    unblinding request (``sbnd_bins_outside_cvband{,_zoomout}.pdf``).

Inputs and a known path mismatch
--------------------------------
The notebook reads the CV envelope from ``sbnd_scanmatch_v1_PROplot.root``
(``ErrorBand/Var0/nu_SBN_sbnd<channel>_preerrband`` graphs). That file
does not exist in ``runs/legacy_root`` (recorded in
``docs/RESIDUAL_TASKS.md``). The nearest existing file,
``sbnd_scanmatch_v1_PROplot_histograms.root``, carries the same
information in its ``errorband`` tree (``y_value``, ``error_y_low``,
``error_y_high`` for the CV, 16 + 14 = 30 points), so the envelope is
read from that tree via ``ENVELOPE_FILE``. The values are those the
graphs would hold, but this was not compared against the missing file.
The envelope is at the exposure of ``xml/contours/sbnd_scanmatch.xml``
(POT 4.196270e+18); the script asserts that envelope and scan CV agree
to 5% in total normalisation.

The notebook wrote the ``outside`` figures to the sibling
``spineosc_sbnd_unblinding_request/figures`` directory; here they go to
``--output`` (default ``figures/``), so copy them over when final.

Run::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/sbnd_scan.py --figure outside
"""

from __future__ import annotations

import re
import sys
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Directory with one ``dmsq*_sinsq*_grid*_PROplot.root`` per grid point.
INPUT = common.LEGACY_ROOT / "parameter_scan"

#: File name pattern of a grid point: dmsq, sinsq, grid index.
SCAN_PATTERN = re.compile(r"dmsq([-\d.eE+]+)_sinsq([-\d.eE+]+)_grid(\d+)_PROplot\.root$")

#: CV and systematic envelope of SBND in the scan's bins. The notebook's
#: ``sbnd_scanmatch_v1_PROplot.root`` is missing; this is the nearest file.
ENVELOPE_FILE = common.LEGACY_ROOT / "sbnd_scanmatch_v1_PROplot_histograms.root"

#: Detectors and channels as named in the scan files. Each subchannel also
#: has a matching ``*slc`` histogram, but those are empty, so the six
#: subchannels are the whole prediction. ``nu`` is the mode and ``_numu``
#: is part of the channel name.
DETECTORS = ("SBND", "ICARUS")
CHANNELS = ("ccqelike", "nproton")
SUBCHANNELS = ("cos", "intime", "nc", "nuecc", "nuothercc", "nutarget")

#: Channel shown in the spectra, ratio and deficit figures.
CHANNEL = "ccqelike"

#: Fixed series order, checked for colourblind separation.
SERIES_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")

#: dmsq values (eV^2) for the spectra; the nearest grid points are used.
DMSQ_TARGETS = (0.3, 1.0, 3.0, 30.0)

#: Order in which the channels are concatenated, on BOTH sides (injected
#: spectra and envelope), so that column j is the same bin in each.
BAND_ORDER = CHANNELS

#: Largest tolerated envelope/scan total-rate mismatch. About 1% is the
#: expected MC-version drift; more means the two are not comparable.
MAX_NORM_MISMATCH = 0.05

#: Figures made when ``--figure`` is not given.
FIGURES_ALL = ("spectra", "ratio", "deficit", "outside")


def find_scan(directory: Path) -> pd.DataFrame:
    """
    List the grid points in a scan directory.

    Parameters
    ----------
    directory : Path
        Directory holding the ``*_PROplot.root`` scan files.

    Returns
    -------
    pandas.DataFrame
        Columns ``dmsq``, ``sinsq``, ``grid`` and ``path``, sorted by
        dmsq then sinsq.

    Raises
    ------
    FileNotFoundError
        If no file matches the scan naming pattern.
    """
    rows = []
    for path in sorted(directory.glob("*_PROplot.root")):
        match = SCAN_PATTERN.search(path.name)
        if match is None:
            continue
        rows.append(
            {
                "dmsq": float(match.group(1)),
                "sinsq": float(match.group(2)),
                "grid": int(match.group(3)),
                "path": path,
            }
        )
    if not rows:
        raise FileNotFoundError(f"no scan files matching {SCAN_PATTERN.pattern} in {directory}")
    return pd.DataFrame(rows).sort_values(["dmsq", "sinsq"], ignore_index=True)


def oscillated(f, detector: str, channel: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Oscillated spectrum at one grid point: the subchannels, summed.

    Parameters
    ----------
    f : uproot.ReadOnlyDirectory
        An open scan file.
    detector, channel : str
        Names as they appear in the file.

    Returns
    -------
    edges, total : np.ndarray
        Bin edges and summed bin contents.
    """
    total = None
    for sub in SUBCHANNELS:
        hist = f[f"Osc_hists/nu_{detector}_{channel}_numu_{sub}"]
        total = hist.values() if total is None else total + hist.values()
    return hist.axis().edges(), total


def unoscillated(f, detector: str, channel: str) -> tuple[np.ndarray, np.ndarray]:
    """
    The no-oscillation CV spectrum, the same in every grid file.

    Parameters
    ----------
    f : uproot.ReadOnlyDirectory
        An open scan file.
    detector, channel : str
        Names as they appear in the file.

    Returns
    -------
    edges, values : np.ndarray
        Bin edges and CV bin contents.
    """
    hist = f[f"ErrorBand/Var0/nu_{detector}_{channel}_numu_cv1d"]
    return hist.axis().edges(), hist.values()


def load_grid(scan: pd.DataFrame) -> tuple[dict, dict, dict]:
    """
    Load the CV and every grid point's oscillated spectrum (~8 s for 900 files).

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.

    Returns
    -------
    cv_ref : dict
        ``(detector, channel) -> (edges, values)`` of the unoscillated CV.
    osc_spectra : dict
        ``(detector, channel) -> array (n_points, n_bins)``.
    osc_rates : dict
        ``(detector, channel) -> array (n_points,)`` of oscillated over
        unoscillated total rate.
    """
    import uproot

    cv_ref, osc_spectra = {}, {}
    with uproot.open(scan.path.iloc[0]) as f:
        for detector in DETECTORS:
            for channel in CHANNELS:
                cv_ref[detector, channel] = unoscillated(f, detector, channel)
                osc_spectra[detector, channel] = np.zeros(
                    (len(scan), len(cv_ref[detector, channel][1]))
                )
    for i, path in enumerate(scan.path):
        with uproot.open(path) as f:
            for detector in DETECTORS:
                for channel in CHANNELS:
                    osc_spectra[detector, channel][i] = oscillated(f, detector, channel)[1]
    osc_rates = {
        key: spectra.sum(axis=1) / cv_ref[key][1].sum() for key, spectra in osc_spectra.items()
    }
    return cv_ref, osc_spectra, osc_rates


def pick_dmsq(scan: pd.DataFrame, targets: Sequence[float] = DMSQ_TARGETS) -> list[float]:
    """
    Grid dmsq values nearest to each target.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.
    targets : sequence of float, optional
        Wanted dmsq values in eV^2.

    Returns
    -------
    list of float
        The nearest grid value for each target.
    """
    values = scan.dmsq.unique()
    return [values[np.argmin(np.abs(values - t))] for t in targets]


def spectrum_at(scan, osc_spectra, detector, channel, dmsq, sinsq) -> np.ndarray:
    """
    Oscillated spectrum at one grid point.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.
    osc_spectra : dict
        Output of `load_grid`.
    detector, channel : str
        Which spectrum.
    dmsq, sinsq : float
        Grid coordinates; must be exact grid values.

    Returns
    -------
    np.ndarray
        Bin contents at that point.
    """
    row = (scan.dmsq == dmsq) & (scan.sinsq == sinsq)
    return osc_spectra[detector, channel][row.to_numpy()][0]


def read_envelope(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    CV, lower and upper systematic errors, concatenated over `BAND_ORDER`.

    Reads the ``errorband`` tree of a ``*_PROplot_histograms.root`` file
    (detector 0 = SBND, total subchannel, prefix ``CV``), one row per
    bin, with channel indices following `CHANNELS`.

    Parameters
    ----------
    path : Path
        The envelope file (`ENVELOPE_FILE`).

    Returns
    -------
    cv, low, high : np.ndarray
        CV and the lower and upper errors, 30 bins in total.

    Raises
    ------
    FileNotFoundError
        If `path` does not exist.
    KeyError
        If a channel has no rows in the tree.
    """
    import uproot

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. The caches exist but the plot step has not run:\n"
            "    scripts/profit_run.sh jobs/contours/sbnd_scanmatch.ini plot"
        )
    with uproot.open(path) as f:
        tree = pd.DataFrame(f["errorband"].arrays(library="np"))
    cv, lo, hi = [], [], []
    for channel_index in range(len(BAND_ORDER)):
        rows = tree[
            (tree["variable"] == 0)
            & (tree["detector"] == 0)
            & (tree["channel"] == channel_index)
            & (tree["prefix"] == "CV")
        ].sort_values("point_index")
        if rows.empty:
            raise KeyError(f"no CV error band for channel {channel_index} in {path.name}")
        cv.append(rows["y_value"].to_numpy())
        lo.append(rows["error_y_low"].to_numpy())
        hi.append(rows["error_y_high"].to_numpy())
    return np.concatenate(cv), np.concatenate(lo), np.concatenate(hi)


def bins_outside(
    scan: pd.DataFrame, cv_ref: dict, osc_spectra: dict, envelope_file: Path = ENVELOPE_FILE
) -> np.ndarray:
    """
    Number of SBND bins outside the CV systematic envelope, per grid point.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.
    cv_ref, osc_spectra : dict
        Output of `load_grid`.
    envelope_file : Path, optional
        File read by `read_envelope`.

    Returns
    -------
    np.ndarray
        Count of bins (of 30) with injected content below CV - sigma_low
        or above CV + sigma_high, one entry per grid point.

    Raises
    ------
    AssertionError
        If envelope and scan CV disagree by more than `MAX_NORM_MISMATCH`.
    """
    cv_band, err_low, err_high = read_envelope(envelope_file)
    scan_cv = np.concatenate([cv_ref["SBND", c][1] for c in BAND_ORDER])
    start = 0
    for channel in BAND_ORDER:
        stop = start + len(cv_ref["SBND", channel][1])
        ratio = cv_band[start:stop].sum() / cv_ref["SBND", channel][1].sum()
        print(f"{channel:9s} envelope/scan CV = {ratio:.4f}")
        start = stop
    total_ratio = cv_band.sum() / scan_cv.sum()
    print(f"combined  envelope/scan CV = {total_ratio:.4f}")
    assert abs(total_ratio - 1) < MAX_NORM_MISMATCH, (
        f"envelope and scan disagree on normalisation by {100 * (total_ratio - 1):.1f}%. "
        "Check the SBND POT in xml/contours/sbnd_scanmatch.xml against the scan."
    )
    injected = np.concatenate([osc_spectra["SBND", c] for c in BAND_ORDER], axis=1)
    outside = ((injected < cv_band - err_low) | (injected > cv_band + err_high)).sum(axis=1)
    print(f"{injected.shape[1]} bins per grid point")
    print(
        f"bins outside: min {outside.min()}, max {outside.max()}, mean {outside.mean():.1f}"
    )
    return outside


def plot_spectra(scan, cv_ref, osc_spectra, picks, ratio: bool = False) -> plt.Figure:
    """
    SBND and ICARUS spectra at maximal mixing for a few dmsq values.

    Maximal mixing means the largest sinsq in the grid. SBND (110 m) and
    ICARUS (600 m) see the same oscillation at different L/E, so the
    same dmsq does very different things to the two spectra.

    Parameters
    ----------
    scan, cv_ref, osc_spectra : see `load_grid`
    picks : sequence of float
        Grid dmsq values, one per series.
    ratio : bool, optional
        Draw the ratio to the unoscillated CV instead of the spectra.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from mutools.plotting.helpers import mark_axis

    sinsq = scan.sinsq.max()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, detector in zip(axes, DETECTORS, strict=True):
        edges, cv = cv_ref[detector, CHANNEL]
        centers = 0.5 * (edges[:-1] + edges[1:])
        if ratio:
            ax.axhline(1.0, color="0.35", lw=1.2, zorder=3)
        else:
            ax.stairs(cv, edges, color="0.35", lw=1.6, label="No oscillation", zorder=3)
        for color, dmsq in zip(SERIES_COLORS, picks, strict=False):
            spec = spectrum_at(scan, osc_spectra, detector, CHANNEL, dmsq, sinsq)
            label = rf"$\Delta m^2 = {dmsq:.3g}$"
            if ratio:
                ax.plot(
                    centers, spec / cv, color=color, lw=1.8, marker="o", ms=3.5, label=label
                )
            else:
                ax.stairs(spec, edges, color=color, lw=1.6, label=label)
        ax.set_xlabel("Reconstructed neutrino energy [GeV]")
        ax.set_ylabel("Oscillated / unoscillated" if ratio else "Events")
        ax.set_title(f"{detector} {CHANNEL}")
        if ratio:
            ax.set_ylim(0, 1.15)
        else:
            ax.set_ylim(bottom=0)
    axes[0].legend(
        title=rf"$\sin^2 2\theta = {sinsq:g}$",
        frameon=False,
        fontsize=8,
        loc="lower left" if ratio else "best",
    )
    mark_axis(axes[0], "SBN Simulation")
    return fig


def plot_deficit(scan, osc_rates) -> plt.Figure:
    """
    Total rate deficit (percent of the unoscillated CV) over the grid.

    Prints the deepest deficit and where it occurs for each detector.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.
    osc_rates : dict
        Output of `load_grid`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, detector in zip(axes, DETECTORS, strict=True):
        table = pd.DataFrame(
            {
                "dmsq": scan.dmsq,
                "sinsq": scan.sinsq,
                "deficit": 100 * (1 - osc_rates[detector, CHANNEL]),
            }
        ).pivot_table(index="dmsq", columns="sinsq", values="deficit")
        mesh = ax.pcolormesh(
            table.columns.to_numpy(),
            table.index.to_numpy(),
            table.to_numpy(),
            shading="nearest",
            cmap="Blues",
            vmin=0,
            vmax=80,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\sin^2 2\theta$")
        ax.set_ylabel(r"$\Delta m^2$ [eV$^2$]")
        ax.set_title(f"{detector} {CHANNEL}")
        fig.colorbar(mesh, ax=ax, label="Rate deficit [%]")
    for detector in DETECTORS:
        deepest = int(np.argmin(osc_rates[detector, CHANNEL]))
        print(
            f"{detector:7s} deepest deficit "
            f"{100 * (1 - osc_rates[detector, CHANNEL][deepest]):.1f}% at "
            f"dmsq={scan.dmsq[deepest]:.3g}, sinsq={scan.sinsq[deepest]:.3g}"
        )
    return fig


def plot_outside(scan, outside, n_bins: int, zoom_out: bool = False) -> plt.Figure:
    """
    Map of the number of SBND bins outside the CV band over the grid.

    Parameters
    ----------
    scan : pandas.DataFrame
        Output of `find_scan`.
    outside : np.ndarray
        Output of `bins_outside`.
    n_bins : int
        Total bins per grid point (the colour scale maximum).
    zoom_out : bool, optional
        Fix the axes to sinsq 1e-3..1 and dmsq 1e-2..1e2.

    Returns
    -------
    matplotlib.figure.Figure
    """
    table = pd.DataFrame({"dmsq": scan.dmsq, "sinsq": scan.sinsq, "outside": outside}).pivot_table(
        index="dmsq", columns="sinsq", values="outside"
    )
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    mesh = ax.pcolormesh(
        table.columns.to_numpy(),
        table.index.to_numpy(),
        table.to_numpy(),
        shading="nearest",
        cmap="Blues",
        vmin=0,
        vmax=n_bins,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    if zoom_out:
        ax.set_xlim(1e-3, 1)
        ax.set_ylim(1e-2, 1e2)
    ax.set_xlabel(r"$\sin^2 2\theta_{\mu\mu}$")
    ax.set_ylabel(r"$\Delta m_{14}^2$ [eV$^2$]")
    fig.colorbar(mesh, ax=ax, label="SBND bins outside CV band")
    return fig


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
        "--figure",
        choices=(*FIGURES_ALL, "all"),
        default="all",
        help="which figure to make (default: all)",
    )
    parser.add_argument("--input", type=Path, default=INPUT, help="scan directory")
    parser.add_argument("--envelope", type=Path, default=ENVELOPE_FILE, help="CV envelope file")
    args = parser.parse_args(argv)

    common.configure_style(args.format)
    scan = find_scan(args.input)
    print(
        f"{len(scan)} grid points: {scan.dmsq.nunique()} dmsq x {scan.sinsq.nunique()} sinsq"
    )
    print(f"  dmsq  {scan.dmsq.min():g} to {scan.dmsq.max():g} eV^2")
    print(f"  sinsq {scan.sinsq.min():g} to {scan.sinsq.max():g}")
    cv_ref, osc_spectra, osc_rates = load_grid(scan)

    wanted = FIGURES_ALL if args.figure == "all" else (args.figure,)
    ext = args.format
    out = args.output
    if "spectra" in wanted:
        fig = plot_spectra(scan, cv_ref, osc_spectra, pick_dmsq(scan))
        print(common.save_figure(fig, out / f"scan_spectra.{ext}"))
    if "ratio" in wanted:
        fig = plot_spectra(scan, cv_ref, osc_spectra, pick_dmsq(scan), ratio=True)
        print(common.save_figure(fig, out / f"scan_spectra_ratio.{ext}"))
    if "deficit" in wanted:
        print(common.save_figure(plot_deficit(scan, osc_rates), out / f"scan_rate_deficit.{ext}"))
    if "outside" in wanted:
        outside = bins_outside(scan, cv_ref, osc_spectra, args.envelope)
        n_bins = sum(len(cv_ref["SBND", c][1]) for c in BAND_ORDER)
        fig = plot_outside(scan, outside, n_bins)
        print(common.save_figure(fig, out / f"sbnd_bins_outside_cvband.{ext}"))
        fig = plot_outside(scan, outside, n_bins, zoom_out=True)
        print(common.save_figure(fig, out / f"sbnd_bins_outside_cvband_zoomout.{ext}"))


if __name__ == "__main__":
    main()
