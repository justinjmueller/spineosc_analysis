#!/usr/bin/env python3
"""
Compare the AMR sensitivity contours (Asimov data, delta-chi2 = 1.6) of different reconstructed-energy binnings.

Reads ``runs/amr_binning/<config>_<variant>/ab_<config>_<variant>_surf_surface_amr.txt`` (written by
``scripts/amr_binning/run_variant.sh``) and, as the baseline of each configuration, the existing
``runs/amr_compare/amr{2,3}_surf_surface_amr.txt`` (1D exclusive = amr2, 1D inclusive = amr3). For every configuration it writes

* an overlay of the delta-chi2 = 1.6 contours of all binnings,
* the limit ``sin^2 2theta`` at fixed ``dmsq`` relative to the baseline binning,
* a table (CSV) with the number of bins, the excluded area, the limit at fixed dmsq, the AMR fit count and the timings.

The AMR text files hold ``log10(sin^2 2theta)  log10(dmsq)  delta_chi2  chi2_abs`` for the adaptively chosen points; the surface
between them is linearly interpolated (Delaunay), the same as ``runs/amr_compare/overlay.py``.

Usage::

    python scripts/amr_binning/compare_binnings.py --config incl1D
    python scripts/amr_binning/compare_binnings.py --config excl1D --output docs/amr_binning
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RUNS = REPO / "runs" / "amr_binning"
BASELINE = {"incl1D": REPO / "runs" / "amr_compare" / "amr3_surf_surface_amr.txt",
            "excl1D": REPO / "runs" / "amr_compare" / "amr2_surf_surface_amr.txt"}
BASELINE_LOGS = {"incl1D": ("process_amr3.time", "surface_amr3.time"), "excl1D": ("process_amr2.time", "surface_amr2.time")}

#: Contour level and the grid on which the interpolated surfaces are evaluated (log10 sin^2 2theta, log10 dmsq).
LEVEL = 1.6
GX, GY = np.meshgrid(np.linspace(-3, 0, 400), np.linspace(-2, 2, 400))

#: dmsq values [eV^2] at which the limit in sin^2 2theta is tabulated.
DMSQ_POINTS = (0.3, 1.0, 1.5, 3.0, 10.0, 30.0)


def load_surface(path: Path) -> np.ndarray:
    """Interpolated delta-chi2 on the common grid, from an AMR text file."""
    from scipy.interpolate import LinearNDInterpolator

    data = np.loadtxt(path)
    return LinearNDInterpolator(np.c_[data[:, 0], data[:, 1]], data[:, 2])(GX, GY), len(data)


def limit_curve(surface: np.ndarray) -> np.ndarray:
    """sin^2 2theta at the contour level, for each dmsq row of the grid (NaN where the level is not reached)."""
    out = np.full(GY.shape[0], np.nan)
    for k in range(GY.shape[0]):
        above = np.nan_to_num(surface[k], nan=0.0) > LEVEL
        if above.any():
            out[k] = 10 ** GX[k, np.argmax(above)]
    return out


def seconds(path: Path, kind: str) -> float:
    """Wall time from a ``time`` log (``-p``: 'real N'; ``-l``: 'N real')."""
    if not path.exists():
        return float("nan")
    text = path.read_text()
    match = re.search(r"real\s+([\d.]+)", text) if kind == "process" else re.search(r"([\d.]+)\s+real", text)
    return float(match.group(1)) if match else float("nan")


def nbins(config: str, variant: str) -> str:
    """Number of bins per channel of a variant, from its XML."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import make_binnings as mb

    xml = RUNS / "xml" / f"{config}_{variant}.xml"
    text = xml.read_text()
    return "+".join(str(len(mb.fit_edges(text[s:e])[1]) - 1) for _, s, e in mb.channel_blocks(text))


def main() -> None:
    """Command-line entry point."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", choices=sorted(BASELINE), required=True)
    parser.add_argument("--output", type=Path, default=REPO / "docs" / "amr_binning")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    entries = {"base": (BASELINE[args.config], REPO / "runs" / "amr_compare" / "logs", *BASELINE_LOGS[args.config])}
    for directory in sorted(RUNS.glob(f"{args.config}_*")):
        variant = directory.name.split("_", 1)[1]
        text = directory / f"ab_{args.config}_{variant}_surf_surface_amr.txt"
        if text.exists():
            entries[variant] = (text, directory / "logs", "process.time", "surface.time")
    surfaces, rows = {}, []
    for variant, (path, logs, process_time, surface_time) in entries.items():
        surface, fits = load_surface(path)
        surfaces[variant] = surface
        rows.append({"variant": variant, "bins": nbins(args.config, variant) if variant != "base" else
                     ("17" if args.config == "incl1D" else "16+14"), "fits": fits,
                     "excluded_area": float(np.nan_to_num(surface, nan=0.0).__gt__(LEVEL).mean() * 12),
                     "process_s": seconds(logs / process_time, "process"), "surface_s": seconds(logs / surface_time, "surface"),
                     **{f"limit_dm{dm:g}": float(np.interp(np.log10(dm), GY[:, 0], limit_curve(surface))) for dm in DMSQ_POINTS}})
    base = {k: v for k, v in rows[0].items() if k.startswith("limit_")}
    for row in rows:
        for key, value in base.items():
            row[key.replace("limit_", "ratio_")] = row[key] / value
    header = list(rows[0])
    with open(args.output / f"binning_{args.config}.csv", "w") as handle:
        handle.write(",".join(header) + "\n")
        for row in rows:
            handle.write(",".join(f"{row[h]:.4g}" if isinstance(row[h], float) else str(row[h]) for h in header) + "\n")

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5))
    for i, (variant, surface) in enumerate(surfaces.items()):
        left.contour(10 ** GX, 10 ** GY, surface, [LEVEL], colors=f"C{i}", linestyles="-" if variant == "base" else "--")
        left.plot([], [], color=f"C{i}", ls="-" if variant == "base" else "--", label=f"{variant} ({rows[i]['bins']} bins)")
        right.plot(10 ** GY[:, 0], limit_curve(surface) / limit_curve(surfaces["base"]), color=f"C{i}", label=variant)
    left.set(xscale="log", yscale="log", xlabel=r"$\sin^22\theta_{\mu\mu}$", ylabel=r"$\Delta m^2_{41}$ [eV$^2$]",
             title=f"{args.config}: $\\Delta\\chi^2={LEVEL}$ (Asimov)")
    left.legend(fontsize=8)
    right.axhline(1.0, color="0.5", lw=0.8)
    right.set(xscale="log", xlabel=r"$\Delta m^2_{41}$ [eV$^2$]", ylabel="limit / baseline limit", ylim=(0.6, 1.8))
    right.grid(alpha=0.3)
    right.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / f"binning_{args.config}.png", dpi=140)
    for row in rows:
        print(f"{row['variant']:14s} bins {row['bins']:>6s} fits {row['fits']:5d} area {row['excluded_area']:.3f} "
              f"process {row['process_s']:6.0f}s surface {row['surface_s']:6.0f}s  "
              + " ".join(f"dm{dm:g}:{row[f'ratio_dm{dm:g}']:.2f}" for dm in DMSQ_POINTS))


if __name__ == "__main__":
    main()
