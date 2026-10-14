"""
Shared paths, figure output and styling for the spineosc plotting studies.

Every script in this directory (``plotting/studies/*.py``) is a stand-in
for one analysis notebook: a documented module with a ``main()`` and a
command-line interface. They share the helpers defined here so that
paths, figure formats and stamps are set in one place.

Run a script with the mutools environment, from anywhere::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/<name>.py --help
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

#: Repository root (this file lives in ``plotting/studies``).
REPO = Path(__file__).resolve().parents[2]

#: Per-study run directories written by PROfit (git-ignored).
RUNS = REPO / "runs"

#: Root-level PROfit outputs moved out of the repository root.
LEGACY_ROOT = RUNS / "legacy_root"

#: Default destination for generated figures (git-ignored).
FIGURES = REPO / "figures"

#: Staging area for reviewed figures, one sub-directory per study.
DOCS = REPO / "docs"

#: Sibling technote working tree that the staged figures are copied into.
TECHNOTE = REPO.parent / "spineosc_technote"

#: Selection and code versions stamped onto figures.
SELECTION_VERSION = "v1.1.0"
CODE_VERSION = "v3.0.4"

#: Watermark used on analysis-in-progress figures.
WATERMARK = r"$\bf{SBN}$ Analysis-in-Progress"


def configure_style(fmt: str = "pdf", dpi: int = 300, style: str = "rootlike") -> None:
    """
    Apply the mutools figure style and saver settings.

    mutools is imported here rather than at module level so that a
    script's ``--help`` works outside the mutools environment.

    Parameters
    ----------
    fmt : str, optional
        File format used by the mutools figure saver.
    dpi : int, optional
        Resolution used when rasterizing.
    style : str, optional
        Name of a bundled mutools matplotlib style sheet.
    """
    import mutools.plotting as mp

    mp.saver.configure(fmt=fmt, dpi=dpi, rasterized=True, ytick_precision=2)
    mp.use_style(style)


def save_figure(fig: plt.Figure, path: Path, dpi: int = 200, close: bool = True) -> Path:
    """
    Write a figure, creating parent directories, and optionally close it.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to write.
    path : Path
        Destination file; the suffix selects the format.
    dpi : int, optional
        Resolution for raster formats.
    close : bool, optional
        Close the figure after saving to free memory in batch runs.

    Returns
    -------
    Path
        The path written.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(fig)
    return path


def base_parser(description: str) -> argparse.ArgumentParser:
    """
    Argument parser with the options common to every study script.

    Parameters
    ----------
    description : str
        Text shown by ``--help``; normally the script's docstring.

    Returns
    -------
    argparse.ArgumentParser
        A parser with ``--output`` and ``--format`` already defined.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--output", type=Path, default=FIGURES, help="figure directory")
    parser.add_argument("--format", default="pdf", help="figure file format")
    return parser
