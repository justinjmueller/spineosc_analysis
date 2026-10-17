"""
Equivalent-POT summary of the v1.1.0 samples.

Tabulates, for each detector (SBND, ICARUS Run 2, ICARUS Run 4) and
sample (on-beam bnblight, dirt, simulation, off-beam), the equivalent
POT and its ratio to the on-beam sample. Replaces
``plotting/summarize_exposure.ipynb``. It makes no figures; it prints
the table (and optionally writes it as CSV).

Inputs are the ``spineosc_numu_<detector>_<tag>_*.root`` selection files
under ``--data-dir`` (default ``$SPINEOSC_DATA/<tag>``, ``~/data/spineosc/<tag>`` if unset), read from the
``events/<key>/events`` tree.

Equivalent POT:

* 'mc' (nominal, dirt, ...) and 'onbeam' (bnblight) samples carry a
  real, directly meaningful POT in ``event_pot``, so their equivalent POT
  is just that sum.
* 'offbeam' (EXT) data has no beam and so no POT of its own (``event_pot``
  sums to 0); it is normalized by gate counting instead: the on-beam
  sample's POT is scaled by the ratio of off-beam gates
  (``event_noffbeambnb``) to the on-beam sample's true BNB gates
  (``event_nbnb``), i.e. the off-beam POT "as if" it had been collected
  at the on-beam POT-per-gate rate. Not every production populates those
  columns (ICARUS's off-beam sample has both sum to 0); then
  ``event_gate_delta`` (off-beam) over ``event_unfolded_nbnb`` (on-beam)
  carries the same information and is used instead.

ICARUS's on-beam/off-beam samples are split by run within their files;
each ICARUS run is summarized under its own key.

Example::

    python plotting/studies/exposure.py --csv exposure_v1.1.0.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

import common
from common import SELECTION_VERSION

#: Default data root; the selection files live in ``<root>/<tag>``.
DATA_ROOT = common.DATA_ROOT

#: Sample name -> how it is normalized ('mc', 'onbeam' or 'offbeam').
SAMPLE_TYPES = {
    "bnblight": "onbeam",
    "dirt": "mc",
    "sim": "mc",
    "offbeam": "offbeam",
}

#: File-name suffix (after ``<detector>_<tag>``) of each sample.
FILE_SUFFIX = {
    "bnblight": "bnblight",
    "dirt": "dirt_wsyst",
    "sim": "wsyst",
    "offbeam": "nosyst",
}

#: Detectors: name -> (file stem, top-level key of each sample within its
#: file). SBND has a single run; ICARUS samples are split by run.
DETECTORS = {
    "sbnd": (
        "sbnd",
        {
            "bnblight": "sbnd_bnblight",
            "dirt": "sbnd_dirt",
            "sim": "sbnd_large",
            "offbeam": "sbnd_offbeam",
        },
    ),
    "icarus_run2": (
        "icarus_run2",
        {
            "bnblight": "icarus_bnblight_run2",
            "dirt": "icarus_dirt_run2",
            "sim": "icarus_run2",
            "offbeam": "icarus_offbeam_run2",
        },
    ),
    "icarus_run4": (
        "icarus_run4",
        {
            "bnblight": "icarus_bnblight_run4",
            "dirt": "icarus_dirt_run4",
            "sim": "icarus_run4",
            "offbeam": "icarus_offbeam_run4",
        },
    ),
}


def equiv_pot(
    events: pd.DataFrame,
    sample_type: str,
    onbeam_events: pd.DataFrame | None = None,
) -> float:
    """
    Equivalent POT of a sample's ``events`` tree.

    Parameters
    ----------
    events : pandas.DataFrame
        The sample's events tree.
    sample_type : {"mc", "onbeam", "offbeam"}
        How the sample is normalized; see the module docstring.
    onbeam_events : pandas.DataFrame, optional
        The bnblight events tree; required for ``"offbeam"``.

    Returns
    -------
    float
        The equivalent POT.

    Raises
    ------
    ValueError
        If ``sample_type`` is unknown, or an off-beam sample is given
        without ``onbeam_events``.
    """
    sample_type = sample_type.lower()

    if sample_type in ("mc", "onbeam"):
        return events["event_pot"].sum()
    if sample_type == "offbeam":
        if onbeam_events is None:
            raise ValueError(
                "onbeam_events (the bnblight sample's events tree) is "
                "required to compute the equivalent POT for an offbeam sample"
            )
        onbeam_pot = onbeam_events["event_pot"].sum()
        onbeam_nbnb = onbeam_events["event_nbnb"].sum()
        offbeam_gates = events["event_noffbeambnb"].sum()
        if offbeam_gates == 0:
            onbeam_nbnb = onbeam_events["event_unfolded_nbnb"].sum()
            offbeam_gates = events["event_gate_delta"].sum()
        return onbeam_pot * (offbeam_gates / onbeam_nbnb)
    raise ValueError(
        f"Unknown sample_type: {sample_type!r} (expected 'mc', 'onbeam', or 'offbeam')"
    )


def load_events(data_dir: Path, tag: str) -> dict[str, dict[str, pd.DataFrame]]:
    """
    Read the events tree of every (detector, sample) pair.

    Parameters
    ----------
    data_dir : Path
        Directory holding the selection files.
    tag : str
        Selection version used in the file names.

    Returns
    -------
    dict
        ``events[detector][sample]`` as a DataFrame.
    """
    import uproot

    events: dict[str, dict[str, pd.DataFrame]] = {}
    for detector, (stem, keys) in DETECTORS.items():
        events[detector] = {}
        for name, key in keys.items():
            path = data_dir / f"spineosc_numu_{stem}_{tag}_{FILE_SUFFIX[name]}.root"
            with uproot.open(path) as rf:
                events[detector][name] = rf[f"events/{key}/events"].arrays(library="pd")
    return events


def summarize(events: dict[str, dict[str, pd.DataFrame]]) -> pd.DataFrame:
    """
    Equivalent-POT table.

    Parameters
    ----------
    events : dict
        Output of ``load_events``.

    Returns
    -------
    pandas.DataFrame
        One row per (detector, sample) with columns ``detector``,
        ``sample``, ``type``, ``equiv_pot`` and ``pot_scale_to_bnblight``
        (the equivalent POT relative to the same detector's bnblight).
    """
    summary = pd.DataFrame(
        [
            {
                "detector": detector,
                "sample": name,
                "type": SAMPLE_TYPES[name],
                "equiv_pot": equiv_pot(
                    df,
                    SAMPLE_TYPES[name],
                    onbeam_events=(
                        det_events["bnblight"] if SAMPLE_TYPES[name] == "offbeam" else None
                    ),
                ),
            }
            for detector, det_events in events.items()
            for name, df in det_events.items()
        ]
    )
    bnblight_pot = summary[summary["sample"] == "bnblight"].set_index("detector")["equiv_pot"]
    summary["pot_scale_to_bnblight"] = summary.apply(
        lambda row: row["equiv_pot"] / bnblight_pot[row["detector"]], axis=1
    )
    return summary


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str], optional
        Argument list; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.add_argument("--tag", default=SELECTION_VERSION, help="selection version")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help="selection-file directory (default $SPINEOSC_DATA/<tag>, else ~/data/spineosc/<tag>)")
    parser.add_argument("--csv", type=Path, default=None, help="also write the table as CSV")
    args = parser.parse_args(argv)

    data_dir = args.data_dir or DATA_ROOT / args.tag
    summary = summarize(load_events(data_dir, args.tag))
    print(summary.to_string(index=False))
    if args.csv is not None:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.csv, index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
