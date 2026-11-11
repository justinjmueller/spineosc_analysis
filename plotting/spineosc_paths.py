"""
Machine-specific locations shared by the plotting and study scripts.

The release ntuples live outside the repository. Point ``SPINEOSC_DATA`` at the directory that
holds the release folders (default ``~/data/spineosc``); the release itself is ``RELEASE``::

    export SPINEOSC_DATA=/path/to/data/spineosc     # contains v1.1.0/spineosc_numu_*.root

``REPO`` is the repository root, found from this file, so scripts never hard-code it.
"""

import os
from pathlib import Path

#: Repository root (this file lives in ``plotting/``).
REPO = Path(__file__).resolve().parents[1]

#: Selection release whose ntuples the scripts read.
RELEASE = "v1.1.0"

#: Directory holding the release folders (override with ``SPINEOSC_DATA``).
DATA_ROOT = Path(os.environ.get("SPINEOSC_DATA", Path.home() / "data" / "spineosc")).expanduser()

#: Directory of the ``RELEASE`` ntuples, e.g. ``spineosc_numu_sbnd_v1.1.0_wsyst.root``.
DATA_DIR = DATA_ROOT / RELEASE
