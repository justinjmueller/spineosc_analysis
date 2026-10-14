# Conventions for plotting/studies scripts

Each script replaces one or more notebooks in `plotting/` with a documented module.

* Style follows mutools (`/Users/mueller/Projects/mutools/src/mutools/plotting/*.py`): module docstring
  that says what the script makes, which inputs it reads and which technote figures it feeds;
  `from __future__ import annotations`; numpy-style docstrings (Parameters/Returns/Raises) on every
  function; type hints; module constants documented with `#:` comments; double quotes; line length
  100; ruff rules E,F,I,B,UP (run `/Users/mueller/Projects/mutools/.venv/bin/ruff check --line-length
  100 --select E,F,I,B,UP <file>` and fix all findings).
* Structure: constants/config at the top, small pure functions (load -> compute -> draw), `main(argv=None)`,
  `if __name__ == "__main__": main()`. Use `argparse` via `common.base_parser`. Every study that has
  hand-set switches in the notebook (tiers, y-limits, which detector, ...) exposes them as CLI options
  or module constants, never as hidden globals. No work at import time. No `display()`/`show`.
* Use `common.py` for REPO/RUNS/LEGACY_ROOT/FIGURES paths, versions and `configure_style`/`save_figure`.
  No absolute `/Users/...` paths in the scripts except through `common` (derive from `REPO`).
* Markdown cells in the notebook become docstrings or comments: keep the physics explanation
  (what a plot shows, caveats, how to read tables), do not keep run-log chatter. Preserve the numbers,
  cuts, binning, labels and file names exactly; this is a refactor, not a re-analysis.
* Do not change the notebook, delete anything, or commit. Do not launch long jobs (PROfit, condor,
  anything > ~2 min). Verify with ruff, `python <script> --help`, `python -c "import ast"`-level parse,
  and, where the input files exist, one quick real invocation (limit to one panel/study via an option).
  Run python with `/Users/mueller/Projects/mutools/.venv/bin/python`.
* Report back: files created, CLI summary, what was verified by running vs only linted, and anything
  in the notebook you could not port or found broken (e.g. missing input paths). Keep the report short.
