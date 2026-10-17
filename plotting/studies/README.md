# plotting/studies

Documented, command-line versions of the analysis notebooks in `plotting/`. Each script has a module
docstring (what it makes, which inputs it reads, which technote figures it feeds), a `main()` and
`--help`. Shared paths, versions and figure helpers live in `common.py`; style rules are in
`CONVENTIONS.md` (mutools style, ruff E/F/I/B/UP at line length 100).

Run with the mutools environment from the repository root:

```
/Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/<script>.py --help
```

| Script | Replaces | Makes |
|---|---|---|
| `signal_selection.py` | `signal.ipynb`, efficiency cells of `notebook.ipynb` | signal/selected histograms, efficiency and purity |
| `resolution.py` | `resolution.ipynb` | resolution stacks and overlays |
| `systematics.py` | `systematics.ipynb` | per-category uncertainty breakdowns, near/far ratio |
| `surface.py` | `surface.ipynb` | sensitivity contour surfaces |
| `exposure.py` | `summarize_exposure.ipynb` | equivalent-POT table |
| `data_access.py`, `data_access_phases.py` | `data_access_results.ipynb` | data/MC panels, p-values, conditional overlays (`--phase 1-5`) |
| `sbnd_scan.py` | `sbnd_data_access.ipynb` | oscillation-scan spectra, deficit, bins outside the CV band |
| `conditional_constraint.py` | `conditional_constraint.ipynb` | SBND -> ICARUS conditional constraint |
| `dent_study.py` | `dent_study.ipynb` | DENT data/MC study tables and figures |

Figures default to `figures/` (or the
study's own folder), never to the technote; copy reviewed ones across by hand.

Checked against their notebooks (2026-10-17, headless runs compared pixel by pixel): `exposure`, `signal_selection`,
`systematics`, `surface` and `resolution` reproduce the notebooks' figures and tables exactly, so those notebooks are
retired (recoverable from commit `d888bdd`). `resolution.py` uses the house saver settings of `common.configure_style`
(300 dpi, rasterized, `ytick_precision=2`), which differ from the resolution notebook's defaults; its data are identical.
Still to check: `sbnd_scan`, `conditional_constraint` (missing PROfit inputs), `data_access*` and `dent_study`, whose
notebooks remain in `plotting/`.
