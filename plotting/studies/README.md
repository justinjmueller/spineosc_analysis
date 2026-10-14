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
| `fakedata_injection.py` | `fakedata_injection.ipynb` | injection panels, closure, detectability, pulls |
| `fakedata_channel.py` | `fakedata_channel_studies.ipynb` | single-channel fake-data questions |
| `fakedata_2d.py` | `fakedata_2d_studies.ipynb` | 2D-configuration fake-data fits |

The fake-data scripts take `--section` and `--list-sections`. Figures default to `figures/` (or the
study's own folder), never to the technote; copy reviewed ones across by hand.

The original notebooks are still in `plotting/` until the scripts' output has been checked against them.
