# Removed: results that used disallowed detector-data-access variables

2026-10-14, on the analyst's instruction. The detector data access does not allow muon length, energy or any energy proxy, or the
angle to the beam axis. The first DENT study had run data/MC comparisons for four observables that fall outside the allowed list:

* `muon_len3d` (straight-line muon length from start and end points)
* `muon_dir_x`, `muon_dir_y`, `muon_dir_z` (direction components; `dir_z` is the polar angle to the beam axis)

None was in the primary set, and no decision in the study depended on them.

## What was removed (308 files; no result values are recorded here)

| Location | Files | Notes |
|---|---|---|
| `runs/dent_study/` | 228 | PROfit outputs, `.bin` caches, results JSON, queue files for the four observables |
| `detsys/profit_systematics_dd/` | 44 | per-systematic ratio files for the four observables |
| `xml/dent_study/sbnd_run1/dd/` | 28 | study XML and ini files, including the Tier A variants (git rm) |
| `detvar/` | 4 | `detsys_dd_*` ratio files (git rm) |
| `work/` | 4 | symlinks to the files above |
| `plotting/dent_validation/` | 2 figures | `fig2_cathode_features.png` (length panels), `fig4_other_variables.png` |
| `docs/dent_study/figures/` | 1 figure | `mechanism_cathode_features.png` (copy of the former) |

## What was edited

* `scripts/dent_study/observables.json` and `make_observables.py`: 21 -> 17 observables. `primary_selection.json`: rows removed.
* `plotting/dent_validation/`: `dentval.py` no longer computes a muon length; `fig134.py` labels removed; `fig2_cathode.py`
  length panels (a) and (b) removed; `scan_all_variables.csv` rows for the four variables removed.
* `plotting/dent_study/report_cells.py`: the cathode-features figure removed. `report.html` and `plotting/dent_study.ipynb`
  were rebuilt (0 failed cells).
* `input/sources.list` and `input/sources.list.example`: the four `detsys_dd_*` entries removed.
* `PREREGISTRATION.md`: the note about `muon_len3d` removed and a deviation entry added.

## Not removed, for the analyst to decide

* **Git history** still contains the 32 tracked files (removals are staged, not committed). Removing them from history needs a
  history rewrite, which was not done.
* **`archive/snapshot_2026-10-08/`** (git-ignored) holds 32 copies of the same files.
* **`mechanism_missing_length.png`** is a simulation-only (CV against DENT) figure about lost track length and uses no data.
* **Azimuth-slice studies** (`xml/dent_study/sbnd_run1/*azimuth*`) build the azimuth from `dir_x` and `dir_y`; they never use `dir_z`.
  They should switch to the allowed `reco_leading_muon_azimuthal_angle` branch.
* **`event_largest_flash_time`** is in the validation scan (`scan_all_variables.csv`). It is neither length, energy nor direction,
  but it is not on the allowed list either.
* Two worker logs in `runs/dent_study/logs/` contain only job start/end lines naming the removed studies.
