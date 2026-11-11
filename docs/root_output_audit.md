# Root-level ROOT output audit (2026-10-08)
Method: each `*.root` in the repo root was matched (1) against the file-name templates read by `plotting/*.ipynb` and `detsys/plotting/*.ipynb`, and (2) against the `tag =` of every ini under `xml/` and `jobs/`, then against the status in `xml/MANIFEST.toml`.
Nothing was moved or deleted. The per-file lists are in `root_output_audit.json`.

| Class | Files | MB |
|---|---|---|
| A keep: read by a technote notebook | 31 | 13 |
| B keep (frozen): source of the frozen SBND data-access section | 375 | 1958 |
| C candidate: output of an archived coverage XML | 11 | 26 |
| D candidate: superseded dataaccess_* tags (replaced by da2_/dent_ tags in runs/) | 117 | 17 |
| F candidate (lower confidence): no XML or notebook in the repo produces or reads this tag | 63 | 18 |

Notes
* Class F is lower confidence: it contains the `*_test0913/0914*` development runs, A/B control runs (`muon_score_ab*`, `gateflip`), `sbndaccess_*` runs from earlier phases of the frozen section, `chi2test`, `PROfit_v1`. Some of the `sbndaccess_*` files may still be needed to reproduce the frozen section; check before deleting.
* The ROOT outputs total about 2 GB and are not where the disk goes. The cached `.bin` files at the root are about 62 GB (about 32 GB under tags no ini produces, 26 GB with no matching ROOT tag) and `runs/` is about 120 GB (`fakedata_full` 44 GB, `cov_v2` 35 GB, `da_v2` 27 GB).

Update 2026-10-08: classes C and D (128 root `.root` files) and the 99 `.bin` caches of their tags (34 GB) were deleted at the user's request. The list is in `root_output_deleted_2026-10-08.txt`. Classes F and A/B were left untouched.

Update 2026-10-08 (second pass): class F (63 root files) and the tier 1 and 2 `.bin` caches (9 files, 3.9 GB) were deleted at the user's request; see `root_output_deleted_2026-10-08.txt` and `root_bin_candidates.md`. Remaining: classes A and B, tiers 3 to 5 of the `.bin` caches.
