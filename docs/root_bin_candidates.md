# Remaining root-level .bin caches (102 files, 31 GB) and deletion candidates (2026-10-08)
`.bin` files are PROfit caches (`*_prop.bin`, `*_syst.bin`, `*_data.bin`). They let `plot` rerun without reprocessing and can be regenerated with PROfit `process`.
Nothing here has been deleted.

| Tier | Files | GB | What | Recommendation |
|---|---|---|---|---|
| 1 | 3 | 2.5 | `dp25_test09242026_{data,prop,syst}.bin`: cache of an archived coverage XML whose ROOT outputs were already deleted (tag differs by the `_vanilla`/`_shifted` suffix, so the earlier pass missed it) | delete |
| 2 | 6 | 1.4 | `sbndaccess_dp_be_*.bin`, `sbndaccess_dp_var0_*.bin`: no XML produces them; their ROOT outputs are class F | delete if class F is deleted |
| 3 | 4+4+4 | 4.3 | caches of `selected`, `signal`, `resolution` (current, active XMLs; minutes to regenerate) | optional; keep if you replot often |
| 4 | 2 | 2.3 | `contours_recoe_{prop,syst}.bin` (read by the conditional-constraint notebook) | keep |
| 5 | 79 | 22.5 | caches of the frozen SBND data-access section (`p4_*` 25 files 10.8 GB, `t1_*` 27 files 6.2 GB, `t2_*` 12 files 2.8 GB, `t3_*` 12 files 2.0 GB, `p5_*` 3 files 0.7 GB) | optional: the ROOT outputs that make the figures are kept, so only quick replotting is lost; hours to regenerate |

Tier 1 and 2 total 3.9 GB. Tiers 1-3 total about 8 GB. Adding tier 5 frees about 31 GB.

Update 2026-10-08: tiers 1 and 2 deleted.
