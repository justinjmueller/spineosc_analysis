# Duplicate-event study (all v1.1.0 samples), 2026-10-10

Scripts: `plotting/syst_study/dup_study.py` (counts), `dup_ratio_impact.py` (effect on the 20 production ratios), `dup_report.py` (tables, figure). Numbers: `duplicates_all.json`, `duplicates_subruns.json`, `duplicates_ratio_impact.json`, `duplicates_data.json`. Figure: `dup_events_excess.png`.
Samples: every tree (`events`, `selected`, `signal`) of the v1.1.0 `nosyst`, `wsyst`, `wsyst_cvrw`, `detvar_nosyst`, `detvar_wsyst`, `dirt_nosyst`, `dirt_wsyst` files (157 trees; the `wsyst` copies have the same duplicate counts as `nosyst`), plus the ICARUS data files and the SBND release detector-access file (SBND data: index columns only, no energy read).

## Definitions
* K1 = (Run, Subrun, Evt). K2 = K1 plus round(1e6 &times; event_max_neutrino_energy), the key `detsys/calculator` uses.
* Excess rows = rows &minus; unique keys. "Energy breaks the degeneracy" = (K1 excess &minus; K2 excess) / K1 excess.
* Full-row duplicate = every stored branch identical (events trees: all 19 columns; selected trees: all reco_/true_/event_ columns).
* `events` trees have one row per event; `selected` and `signal` trees have one row per interaction, so several rows per event are legitimate there. `event_max_neutrino_energy` is an event-level quantity, so K2 cannot separate interactions within one event.

## 1. events trees (one row per event)
| file | tree | rows | K1 excess | K2 excess | K1 degeneracy broken by energy | full-row identical excess | K2-duplicate rows identical |
|---|---|---|---|---|---|---|---|
| icarus_run2_detvar_nosyst | icarus_detvar_charge_scale_run2 | 1,295,142 | 300 (0.023%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run2_detvar_nosyst | icarus_detvar_sce0x_run2 | 1,453,190 | 300 (0.021%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run2_detvar_nosyst | icarus_detvar_sce2x_run2 | 1,264,736 | 250 (0.020%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run2_detvar_nosyst | icarus_wiremod_xthetaxw_run2 | 1,809,767 | 550 (0.030%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run2_detvar_nosyst | icarus_wiremod_yz_run2 | 1,786,548 | 450 (0.025%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run2_dirt_nosyst | icarus_dirt_run2 | 490,534 | 300 (0.061%) | 300 | 0.0% | 300 | 600 of 600 |
| icarus_run2_nosyst | icarus_offbeam_run2 | 1,282,878 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| icarus_run2_nosyst | icarus_run2 | 1,649,739 | 400 (0.024%) | 0 | 100.0% | 0 | 0 of 0 |
| icarus_run4_detvar_nosyst | icarus_detvar_charge_scale_run4 | 2,967,001 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| icarus_run4_detvar_nosyst | icarus_detvar_sce0x_run4 | 2,659,861 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| icarus_run4_detvar_nosyst | icarus_detvar_sce2x_run4 | 2,690,681 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| icarus_run4_detvar_nosyst | icarus_wiremod_xthetaxw_run4 | 2,625,102 | 13,878 (0.529%) | 13,878 | 0.0% | 13,878 | 27,756 of 27,756 |
| icarus_run4_detvar_nosyst | icarus_wiremod_yz_run4 | 3,481,280 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| icarus_run4_dirt_nosyst | icarus_dirt_run4 | 1,389,297 | 181,780 (13.084%) | 181,780 | 0.0% | 181,780 | 363,560 of 363,560 |
| icarus_run4_nosyst | icarus_offbeam_run4 | 208,989 | 220 (0.105%) | 220 | 0.0% | 138 | 266 of 435 |
| icarus_run4_nosyst | icarus_run4 | 2,806,621 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_charge_scale | 6,306,804 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_cv | 1,627,100 | 98,210 (6.036%) | 0 | 100.0% | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_dent | 1,579,145 | 90,356 (5.722%) | 0 | 100.0% | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_m5_gain | 6,320,856 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_mattmod_yz | 6,106,384 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_newcv | 4,006,250 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_newdent | 3,800,142 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_p5_gain | 6,078,117 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_sce0x | 6,303,240 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_sce2x | 6,257,494 | 164 (0.003%) | 164 | 0.0% | 164 | 328 of 328 |
| sbnd_detvar_nosyst | sbnd_detvar_smear15 | 6,309,654 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_detvar_nosyst | sbnd_detvar_wiremod_xthetaxw | 5,556,318 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_detvar_wiremod_yz | 6,346,087 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_detvar_nosyst | sbnd_intime | 2,384,909 | 4,665 (0.196%) | 4,665 | 0.0% | 0 | 0 of 9,330 |
| sbnd_dirt_nosyst | sbnd_dirt | 1,330,763 | 1,250,763 (93.988%) | 33,924 | 97.3% | 1,385 | 1,442 of 67,086 |
| sbnd_nosyst | sbnd_large | 6,341,202 | 208 (0.003%) | 208 | 0.0% | 208 | 416 of 416 |
| sbnd_nosyst | sbnd_offbeam | 1,643,131 | 0 (0.000%) | 0 | n/a | 0 | 0 of 0 |
| sbnd_nosyst | sbnd_small | 157,812 | 114 (0.072%) | 0 | 100.0% | 0 | 0 of 0 |

Four distinct situations (the table above is for the `nosyst`/`detvar_nosyst`/`dirt_nosyst` files):

1. **Index collisions fully resolved by the energy** (K2 excess 0, 100% broken): ICARUS Run 2 nominal (400 of 1.65 M events, 0.024%) and its detvars (250&ndash;550), SBND `sbnd_small` (114), and the older SBND `sbnd_detvar_cv` / `sbnd_detvar_dent` pair (98,210 = 6.0% and 90,356 = 5.7%). These are different events that share a Run/Subrun/Evt, so adding the energy is the correct fix. The current SBND pair `newcv` / `newdent` and the ICARUS Run 4 nominal and all but one Run 4 detvar have no collisions at all.
2. **Genuine duplicates, not separable by energy** (K2 excess = K1 excess, rows identical in every branch):
   * ICARUS Run 4 `wiremod_xthetaxw`: 13,878 duplicated events (0.529%), all 27,756 rows identical. Five subruns: 12968 (7,065 of 21,804 events duplicated), 12966 (2,893 of 10,573), 12970 (2,223 of 23,545), 12965 (1,637 of 8,357), 12967 (60 of 120). No other Run 4 detvar has any.
   * ICARUS Run 4 dirt: **181,780 duplicated events (13.08%)**, all identical, in 41 of 191 subruns; in 27 of them exactly half of the rows are copies, i.e. every event appears twice.
   * ICARUS Run 2 dirt: 300 (0.061%), identical, 6 subruns with exactly 50 copies each.
   * SBND `sbnd_large` and every SBND variation sharing its generation: 208 duplicated events (0.0033%), identical, in 10 subruns where every event appears twice (164 in SCE2x). They sit in the nominal and the variations alike, so they cancel in ratios.
3. **Not resolved because there is no neutrino** (off-beam and in-time cosmics have no neutrino energy, so K2 = K1): SBND in-time 4,665 collisions (0.196%), all distinct events; ICARUS Run 4 off-beam 435 rows in duplicated groups (266 identical, 169 differing only in `event_nreco`); SBND off-beam and ICARUS Run 2 off-beam: none.
4. **SBND dirt:** 93.99% of events share their index with another event (all Run &isin; {1, 2, 3, 300, 301}); the energy resolves 97.3% of that, leaving 33,924 collisions of which only 1,442 of 67,086 rows are identical. The rest differ in `event_nreco` / `event_ntrue`: the same generated neutrino reconstructed more than once, not duplicated events.

**Answer to "how much does the energy break the degeneracy":** completely when the collision is between different events (100% in ICARUS Run 2, SBND small/cv/dent), 97.3% in SBND dirt, and not at all for genuine duplicates (0%) or neutrino-free events (0%).

## 2. selected trees (one row per reco interaction)
| file | tree | rows | K2 excess (distinct interactions in one event + identical copies) | identical-row excess | share of K2 excess that is identical |
|---|---|---|---|---|---|
| icarus_run2_detvar_nosyst | icarus_detvar_charge_scale_run2 | 95,839 | 61 (0.064%) | 0 | 0.0% |
| icarus_run2_detvar_nosyst | icarus_detvar_sce0x_run2 | 107,966 | 60 (0.056%) | 0 | 0.0% |
| icarus_run2_detvar_nosyst | icarus_detvar_sce2x_run2 | 95,196 | 61 (0.064%) | 0 | 0.0% |
| icarus_run2_detvar_nosyst | icarus_wiremod_xthetaxw_run2 | 136,009 | 86 (0.063%) | 0 | 0.0% |
| icarus_run2_detvar_nosyst | icarus_wiremod_yz_run2 | 133,941 | 87 (0.065%) | 0 | 0.0% |
| icarus_run2_dirt_nosyst | icarus_dirt_run2 | 69 | 0 (0.000%) | 0 | 0.0% |
| icarus_run2_nosyst | icarus_offbeam_run2 | 190 | 0 (0.000%) | 0 | 0.0% |
| icarus_run2_nosyst | icarus_run2 | 123,223 | 71 (0.058%) | 0 | 0.0% |
| icarus_run2_wsyst_cvrw | icarus_run2 | 122,995 | 56 (0.046%) | 0 | 0.0% |
| icarus_run4_detvar_nosyst | icarus_detvar_charge_scale_run4 | 218,010 | 133 (0.061%) | 0 | 0.0% |
| icarus_run4_detvar_nosyst | icarus_detvar_sce0x_run4 | 196,256 | 106 (0.054%) | 0 | 0.0% |
| icarus_run4_detvar_nosyst | icarus_detvar_sce2x_run4 | 200,873 | 105 (0.052%) | 0 | 0.0% |
| icarus_run4_detvar_nosyst | icarus_wiremod_xthetaxw_run4 | 195,608 | 1,107 (0.566%) | 988 | 89.3% |
| icarus_run4_detvar_nosyst | icarus_wiremod_yz_run4 | 257,865 | 159 (0.062%) | 0 | 0.0% |
| icarus_run4_dirt_nosyst | icarus_dirt_run4 | 180 | 1 (0.556%) | 0 | 0.0% |
| icarus_run4_nosyst | icarus_offbeam_run4 | 35 | 0 (0.000%) | 0 | 0.0% |
| icarus_run4_nosyst | icarus_run4 | 207,974 | 114 (0.055%) | 0 | 0.0% |
| icarus_run4_wsyst_cvrw | icarus_run4 | 207,570 | 84 (0.040%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_charge_scale | 546,065 | 2,554 (0.468%) | 12 | 0.5% |
| sbnd_detvar_nosyst | sbnd_detvar_cv | 141,252 | 707 (0.501%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_dent | 132,870 | 650 (0.489%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_m5_gain | 538,981 | 2,516 (0.467%) | 12 | 0.5% |
| sbnd_detvar_nosyst | sbnd_detvar_mattmod_yz | 530,526 | 2,474 (0.466%) | 12 | 0.5% |
| sbnd_detvar_nosyst | sbnd_detvar_newcv | 292,518 | 1,433 (0.490%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_newdent | 268,746 | 1,272 (0.473%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_p5_gain | 515,206 | 2,367 (0.459%) | 13 | 0.5% |
| sbnd_detvar_nosyst | sbnd_detvar_sce0x | 539,123 | 2,501 (0.464%) | 12 | 0.5% |
| sbnd_detvar_nosyst | sbnd_detvar_sce2x | 555,026 | 2,643 (0.476%) | 9 | 0.3% |
| sbnd_detvar_nosyst | sbnd_detvar_smear15 | 545,163 | 2,564 (0.470%) | 10 | 0.4% |
| sbnd_detvar_nosyst | sbnd_detvar_wiremod_xthetaxw | 478,125 | 2,198 (0.460%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_detvar_wiremod_yz | 541,289 | 2,530 (0.467%) | 0 | 0.0% |
| sbnd_detvar_nosyst | sbnd_intime | 1,356 | 0 (0.000%) | 0 | 0.0% |
| sbnd_dirt_nosyst | sbnd_dirt | 624 | 1 (0.160%) | 0 | 0.0% |
| sbnd_nosyst | sbnd_large | 551,732 | 2,564 (0.465%) | 12 | 0.5% |
| sbnd_nosyst | sbnd_offbeam | 619 | 2 (0.323%) | 0 | 0.0% |
| sbnd_nosyst | sbnd_small | 12,555 | 53 (0.422%) | 0 | 0.0% |
| sbnd_wsyst_cvrw | sbnd_large | 548,660 | 2,315 (0.422%) | 12 | 0.5% |

* In every MC sample 0.05&ndash;0.5% of selected rows share K2 with another row (SBND 0.46%, ICARUS Run 2 0.06%, Run 4 0.055%). These are **not duplicates**: the rows differ (several interactions in one event). The calculator counts all of them.
* Fully identical rows are rare: 9&ndash;13 per SBND tree (from the 208 duplicated events), none in most ICARUS trees.
* **The exception is ICARUS Run 4 `wiremod_xthetaxw`**: 988 identical excess rows (0.51% of 195,608), 1,968 of 2,210 rows in K2-duplicated groups identical (89%); the rest (242) is at the normal multi-interaction level of its siblings (200&ndash;320). The `signal` tree is affected the same way (1,115 identical excess rows).
* The duplicated events of the ICARUS Run 4 and Run 2 **dirt** samples are not duplicated in their `selected` trees (2 identical of 180 rows; 30 selected rows belong to duplicated events).

## 3. Effect on the 20 production ratio histograms
Counting every row versus removing rows identical in every branch (after the pool restriction and the calculator cut):

* 19 of 20 pairs change by at most 0.01% in any bin and 0.00&ndash;0.01% integrated (SBND 7 identical rows on each side, which cancel).
* **ICARUS Run 4 WireModxThetaXW: integrated +1.07% becomes +0.52%, maximum single-bin shift 0.67%** (816 identical rows removed from the varied side, 0 from the nominal).

## 4. Where else the duplicates matter (bookkeeping, not ratios)
* **ICARUS Run 4 dirt normalization.** The POT of the dirt sample (histogram and `pot` in the XML, 1.774677e21) equals the sum of `event_pot` over all event rows including the copies; without the copies it is 1.5424e21 (a factor 1.1506 smaller). The `selected` tree is not doubled, so the dirt prediction would be scaled with too much POT and come out low by about 13%. Dirt is 180 of about 208,000 selected Run 4 MC rows, so the effect on the total is of order 0.01%, but the bookkeeping is wrong.
* **ICARUS Run 4 data exposure** (`icarus_run4_v1.1.0_bnblight.root`): 891 identical event rows (0.57% of 155,321; also 769 collisions of distinct events). The file's POT (1.8873e19) and livetime (4.1722e6) include them; without them 1.8783e19 and 4.1525e6 (&minus;0.48%). The `selected` tree has no identical rows.
* ICARUS Run 4 varied-sample POT for `wiremod_xthetaxw` is overcounted by 0.5% as well (1.6593e21 vs 1.6506e21).
* Other data: SBND release detector-access file: events tree 0 excess; selected 174 (0.27%) rows sharing an index, all distinct rows (several interactions). ICARUS release Run 2: 12 shared indices, distinct. ICARUS Run 2 `bnblight`: 1.

## 5. Conclusions
1. The calculator key (K2) is adequate for pairing CV and variation events wherever the collisions are between different events, which is every MC sample except those with true duplicates and the neutrino-free samples.
2. True duplicates (identical rows, K2 cannot help) exist in exactly these samples: ICARUS Run 4 `wiremod_xthetaxw` (0.53%), ICARUS Run 4 dirt (13.1%), ICARUS Run 2 dirt (0.06%), SBND `sbnd_large` and its variations (0.003%), ICARUS Run 4 off-beam (138 excess rows, 0.07%), ICARUS Run 4 data (0.57%).
3. They come in whole-subrun blocks (every event of a subrun twice, or part of one), which points at input files listed or processed twice, not at random event-level duplication.
4. Only one of them changes an analysis input in a measurable way through the ratios (ICARUS Run 4 WireModxThetaXW, +0.55% integrated, up to 0.67% per bin); two change exposure/normalization bookkeeping (ICARUS Run 4 dirt POT +15%, ICARUS Run 4 data POT/livetime +0.48%).
5. Recommended: find why those subruns are doubled (file lists), regenerate with de-duplication, recompute the Run 4 dirt and data exposures, and add an identical-row check (hash of all branches) to production; in the calculator, drop identical rows and, for pairing, add an interaction-level key to K2.

Caveats: full-row identity is evaluated on the branches stored in these files; "identical" does not prove that the underlying input files were duplicated. The dirt POT conclusion assumes that PROfit scales the dirt sample by the XML `pot` attribute, which I did not run; the selected tree being single-counted is checked on 180 rows (30 in duplicated events, 2 identical rows).

## 6. Do the duplicates occur in streaks? (follow-up)
Script `plotting/syst_study/dup_streaks.py`, `dup_streak_fig.py`; numbers `duplicates_streaks.json`; figure `dup_streaks.png`. Method: in storage order, every row identical to an earlier row is mapped to its first occurrence; a streak is a maximal run of consecutive repeated rows whose originals are also consecutive (constant offset), i.e. a contiguous block of the tree copied verbatim.

**Yes. In every sample with true copies the repeated rows are almost entirely in streaks, and the streaks are whole input files.**

| tree | repeated rows | streaks | rows in streaks of >= 10 | streak length | where the copies sit |
|---|---|---|---|---|---|
| ICARUS Run 4 dirt | 181,780 | 3,431 | 99.9% | mostly exactly 50 (183 of the first 200), the rest 100, 150, ... 350; 96.5% multiples of 50 | in bands at about 24%, 64&ndash;72%, 80&ndash;82% and 87&ndash;100% of the tree (89% of the last 5% slice is repeats); each band is an in-order copy of a contiguous region of the first 24% of the tree |
| ICARUS Run 4 WireModxThetaXW | 13,878 | 260 | 99.9% | 178 of the first 200 exactly 50, 13 exactly 100, one 250, a few partial (6&ndash;43, 70); 96% multiples of 50 | entirely inside the first 1.8% of the tree (entries 1,725&ndash;46,499), runs 12965&ndash;12970 |
| ICARUS Run 2 dirt | 300 | 6 | 100% | exactly 50, all starting on the 50-event grid of their subrun | last 2% of the tree |
| SBND large | 208 | 10 | 100% | 15&ndash;26, each a whole (Run, Subrun) segment (e.g. Run 3554 Subrun 81, Evt 1&ndash;98, twice) | 80&ndash;95% of the tree, adjacent to (offset 24) or within about 1,200 entries of the original |
| ICARUS Run 4 off-beam | 138 | 27 | 0% (all streaks of 3&ndash;9 rows) | segment-sized (median segment 5 rows) | scattered |
| ICARUS Run 4 data (bnblight) | 891 | 186 | 0% (all streaks of 5 or fewer rows) | segment-sized (median segment 5 rows) | scattered |

* The block size matches the file size: the dirt and detvar trees are concatenations of 50-event segments (median (Run, Subrun) segment length 50), and the streaks are 50, 100, 150, ... rows long with identical Evt ranges. Streak starts coincide with the 50-event grid of their subrun for 61&ndash;100% of the copies and 53&ndash;100% of the originals (it cannot be 100% when partial files shift the grid).
* Copies always come after their originals in the tree, always carry the same Run and Subrun as the original (100%), and are bit-identical in every branch.
* In the Run 4 dirt sample the trees are not sorted by run: one (Run, Subrun) is spread over many separated 50-event chunks, and the copies form several contiguous bands that each replay a stretch of the first 24% of the file order. That is what a file list that contains the same stretch of files two or three more times looks like.
* The Run 4 WireModxThetaXW duplicates are confined to the first 1.8% of the tree, i.e. the first stretch of files was listed or processed twice.
* The off-beam and data duplicates (138 and 891 rows) are not streaky at the file scale (streaks of at most 9 and 5 rows) but their streak length equals their segment length (5), so they are consistent with whole-segment copies of 5-row files; with five rows per file that is weaker evidence.

**Interpretation (supported, not proven):** the duplicates are copies of whole input files (50 events for ICARUS MC, 15&ndash;26 for SBND), placed after the originals, which points at repeated entries in the processing file lists (for example resubmitted or doubly merged jobs) and not at event-level corruption. Confirming needs the file list or the job logs of those productions.
