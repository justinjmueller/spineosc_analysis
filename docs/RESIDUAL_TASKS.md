# Residual tasks

Open items that are not blocking the current work. Add new entries at the bottom; move finished ones to "Done" with the date.

## Open

### Technote plots to update to v1.1.0
Opened 2026-10-07. Technote: `/Users/mueller/Projects/spineosc_technote` (figures under `figures/`, referenced from `sections/*.tex`).

**Method and caveat.** The list comes from the `\includegraphics` calls in the `.tex` files plus the modification dates of the
figure files. The figure files carry no version tag, so "stale" below is inferred from dates older than the v1.1.0 work
(2026-10-06/07). Confirm each group against its producing notebook before regenerating.

**Likely stale (figures last written 2026-05 to 2026-09-30, before the v1.1.0 inputs, cvwgt, crosser cut and new exposures):**

| Section (`.tex`) | Figure dir | Used | Date | Producing notebook | Needs |
|---|---|---|---|---|---|
| `systematics.tex` | `figures/systematics/` | 104 | 09-14 to 09-30 | `plotting/systematics.ipynb` (+ detvar/WireMod/SCE studies) | New staged set in `docs/systematics_v1.1.0/` (60 figures) replaces the PROfit fractional-uncertainty plots. The WireMod, SCE, charge-scale and PID individual-variable studies are separate inputs and are not covered by it. |
| `systematics.tex` | `figures/data_sim_comparisons/datamc_*` | 24 | 09-14 to 09-17 | `plotting/notebook.ipynb` | Regenerate with v1.1.0 XMLs and exposures. |
| `data_sim_comparisons.tex` | `figures/data_sim_comparisons/` | 72 | 09-14 to 09-17 | `plotting/notebook.ipynb` | Same: SBND, ICARUS Run 2 and Run 4 data/MC, single- and multi-proton. |
| `signal_and_selection.tex` | `figures/signal_and_selection/` | 60 | 09-14 to 09-17 | `plotting/signal.ipynb` | Efficiency and purity plots. The SBND selection now includes the crosser veto, so SBND at least must change. The folder also holds an apparent stray copy, `efficiency_icarus_run*_muon_azimuthal_angle (1).pdf`, to remove. |
| `data_quality.tex`, `appendix_data_quality.tex` | `figures/data_quality/`, `figures/beam_quality/` | 6 + 46 | 09-14 to 09-30 | not identified | Exposure and POT summary must match `plotting/summarize_exposure.ipynb`; the per-run quality plots follow the v1.1.0 files. |

**Dated 2026-10-06/07 (regenerated recently); verify the version used rather than assuming:**

| Section | Figure dir | Used | Notes |
|---|---|---|---|
| `detector_data_access.tex` | `figures/detector_data_access/` | 62 | Dated 10-07. Check that the XMLs behind them use cvrw+cvwgt, the SBND crosser cut and the high-y/z veto on data (the original `_data.xml` files omitted the veto in 35 of 39 SBND studies). |
| `fake_data_studies.tex` | `figures/fake_data/` | 106 | Dated 10-07. Check that the inputs are v1.1.0. |
| `sensitivity_studies.tex` | `figures/sensitivity_studies/` | 6 | Dated 10-07. Check against the updated systematics (Binding Energy, MattMod Y-Z, DENT with the crosser cut). |
| `sbnd_data_access.tex` | `figures/sbnd_data_access/` | 35 | Dated 10-06. Check the version. |

**Progress (2026-10-07, evening).** Regenerated with v1.1.0 inputs and the SBND crosser + 4 cm end-point cuts, and copied into the technote working tree (unstaged): `signal_and_selection` (efficiency table
and 54 figures), `systematics` (60 PROfit figures), `detector_data_access` (142 figures; the 28 SBND studies in `runs/da_v2`, ICARUS Run 2 from `runs/dent_study`; the end-x zoom uses a central bin of +-5 cm).
New in the technote: DENT subsection (`systematics.tex`) and the cathode-region subsection (`detector_data_access.tex`) with 7 figures from `docs/cathode_motivation/`. In progress: `data_sim_comparisons`
(`xml/coverage_v2`, `runs/cov_v2`, notebook `detsys/plotting/datamc.ipynb`). Still open: data quality, other systematics studies (WireMod/SCE/charge-scale/PID), verifying fake-data and sensitivity inputs.

**Scan of the technote's figures (2026-10-07, by modification time; SBND data access is exempt).** Fresh (regenerated today): `data_sim_comparisons` (96), `detector_data_access` (66 referenced), 36 of the 48 `signal_and_selection`
figures, 75 of the 102 `systematics` figures. Not yet updated:
* `signal_and_selection`, 12 `resolution_*` figures (6 particle KE, 6 energy; SBND, ICARUS Run 2 and Run 4; 1p and Mp). Producer `plotting/resolution.ipynb` with `xml/resolution/{inclusive,exclusive}.xml`, which still use the old SBND trees, no `cvwgt`, and no crosser/end-point cuts. Needs the same migration as the selected XMLs.
* `systematics`, 27 figures from the detector-model studies: charge-scale bias and fitted scale factor (`*_recomb_relbias`, `*_fitted_*`, `*_varied_relbias`), SCE maps (`SCE_*`, `SBND_0x/1x`, `ICARUS_0x/1x`), WireMod maps and 1D slices (`depositionmaps`, `ratiomaps`, `chi2maps`, `1D_*`), trigger uncertainty (3), CORSIKA off-beam comparison (1). Their producers are outside the repo; the selection-dependent ones (trigger, CORSIKA comparison, varied relbias) need checking against the new SBND selection, the map and SCE studies probably do not.
* `fake_data` (106) and `sensitivity_studies` (6): made before the cuts from the old contours inputs; need the contour reruns (long) after the contours/suite XMLs (already updated).
* `data_quality` (51) and `beam_quality` (1): dated September; check against the v1.1.0 files and the exposure summary. `reconstruction` (1) is an architecture drawing and is unaffected.

**Not affected by the version:** `reconstruction.tex` (`spine_architecture.png`).

**Text to revisit alongside the plots.** `datasets_and_infrastructure.tex` has a commented-out paragraph and a
`spine_prod` table (lines ~9 and ~66-88) that still describe the Spring 2025 production and the earlier weights. The dataset
and exposure tables and the systematics table need the v1.1.0 numbers.

**Binding Energy (done 2026-10-07).** The `BindingEnergy_*` histograms in `detvar/detsys-v1.1.0v1.root` were regenerated from the v1.1.0 `_cvrw` files with `cvwgt` and the v1 selection (`detsys/binding_energy/make_be_v1.py`, `inject_be.C`). The 2D, delta-p and QE-energy Binding Energy histograms used by the `sbnd_data_access` XMLs were not regenerated.

**Update (2026-10-17).** The resolution group is done (see the overnight log). The tracked fake-data plotting scripts (`plotting/studies/fakedata_*.py`, `plotting/scripts/build_fakedata_*`, `fd2d.py`) and the `spinefd` study were moved to `archive/fakedata_2026-10/` (git-ignored), so the `fake_data` figures have no tracked producer; regenerate with `fdgen` plus new plotting code, or restore from the archive. `sensitivity_studies` is produced by `plotting/studies/surface.py` (verified).


### Re-add updated fractional uncertainty plots on the near/far ratio
Opened 2026-10-07.

The old near/far ratio figure (`inclusive_fractional_uncertainty_ratio.pdf`) and its paragraph were removed from `sections/systematics.tex`
(with the old inclusive per-detector figure), because they predated the current analysis and the caption was wrong. The ratio view is still
useful for showing how much the systematics cancel between SBND and ICARUS, so add an updated version later. PROfit already writes the ratio
plot for the current suite (`runs/systematics/systs_recoE_ratio_fractional_systematics.pdf`); it needs checking against the v1.1.0 inputs, labelling per
detector pair (SBND with ICARUS Run 2 and Run 4), and new text.

### Cover the SBND APA gap in the detector data access section
Opened 2026-10-07.

**Status (2026-10-07 overnight).** APA-gap figures (vertex z and muon end z, gap at z = 250 cm confirmed) and text were added to `detector_data_access.tex`. Still open: the angular dependence of the position comparisons.


The conveners recommended checking data/simulation agreement in the vertex and muon end point positions near the SBND cathode **and APA gaps**, including
the dependence on track angle. `sections/detector_data_access.tex` currently covers only the cathode (and says in the Results introduction that the
SBND APA gaps are not covered). Add the SBND APA-gap comparisons (vertex and muon end point position near the gaps, and the angular dependence), then
update that sentence, the Summary table, and the conclusions. Also covered by the same recommendation: the angular dependence of the position comparisons,
which the section does not show either.

### Add a quantitative summary table to the detector data access section
Opened 2026-10-07.

A table of shape-only chi-square and p-values (per variable, per channel, SBND and ICARUS Run 2) was drafted for the summary of
`sections/detector_data_access.tex` and then dropped, because the values were not right. The section is qualitative for now.
Before adding a table: regenerate the figures from the current inputs, decide how the chi-square is defined and interpreted (the covariance is dominated by
correlated systematics, so the tail probability assuming ndf = nbins - 1 is only a guide), and make the table, the figure legends, and the text agree. Then
re-check each qualitative statement in the section against the final numbers. The baseline numbers used for the draft are in `runs/dent_study/results/`
(tier `baseline`, variant `v1`), but they should not be reused without checking the points above.

### Add a larger azimuthal-angle scan for the PID and primary scores
Opened 2026-10-07.

`sections/detector_data_access.tex` splits the muon and leading-proton PID and primary scores into only three azimuthal bins. Add a finer scan in the muon
azimuthal angle for SBND and ICARUS Run 2, for all four scores (muon and leading-proton PID and primary), so that angle-dependent mis-modeling is localized
better. For example, SBND and ICARUS studies with 8 azimuthal bins already exist in the detector data access set (`*_azimuth8` for SBND and
`*_coarse_azimuth8` for ICARUS Run 2) and were run for the DENT study; they could be the starting point. Decide the binning and how to present the many panels (for example
a summary of the level of agreement against angle bin in place of one plot per bin), then update the section text and conclusions. Keep the SBND azimuth convention
(the DENT definition, `atan2(start_dir_x, start_dir_y)`) consistent with the existing three-bin figures.

### Write a solid "impact to the analysis" narrative for the detector data access section
Opened 2026-10-07.

The end of `sections/detector_data_access.tex` has only a short "Implications for the analysis" paragraph, which says the primary-score shape differences lie largely away
from the 0.1 selection threshold and that their effect on the selection efficiency is likely limited but not quantified. Replace it with a complete narrative that connects each finding
to its effect on the analysis: for each disagreement (SBND muon end point near the cathode, SBND muon primary score, leading-proton primary score in both detectors, modest SBND
leading-proton PID tension), say what the observable feeds in the selection, how much of the population sits near the relevant threshold, and what the consequence is for the
selection efficiency, the detector systematic uncertainties (including whether the current variations cover it), and the sensitivity. Support the statements with numbers where they can be
made, for example the fraction of events within a given distance of the 0.1 primary-score threshold, the efficiency change from reweighting simulation to data in the
disagreeing variables, and the effect of the SBND crosser and end-point cuts. State any residual uncertainty or follow-up (including the larger azimuthal scan and the APA gap).

### Ask Nupur to update the labels on the N-1 studies
Opened 2026-10-07.

Request Nupur to update the labels on the N-1 studies so that they match the labels used in the other plots (cut names, category and detector labels, and the
watermark and version stamps). I did not find "N-1" in the technote's `.tex` files, so the figures may be absent from the current text or live elsewhere; confirm where they
are and which labels differ before sending the request.

### Review the remaining root ROOT outputs
Opened 2026-10-13. Updated 2026-10-17.

Done 2026-10-17: removed `test_detsysonly.{ini,xml}` and `runs/xml_backup_20261008.tgz`; archived the one-time migration scripts `build_cov_v2.py` / `build_da_v2.py` (`archive/cleanup_2026-10/`); the tier 1 and 2 caches of `docs/root_bin_candidates.md` were already deleted. Kept: tier 3 (active XMLs, minutes to regenerate), tier 4 (read by the conditional-constraint notebook) and tier 5 (frozen SBND data-access caches, hours to regenerate); delete tier 3 or 5 only if the space is needed.

### Verify the remaining plotting/studies scripts and retire their notebooks
Opened 2026-10-13. Updated 2026-10-17.

Verified 2026-10-17 by running the notebooks headless on scratch copies and comparing pixel by pixel and number by number (see `plotting/studies/README.md`): `exposure`, `signal_selection`, `systematics`, `surface`, `resolution`, `data_access*` (80 figures, all phases) and `dent_study` (26 figures, 541 printed numbers) reproduce their notebooks. Those notebooks are retired, except `plotting/dent_study.ipynb`, which `build_report.py` regenerates together with `docs/dent_study/report.html`. Differences that are deliberate: `resolution.py` saver settings; phase-4 `constrained_*` at 200 dpi (identical at the notebook's 100 dpi). The comparison tool was a scratch script (not in the repo).

Still to check, with the notebooks kept in `plotting/` until then (deferred 2026-10-17: needs a PROfit rerun and a decision on which build is v3.0.4; for now both scripts count as linted only):
- `conditional_constraint.py`: the notebook's input `runs/legacy_root/contours_recoe_covar_PROplot.root` does not exist. It can be regenerated from the kept caches `contours_recoe_{prop,syst}.bin` with `plot --with-covar` (`jobs/contours/contours_recoe.ini`); needs a PROfit v3.0.4 build.
- `sbnd_scan.py`: the notebook fails 3 cells because `sbnd_scanmatch_v1_PROplot.root` is missing, and there are no `sbnd_scanmatch` caches, so regenerating it needs `process` then `plot` (`jobs/contours/sbnd_scanmatch.ini`). The spectra, ratio and deficit figures were not compared (the notebook only displays them).
- Known leftovers: `plotting/dent_study/replica.py`/`dentstudy.py` paths fixed; output destinations changed (scripts write to `--output`, not the technote or unblinding-request directories); `systematics.py` replaced the `SYSTS_*` environment variables with options.

### Confirm the proposed manifest families
Opened 2026-10-17.

`xml/MANIFEST.toml` classifies all 837 XML and ini files (0 unclassified, no `needs-review` left after the fake-data and DENT-improvement XMLs were archived), but 10 families are still marked "(proposed)" in `reviewed` and need the owner's confirmation. List them with `python scripts/figure_provenance.py manifest <dirs>` or `grep proposed xml/MANIFEST.toml`.

## Done

* 2026-10-17: `scripts/detvar_smoothing/` (closure, held-out, smoother bias/variance, split closure) archived to `archive/cleanup_2026-10/scripts_detvar_smoothing/` together with the DENT-improvement XMLs and scripts it depended on; the smoothing reports in `docs/detvar_smoothing/` and `plotting/studies/detvar_smoothing_figures.py` stay. Restore all three archives to rerun the study.
* 2026-10-07: Binding Energy histograms regenerated from the v1.1.0 `_cvrw` files (details in the technote item above).
* 2026-10-07: Resolution XMLs migrated and the 12 resolution figures regenerated (overnight log).
* 2026-10-13: Root-level ROOT outputs moved to `runs/legacy_root/`; classes C, D, F outputs deleted (`docs/root_output_audit.md`).
* 2026-10-14: Disallowed detector-data-access results (muon length and direction) removed (`docs/dent_study/REMOVED_DISALLOWED.md`).
* 2026-10-17: Fake-data work (except `fdgen/`) archived to `archive/fakedata_2026-10/`; `xml/dent_improvement` and `tmp_debug` archived; manifest and README updated.
* 2026-10-17: Notebooks `summarize_exposure`, `signal`, `resolution`, `systematics`, `surface` and `notebook` verified against their `plotting/studies` scripts and retired.

### Follow-up on why the charge-scale / resolution systematic jumps at high energies (closed)
Opened 2026-10-07.

**Observation.** In the SBND multi-proton reconstructed-energy uncorrelated-systematics figure
(`docs/systematics_v1.1.0/systematics_sbnd_run1_uncorrelated_recoe_multip.pdf`), the "Charge Scale / Resolution" curve rises
smoothly to about 0.05 at 1.4-1.5 GeV, then jumps to about 0.14 in the last bins (1.5-1.7 and 1.7-2.0 GeV). That step
dominates the total uncertainty (0.148) and is what set the new `YLIM['uncorrelated']` maximum.

**Why it matters.** A step this sharp, aligned with the bin edge at 1.5 GeV, looks more like a binning or statistics
artefact than a physical energy dependence.

**Things to check.**
- MC statistics in the 1.5-1.7 and 1.7-2.0 GeV bins of the detvar samples (CV and varied) behind the ratio histograms.
- Whether the bins were merged at 1.5 GeV in the `contours_recoe.xml` / suite binning, and whether the detsys ratio
  histograms use the same merging.
- Whether the effect is also present in the single-proton channel and in true-energy binning.
- Whether ICARUS shows the same feature.
- Sensitivity to the cathode-crosser and high-y/z cuts, since the detvar ratios were rebuilt with them
  (`detvar/detsys-v1.1.0v1.root`).
- The origin of the variation (an overall energy-scale shift versus a resolution smearing migrating events across the
  last bin edge).

**Findings (2026-10-07).** SBND uses `ChargeScale15` (`sbnd_detvar_smear15` against `sbnd_large`). The ratio histogram has one 1.5-2.0 GeV bin
(0.867), so the two analysis bins above 1.5 GeV (1.5-1.7 and 1.7-2.0) get the same value. The ratio was reproduced from the ntuples and the
drop is a selection effect, not an energy shift: of the 624 CV events in that bin, 102 are selected only in CV and 29 only in the varied sample
(net -11.7%), while energy migration of events selected in both samples is only about -1.3%. The fraction of CV-selected events missing from the
varied selected tree rises with energy (about 7% below 1.3 GeV, 10% at 1.3-1.5, 17% at 1.5-2.0, against about 5% in the reverse direction), and
the lost events have more protons (25% with five or more, against 15% for the kept events). The cut is the pion veto (`reco_pion_multiplicity == 0`): 12 of the 21 lost events that can be traced
through the varied sample's signal tree fail it, and in truth terms about 30-38% of selected events with a true pion leave the selection at
every energy, against 4-7% of events without one. The share of true-pion events in the selected sample rises steeply with energy (about 5% at
0.7-1.0 GeV to 32% of non-crossers at 1.5-2.0 GeV), and the crosser cut removes crossers, which are less pion-rich (20%). Still open: whether
to split or smooth the high-energy ratio bin, and why the varied reconstruction identifies pions the CV missed.

**Resolution (2026-10-07, `plotting/dent_study/chargescale_highE.py`).** With paired-bootstrap errors on the current selection (crosser and 4 cm end-point cuts), splitting the 1.5-2.0 GeV
bin gives 0.870 +- 0.027 (1.5-1.7 GeV) and 0.856 +- 0.068 (1.7-2.0 GeV), consistent with the single-bin value 0.868 +- 0.022, so a split adds no information. The step from the 1.3-1.5 GeV bin
(0.953 +- 0.009) to 1.5-2.0 GeV is significant (about 3.5 sigma) and lies about 3 sigma below a linear trend fitted above 0.7 GeV, so a smooth fit would understate it. Recommendation: keep the
single bin, document the pion-veto origin in the note, and do not smooth.

**Decision (2026-10-13).** Keep the single, unsmoothed 1.5-2.0 GeV bin (smoothing is not obviously more physical). Still to do in the technote: a sentence on the pion-veto origin of the step.

**Where to look.** Detvar calculator `detsys/calculator/configuration.h`; ratio histograms in `detvar/`; frac_syst trees
in `runs/systematics/`; plotting in `plotting/systematics.ipynb`.

## Overnight log (2026-10-07)
* Azimuth scan redone (`plotting/dent_study/azimuth_scan.py`, `azimuth_systematics.py`): short labels, proton scores vs proton azimuth, only active SBND variations (charge scale, WireMod, non-uniformity, DENT). Muon PID azimuth pattern follows the charge-scale shape; proton PID dips near +-90 deg exceed every active variation.
* theta_xW scan (`plotting/dent_study/thetaxw_scan.py`, `docs/azimuth_scan/thetaxw_sbnd.png`): data-MC mean score vs theta_xW per wire plane (collection vertical, induction +-60 deg). WireMod X-theta_xW is nearly flat in theta_xW and does not improve the shape agreement (offset-removed chi2 mostly unchanged; worse for collection muon scores, where the detvar shift has spikes at +-35 deg that look like sample-statistics or spline artefacts). The charge-scale variation reproduces the data shapes (collection proton PID dip near theta=0). To do: ICARUS (collection horizontal, induction 1 horizontal, induction 2 +-30 deg from the vertical, alternating by TPC; use TPC-dependent angles), check the WireMod detvar sample statistics at +-35 deg.
* Stretch goals queued, in order: (1) per-systematic magnitude deep dive; (2) DENT narrative audit; (3) events+selected intermediate table; (4) ternary PID plots (`spineosc_pid`), 2x muon boost; (5) APA gap figures.
* APA gap figures rendered (`docs/detector_access_v1.1.0/dataaccess_sbnd_run1_{vertexz,muon_endz}_apagap_incl.pdf`, 4 cm bins over z = 210-290 cm): shape chi2/ndf 24.5/19 (vertex z), 31.4/19 (muon end z); only a single low bin at 250-254 cm in vertex z (ratio ~0.8), no systematic hole. NOT yet in the technote: gap z position (assumed ~250 cm) needs confirmation.
* Ternary PID plots (`plotting/dent_study/ternary_pid.py`, `docs/ternary_pid/`): the selected muon candidates all lie in the 2x-boost muon region (boost is in the selection). Share of candidates in the boost band (2x wins, argmax does not): nominal 4.4% (pair pool) -> charge scale 2.2%, WireMod 3.8->4.0%, non-uniformity 4.2->5.1%, DENT 3.8->3.6%. Signal share in the band 73% vs 89% outside (signal = Target subchannel). Caveat: displacement fields are conditioned on being selected in both samples; regression toward the muon vertex is expected in any re-simulation.
* DENT audit written: `docs/dent_study/AUDIT.md` (pile-up robust; crosser fraction depends on cvwgt, 29.0% vs data 27.5%; text caveat proposed, technote not yet changed). APA gap figures and text added to `detector_data_access.tex` (gap at z = 250 cm confirmed). Azimuth study excluded from the technote pending more investigation.
* Per-systematic deep dive, first pass: `docs/syst_deepdive/README.md` (large symmetric churn at the selection boundary; net ratios are small imbalances; need a same-seed nominal rerun or the idea-1 table to separate noise from physics).
* XML manifest drafted: xml/MANIFEST.toml + scripts/xml_manifest.py (--check/--list); all 1138 xml/ini files classified (170 active, 314 frozen, 475 archived, 175 needs-review, 3 deprecated, 1 reference); entries marked 'proposed' until confirmed.
* XML consolidation (2026-10-08, backup `runs/xml_backup_20261008.tgz`): ICARUS detector-access XMLs now live only in `xml/detector_data_access/icarus_run2` (ICARUS has no DENT); SBND detector access v2 replaced `xml/detector_data_access/sbnd_run1`; coverage v2 replaced `xml/coverage`, 19 un-migrated old coverage files moved to `xml/archive/coverage`; `sbnd_scanmatch.xml` reverted to HEAD and `contours_pn`, `contours_sbnd2d`, `sbnd_scanmatch` marked frozen. `build_cov_v2.py` / `build_da_v2.py` are one-time migrations and refuse to rerun without --force. Manifest: `xml/MANIFEST.toml`, lint `scripts/xml_manifest.py --check`.
* Resolution XMLs migrated (`xml/resolution/*` via `fix_selected`: v1.1.0 trees, cvrw + cvwgt, veto, SBND crosser and end-point cuts); PROfit plot rerun in `runs/resolution`, notebook stamp updated to selection v1.1.0.
* Resolution figures regenerated (12 files in technote figures/signal_and_selection, selection stamp v1.1.0); text quotes no numbers so no text change needed. Stale-figure inventory: resolution group now done.
* Deleted root-level ROOT outputs of classes C and D plus their .bin caches (128 + 99 files, 34 GB); see docs/root_output_audit.md.
* Deleted class F root outputs (63 files) and tier 1+2 .bin caches (9 files, 3.9 GB). Remaining root outputs: classes A and B plus .bin tiers 3-5 (see docs/root_bin_candidates.md).
* Root outputs moved to runs/legacy_root/ (590 files, 31 GB, plus parameter_scan/); notebooks signal, resolution, surface, conditional_constraint, sbnd_data_access, data_access_results now read from there. Pre-existing path mismatches found (not caused by the move): conditional_constraint.ipynb reads contours_recoe_covar_PROplot.root and sbnd_data_access.ipynb reads sbnd_scanmatch_v1_PROplot.root, neither exists (nearest: contours_recoe_v1_PROplot.root, sbnd_scanmatch_v1_PROplot_histograms.root); data_access_results.ipynb reads dp_test09142026_v1_PROplot.root, deleted with class F (frozen, earlier phase).
* Overnight detector-systematics study delivered: docs/syst_study/REPORT.md and report/index.html (see summary there).
