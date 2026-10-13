# DENT data/MC study: analysis plan (written before the new results were produced)

Status: written 2026-10-07 (early hours), before any Tier A/B/C result for the wide-net observables existed. Changes after this
point are recorded in the "Deviations" section at the bottom, never silently.

## 1. Questions
1. Does the SBND DENT detector-variation simulation improve the **shape** agreement between simulation and the SBND
   detector-access data, observable by observable?
2. Where does it match the data and where does it not? In particular, are the features that make the DENT variation large
   (crossing-track migration at the cathode, the near-cathode endpoint pile-up, the muon calorimetric shift) present in data?
3. Is the DENT systematic as assessed (a one-sided per-bin ratio from the detvar CV/DENT samples, applied to the production CV)
   consistent with using the DENT detvar sample directly?

## 2. Data, samples, selection
* Data: `spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root`, tree `events/sbnd_bnblight_full/selected`
  (event count not quoted). This file has no energy or momentum branch. **No energy or momentum quantity is used anywhere.**
* Production MC (Tier B): `spineosc_numu_sbnd_v1.1.0_wsyst_cvrw.root` (CV-reweighted, `cvwgt` multiplied into every beam-MC
  branch), with the full interaction-model, flux, geant4 and detector systematic suite of the existing data-access XMLs.
* Detvar MC (Tier A): `sbnd_detvar_newcv` and `sbnd_detvar_newdent` from `spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root`
  (no weights, so no interaction-model systematics and no `cvwgt`; MC statistics only).
* Selection (identical for MC and data except as noted): `reco_proton_multiplicity >= 1`, `reco_fiducialize_cathode == 1`,
  `reco_veto_sbnd_highy_highz == 1`. The MC additionally applies the flash-PE threshold (`0.647 * event_largest_flash_pe > 2000`),
  which is **trigger emulation and must not be applied to data**.

## 3. Rules that are not negotiable
1. **No fitting.** DENT is injected only at fixed strengths 0 (CV), 0.5 and 1.0 (full DENT). No minimisation over its strength,
   no profile, no scan used to pick a best value. (Earlier exploratory interpolation estimates of a "DENT strength" were of
   exactly this kind; they are not used and not repeated.)
2. Every PROfit run uses `area-norm = true` **and** `shapeonly = true`. Absolute comparisons are not made.
3. Shape-only chi-square and its p-value (`scipy.stats.chi2.sf(chi2, ndf)`, `ndf = nbins - 1`) are the quantitative metrics.
4. Systematics are assessed in the observable's own binning. Slicing into channels is allowed.
5. Plots follow `plotting/data_access_results.ipynb` (mutools, `use_style('rootlike')`).

## 4. Observables
* Wide-net registry: `scripts/dent_study/observables.json` (24 energy-free observables present in the data tree), binned with
  equal-population quantiles of the **detvar CV MC only**.
* **Primary set** (pre-specified criterion, MC-only): observables with `sep/ndf >= 1`, where `sep` is the chi-square between
  the CV and DENT shape at data-sized statistics. Result (`scripts/dent_study/primary_selection.json`):
  `muon_calo_csda, crosser_flag, muon_proton_score, muon_abs_endx, muon_endx_full, muon_pion_score, proton_mip_score,
  proton_calo_csda, proton_muon_score, proton_pion_score`.
* Everything else (remaining registry observables and the regenerated detector-access studies) is **diagnostic**: reported,
  not used for any decision.
* Note: `muon_len3d` has `sep/ndf = 0.07`. The muon straight-line length is only weakly sensitive to DENT even though DENT
  truncates some tracks at the cathode; most flipped events lose a few cm.

## 5. Variants
| id | prediction | covariance | purpose |
|---|---|---|---|
| A_cv, A_dent | detvar CV / detvar DENT sample as the MC | MC stat only | direct comparison; no ratio assumption |
| B v1 | production CV (cvrw+cvwgt) | full, DENT included | the standard analysis |
| B exDENT | production CV | full, DENT removed | how much DENT as an uncertainty absorbs |
| B inj0p5, inj1p0 | CV shifted by 0.5 / 1.0 x DENT | full, DENT included | effect of the DENT central-value shift |
| C | Poisson fake data from CV or from CV+DENT, 50 seeds | as B | expected distributions |

`--exclude-systs` and `--inject-systs-cv` cannot be combined (the injection is silently ignored), so "shifted central value with
DENT removed from the covariance" is not available; it is not approximated.

## 6. Decision quantities
* `dchi2 = chi2(theta=0) - chi2(theta=1)` in the B variants with DENT in the covariance; positive favours DENT.
* Tier C gives the distribution of `dchi2` for fake data drawn from CV and from CV+DENT (same seed, paired). The data value is
  placed in both: `p_CV = P(dchi2 >= obs | truth CV)` and `p_DENT = P(dchi2 <= obs | truth DENT)`. Only those two numbers, and the
  plain `chi2`/`p` of each variant, are quoted as evidence. No best-fit strength is quoted.
* Multiplicity: the ten primary observables are tested; Holm-Bonferroni adjusted values are reported alongside the raw ones. The
  observables share events and are correlated, so the adjustment is conservative and the number of independent tests is smaller.
* Tier C with Poisson-only fake data does not throw the systematic nuisance parameters; its `dchi2` distributions are therefore
  narrower than reality. A subset also gets `--pseudo-experiment` throws if time permits; the difference is reported.

## 7. Quality controls (all recorded per study in `runs/dent_study/results/*.json`)
* `shapeonly` is shown to be live per study by re-running without it: chi2 must differ. (`ndf` cannot prove it.)
* Seed pinned (`seed = 20261007`); Tier C seeds 1..50.
* `cvwgt` present on every beam-MC branch; main-MC file is the `_cvrw` file.
* Data and MC selections identical (see section 2). The original `_data.xml` files omitted the high-y/z veto in 35 of 39 SBND
  studies; the regenerated ones add it.

## 8. What I had already seen (honesty section)
Before this plan, in exploratory work: (a) I compared data to the CV and DENT detvar samples for all 46 shared energy-free
variables, using shape chi-squares; (b) I identified the crossing-flag flips, the near-cathode pile-up and the muon calorimetric shift
in MC, and saw that the pile-up was absent in data; (c) those exploratory data comparisons applied the flash-PE cut to data too (0.8% of
events), which was a mistake. Consequences: the *selection* of the primary set uses only MC quantities, but the analyst was not blind
to the outcome for the cathode observables. The cathode-feature results are therefore **hypothesis-prompted, not independent
confirmations**, and the report labels them as such.

## Deviations (appended as they happen)
* 2026-10-07 03:35: removed `pion_mult`, `photon_mult`, `electron_mult` from the registry (24 -> 21 observables). The selection forces all three to exactly 0, so they carry no information and would be rejected by PROfit's empty-bin check. None was in the primary set.
* 2026-10-07 03:30: exploratory data dumps had wrongly applied the flash-PE cut to data; re-made without it. The headline cathode numbers changed by <0.1 sigma-equivalents (see key_numbers.csv); the corrected files are the ones used from here on.
* 2026-10-07 03:59 (AFTER seeing the first Tier C result, `crosser_flag`): the data's dchi2 fell between the truth-CV and truth-DENT distributions, outside both. Added a third fixed fake-data truth, DENT x 0.5, for all Tier C jobs started after this time, and switched the reported Tier C p-values to two-sided empirical ones (floor 2/(n+1)). This is a third pre-specified point, not a fit; no strength is estimated. `crosser_flag` will be re-run with the new truth.
* 2026-10-07 04:30: found that the two crossing-split studies (`muon_endx_zoom_crosser`, `dent_crosser/muon_score`) register DENT under per-population names (`_stay`, `_xing`, `_xing_plusx`, `_xing_minusx`), so the standard `DENT_sbnd_Run1` exclusion/injection was ignored by PROfit and the driver recorded the injected variants as 'ok' with chi2 identical to v1 (invalid). The driver now discovers the active DENT names per study, injects/excludes all of them, and marks a variant failed if PROfit logs 'Unrecognized spline' / 'matches no systematic' or if chi2 equals v1 exactly. Both studies were re-queued; the report prints an audit of any remaining inert injection.
* 2026-10-07 06:15: Tier C switched from Poisson-only throws to PROfit `--pseudo-experiment` (Poisson + Gaussian nuisance throws through the total covariance). A six-seed test showed the Poisson-only dchi2 spread (sd 0.2-0.3) is several times narrower than the pseudo-experiment spread (about +-1-3), so Poisson-only p-values cannot calibrate the data. Poisson-only results are retained as a supplementary column; the pseudo-experiment version is the one the decision rule uses. All ten primaries are re-run with it.
* 2026-10-07 12:45: found that `--inject-systs-cv` re-derives the systematic covariance around the SHIFTED prediction (fractional DENT covariance first bin 0.84 -> 0.51 for `muon_calo_csda`), so PROfit's chi2(theta=1) mixes the central shift with a covariance change. For `muon_calo_csda` PROfit reports dchi2 = -18.8 (the shift 'worsens' agreement) while the same covariance held at its theta=0 form gives +0.2, and the stat-only chi2 falls from 14,606 to 3,556. A numpy replica of PROfit's chi2 (Neyman diag(n) + F o p p^T with F = collapsed_total_frac_cov) reproduces PROfit to <0.1 in chi2 and is used to (a) report the fixed-covariance dchi2 and (b) run 20,000 pseudo-experiments per truth (Poisson + Gaussian systematic throws through F) for calibrated p-values. The PROfit-native dchi2 and its 50-throw Tier C are kept and shown alongside, clearly labelled; the replica-based quantities are the clearer measure of whether the DENT shift improves the agreement.
