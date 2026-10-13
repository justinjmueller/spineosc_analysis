# Does the SBND DENT simulation improve data/MC agreement?

*A shape-only, no-fit study with the SBND release detector-access data. Generated 2026-10-07 from `plotting/dent_study/report_cells.py`.*

## Summary

**Question.** DENT is the detector variation that models the cathode-crossing (and related) effects in SBND. After the
cathode-crosser cut it is one of the larger uncorrelated systematics. Does simulating it (the DENT detvar sample, or the
production CV shifted by the DENT ratio) describe the data better than the CV does, and where does it match or fail to match?

**Method in one paragraph.** SBND release detector-access data (the file has no energy or momentum
column, and none is used) are compared with the SBND simulation, always area-normalised and shape-only, in each observable's own
binning, with shape chi-square and p-value as the metrics. Nothing is fitted: DENT enters only at fixed strengths 0, 0.5 and 1.0.
Tier A uses the detvar CV and DENT samples directly (MC statistics only). Tier B uses the production CV (`cvrw` plus `cvwgt`)
with the full systematic covariance. Tier C asks how well one observable could tell the two apart, using fake data.

**Findings.**

1. **The data do not show that DENT applied to the central value improves the agreement, and they do not reject it.** With the
   production covariance, shifting the CV by the full DENT effect changes the shape chi-square by less than 2 units for 9 of
   the 10 primary observables (5 improve, 5 worsen). The tenth, `muon_calo_csda`, is dominated by a covariance re-centering artefact
   (section 7). Across the 111 SBND regenerated detector-access channels the largest changes are +5.8 and -3.1 units.
2. **Tier A (MC statistics only) splits by physics.** The DENT sample is closer to the data for the seven calorimetry and
   particle-ID observables (largest for `muon_calo_csda`, 10,534 to 2,188 on 19 degrees of freedom), and farther for the three
   cathode-geometry observables (`crosser_flag`, `muon_abs_endx`, `muon_endx_full`; for `muon_abs_endx` 51 to 338 on 19 degrees of freedom).
   The calorimetry/PID improvement disappears as the other systematics are added to the covariance (section 6), so it is not
   evidence that DENT is the cause of those data/MC differences.
3. **The cathode features that make DENT large are only partly confirmed by the data.** The fraction of events flagged as crossers
   is 27.46 +- 0.18 % in data against 27.92 % (CV) and 25.36 % (DENT). Non-crossing muons ending within 5 cm of the cathode
   (the "pile-up") are 1.95 % of events in data against 2.07 % (CV) and 3.78 % (DENT): the data follow the CV, and the large DENT
   pile-up is not present. The deficit of crossing muons that end within 4 cm of the cathode is, however, in the direction DENT
   predicts (data 0.68 %, CV 1.14 %, DENT 0.55 %). These cathode comparisons were seen before the study plan was written, so they are
   hypothesis-prompted, not independent confirmations.
4. **DENT as an uncertainty is not redundant.** Removing DENT from the covariance makes the muon end-x observables fail
   (`muon_abs_endx`: p = 0.46 with DENT, 0.003 without; `muon_endx_full`: 0.15 against 2e-5). It absorbs a real data/MC difference.
5. **No single observable has the power to decide.** With pseudo-experiments through the fixed covariance, the separation between
   the "CV" and "full DENT" expectations is d' = 0.5 to 1.0 per observable, and all ten primaries have two-sided p >= 0.15 under all
   three truths (CV, half DENT, full DENT). The observables are correlated, so their p-values cannot be multiplied together.

**Recommendation.** Keep DENT as an uncertainty, and do not shift the central value by it on the strength of this study. The
mechanism that most plausibly makes DENT large at high energy is the near-cathode endpoint pile-up, which the data lack. A
near-cathode endpoint treatment (an endpoint cut, or a dedicated near-cathode observable) is the natural follow-up and is not part of this study.

**What this study cannot say.** It uses no energy or momentum quantity, so it cannot test the effect of DENT on the
reconstructed-energy distributions that drive the analysis. The detvar samples are neutrino-only and their CV is not the
production nominal. ICARUS has no DENT variation; its regenerated studies are reported only as a baseline.
