# Smoothing the detector-systematic ratios: closure first, then smoothing

Status: plan written 2026-10-14, before any run. Simulation only (no data are read). All fits float **all** systematics (analyst's decision).
Record every later change under "Deviations".

## Question
The detector systematics enter as per-bin ratios of a varied sample to a nominal one. Bin-to-bin fluctuations in those ratios are partly MC-statistical
noise (paired bootstrap, `docs/detvar_stat/`: SCE2x and the DENT ratio after the cuts are consistent with noise; ChargeScale at 1.5-2 GeV is real structure). Can a
smoothing method remove the noise without erasing real structure, and how well does the smoothed ratio reproduce the variation?

## Two closure tests, global fit with every systematic floating
1. **Injection closure (machinery).** Fake data = the CV prediction with one systematic injected at +1 sigma (`--inject-systs <name> 1`), global fit.
   Expect best-fit pull = 1 for the injected parameter and chi2 near zero (about 0 for an Asimov-like fake dataset; compare with the chi2 of the same fit with no injection).
   This checks the ratio histogram, the spline interpolation, the covariance and the fit, not the ratio's fidelity.
2. **Sample closure (representation).** Fake data = the variation sample itself, the MC = the nominal sample with the ratio systematic. The chi2 and the
   best-fit pull of the corresponding parameter measure how well the ratio represents the variation. Caveats built into the design:
   * The ratio is made from the same events as the injected sample, so closure is trivially perfect in the binning the ratio was made in. The test is informative
     for what is **not** identical: independent events (half-split below), a different observable or binning, and the production reweighting.
   * The nominal MC must be the **unweighted** nominal (`wsyst` file, no `cvwgt`) so that the only difference from the variation sample is the detector effect.
   * The expected chi2 is not 0: the injected sample carries its own MC-statistical fluctuation, so the reference distribution is built from the paired bootstrap.

## Half-split (cross-validation) so noise and bias can be separated
Split the variation sample's events at random into halves A and B. Derive the ratio from half A, inject half B (and the reverse). An unsmoothed ratio carries the noise
of A into the prediction; a smoothed ratio carries less noise but may carry bias. The chi2 and pull on the held-out half is the figure of merit. Repeat over many splits.

## Smoothers and metrics (phase 2, after the baselines)
Candidates: low-order polynomial or Legendre fit, cubic spline or LOESS, Gaussian kernel, all weighted by the bootstrap sigma; optionally on log(ratio); with and without
a normalization-preserving constraint; the bandwidth or order chosen by cross-validation, not by eye.
Metrics: bias^2 and variance per bin from smoothing each bootstrap replica and comparing with a high-statistics reference (the full sample pooled); held-out closure chi2 and
pull; preservation of real structure (the ChargeScale 1.5-2 GeV bin, z = -5.2); behaviour on noise-dominated ratios (SCE2x, DENT after the cuts).

## Test systematics and observables
ChargeScale15 (real structure), SCE2x (noise dominated), DENT at S2 (small, mostly noise), and WireModxThetaXW (large, smooth). Observable: start with the SBND-only Tier B
configurations already built (`xml/dent_improvement/sbnd_run1/`, the PID scores and calorimetric differences, S0/S2, full systematic suite, per-observable ratio files);
then the reco-energy ratios of the analysis (`detsys-v1.1.0v2.root`) in a closure configuration.

## Gates
G1: test 1 passes on the pilot (pull 1, chi2 near 0) for each test systematic. G2: the half-split baseline (unsmoothed) is established with its spread. Only then
smoothers are compared.

## Deviations
* 2026-10-14 (before the closure runs): SCE is not a floated systematic in the SBND analysis suite (the SBND detector systematics are DENT, ChargeScale15, WireModxThetaXW,
  MattModYZ, ProtonMiss, ToroidCalibration and TriggerEmulation), so SCE2x cannot be a closure target. The noise-dominated test case becomes DENT at S2, with MattModYZ as a
  second small-shift case. Pilot of test 1 (muon score, S2, ChargeScale15 injected at +1 sigma, all systematics floating, absolute normalization): global best-fit chi2 0.9977,
  injected parameter fitted at 0.9977, all other nuisances within 0.05 of 0; 17 seconds per fit. A total chi2 of about pull^2 is expected, because PROfit's chi2 includes the
  nuisance penalty; the data-only part is close to 0.
* 2026-10-14 G1 (closure test 1), muon_score_inclusive_s2, +1 sigma injected, all systematics floating, absolute normalization (`scripts/detvar_smoothing/closure.py test1`,
  `test1_muon_score_inclusive_s2.csv`): fitted injected parameter 0.959 (DENT), 0.998 (ChargeScale15), 0.847 (WireModxThetaXW), 0.978 (MattModYZ); total chi2 equal to the pull
  squared within 0.01 apart from WireModxThetaXW (0.847 at a pull of 0.847, so a data-only chi2 of 0.13); largest other pull 0.015-0.055; no-injection baseline chi2 = 0.
  The WireModxThetaXW shrinkage is the expected result of a partial degeneracy with other floating parameters (MaNCRES) and does not indicate a representation problem. G1 passes
  for this observable; it should be repeated for the other studies before it is relied on.

## Synthetic-truth bias and variance of the smoothers (simulation only; `scripts/detvar_smoothing/smoother_bias_variance.py`, `bias_variance.csv`; 2026-10-14)

Design: six truth ratio shapes (null, linear, quadratic, rising, wave, and a linear shape with a -10% step in the last bin, "edge") in the 17 analysis reco-energy bins; Gaussian
noise per bin at the paired-bootstrap sigma of each of four systematics (DENT after the cuts, ChargeScale15, WireModxThetaXW, MattModYZ) scaled to a one-third sample (the part size of the
sample closure); 200 noisy draws per scenario; every smoother applied to each draw. Relative MSE = (bias^2 + variance summed over bins) / (sum of sigma^2, i.e. the unsmoothed error).
Below 1 is a gain. (A first version of the cross-validated local-linear smoother always chose the narrowest bandwidth because its leave-one-out error was trivially zero; that was a bug,
fixed before the numbers below.)

**Mean over the six truth shapes**

| method | ChargeScale15 | DENT | MattModYZ | WireModxThetaXW |
|---|---|---|---|---|
| constant | 4.13 | 0.94 | 7.89 | 2.21 |
| gauss1.5 | 0.81 | 0.27 | 1.46 | 0.48 |
| gauss3 | 1.58 | 0.38 | 3.0 | 0.84 |
| gauss_cv | 0.69 | 0.28 | 1.18 | 0.46 |
| local1.5 | 0.63 | 0.4 | 0.91 | 0.5 |
| local3 | 0.6 | 0.46 | 0.77 | 0.51 |
| local_cv | 0.65 | 0.38 | 1.0 | 0.5 |
| none | 1.05 | 1.03 | 1.03 | 0.97 |
| poly1 | 2.04 | 0.51 | 3.95 | 1.12 |
| poly2 | 1.29 | 0.39 | 2.48 | 0.76 |
| poly3 | 0.67 | 0.3 | 1.13 | 0.46 |
| poly_cv | 0.74 | 0.28 | 1.29 | 0.46 |
| spline_gcv | 0.72 | 0.58 | 0.76 | 0.62 |

**Worst case over the six truth shapes**

| method | ChargeScale15 | DENT | MattModYZ | WireModxThetaXW |
|---|---|---|---|---|
| constant | 13.3 | 3.02 | 25.61 | 7.15 |
| gauss1.5 | 3.02 | 0.76 | 5.79 | 1.73 |
| gauss3 | 3.22 | 0.75 | 6.16 | 1.77 |
| gauss_cv | 2.9 | 0.77 | 5.42 | 1.71 |
| local1.5 | 1.95 | 0.68 | 3.5 | 1.23 |
| local3 | 1.35 | 0.61 | 2.34 | 0.91 |
| local_cv | 2.25 | 0.71 | 4.26 | 1.42 |
| none | 1.15 | 1.11 | 1.12 | 1.07 |
| poly1 | 5.25 | 1.3 | 10.59 | 2.84 |
| poly2 | 3.76 | 0.96 | 7.59 | 2.04 |
| poly3 | 2.98 | 0.8 | 5.74 | 1.77 |
| poly_cv | 3.6 | 0.87 | 6.95 | 1.98 |
| spline_gcv | 1.03 | 0.71 | 1.28 | 0.84 |

**Bias in the last bin for the edge truth, in % (true value -10%; 0 means the step is preserved)**

| method | ChargeScale15 | DENT | MattModYZ | WireModxThetaXW |
|---|---|---|---|---|
| constant | 7.3 | 7.3 | 7.3 | 7.3 |
| gauss1.5 | 8.5 | 8.5 | 8.7 | 8.7 |
| gauss3 | 8.7 | 8.8 | 8.8 | 8.8 |
| gauss_cv | 8.3 | 8.4 | 8.3 | 8.5 |
| local1.5 | 6.4 | 6.2 | 6.4 | 6.5 |
| local3 | 4.5 | 4.2 | 4.8 | 4.7 |
| local_cv | 6.9 | 6.8 | 7.2 | 7.2 |
| none | -0.1 | -0.4 | 0.1 | 0.6 |
| poly1 | 9.8 | 9.9 | 9.8 | 9.9 |
| poly2 | 9.3 | 9.4 | 9.5 | 9.5 |
| poly3 | 8.3 | 8.2 | 8.5 | 8.5 |
| poly_cv | 9.4 | 9.2 | 9.5 | 9.3 |
| spline_gcv | 2.0 | 2.1 | 2.0 | 2.6 |

Reading (the conclusions depend on this truth set and on the noise profiles, which are real):
* **The noise-dominated systematic (DENT after the cuts) gains from almost any smoother** (relative MSE 0.27-0.58, worst case below 1.1 for all but the low-order polynomials).
* **Systematics with real structure gain little or lose.** ChargeScale15, MattModYZ and WireModxThetaXW have small per-bin noise, so a smoother mostly adds bias; every smoother except the
  smoothing spline and the widest local-linear has a worst case above 1.2 (up to 25x for a constant fit).
* **The edge case separates the methods.** Polynomials, Gaussian kernels and the cross-validated variants erase 7-10% of a 10% step in the last bin. The smoothing spline with the
  GCV penalty keeps it (2% residual bias), the wide local-linear fit is intermediate (4.5%).
* **The smoothing spline with a GCV penalty is the only smoother whose worst case stays near 1 for all four systematics** (0.71-1.28) while gaining on average (0.58-0.76). It is the leading
  candidate for the PROfit closure comparison. Polynomials of order 3 and Gaussian kernels gain more on the noisy systematics but destroy real structure.

## G2: sample closure baseline, unsmoothed ratios (PROfit, all systematics floating; `split_closure.py --all`, `closure2_muon_score_inclusive_s2.csv`; 2026-10-14)

Muon score, S2, three-way event split (ratio part / MC part / pseudo-data part), absolute normalization at 1e20 POT, three role permutations. Data = the variation sample's
events of the data part; null = the reference sample's events of the same part (no variation). Pull = the fitted value of the parameter under test (1 would mean the ratio represents
the variation at one sigma); excess chi2 = total chi2 with the variation data minus the null.

| Systematic | Fitted pull, varied data (3 splits) | Fitted pull, null | Excess chi2 (mean, sd) |
|---|---|---|---|
| ChargeScale15 | 0.89, 0.80, 0.89 (mean 0.86, sd 0.05) | 0.00, 0.06, 0.00 | 5.7, 7.2 (14.0, 0.9, 2.3) |
| MattModYZ | 1.00, 1.00, 0.74 (mean 0.91, sd 0.15) | 0.17, 0.40, 0.06 | 8.2, 1.7 |
| DENT | 0.08, 0.00, 0.05 (mean 0.04, sd 0.04) | 0.05, 0.00, 0.00 | -0.3, 4.0 |
| WireModxThetaXW | 0.37, 0.00, 0.00 (mean 0.12, sd 0.21) | 0.26, 0.01, 0.00 | 0.1, 1.2 |

Reading: ChargeScale15 and MattModYZ are represented at about 85-90% of full strength. The DENT effect at S2 on the muon score is not distinguishable from the null (as the bootstrap
predicted), and WireModxThetaXW is not pulled because its muon-score effect is degenerate with MattModYZ ("Non-uniformity SBND" takes the pull instead, 0.36-0.49). The excess chi2 of
ChargeScale15 varies from 0.9 to 14 across the three splits, a direct measure of the noise of a one-third-sample ratio. Three splits are too few for a strong comparison; the smoothed variants
below use the same splits so that the comparison is paired.

## G3: smoothed ratios in the PROfit sample closure (`split_closure.py --multi`, `closure_summary.py`, `closure2m_muon_score_inclusive_s2.csv`, `closure_by_split.csv`; 2026-10-14)

Same design as G2 (muon score, S2, three splits, all systematics floating, 1e20 POT), but one processed MC per split serves all four systematics (all four ratios replaced by their
part-made versions at once; the unsmoothed multi baseline reproduces the per-systematic G2 baseline: ChargeScale15 0.86 +- 0.03 and excess chi2 5.2, MattModYZ 0.91 +- 0.16 and 8.2).
Smoothers act on the ratio bins inside the fit window (0.3-1.0) with weights 1 / sigma^2 from the part sample; paired changes are relative to the unsmoothed ratios of the same split (mean, sd over 3 splits).

| Systematic | Variant | Pull, varied data | Excess chi2 | Paired change in excess chi2 |
|---|---|---|---|---|
| ChargeScale15 | unsmoothed | 0.86 | 5.2 | |
| | smoothing spline (GCV) | 0.91 | 2.8 | -2.4 +- 2.6 |
| | local-linear, bandwidth 3 | 0.93 | 6.5 | +1.3 +- 2.7 |
| | polynomial, CV order | 0.89 | 13.2 | +8.0 +- 1.4 |
| MattModYZ | unsmoothed | 0.91 | 8.2 | |
| | smoothing spline | 1.00 | 4.0 | -4.3 +- 1.2 |
| | local-linear | 1.00 | 3.9 | -4.3 +- 2.2 |
| | polynomial, CV order | 1.00 | 4.9 | -3.4 +- 4.7 |
| DENT | unsmoothed | 0.06 | -0.8 | |
| | smoothing spline | 0.10 | 0.0 | +0.8 +- 0.1 |
| | local-linear | 0.09 (null pull 0.16) | 0.4 | +1.2 +- 0.1 |
| | polynomial, CV order | 0.07 | 1.0 | +1.8 +- 0.4 |
| WireModxThetaXW | all variants | 0.00-0.12 | -0.1 to 0.5 | -0.5 to -1.1 |

The null chi2 is unchanged by smoothing (-0.15 to -0.87, within noise).

Reading:
* **The smoothing spline with a GCV penalty is the best and safest of the three.** It moves the pulls of the real effects closer to 1 and lowers their excess chi2 (ChargeScale15 5.2 -> 2.8 and its spread across splits 7.1 -> 4.6; MattModYZ 8.2 -> 4.0), and does not harm the noise-dominated cases.
* **A polynomial with a cross-validated order damages ChargeScale15** (excess chi2 5.2 -> 13.2, paired +8.0 +- 1.4): it removes real structure, as the synthetic bias-variance study predicted for polynomials and kernels.
  The local-linear fit is intermediate. The PROfit closure therefore agrees with the synthetic-truth study on the ranking of the methods.
* DENT after the cuts (and WireModxThetaXW) show no benefit from any smoother here: for DENT the excess chi2 rises by 0.8 to 1.8, small compared with its split-to-split sd of 5, and the fitted pull does not change. This observable is not
  a sensitive test of noise reduction for those two.

Limits: three splits only (the paired changes have sds of 1-5 chi2 units); one observable (the muon score at S2) whose real structure is smooth; the reco-energy ratios the analysis uses, and their sharp last-bin feature, are not tested in PROfit
(only in the synthetic study); the ratio is smoothed in the observable it was made in, not across channels or subchannels.


## Update 2026-10-15: rolling windows, reco energy, and a POT-normalization bug

**Bug (found, fixed, reruns in progress).** `split_closure.ratio_from_part` divided raw event counts of the varied and
nominal samples without scaling by POT. The samples differ: nominal 5.59e20, ChargeScale15 5.57e20 (-0.5%), MattModYZ
5.39e20 (-3.7%), WireModxThetaXW 4.90e20 (-12.4%); DENT is restricted to common events and was unaffected. The part-made
ratios of MattModYZ and WireModxThetaXW therefore sat 4% / 12% low in every bin, and every PROfit closure result made before
this fix for those two (and the "WireMod degenerate with MattModYZ" remark) is unreliable. The fixed ratios agree with the
production ratios (mean over bins 2-9: 0.992/0.999/0.997 for ChargeScale15/WireModxThetaXW/MattModYZ). Old CSVs are kept as
`*_prePOTfix.csv`. The synthetic bias-variance study and the held-out check below do not depend on the normalization
(held-out numbers are unchanged to two decimals after the fix).

**Rolling window** (`roll3`, `roll5`: inverse-variance-weighted mean over 3 or 5 bins, truncated at the edges; `median3`).
Synthetic truth, relative MSE (1 = unsmoothed; mean over smooth truths / edge truth): roll3 0.2-0.45 / 2.9-5.5, spline 0.55-0.7 / 0.95-1.2.
Last-bin bias for a -10% step: roll3 8.2%, roll5 8.8%, local3 4.4%, median3 4.7%, spline_gcv 1.7-2.6%.

**Held-out predictive error** (`heldout.py`; smoothed ratio of one third vs raw ratio of another, paired over six ordered pairs,
chi2 per bin, change vs unsmoothed): reco energy: roll3 -0.11, poly1 -0.12, gauss1.5 -0.11, spline -0.07, constant +0.11;
muon score: roll3 -0.38 +- 0.15, median3 -0.42, local1.5 -0.42, spline -0.39 +- 0.09, roll5 +0.74, poly2/3 0, gauss3 +4.5.
The score is dominated by well-measured bins, so it says little about the last bin.

**Reco-energy PROfit closure (pre-fix ratios, 3 splits):** not discriminating: excess chi2 within +-0.4 for all variants
(DENT +3.5 +- 4 for every variant, including unsmoothed), pulls 0-0.35; the energy shape is insensitive to these variations
and other nuisances absorb them. Rerun after the fix is in progress (energy, then muon score with roll3/roll5).
