# Q2: statistical content of the detector-systematic ratios

Question: "Given that only small stochastic effects remain (identical events in CV and varied samples), is there any evidence of statistical effects being assessed in the ratio/splines?"

Short answer: **Yes, and it is quantitatively important for about half of the variations.** Using identical events removes the independent-sample MC error (a naive independent-Poisson error would be 1.9-5.2x, median 2.9x, larger), but it does not remove the *churn noise*: each variation moves 2-13% of selected events in or out of the selection, and the net ratio in a bin is a small difference of two large, nearly cancelling flows. What remains is a per-bin error of 0.4-0.6% (bulk bins, median) with 0.1-2% extremes, and 0.5-0.7% (ICARUS) or 1.4-4.7% (SBND) in the last bin. PROfit uses each bin raw (hist1d, no error), so this noise enters the covariance as if it were systematic. The noise is well understood (calibrated below, no inflation factor needed) and can be removed by a validated smoothing.

All numbers: pair tables of `runs/syst_study/pairs/` (reproduce the production histograms to ~1e-3; exception: ICARUS Run4 WireModxThetaXW, Section 5), cvwgt not used, no event weights. Scripts `plotting/syst_study/q2_*.py`; JSON in `docs/syst_study/q2/`. Summary JSON: `q2/q2_summary.json`.

## 0. Method and assumptions

* Ratio r_b = V_b/N_b, events are the independent units. Per-event contribution d_ib = 1[bin_v=b] - r_b 1[bin_n=b]. Delta method: Cov(r_b,r_b') = sum_i d_ib d_ib' / (N_b N_b'). In terms of the cross-classified table, sigma_b^2 N_b^2 = r_b^2 (L_b + O_b) + (G_b + I_b) + S_b (1-r_b)^2 with L = lost (only nominal selected), G = gained (only varied), O/I = both selected, migrating out of / into bin b, S = both selected, same bin. These are the counts the calculator could output per bin to get the error analytically. Neighbouring-bin noise correlation is -0.17 (-0.22..-0.11) because migrating events leave one bin and enter the next.
* Joint bootstrap: Poisson(1) weights generated from a hash of the event key and replicate index, so the same event gets the same weight in every pair table (needed for cross-variation covariances; SBND variations share the nominal sample).
* Smooth models: GLS Legendre polynomials in ln E (x scaled to [-1,1]) and a penalised (2nd divided difference) Whittaker smoother with leave-one-out CV; the polynomial basis in ln E is an assumption (results are similar for Whittaker, Section 2).
* Everything here measures *event-sampling + re-simulation* noise (the per-event flips are random draws, so both are inside the per-event variance). It cannot detect processing defects that are not random across events (Section 5 finds one).

## 1. Calibration of the per-bin uncertainty

Figure `q2_fig1_calibration.png`.

{{T1}}

* Bootstrap vs analytic: median ratio 1.002 over the 20 pairs (range 0.991-1.008), per-bin range 0.97-1.04; bootstrap bias < 0.07 sigma. The two agree because both are the same first-order quantity; this is a check of the implementation, not independent validation.
* Split-half (independent test): 400 random event splits; E[(r_A - r_B)^2]/E[sigma_A^2+sigma_B^2] = 0.99 (per pair 0.89-1.04; per-bin sd 0.07) and chi2/ndf of (r_A - r_B) = 1.00 (0.98-1.03). So the analytic variance describes the scatter between independent subsamples. Caveat: a split of a fixed sample tests the per-event variance and the independence of events, not the realised noise of the full sample; the latter is tested in Section 2 (residual pulls).
* Event correlations: splits by Run (SBND, 4400 runs, ~80 events per run) or by (Run,Subrun) (ICARUS, 171-321 blocks): chi2/ndf 1.007 (0.96-1.12); delete-one-block jackknife/event-level sigma = 1.005 (0.975-1.043), block bootstrap 1.000; SBND Run//25 super-blocks (181 blocks) 1.002. **No inflation factor is needed** (statistical noise of the jackknife itself is ~5% for G=171). Per-run lost/gained fractions are binomial (sd of z = 1.0-1.1) for all pairs except one.
* The exception: ICARUS Run4 WireModxThetaXW has block-jackknife/event = 1.16 (split chi2/ndf 1.36), caused by 4-6 runs (12965-12970) where the varied sample contains duplicated events (events tree: 27756 duplicated rows). After removing those runs it is 1.01. See Section 5.
* The paired error is 2-5x smaller than an independent-sample error: pairing is doing its job, and the remaining error is exactly the churn noise.

## 2. Signal versus noise in the observed shapes

Figures `q2_fig2_shapes_*.png`, `q2_fig3_signal_noise.png`. "Class" rules (fixed in advance; P=0.01): (iii) chi2 of r-1 against zero (17 bins, full covariance) has p > 0.01, or BIC prefers the zero model; (ii) otherwise, BIC prefers a constant or the extra polynomial terms are not significant (nested chi2 test p > 0.01); (i) otherwise. A "*" marks a residual chi2 about the BIC model with p < 0.05 (excess non-smooth structure).

{{T2}}

Findings:

* **Noise model validated at the level of the realised data.** Residual chi2 about fixed-degree GLS fits pooled over the 20 pairs: degree 2: {{F2}}, degree 3: {{F3}}, degree 4: {{F4}} per dof; the p-values of the 20 residual chi2 are consistent with uniform (KS p = {{F3KS}} at degree 3). Pooled bulk-bin pulls about the smooth fit have sd 0.93 (expected ~0.89 given the fitted parameters). So the bin-to-bin scatter about the smooth trends is statistical, with the analytic sigma. Individual outliers: SBND WireModYZ residual chi2/ndf = 0.24 (p-lower-tail 0.002, one in 20 is not alarming) and SBND ChargeScale 1.76 and ICARUS-R4 SCE0x 1.82 (p = 0.04, 0.03).
* Classification: of the 20 production pairs 11 are (i), 5 are (ii), 4 are (iii) (the Run4 WireModxThetaXW changes from (iii) to (i) when cleaned). (iii): SBND WireModxThetaXW (chi2(0)/17 = 0.5; all 17 bins within 2 sigma; mean |r-1| 0.37% *below* the 0.59% expected from pure noise), ICARUS-R2 WireModYZ (p = 0.06), ICARUS-R4 WireModYZ (p = 0.77), ICARUS-R4 WireModxThetaXW (raw, dup-biased). (ii) constant offsets: SBND SCE2x (+1.5%), ICARUS-R2 ChargeScale (-1.1%) and SCE2x, ICARUS-R4 ChargeScale (-0.9%) and SCE2x.
* One-sidedness: for pure noise E|r-1| = 0.80 sigma per bin (rms 1 sigma); e.g. sigma = 0.6% gives a "systematic" with mean |r-1| ~0.5% from nothing. The observed mean |r-1| exceeds that noise floor by 2-15x for all large (i) variations (GainHigh 5.7% vs 0.48%, ChargeScale 2.2% vs 0.39%, WireModYZ 2.7% vs 0.60%) but is at or below the floor for the (iii) pairs. The number of bins with |r-1| < 2 sigma is 17/17 for the (iii) pairs but 11-15/17 for noise-compatible ICARUS shapes (pure noise: 16.2/17), and 1-5/17 for the large SBND shapes.
* Fraction of bin-to-bin variance consistent with noise (first/second differences; column "d1/d2"): 1.0 (i.e. all) for SBND WireModxThetaXW, DENT, SCE2x, ICARUS-R2 WireModYZ/SCE2x/WireModxThetaXW, ICARUS-R4 ChargeScale/SCE2x/WireModYZ; 0.4-0.5 for SBND ChargeScale/GainLow; ~0.8 SCE0x. Second-difference chi2/ndf (expected 1 + curvature) pooled mean 0.92.
* Signal amplitudes versus noise: debiased signal rms/noise rms is 5-8 (SBND ChargeScale, ChargeScale15, GainHigh/Low), 2-3 (MattModYZ, WireModYZ, ICARUS ChargeScale), and ~1 or less for SCE/WireModxThetaXW/WireModYZ in ICARUS and SBND WireModxThetaXW/DENT. Note the rms noise includes the large last-bin error (SBND), so the per-bin comparison is more favourable for the bulk bins.
* Out-of-sample test (`q2_splitcv.json`, 100 splits, both directions; T = sum_b (r_B,b - f_A,b)^2/sigma_b^2, floor = irreducible noise of half B = 2 per bin; the raw half-A ratio gives T/floor = 2.0 by construction):

{{T7}}
 The smooth fit predicts the held-out raw ratio of the other half better than the raw half does, for all 20 pairs; it removes ~88% (poly-BIC) / ~80% (Whittaker) of the removable excess. Setting everything to zero is catastrophic for most (T/floor up to 114) and good only for the (iii) pairs. Last bin: smoothing it helps on average (1.42 vs 1.97) but hurts for SBND ChargeScale, ChargeScale15 and GainLow (T/floor 2.8-3.4 vs ~2.0 raw): the extrapolated deg-2 prediction misses the last bin by -2.1, -2.7, -3.0 sigma for exactly those three (ICARUS pulls up to +2.3). A pooled empirical-Bayes edge variance is not significant (tau = {{TAU}}, 68% UL {{TAUUL}}) because it is heterogeneous (3/20 edge-affected, all energy-scale-like, all negative, as expected from overflow migration out of the 2.0 GeV edge).

## 3. Noise in the covariance PROfit builds

C = (r-1)(r-1)^T has E[C_raw] = f f^T + Sigma_noise: the noise adds the full (diagonal-dominated) noise covariance, trace += sum sigma_b^2, and is not removable by averaging because the one-sided (r-1)^2 never averages to zero (mirrored/symmetrised use gives the same, since the square is sign-blind). Smooth fit = recommended algorithm (Section 6). Figure `q2_fig4_covariance.png`.

Per variation (trace in 1e-4; noise/raw is expected-noise to observed trace; a ratio near or above 1 means the "systematic" is noise):

{{T3}}

Per detector (variations summed as independent covariances; "xml-active" = ChargeScale, SCE (ICARUS), WireModxThetaXW, WireModYZ, DENT as in standard.xml; Gain is commented out there; the ICARUS-R4 total uses the clean WireModxThetaXW). MC: truth = smooth fit, noise = joint bootstrap deviations of the same events.

{{T4}}

* SBND: noise is 6.8% of the trace of all ten variations but 21% for the xml-active set (noise trace = 0.47 x the ChargeScale trace). The first eigenvector is stable (overlap 0.997); the noise-only matrix would have 58% of its trace in the first eigenvector (dominated by the last bin and by the common-mode noise from the shared nominal sample), so noise is mostly an incoherent diagonal addition. In SBND 56-65% of each variation's noise trace sits in the last bin.
* ICARUS: noise is 40-43% of the raw trace; the noise trace is 1.7-2.0x that of the dominant (ChargeScale) physical covariance. Smoothing drops the effective rank from 2.7 (2.1) to 1.5 (1.4), moves the first-eigenvalue fraction from 0.50 (0.64) to 0.81 (0.83), and the adjacent-bin correlation of the summed matrix from 0.43 (0.65) to 0.99. Raw single-variation matrices are rank-1 with correlation +-1; for pairs consistent with noise ~45% of bin-pair signs flip between raw and the fit.
* Interpretation: for ICARUS the raw covariance contains a bin-to-bin uncorrelated "systematic" of order 0.7% per bin that has no physical origin; for SBND it matters mainly for SCE0x/2x, WireModxThetaXW, DENT (noise/raw 0.3-4.2) and for the last bin.

## 4. Consistency tests as noise probes

Covariances between ratio vectors use the exact analytic cross-covariance over shared events (SBND variations share the nominal and the pool; verified against the joint bootstrap, diagonal ratio 0.88-1.08). Mean noise correlation between different SBND variations is 0.26 (up to 0.62): the shared nominal makes noise positively correlated across variations, so sums of SBND covariances have a common-mode noise term. ICARUS Run2 and Run4 share no events (0 common keys).

{{T5}}

Figure `q2_fig5_consistency.png`.

* ICARUS Run2 vs Run4 (independent samples): pooled chi2/ndf over the five variations {{ICC}}. Differences are consistent with pure statistical noise; no evidence for real Run2-Run4 differences or for under-estimated errors at the present sensitivity (SCE0x 1.59, p = 0.06 is the largest). They could therefore be fitted jointly (reduces noise by ~sqrt(2)) if the physics expectation is the same shape; this was not applied because the detector conditions of the runs may legitimately differ.
* SBND ChargeScale vs ChargeScale15: not equal (different amplitude; chi2/ndf 4.6), but the proportional model y(CS) = 0.78 y(CS15) fits with chi2/16 = 1.32 (p = 0.18) using the exact cross-covariance (rho_noise = +0.62). The differences are a smooth scale; the noise structure is as predicted.
* SCE0x vs SCE2x: antisymmetry y(2x) = -y(0x) holds within noise in ICARUS R2/R4 (p = 0.06/0.10) and marginally in SBND (chi2/ndf 1.93, p = 0.012); the quadratic part (0x+2x)/2 has rms 0.71%/0.79%/0.50% against expected noise 0.79%/0.60%/0.47%, i.e. compatible with zero. The response is linear in the SCE scaling within errors, and the linear part is real in SBND (1.51% vs 0.50% noise) but ~1.5-2 sigma-level per bin in ICARUS. Consequence: a two-point spline through (0x, nominal, 2x) interpolates two noisy points per bin; the half-difference (y2-y0)/2 has half the noise variance of either point.
* SBND WireModYZ vs MattModYZ: equality rejected (amplitudes 3.2% vs 1.1%), but the proportional model (alpha = 6.1) has chi2/16 = 2.2 (p = 0.004): MattModYZ (about 6x smaller) is not purely a scaled WireModYZ at the 2-3 sigma level; with MattModYZ noise (0.44%) scaled by 6 the comparison has little power. GainHigh+GainLow: strongly non-antisymmetric (quadratic part 2.5% vs noise 0.6%): real nonlinear response, not noise.
* Reading: chi2 beyond statistical in a comparison means either real physics (smooth differences, which are also seen as significant fits of the difference with 1-2 polynomial terms) or under-estimated error; the clean outcome of all tests where the same physics is expected (Run2/Run4, CS/CS15 scaled, SCE linearity) is chi2/ndf ~ 1-1.6, which supports the error model.

## 5. Pool and duplicate-key checks

Figure `q2_fig6_pool.png`.

{{T6}}

* Pool sizes differ between variations by up to 25% for the same CV (SBND: 5.1M (WireModxThetaXW) to 6.3M events; DENT has its own, disjoint CV sample of 3.7M events and no common pool; ICARUS-R2 1.25-1.6M, R4 2.1-2.8M). The common pool is 47-80% of the rows. Restricting to the common pool leaves the ratios unchanged: chi2 {{CPC}} summed over the 19 pairs, z sd {{CPZ}}, largest single bin 2.5 sigma; the error grows by 1.1-1.5x as expected from the lost statistics. So the pool definition does not bias the ratios, but differing pools make the noise of different variations differ (~20-30%), and missing events (failed jobs) cost precision for free: re-running them is the cheapest noise reduction.
* Duplicated keys (0.6% of selected rows in SBND, 0.1% in ICARUS) are present in both nominal and varied samples; the production ratio (all rows) differs from the unique-key ratio by a median rms of 0.08 sigma per bin (max 0.31 sigma): negligible.
* **Defect found: ICARUS Run4 WireModxThetaXW.** In runs 12965-12970 the varied sample contains duplicated events (27756 duplicated rows in the events tree; 1786 in the selected tree vs ~160 expected). The production ratio counts them twice, giving a coherent shift: production - clean = {{PRODCLEAN}} % (weighted over 17 bins; chi2 vs the clean ratio {{PRODCHI}}/17), and the pair table, which drops duplicated keys, instead manufactures 456+178+117 'lost' events (z = 52 in the per-run lost fraction) and a -0.65% shift in the opposite direction. The ratio with those runs removed is the trustworthy one (+0.56 +- 0.14% constant, smooth class (i)); the raw production histogram is in the noise class (iii)/(ii) only because the duplicates have cancelled the signal. No other pair has a run-level anomaly (|z| < 4.3 on all others).

## 6. Recommendation

Should the ratios be smoothed before use in PROfit? **Yes**, for every variation (including those with large signals, where smoothing is harmless: GainHigh, SBND ChargeScale have noise/raw 1-4% and change by < 8% in trace), and the largest gain is for ICARUS and for SBND SCE/WireModxThetaXW/DENT.

Validated algorithm ("poly-BIC, last-bin pre-test"):
1. Compute y_b = r_b - 1 and the covariance Sigma (analytic formula of Section 0, or a hashed-Poisson bootstrap; both agree to 0.3%).
2. Bulk bins (1-16): GLS fit of Legendre polynomials in x = ln E rescaled to [-1,1] with Sigma, candidates {zero, degree 0..4}; choose by BIC = chi2 + k ln 16. (AIC/LOO give the same degree for most pairs; BIC is the more conservative.) The zero model makes (iii) shapes vanish and the constant model makes (ii) shapes flat.
3. Last bin (1.5-2.0 GeV): predict it by extrapolating the bulk fit (degree <= 2); if |y_17 - prediction|/sqrt(sigma_17^2 + var_pred) < 2, use the precision-weighted combination of the raw value and the prediction with an extra edge variance tau^2 = (0.95%)^2; otherwise keep the raw value (overflow-migration edge effect; seen in 3/20, all energy-scale-like). Conservative alternative: always keep the last bin raw.
4. Build the covariance from r_hat = 1 + f_hat. Optionally add the (small) fit-parameter uncertainty H Sigma H^T.
5. Cross-check with the Whittaker-CV smoother; a bin difference > 1 sigma between the two flags a problem. Do not use the "debiased" covariance y y^T - Sigma (not positive semi-definite, high-variance).

Validation (closure with the data's own noise; truth scenarios from smooth fits of the 20 shapes, with +-1, +-2 sigma last-bin edge offsets, zero, and constant offsets of 0.5 and 1 rms sigma; MSE relative to the raw ratio = 1.0; mean over 20 variations; `q2_fig7_validation.png`):

{{T8}}

* Smooth truths: noise variance removed to 0.34-0.52 of raw (all-bin variants 0.21-0.68, last-bin raw 0.42-0.49); zero/constant truths: 0.03-0.35 (the residual is the raw last bin). Last-bin edge features of 1-2 sigma bias the pre-test variant by 0.2-0.3 sigma_last (MSE 0.5-0.7 of raw), the all-bin smoothers by up to 1.2 sigma and make the last bin worse than raw (MSE up to 3x); "last raw" is unbiased in every scenario. Out-of-sample on the data: Section 2 table.
* Risks: (a) the high-energy edge bin: both smoothness and polynomial extrapolation fail when the overflow-migration effect is real; (b) a high-degree polynomial can oscillate (cap at degree 4); (c) smoothness in ln E is assumed: genuine sharp PID-boundary structure (e.g. SBND ChargeScale residual p = 0.04, ICARUS-R4 SCE0x p = 0.03) is partially removed (the excess is ~0.6 sigma per bin, below the noise); (d) the first bin (0.3-0.4 GeV) is threshold-affected for GainHigh (degree 1 misses it by 10 sigma; degree >= 2 is required, which BIC selects); (e) selecting by BIC on the same data slightly underestimates the uncertainty of the selected curve.

Per-variation recommendation (constant +- err is the GLS constant over all 17 bins; for (ii)/(iii) this is the quantity to retain, e.g. as a normalisation-type uncertainty with the stated error rather than 17 independent noisy bins):

{{T9}}

Noise-dominated systematics (noise trace >= 0.5 of the raw trace): SBND WireModxThetaXW (4.2), DENT (2.1), ICARUS-R2 SCE2x (0.54), WireModxThetaXW (0.56), WireModYZ (0.92), ICARUS-R4 SCE2x (0.81), WireModYZ (1.76). For these, fitting a constant or zero is statistically preferred and the 17-bin raw covariance is not meaningful.

Other practical actions: (1) fix/deduplicate the ICARUS-R4 WireModxThetaXW varied sample (and add a duplicate-event check to the calculator); (2) increase statistics where the pool is smallest (failed jobs), cheaper than any modelling; (3) widen the last bin or add an overflow bin if the 1.5-2.0 GeV bin matters (its SBND error is 1.4-4.7%); (4) output category counts (L, G, in, out, stay) per bin in the calculator so that the analytic Sigma is available to PROfit/smoothing without a bootstrap.

## 7. Limitations and honest caveats

* The delta-method variance is first order; it is validated by the bootstrap (includes non-linearity, 0.97-1.04 per bin) and splits, but the last bin of SBND (few hundred events) is non-Gaussian.
* The BIC degree and the class labels depend on the thresholds (P = 0.01); classifications near the boundary (ICARUS-R2 WireModYZ p = 0.06, ICARUS-R4 WireModxThetaXW) are not robust. Under-dispersion (WireModYZ SBND, residual chi2/ndf 0.24) is not significant after look-elsewhere, but if real it would mean the analytic error is conservative for that pair.
* Combinations (SCE as a two-point spline, "SCE_icarus_RunX", Gain, trigger emulation, ProtonMiss etc.) were not analysed except SCE0x/2x as separate pairs; `Gain` is commented out in standard.xml.
* Systematic non-statistical issues in the varied samples (like the duplicates) are invisible to a statistical error model; the run-level binomial scan is the only tool used to find them.
* The recommended last-bin treatment trades a small bias risk for noise reduction; the conservative choice (last bin raw) is the safe default if the edge bin is physics-critical.

## Files

Scripts: `plotting/syst_study/q2_{cache,lib,boot,fit,smooth,poolkeys,clean,s1_calibration,s1b_clean_block,s2_signal_noise,s2b_splitcv,s3_covariance,s4_consistency,s5_pool,s6_validate,s7_summary,figs,report}.py` (+ `q2_report_text.md`, the text template). Results: `docs/syst_study/q2/*.json`, figures `q2_fig1..7_*.png`.
