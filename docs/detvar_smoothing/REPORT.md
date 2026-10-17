# Smoothing the detector-systematic ratios: results so far

Written 2026-10-14. Simulation only. Details, designs and every intermediate result are in `PLAN.md`; scripts are in `scripts/detvar_smoothing/`.

## Summary
* **Closure test 1 (inject 1 sigma, global fit, all systematics floating) passes** for the muon score at S2 for DENT, ChargeScale15, WireModxThetaXW and MattModYZ (fitted injected parameter 0.96, 1.00, 0.85, 0.98; total chi2 equals the
  pull squared because PROfit's chi2 includes the nuisance penalty). SCE is not a floated systematic in the SBND suite, so it is not a test case.
* **Closure test 2 (variation sample as pseudo-data)** uses a three-way event split (ratio part, MC part, pseudo-data part) so the ratio noise is independent of everything else. With unsmoothed ratios, ChargeScale15 is pulled to
  0.86 +- 0.03 and MattModYZ to 0.91 +- 0.16; DENT at S2 is indistinguishable from the null; WireModxThetaXW is degenerate with MattModYZ. The excess chi2 of ChargeScale15 varies from 0.5 to 13 across the three splits: that is the noise of a one-third-sample ratio.
* **Smoothers against a known truth** (six truth shapes, noise at the real paired-bootstrap level): noise-dominated ratios gain from almost any smoother; ratios with real structure gain little or lose; polynomials, Gaussian kernels and cross-validated variants erase most of a 10% last-bin step,
  while the smoothing spline with a GCV penalty keeps it (2% residual bias) and has a worst case near 1 (0.71-1.28).
* **In the PROfit closure the smoothing spline is best:** it moves the pulls of real effects towards 1 and lowers their excess chi2 (ChargeScale15 5.2 -> 2.8, MattModYZ 8.2 -> 4.0) without harming the noise-dominated cases; a cross-validated polynomial makes ChargeScale15 worse (5.2 -> 13.2).

## Limits
Three splits; one observable (muon score, S2); the analysis's reco-energy ratios and their sharp last-bin structure are tested only in the synthetic study; smoothing is applied in the observable the ratio was made in. A conclusion about smoothing the analysis ratios needs the reco-energy closure and more splits.

## Open question for the analyst
Whether to extend the PROfit closure to a reco-energy configuration (the ratios used by the analysis) before choosing a smoother; the synthetic study says the last-bin feature is where methods differ most.
