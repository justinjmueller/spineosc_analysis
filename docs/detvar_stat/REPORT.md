# Stochastic (finite-MC-statistics) contribution to the detector-model uncertainties

Script: `plotting/studies/detvar_stat_bootstrap.py`. Simulation only. SBND, S2 selection, the detsys energy binning (17 bins).

## Method
Each detector systematic is a per-bin ratio varied / nominal. Paired Poisson bootstrap: every event of the union of the two matched samples
gets one Poisson(1) weight per replica, the same in both samples, so the pairing is kept; events in only one sample are resampled
independently. 500 replicas. Per bin: `delta = r - 1`, `sigma_stat` (replica spread), `z = delta / sigma_stat`, and the noise-subtracted size
`sqrt(max(delta^2 - sigma^2, 0))`. Nothing is smoothed or modified.

## Calibration (`--null-test`)
Random half-splits of one sample (true ratio 1), 300 repeats: z has rms 0.95 (1.3-1.5 GeV) and 1.02 (1.5-2.0 GeV), with 2.7% and 5.0% of
values above |z| = 2 (ideal 1.00 and 4.6%). Over all bins rms 1.00. The bootstrap sigma is the right size, including in the low-count bins.

## Results (rms over the 17 bins)
| Systematic | rms delta | rms sigma_stat | noise fraction | bins with abs(z) < 2 |
|---|---|---|---|---|
| DENT | 1.0% | 1.5% | 208% | 12 |
| ChargeScale | 3.6% | 0.7% | 3.5% | 3 |
| ChargeScale15 | 4.3% | 0.7% | 2.7% | 2 |
| SCE0x | 2.3% | 0.9% | 15.7% | 5 |
| SCE2x | 0.5% | 0.9% | 296% | 16 |
| WireModxThetaXW | 12.6% | 1.0% | 0.6% | 0 |
| WireModYZ | 3.0% | 1.1% | 13.2% | 6 |
| MattModYZ | 3.5% | 0.5% | 1.9% | 1 |
| GainHigh | 7.6% | 0.9% | 1.3% | 2 |
| GainLow | 5.4% | 0.7% | 1.5% | 0 |

"Noise fraction" is the sum of sigma^2 over the sum of delta^2. For most variations the shifts are real and the noise is a few percent of
the squared shift. **SCE2x is consistent with zero in 16 of 17 bins; its current treatment as a 1 sigma systematic is almost all noise.**
SCE0x and WireModYZ are 13-17% noise. **DENT is also consistent with noise: rms shift 1.0% against a statistical rms of 1.5%, 12 of 17 bins within 2 sigma.**

## The 1.5-2.0 GeV bin (multi-proton tail)
ChargeScale -10.6 +- 2.1% (z = -5.2), ChargeScale15 -13.7 +- 2.2%, WireModxThetaXW -11.3 +- 3.0%, GainHigh +14.2 +- 2.7%, GainLow -16.2 +- 2.1%
(all significant), DENT 0.0 +- 4.6%, SCE0x -5.2 +- 2.8%, WireModYZ +6.5 +- 3.4%, SCE2x -0.3 +- 3.0%, MattModYZ +0.2 +- 1.6%.
The ChargeScale value is a real effect, not a fluctuation (see `docs/detvar_churn/REPORT.md` for where it comes from).

## Limits
* It quantifies noise from the finite event sample given the pairing. It does not say what the correct systematic is, and it does not
  replace the ratios; a noise-subtracted shift is a diagnostic.
* Systematics other than DENT use the production nominal sample, of which the variations are only ~94% paired; the unpaired events add noise,
  and the bootstrap includes it.
* The weights are uniform; there is no sample reweighting in these ratios.

## Correction (2026-10-14)
The first version of this report gave DENT an rms shift of 5.8% and a 6.5% noise fraction. That was wrong: the DENT pair (`sbnd_detvar_newcv`, 2.89e20 POT, and
`sbnd_detvar_newdent`, 2.74e20 POT) differs in POT by 5.4%, and the bootstrap had not been restricted to the events common to both samples, so the POT ratio
appeared as an offset of about -5% in every bin. The analysis (and the detsys calculator) uses only the common events, and the script now does the same for DENT
(`COMMON_ONLY`). The other variations use the production nominal sample, whose POT matches the detvar samples to within 1% except GainHigh (5.36e20 against 5.59e20,
-4%), which may carry a similar small offset; this was not corrected.
