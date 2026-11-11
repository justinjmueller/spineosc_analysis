# DENT narrative audit (2026-10-08)

Script: `plotting/dent_study/dent_audit.py` (numbers in `audit.json`); the truth-mechanism check was run inline (see below).
All SBND, fiducial + high-y/z veto + flash threshold on MC, >= 1 proton, detector-access data (no energy used for data).

## What holds
1. **The end-point pile-up is the robust discriminator.** Fraction of non-crossing selected events with |x_end| < 4 cm:
   data 2.10 +- 0.07 %, production CV (cvwgt) 2.35 %, detvar CV 2.31 %, DENT 4.54 %. Data is *below* the CV by ~10 % and a factor ~2 below DENT.
   cvwgt barely moves it (2.32 -> 2.35). East and west TPC sides agree within each sample (CV 2.33/2.30, DENT 4.52/4.56, data 1.95/2.25 +- 0.1, ~2 sigma apart).
2. **Mechanism.** Events that are in the DENT pile-up but not in the CV pile-up (n = 2441 of 207k common): 70 % were flagged crossers in the CV sample, median CV |x_end| = 4.0 cm (= just outside the window);
   by truth 70 % end within 4 cm of the cathode and 31 % are true crossers. The CV pile-up itself is 96 % true |x_end| < 4 cm and 2 % true crossers. So DENT converts crossers/near-cathode tracks into
   truncated tracks ending at the cathode.
3. **Factorial of the cuts (DENT/CV - 1, integrated, 1p / Mp):** no cut -1.9 / -2.4 %; crosser cut only +1.5 / +1.3 % (high-energy bins +2.5 / +4.6 %); end-point cut only -3.1 / -3.3 %; both -0.95 +- 0.16 / -0.55 +- 0.31 %.
   The two cuts act with opposite signs and are non-additive; neither alone flattens the ratio. Residual with both: -1 % in 1p (6 sigma), <= 1.6 % in the highest Mp bin.

## What needs a caveat
4. **The crosser fraction is not a clean DENT discriminator.** Detvar CV 27.9 %, data 27.46 +- 0.18 %, DENT 25.4 % (as quoted), but the *production* CV with cvwgt is **29.0 %** (27.8 % without cvwgt).
   The interaction-model reweight alone moves the crosser fraction by +1.2 points. Against the production CV the data sit 1.5 points low, i.e. between CV and the DENT-shifted CV (~26.4 %).
   The crosser fraction therefore depends on the generator model and should not be presented as evidence that the data prefer the CV over DENT. The pile-up (point 1) is unaffected.
   Suggested rewording of `detector_data_access.tex` line ~54 and `systematics.tex` line ~218: quote the crosser fraction as a description of the mechanism (DENT lowers it) and note that its absolute value depends on the interaction-model reweighting; base the data comparison on the pile-up only.
5. A true null test of the DENT/CV pairing is not possible: the detvar CV shares no events with `sbnd_detvar_cv` or with the production nominal (different generation), so the baseline noise of a nominal-vs-nominal ratio cannot be measured. The bootstrap errors are the only handle.

## Not yet done
* x = 0 calibration offset (the pile-up window is symmetric; the east/west difference in data is 2 sigma and consistent with zero offset).
* Whether the pile-up excess survives in ICARUS-like observables (not applicable) and the PID impact for non-crossers (previously ~flat).
