# Detector-systematics investigation: plan (written 2026-10-10, before the analyses)

Goal: a critical, quantitative account of what drives each detector systematic (SBND: ChargeScale / ChargeScale15, SCE0x/2x, WireModxThetaXW, WireModYZ, MattModYZ, Gain, DENT, trigger emulation; ICARUS Run 2/4: ChargeScale, SCE0x/2x, WireModxThetaXW, WireModYZ, trigger emulation).
Inputs: v1.1.0 detvar files and the production ratio histograms `detvar/detsys-v1.1.0v2.root`. cvwgt is not used anywhere.
Infrastructure: `plotting/syst_study/{common,build_pairs,build_pairs2}.py` build paired event tables (runs/syst_study/pairs*/). The pair-table ratio reproduces the production histogram exactly (ChargeScale15_sbnd_Run1: max |diff| = 0.0).

## Questions and how each will be answered
Q1 Why do some ratios trend with energy?
  Hypotheses: (H1a) energy-scale shift x steep spectrum; (H1b) overflow migration past the last bin edge (2.0 GeV); (H1c) selection churn growing with energy (pion-veto/PID boundary, longer tracks, containment);
  (H1d) background/class mixture changing with energy; (H1e) multiplicity-dependent effects (1p vs Mp), (H1f) statistical.
  Tests: ratio decomposition into lost / gained / migration; energy-only and selection-only counterfactuals; event-level shifts of E_vis, muon KE, proton KE, scores; spectral-slope prediction r-1 ~ -d ln N/d ln E x delta; truth-class composition.
Q2 Statistical content of the ratios/splines.
  Tests: paired Poisson bootstrap error per bin; second-difference (smoothness) statistic against the bootstrap error; split-half cross-validation (does a smooth fit predict the held-out half better than the raw binned ratio?);
  expected |r-1| noise floor; comparison of near-null variations.
Q3 Direct PROfit DetVar reproduction (delegated; see Q3_direct_detvar.md).
Q4 Charge scale deep dive: per-particle (muon, leading proton) kinematic and score shifts, which selection decisions flip, SBND vs ICARUS, ChargeScale vs ChargeScale15 definitions, and the sensitivity impact (surf files).

## Standing caveats to keep in mind
* Only selected events are stored; the reason a non-selected event fails cannot be read for the side where it fails except where the row is in the selected tree but fails a post-cut (pairs2 tables).
* Ratios are computed on events common to both samples; pools differ per variation.
* Production ICARUS ratios are computed WITHOUT reco_user_containment_cut_filterEE although the analysis XMLs apply it (to be quantified).
* The ratio histograms are used per bin without smoothing.
