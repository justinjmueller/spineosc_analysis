## 1. Data, samples and rules

* **Data.** `spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root`, tree `events/sbnd_bnblight_full/selected`
  (area-normalised comparisons only; no event counts are quoted). This file is cleared for detector-level studies and has no energy or momentum branch. SBND data files that
  carry energy are blinded and were not opened.
* **Production MC (Tier B).** `spineosc_numu_sbnd_v1.1.0_wsyst_cvrw.root`, CV-reweighted (`cvwgt` in every beam-MC weight), with the full
  interaction-model, flux, geant4 and detector systematic suite.
* **Detvar MC (Tier A).** `sbnd_detvar_newcv` and `sbnd_detvar_newdent`. They have no `cvwgt` and no interaction-model
  systematics, so only MC statistics enter their covariance.
* **Selection.** At least one reconstructed proton, `reco_fiducialize_cathode == 1` and `reco_veto_sbnd_highy_highz == 1`. The MC also
  applies the flash-PE threshold, which is trigger emulation and is not applied to data.
* **Rules.** No fitting; every PROfit job uses `--area-norm` and `--shapeonly`; shape chi-square and p-value
  (`ndf = nbins - 1`) are the metrics; systematics are assessed in the observable's own binning; plots follow
  `plotting/studies/data_access.py`.
* **Observables.** Twenty-one energy-free observables were registered, binned by equal-population quantiles of the detvar CV only.
  The **primary set** (ten observables) was chosen by an MC-only criterion: the CV-versus-DENT shape chi-square at data-sized
  statistics, per degree of freedom, is at least 1. The analyst had already seen exploratory data comparisons (section 9), so the choice
  was not blind. Everything else is diagnostic.
* **Plan.** The plan and every change to it after the first result are in `docs/dent_study/PREREGISTRATION.md`.
