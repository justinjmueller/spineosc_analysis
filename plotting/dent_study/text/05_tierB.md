## 5. Tier B: production CV with the full systematic covariance
Prediction: `spineosc_numu_sbnd_v1.1.0_wsyst_cvrw.root` (`cvwgt` applied) with every systematic, area-normalised and shape-only.
Four variants of the same prediction: `v1` (standard, DENT in the covariance), `exDENT` (DENT removed from the covariance), and the CV
shifted by 0.5 or 1.0 times the DENT ratio (`inj0p5`, `inj1p0`) with DENT kept in the covariance. DENT is never fitted. The
shifted-CV and DENT-removed variants cannot be combined (PROfit then ignores the injection), so that combination is not shown.

**What the table says.**
* The production CV is already acceptable in six of the ten primaries (p = 0.11 to 0.59 for `crosser_flag`, `muon_proton_score`,
  `muon_abs_endx`, `muon_endx_full`, `proton_calo_csda`, `proton_muon_score`). It fails for `muon_calo_csda` (p = 1.7e-10),
  `proton_pion_score` (1.2e-8), `proton_mip_score` (2.2e-4) and `muon_pion_score` (0.010). Injecting DENT repairs none of them: for
  `proton_pion_score`, `proton_mip_score` and `muon_pion_score` the change is at most 1.4 units, and the `muon_calo_csda` change is an artefact (section 7).
* The change in chi-square from the full DENT shift is below 2 units in nine of ten primaries, with five improving and five worsening.
* **DENT as an uncertainty matters for the endpoint.** Removing it from the covariance raises chi-square for `muon_abs_endx` from 19.0 to 40.3
  (p 0.46 to 0.003), for `muon_endx_full` from 29.9 to 62.0 (p 0.15 to 2e-5), and for `muon_pion_score` from 30.4 to 41.0.
