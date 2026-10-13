## 8. Supplementary: PROfit's own pseudo-experiments (50 throws)
These are PROfit `--pseudo-experiment` throws (Poisson plus Gaussian nuisance throws) with the DENT strength injected into the fake data at 0, 0.5 or 1.0, then
the same chi-square change as in section 5. They are kept for completeness. Because PROfit re-centers the covariance on the injected prediction, they
carry the artefact described in section 7 (for example `muon_calo_csda` has a mean change of 162 with a standard deviation of 268, which is not meaningful).
The calibrated version is the replica analysis in section 7. Poisson-only fake data (an earlier version of this tier, replaced
on 2026-10-07) was too narrow to calibrate anything: its spread was 0.2 to 0.3 units against 1 to 3 for the pseudo-experiments.
