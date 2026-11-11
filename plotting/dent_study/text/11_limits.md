## 10. Limitations and what was not done
* **No energy.** Energy and momentum quantities are blinded and embargoed for SBND data, so nothing here tests DENT on the energy distributions, which is where its uncertainty matters. A study of the high-energy region would need to be run by someone cleared to look.
* **Power.** Single observables have d' of about 0.5 to 1. The set of observables shares events and cannot be combined into one test without a joint covariance, which was not built.
* **Detvar CV is not the production nominal.** Tier A chi-squares mix DENT with other data/MC differences. The Tier B work removes this by using the production CV and the DENT ratio.
* **DENT enters as a one-sided ratio.** The ratio is taken from the detvar CV/DENT pair and applied to the production CV. Whether that is equivalent to using the DENT sample directly was checked
  for the native PROfit detvar mode (unmatched-event mode agrees); a Tier B run with the native detvar mode was not done.
* **Cathode comparisons are hypothesis-prompted** (the analyst had seen them before the plan was written).
* **Trigger emulation.** The flash-PE cut is applied to the MC only. Its efficiency is not part of the comparison.
* **Not examined:** WireMod Y-Z (SBND uses MattMod Y-Z in the systematics), the cathode-region pile-up with an endpoint cut applied to the data comparison, and any ICARUS DENT-like variation.
