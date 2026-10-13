## 11. Deviations from the plan (all dated in `docs/dent_study/PREREGISTRATION.md`)
* Removed three multiplicity observables (pion, photon, electron) that the selection forces to zero (24 to 21 observables). None was primary.
* The first exploratory data dumps wrongly applied the MC flash cut to data. They were re-made without it; the headline cathode numbers moved by less than 0.1 sigma.
* After the first Tier C result (`crosser_flag`), a third truth (half DENT) and two-sided p-values were added. This is a third fixed point, not a fit.
* Found that two crossing-split studies register DENT under per-population names, so injections were silently ignored. The driver now detects this; both were re-run and the report audits for any inert injection (none remain).
* Replaced Poisson-only Tier C by pseudo-experiments, then by the fixed-covariance replica, after finding that PROfit's injection re-derives the covariance (section 7).
* Changes after seeing results are listed so that the p-values above can be read with the right amount of caution.
