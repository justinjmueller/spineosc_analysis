## 7. Is the DENT shift an improvement? Fixed-covariance analysis and pseudo-experiments
**Why PROfit's own change in chi-square is not enough.** With `--inject-systs-cv` PROfit re-derives the systematic covariance
around the *shifted* prediction. For `muon_calo_csda` the first-bin fractional DENT covariance falls from 0.84 to 0.51 and the
statistics-only chi-square from 14,606 to 3,556, so PROfit's chi-square at strength 1.0 mixes the central-value shift with a change of the
covariance. PROfit reports a change of -18.8 for that observable; with the covariance held at its strength-0 form it is +0.20.

**The replica.** A numpy replica of PROfit's shape-only chi-square, `chi2 = (n-p)^T [diag(n) + F o (p p^T)]^-1 (n-p)`, with `F` the
collapsed total fractional covariance and `p` the prediction area-normalised to the data, reproduces PROfit at all three strengths
(maximum difference printed below). Holding `F` at its strength-0 value for both predictions gives the **fixed-covariance change in chi-square**,
`chi2(0) - chi2(1)`; positive favours the DENT-shifted prediction. The same replica runs 20,000 pseudo-experiments per truth (Poisson
fluctuations plus Gaussian systematic throws through `F`), with truth at strength 0, 0.5 or 1.0, to calibrate it. The p-value is two-sided and
empirical, with a floor of 2/(N+1). No strength is estimated; the three truths are fixed points.
