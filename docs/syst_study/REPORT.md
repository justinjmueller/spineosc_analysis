# Detector systematics study (v1.1.0), 2026-10-10

Full report: `docs/syst_study/report/index.html` (published as an artifact; figures in `report/figs/`). Sub-reports: `Q2_statistical_content.md`, `Q3_direct_detvar.md`, `PerSystematic_report.md`. Plan: `PLAN.md`.

## Key findings
1. Established: all 20 production ratio histograms reproduce exactly from event-level paired tables. Each variation selects 2-13% of events in only one sample (each direction); the net is 1-12% of that churn; 95-99.6% of flips are production-selection flips.
2. Established: `reco_visible_energy` is range based (medulla `particle_variables.h:234-250`), 67-78% of events keep identical energy under charge-type variations. Energy trends are selection effects at species boundaries.
3. Q1 (supported): the SBND high-energy sample is multi-proton (92% at 1.5-2 GeV vs 39% in ICARUS); flip probability grows with multiplicity; a PID+multiplicity model without energy reproduces the net up to 1.5 GeV and about half of the last-bin deficit (about 600 events).
4. Q2 (established): noise is in the ratios (bootstrap calibrated, split-half chi2/ndf 1.00). 11 smooth, 5 constant, 4 consistent with zero. Noise is 40-43% of the raw covariance trace for ICARUS, 21% for the SBND active set. Smoothing recommended (GLS polynomial in ln E, last bin raw if pull > 2 sigma).
5. Q3 (established): native PROfit DetVar reproduces the ratios to 1e-7 with an intermediate tree carrying every pool event and a `selected` weight (zero-weight patch or 1e-9 epsilon weight). Matched mode as used today cannot.
6. Q4 (supported/open): SBND ChargeScale is a multi-proton selection effect (1p flat, Mp -15% in the last bin), ICARUS a flat -1%. Event-weighted effect sizes are comparable (1.7% vs 1.0-1.2%). The PROfit "no X" scans have a common chi2 offset (every variant ~1.16x), so they cannot rank systematics. Toy: similar SBND/ICARUS effect only if the cross-detector normalization is tight.
7. Defects: ICARUS Run 4 WireModxThetaXW has duplicated events (shipped ratio biased +0.5 to +1.1%); one ratio is applied to both proton channels; ICARUS ratios omit filterEE; truth missing in SBND WireMod varied trees; `detsys/calculator` lacks the endx cut (use `calculator_endx4`).

8. Clarified 2026-10-10: SBND ChargeScale = nonlinear charge scale + ~5% smearing; ChargeScale15 = same with 15% smearing (analysis uses it; fitted SBND value is 5.1%, ICARUS 14.9%). 15% vs 5% raises the SBND effect by about 28% (last bin -13.2% vs -10.3%, integrated -0.65% vs -0.46%); the pure nonlinear scale carries ~67-79% of the ChargeScale15 effect (two-point extrapolation). The 15% is deliberate (same smearing at both detectors although SBND fitted lower).
