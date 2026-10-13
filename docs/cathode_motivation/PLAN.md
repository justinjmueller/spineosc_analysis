# Cathode-crosser and end-point cut: motivation study (pre-registration of the 4 cm criterion)

Written 2026-10-07, before the cut-value scan was run in its final form. Disclosure: the analyst had already seen exploratory comparisons at 3, 5, 10 and 20 cm
(`plotting/dent_validation/endpoint_cut_scan.C`) and the 4 cm / crosser-cut comparisons in the analysis binning, so the choice is hypothesis-prompted and the scan is a
quantitative justification, not an independent test. The criterion below is fixed before the scan figures are made and is evaluated at every cut value.

## Samples and selection
* MC: detvar CV (`sbnd_detvar_newcv`) and DENT (`sbnd_detvar_newdent`) of `spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root`, events matched by
  (Run, Subrun, Evt, event_max_neutrino_energy) over the generated-event lists of both samples, as in the detsys calculator. Selection: `reco_fiducialize_cathode == 1`,
  `reco_veto_sbnd_highy_highz == 1`, `reco_cathode_crosser == 0`, then `abs(reco_leading_muon_end_x) >= X`. Channels by `reco_proton_multiplicity` (1, >1).
* Data: SBND release detector-access file (energy-free columns only). Only fractions of selected events are shown. No track length or length proxy is used anywhere.
* Signal-efficiency cost: production `sbnd_large` selected tree, true categories 0 and 2, and the signal tree for the truth-level cut.

## Scan
X in {0, 1, 2, 3, 4, 5, 6, 8, 10} cm (X = 0 is the crosser-cut selection).

## Quantities (each with a paired bootstrap uncertainty, 500 resamples of the matched events)
* A (DENT systematic flatness): in the coarse reco-energy bins (1p: 0.3-0.7, 0.7-1.0, 1.0-1.5 GeV; Mp: 0.4-0.8, 0.8-1.2, 1.2-2.0 GeV), the DENT/CV - 1 values; summaries are the
  weighted slope against bin centre and the integrated DENT/CV - 1.
* B (pile-up removed): the DENT/CV ratio of events with |x_end| in the 4 cm band just above X, [X, X+4) cm, among non-crossers; and the data/CV/DENT fractions of selected events with |x_end| < X.
* C (cost): fraction of selected true-signal events (CV, crosser-cut selection) removed by the cut, and the fraction of true signal with |true_leading_muon_end_x| < X.

## Criterion
The chosen X is the smallest scan value for which all hold: (A) the slope is consistent with zero within 1 sigma and |integrated DENT/CV - 1| <= 1.5% in both channels;
(B) the DENT/CV ratio in [X, X+4) cm is within 10% of unity; (C) the signal removed is at most 3%. If no value satisfies all, the figure reports that.

## Constraints on presentation
Area-normalized or fractions only, no event counts, no energy from data, PROfit/mutools style, "SBND Simulation" or "SBND Data Access" labels.

## Outcome (2026-10-07)
Scan results are in `docs/cathode_motivation/scan.json`. Evaluated against the criterion: the smallest passing value is **3 cm** (band ratio 1.04 +- 0.02, signal cost 1.9%); 4 cm also passes all three tests
(band ratio 0.97 +- 0.02, cost 2.3%). Decision (2026-10-07): keep **4 cm**, one step of margin, and state in the note that 3 cm is the smallest value that passes. The criterion was not changed after seeing the result.
