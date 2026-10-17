# DENT cathode track-breaking study: plan (written 2026-10-14, before the first run)

Question: how much of the DENT shift in the leading-muon PID scores comes from (a) tracks that are
**broken** at the cathode and not stitched (reco crosser in the CV, non-crosser in DENT), versus
(b) tracks that remain **stitched** (crosser in both) but acquire a **gap** (missing charge/length)?

Disclosure: the analyst has already seen the aggregate DENT cathode results (`docs/cathode_motivation`),
so the hypothesis is prompted by them; this study adds per-event matching, it is not independent of them.

## Samples (MC only, no data)
Detvar file `spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root`, `signal` and `selected` trees.
* Primary pair: `sbnd_detvar_newcv` -> `sbnd_detvar_newdent` (the pair the production DENT ratio uses).
* Replication: `sbnd_detvar_cv` -> `sbnd_detvar_dent`.
* Control for re-simulation noise: `sbnd_large` (nosyst file) -> `sbnd_detvar_charge_scale` and `smear15`,
  a variation with no cathode-specific physics. Excess of DENT over control is attributed to DENT.
Events are matched on (Run, Subrun, Evt, event_max_neutrino_energy); duplicated keys are dropped.

## Populations and classes
True signal events with a reconstructed leading muon in both samples. `true_cathode_crosser` splits
true crossers from true non-crossers. For true crossers, by the reco crossing flag (CV -> variation):
* **break**: crosser -> non-crosser;
* **stitched**: crosser -> crosser, split into
  * **gap**: start and end points both move by < TOL_POINT (5 cm) but the track length falls by more
    than TOL_LENGTH (3 cm). Same end points with less path length is the signature of missing charge
    inside a still-joined track;
  * **intact**: all other stitched tracks;
* **gain**: non-crosser -> crosser (reported, not decomposed);
and true non-crossers as a cathode-free reference.
The two tolerances are parameters; the report includes a scan so the split is not a single cut.

## Quantities (paired bootstrap, 500 resamples of matched events)
1. Transition matrix of the reco crossing flag among true crossers, DENT and controls.
2. For each class: mean change of the leading-muon PID scores (muon softmax, primary softmax, MIP softmax)
   and of the length and calorimetry/range difference, variation minus nominal.
3. Decomposition of the total mean score shift among true signal events into class contributions
   (class fraction x class mean shift); controls give the baseline for each class.
4. Where the break ends up: |end_x| of the DENT track for broken events versus the nominal.

## Interpretation rule (fixed now)
Call a cause "present" only if its DENT contribution exceeds the same class's control contribution by
more than 2 bootstrap sigma. If the endpoint-based gap definition is not stable under the tolerance scan
(share changes by more than a factor of two across the scan) say so rather than quoting one number.

## Outcome (2026-10-14)
See `REPORT.md`. Break present (about 4.7 % of true crossers, 115-140 times the controls, replicated). Gap not isolated: the
tolerance-scan stability condition failed (0.01 % to 17 %), so no gap number is quoted.
