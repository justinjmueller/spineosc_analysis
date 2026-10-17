# DENT cathode track breaking: per-event results (2026-10-14)

Plan and class definitions: `PLAN.md` (written before the first run). Code: `plotting/studies/dent_cathode_break.py`.
Tables: `transitions.csv`, `decomposition.csv`, `tolerance_scan.csv`, `summary.json`. Figures (git-ignored, regenerate with the
script): `crossing_transitions`, `pid_decomposition`, `gap_signature`, `tolerance_scan`.

Updated with sections 4 (selection fate) and 5 (dependence).

Samples: detvar `signal` trees, true signal with a reconstructed leading muon in both samples, cathode-fiducial and high-y/z vetoes
applied, matched per event, unweighted MC (232,511 matched DENT events; about 63,500 true crossers). Errors are Poisson-bootstrap
(500 resamples). Controls: `sbnd_large` -> ChargeScale and ChargeScale15.

## 1. Breaking at the cathode is real and specific to DENT

Fraction of true cathode crossers that are reco crossers in the CV and non-crossers in the variation:

| pair | break (X -> NX) | gain (NX -> X) |
|---|---|---|
| DENT (newcv -> newdent) | 4.73 +- 0.08 % | 0.109 +- 0.013 % |
| DENT replication (cv -> dent) | 4.70 +- 0.12 % | 0.096 +- 0.017 % |
| ChargeScale | 0.034 +- 0.005 % | 0.044 +- 0.006 % |
| ChargeScale15 | 0.041 +- 0.006 % | 0.048 +- 0.006 % |

DENT breaks about 115 to 140 times more crossers than the controls, and the two independent DENT pairs agree. Caveat: the controls pair
`sbnd_large` with a different sample than the CV used for DENT, so they bound the re-simulation noise only approximately.

## 2. Effect of breaks on the PID shift: opposite sign to the net shift

Class contributions to the mean shift (variation minus nominal) of the leading-muon muon softmax; contributions sum to the total.

| class | fraction | mean shift in class | contribution |
|---|---|---|---|
| break | 1.29 % | -0.070 | -0.0009 +- 0.00005 |
| gap (5 cm / 3 cm) | 0.27 % | -0.002 | -0.000005 +- 0.000014 |
| intact stitched | 25.6 % | +0.0052 | +0.0013 +- 0.0001 |
| non-crosser | 72.8 % | +0.0036 | +0.0026 +- 0.0002 |
| total | | | +0.0031 +- 0.0003 |

The replication gives the same pattern (break -0.0010, total +0.0021). Broken tracks have a lower muon score, so breaking pulls the
average down; it cannot be the source of the net increase. The same sign holds for primary softmax (break -0.0003 of a +0.0005 total)
and MIP softmax (-0.0001 of +0.0034). Most of the net muon-softmax shift is carried by non-crossers (+0.0026 of +0.0031), which have no
cathode dependence by construction; this class contains true crossers reconstructed as non-crossers in both samples (about 0.6 % of crossers).
The controls (+0.024, +0.028) show a much larger shift from a genuine calibration change, so the DENT shift is small in absolute terms.

## 3. The "gap" cause is not isolated by this definition

Gap = stitched in both, start and end points move by less than 5 cm, path length falls by more than 3 cm. Its contribution is
consistent with zero (above), but the class is not a stable population: the gap share of stitched crossers changes from 0.01 % to 17 %
across the tolerance scan (end point 2-12 cm, length loss 1-8 cm), and the controls are non-zero at tight tolerances (0.01-0.7 %).
`gap_signature.pdf` shows why: DENT broadens a continuous core of small length and end-point changes, with no separate cluster.
By the pre-stated rule no single gap number is quoted, and a path-length-only gap is not supported by these data.

A calorimetric signature does exist in stitched tracks. Stitched crossers have a mean change of the calo-minus-CSDA difference of
-0.032 (replication -0.034) against -0.002 for non-crossers, contributing -0.0083 of the -0.0098 total shift of that variable, while
their path length is essentially unchanged (-0.02 cm). Charge missing inside a joined track (a calorimetric rather than a length gap) is
consistent with this, but the sign convention of `reco_leading_muon_calo_csda_diff` was not checked here and the quantity was not
used to define a class, so this is a pointer, not a measurement.

## 4. Selection fate (added 2026-10-14)

Events are classed as before when reconstructed in both samples ("unreco" otherwise) and followed into the `selected` tree of each
sample (`selection_fate.csv`, `selection_fate.pdf`). Shares are of the events selected in either sample.

Selected-tree definition (before the crosser veto and end-point cut), DENT:

| class | selected in both | nominal only | DENT only |
|---|---|---|---|
| break | 78.3 % | 18.0 % | 3.8 % |
| intact stitched | 90.9 % | 5.2 % | 3.9 % |
| non-crosser | 89.2 % | 5.6 % | 5.3 % |
| ChargeScale non-crosser (control) | 94.8 % | 3.0 % | 2.2 % |

Broken tracks are lost from the selection about three times as often as intact ones (18 % against 5 %) and are rarely gained, the only class
with a strongly asymmetric churn. The replication gives the same shares (77.3 / 17.6 / 5.1 %).

With the analysis crosser veto and the 4 cm end-point cut, the stitched classes cannot be selected in either sample, while the
broken tracks are selected in DENT only: 746 events (313 in the replication), because the CV crosser veto removes them and DENT's broken
track escapes it. Events in the non-crosser class churn in both directions (9,117 nominal only, 8,065 DENT only), so the
non-crosser net change (-1,052) is larger than the gain from breaks (+746); the overall selected count is about 0.2 % lower in DENT.
Events with no reconstructed muon in one sample (unreco) churn symmetrically (2,433 against 2,437).
This shows that, under the analysis selection, track breaking is a direct source of events that enter the DENT selection and are absent from the CV.

## 5. Dependence of the break rate (added 2026-10-14)

Break rate among true crossers that are reco crossers in the CV, in bins of nominal-sample variables (`break_rate_dependence.csv`, `.pdf`):

| |start direction x| | DENT | replication | ChargeScale |
|---|---|---|---|
| 0.0-0.2 (parallel to the cathode plane) | 27.7 % | 26.1 % | 0.10 % |
| 0.2-0.4 | 8.0 % | 7.9 % | 0.04 % |
| 0.4-0.6 | 4.1 % | 4.0 % | 0.03 % |
| 0.6-0.8 | 3.2 % | 3.2 % | 0.03 % |
| 0.8-0.9 | 2.5 % | 2.5 % | 0.02 % |
| 0.9-1.0 (along the drift) | 3.0 % | 2.9 % | 0.04 % |

Against |reco vertex x| the rate is 11.4 % at 5-10 cm, falls to about 3.4 % at 40-70 cm and rises to 5.6 % at 120-200 cm (replication 9.9, 3.7 and 5.3 %);
the 0-5 cm bin is empty. The break rate is thus dominated by tracks nearly parallel to the cathode plane, which spend a long path close to it,
and by vertices close to the cathode. A plausible reading (a hypothesis, not tested here) is that DENT displaces charge near the cathode
so that the two halves of a shallow track no longer meet the stitching criteria. The two independent DENT pairs agree bin by bin, and the controls stay two to three orders of
magnitude lower in every bin.

## Limits
* The PID decomposition covers events reconstructed in both samples; selection churn is treated separately in section 4 and is not part of the PID decomposition.
* Classes use reco crossing flags of matched events; they do not isolate where along the track a break occurs.
* Signal trees with the vetoes above, not the full selection including PID thresholds, and no CV weights.
* The decomposition shows association, not mechanism: breaking and PID changes are both consequences of the varied charge response.

## Suggested follow-ups
* Define a gap by the calorimetric deficit (per-hit dE/dx near the cathode) rather than path length; needs per-hit information not in these ntuples.
* Test the shallow-angle hypothesis with the true track angle and true displacement, and check whether a stitching-criterion change recovers the broken tracks.
