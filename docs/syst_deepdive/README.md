# Per-systematic deep dive, first pass (SBND, final selection)
Script `plotting/dent_study/syst_deepdive.py`; figure `deepdive_sbnd.png`; numbers `deepdive.json`.
Events common to nominal and variation; "lost" = selected in the nominal only, "gained" = selected in the variation only, "migration" = selected in both but in a different energy bin.

| Variation | net (var/nom - 1) | lost | gained | signal share of lost (nominal 88.5 %) |
|---|---|---|---|---|
| Charge scale | -0.64 % | 6.3 % | 5.7 % | 0.60 |
| WireMod X-theta_xW | -0.13 % | 7.6 % | 7.4 % | 0.61 (gained: truth not filled in this sample, ignore) |
| Non-uniformity | -0.26 % | 2.1 % | 1.9 % | 0.62 |
| DENT | -0.85 % | 13.1 % | 12.2 % | 0.65 |

Observations
* The net ratio is a small imbalance between two much larger, nearly cancelling flows. Every variation re-selects 2-13 % of the events in each direction.
* The flows sit at the selection boundary: the lost fraction rises to 20-40 % for events with a muon score of 0.2-0.4 (the edge of the 2x-boosted muon region) and for short muons (< 50 cm), and falls to a few % for confident muons.
  Lost events are background-enriched (signal share 0.60-0.65 vs 0.885).
* Charge scale is the only variation with a systematic tilt in the score: lost > gained at high muon score, gained > lost at low muon score (it moves candidates toward the muon vertex, see the ternary migration plot).
  WireMod and DENT are nearly symmetric at every score: their net effect is the small residual of the churn, and is more sensitive to statistical noise in the churn than to a physical shift.
* Non-uniformity has the smallest churn but a clear asymmetry (low muon score lost twice as often as gained), which produces the 4 % net gain in the highest 1p bin.
* The highest-energy bin of every variation has a large energy-migration component (overflow past 2 GeV).

Open questions / next
* A nominal-vs-nominal null for the churn cannot be built from the files (no event overlap between campaigns). The churn size for a re-simulation with zero physics change is therefore unknown;
  if the churn is largely re-reconstruction noise, the one-sided ratio systematics partly measure noise. Needs a production-side check (same-seed rerun of the CV) or a split of the churn by true class.
* Trace each churn population to the failing cut (needs the idea-1 intermediate table; the `selected` trees only contain selected events, so a per-cut failure reason cannot be read for the unselected side).
