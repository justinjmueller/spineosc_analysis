# The SBND APA-gap deficit in the detector data access

Written 2026-10-17. Reproduce with `plotting/studies/apagap.py` (figures and `gap_summary.json` in this directory).

## Scope and rules

* Data: `spineosc_numu_sbnd_release-v1.1.0_detector_access_DONOTLOOK.root`, tree `events/sbnd_bnblight_full/selected`, read through the
  detector-data-access whitelist only (vertex, muon and proton start and end points, calorimetric difference, PID and primary scores, selection
  flags, and `Run` as a label). No energy, length, polar angle or direction component; the script refuses any other branch. Two-dimensional
  views combine two allowed coordinates and never form a difference of them. The data file has no azimuth branch.
* Simulation: `sbnd_large` `selected` tree of the `cvrw` file with `cvwgt` and the flash-PE trigger emulation, neutrino categories only (the
  selected tree holds no cosmic events), the same selection as `xml/detector_data_access/sbnd_run1/*_inclusive.xml`.
* **Area-normalised only, with no sense of the size of the data.** The simulation is scaled to the data of the selection and nothing else is
  normalised or fitted. This directory, the script's output and the provenance sidecars carry no event count, total, rate or file size of the
  data: histograms are fractions of the selected events, the tables hold ratios and pulls, and the data file is recorded by hash only. Statistical
  error bars and pulls are shown because the comparison needs them; like the error bars of any area-normalised plot they reflect how populated a
  bin is.
* Pulls are `(data - MC) / sqrt(MC + MC statistical variance)`: data Poisson plus simulation statistics, **no systematic uncertainty**, except where
  the PROfit systematic band is shown for scale.

## What is seen

| Quantity | Data/MC | Pull |
|---|---|---|
| vertex z in [251.5, 254.5) cm (the gap window) | 0.73 | -4.2 |
| vertex z in [254.5, 256.5) cm (just above) | 1.20 | +2.7 |
| vertex z in [249.5, 251.5) cm (just below) | 0.98 | -0.3 |

* The gap window is the most deficient 3 cm window anywhere in z = 226-280 cm. With 18 independent windows in that range the trial-corrected
  probability is about 2e-4 (statistics only), against 1e-5 locally. The existing PROfit comparison in 4 cm bins sees it as one bin at 0.815 of the
  prediction (250-254 cm) next to 1.075 in 254-258 cm, with a systematic band of only about 2%, so the deficit is far outside the detector-systematic
  band; its shape chi-square of 24.5/19 hides it because the dip is diluted over the bin and spread over 20 bins.
* In 0.5 cm bins (`gap_vertex_z`) the data/MC ratio falls from 0.83 at 251.5 cm to 0.52 at 253.5 cm and recovers to 1.2 above 254.5 cm, with a
  spike of 1.53 at 255.5 cm. The simulation has the same structure to a lesser degree (a dip at 253-254 cm and a spike at 255 cm), so part of the
  feature is reproduced and the data have more of it.
* The vertex, the muon start and the proton start show the same dip in the 2 cm bin at 252-254 cm (data/MC 0.68, 0.70, 0.69; pulls -4.2, -3.9,
  -4.0). They are largely the same coordinate, so this is one observation, not three. The muon end shows no dip (0.98); the proton end shows an
  **excess** in that bin (1.25, +3.4 sigma) (`gap_coordinates`).

## Correlations (all post hoc, statistics only; `gap_splits`, `gap_proton_end_vs_vertex`)

* **Run:** the deficit is in every period (data/MC 0.92, 0.62, 0.71, 0.62 over four Run ranges, each +-0.1): not one bad run.
* **x:** 0.58-0.61 for x = -120 to 60 cm and 0.81-0.92 at x < -120 cm and x > 60 cm (0.84, 0.92 on the +x side beyond 60 cm); **y:** about 0.55-0.72
  except 0.9 at y = 0-50 and 100-170 cm. Each slice has +-0.15 statistical error.
* **PID:** concentrated in events with a saturated leading-proton score (softmax >= 0.99: data/MC 0.58, -5.2 sigma; the other proton score slices
  are consistent with 1, 0.91-1.09). By primary score the deficit is in the 0.7-1.0 slices (0.62, 0.69) and not in 0.3-0.7 (0.98). The muon score is
  low in every slice (0.33, 0.73, 0.78).
* **Where the events go:** for vertices at 252-254 cm with the proton ending at 258-264 cm (across the gap, toward +z) the ratio is 0.45 (-2.7
  sigma); for vertices at 254-256 cm with that same proton end it is 1.48 (+2.6 sigma) (`gap_proton_end_vs_vertex`). Two cells beyond 2.5 sigma
  among about 40 is not strong evidence on its own, but it is the pattern a vertex pushed from inside the gap to its +z edge would give.

## Reading

The data are consistent with events whose vertex lies in the dead region around z = 253 cm being **lost or displaced to the edges of the gap**,
mostly the +z edge, and the simulation under-represents that. The deficit in the window is larger than the excess just above it, so part of the
loss is real inefficiency and part is migration. The PID dependence (clean, saturated proton scores disappearing) fits events that lose hits at the
gap. None of this is an energy effect, and none of it is established as the cause: the data have no truth.

## How long is the gap? (`plotting/studies/apagap_length.py`)

Two probes, same data and rules (area-normalised, allowed variables, no counts; all lengths in cm of z).

**1. Equivalent length (no model).** `integral (1 - data/MC) dz` over a window is the stretch of z the missing events would occupy if the
detector were uniformly efficient: **0.79 +- 0.16 cm** for the deficit side (251.5-254.5 cm), **0.43 +- 0.17 cm** for the pile-up above it
(254.5-256.5 cm), and **0.36 +- 0.23 cm** net. So most of the displaced stretch is compensated just above the gap and only a fraction is lost.

**2. Box scan with the simulation's own reconstruction.** The simulation's true vertex is smooth across the gap and its reconstruction
resolves it to a robust 0.33 cm, so a no-gap pseudo-reconstruction (true vertex plus the residuals measured away from the gap) is a clean reference.
Events whose *true* vertex lies in a box `[a, b)` are lost, or displaced to the +z edge of the box, and the reconstructed vertex histogram is
compared with the data in 0.5 cm bins (Poisson likelihood, prediction area-normalised, `length_scan`, `length_summary.json`).

| Model | a [cm] | b [cm] | Length L [cm] | -2 ln L better than no gap | Fractions |
|---|---|---|---|---|---|
| displaced to the +z edge | 251.75 | 255.5 | **3.75** | 62 | 0.3 displaced |
| lost | 252.5 | 254.0 | 1.5 | 27 | 0.5 lost |
| lost and displaced | 251.75 | 255.5 | 3.75 | 62 | all displaced, none lost |

* **Displacement beats loss by 35 units of -2 ln L**, and giving the model both leaves loss at zero: the gap mainly moves vertices (to the +z
  edge) and does not just remove them. Equivalently, about 30% of the events with a true vertex in the 3.75 cm stretch are displaced.
* **Uncertainty:** Poisson toys of the best-fit model refit the length with a spread of 0.54 cm (0.72 cm after inflating by the 1.78 per bin
  deviance scatter that the rest of the range also shows, from template statistics and broad shape differences); the fitted length is 3.5 +- 0.5
  on average for toys of a true 3.75 cm. **Length 3.75 cm, about 0.7 cm uncertain (statistics only).** The profile alone gives 3.5-4.0 cm; the toys
  are the honest number. The result does not move with the fit range (240-268, 236-272 cm) or the residual pool (margin 5 or 12 cm), but the scan
  grid is 0.25 cm and the edges are quantised to it.
* **Goodness:** inside 248-259 cm the per-bin deviance drops from 4.4 (no gap) to 1.6 (best model), the same as the 1.8 per bin outside; the
  model describes the gap region as well as the rest of the range is described.
* **Closure on the simulation itself:** the same scan applied to the simulation's reconstructed vertices (it has its own built-in gap) gives
  `a = 254.5, b = 255.25`, L = 0.75 cm with 0.9 displaced, but toys of that case spread widely (1.8 cm) because the effect is small. The data's gap is
  about **3 cm longer than the simulation's**, spread over a region where only part of the events move.
* **Independent check with end points:** the proton and muon end z (`length_edges`) pile up inside the gap (ratio about 1.3-1.8 at 253.5-254.5 cm) and
  fall sharply from 255.5 cm (about 0.4-0.5 at 255.5-256.0 cm): tracks that would end beyond the +z edge are cut short, which puts the +z edge at
  255.5 cm, the same place as the vertex fit.

**Reading.** The affected stretch is about 3.75 cm long (251.75 to 255.5 cm, uncertainty about 0.7 cm), against about 0.75 cm in the simulation,
and it acts by displacing vertices to the +z edge more than by removing them. The 3.75 cm is the extent of the region where the data reconstruction is
affected; it is not by itself the physical size of a dead region. Caveats: statistics only; the model assumes the displacement goes to the +z edge
and uses the simulation's resolution; edges are on a 0.25 cm grid.

## Suggested next steps

1. **Use the simulation's truth** (not restricted for MC): how truth vertices in 250-256 cm are reconstructed, to separate loss from migration in the
   simulation and see whether its dip and spike are the same mechanism at a smaller size.
2. Compare the dip position and width with the APA-frame geometry to see whether the effective dead region is wider or shifted in the data.
3. If the effect is to be handled in the analysis: an explicit gap variation (a wider effective dead region, or a z veto over 251.5-256.5 cm, which
   removes a small fraction of the selected events).

## Caveats

Statistical significance only (no systematic variations other than the 2% PROfit band shown for scale); the splits are post hoc and their errors are
large; the simulation's weights are `cvwgt` with no further reweighting; the dirt tree is empty in these inputs, as elsewhere in the detector data
access.
