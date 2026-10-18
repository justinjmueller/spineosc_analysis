# Reconstructed-energy binning and the sensitivity contours

Written 2026-10-18. Asimov (fake) data only, so none of the data-access restrictions apply. Scripts: `scripts/amr_binning/`; tables
`binning_incl1D.csv`, `binning_excl1D.csv`; figures `binning_incl1D.png`, `binning_excl1D.png`.

## What was done

For the two one-dimensional configurations of `runs/amr_compare` (1D inclusive: one `numu` channel, 17 bins, 0.3-2.0 GeV; 1D exclusive:
`ccqelike` 16 bins 0.3-1.5 GeV and `nproton` 14 bins 0.4-2.0 GeV), the reconstructed-energy binning of each channel was changed and the
adaptive-mesh-refinement (AMR) surface scan to **delta-chi2 = 1.6** was rerun, with the same options as the baseline (`--surface-amr --amr-initial 12
--amr-levels 3 --amr-levels-chi2 1.6`, seed 405, Asimov data, about 1,160-1,200 fits each). Only the edges of the fit variable change; samples, selection and
systematics are untouched (the XML variants differ from the baseline in the `edges` attribute alone).

Binnings (per channel, keeping the channel's range):

| Name | Definition | Inclusive bins | Exclusive bins |
|---|---|---|---|
| base | the analysis binning (the existing `amr3`/`amr2` runs) | 17 | 16 + 14 |
| merge2 / merge4 | every second / fourth baseline edge | 9 / 5 | 8+7 / 4+4 |
| uniform 0.10 | uniform 0.10 GeV bins | 17 | 12 + 16 |
| uniform 0.05 | uniform 0.05 GeV bins | 34 | 24 + 32 |
| equalpop N | equal expected population per bin, from a fine-binned (0.025 GeV) population of the Asimov prediction summed over detectors | 5, 9, 17 | 8+8, 15+15 |

**Not feasible:** uniform 0.025 GeV bins. PROfit aborts in its empty-bin check (the reconstructed energy only starts at about 0.33 GeV and the top bins
hold no events), so the fit would be singular. The fine population itself was obtained with `--force` and is only used to place equal-population edges.

**Detector systematics:** PROfit looks up the detector-systematic ratio histograms per event by its energy value, so they keep their native 17-bin
resolution in every variant; finer fit bins see the same piecewise-constant ratios, coarser ones average them.

**Not done:** the two-dimensional configurations (`excl2D`, `incl2D`, the base config), whose baseline scans take 18 to 90 minutes each.

## Results

Limit = sin^2 2theta at delta-chi2 = 1.6, **relative to the baseline binning** (above 1 is a weaker limit), at fixed dmsq. Contour areas are in log10 units
(the scan covers 12). Ratios at dmsq above 3 eV^2 have a few-per-cent scatter from the adaptive contour's resolution.

### 1D inclusive

| Binning | bins | Area | 0.3 eV^2 | 1.5 | 3 | 10 | 30 |
|---|---|---|---|---|---|---|---|
| base | 17 | 3.018 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| uniform 0.05 | 34 | 3.038 | 1.03 | 1.02 | 1.00 | 0.97 | **0.87** |
| uniform 0.10 | 17 | 2.968 | 1.07 | 1.07 | 1.05 | 1.04 | 0.93 |
| equalpop 17 | 17 | 2.940 | 1.07 | 1.07 | 1.13 | 1.04 | 1.17 |
| merge2 | 9 | 2.904 | 1.07 | 1.09 | 1.07 | 1.10 | 1.11 |
| equalpop 9 | 9 | 2.840 | 1.11 | 1.15 | 1.21 | 1.19 | 1.28 |
| merge4 | 5 | 2.735 | 1.21 | 1.23 | 1.15 | 1.23 | 1.26 |
| equalpop 5 | 5 | 2.695 | 1.19 | 1.21 | 1.41 | 1.31 | 1.46 |

### 1D exclusive

| Binning | bins | Area | 0.3 eV^2 | 1.5 | 3 | 10 | 30 |
|---|---|---|---|---|---|---|---|
| base | 16+14 | 3.050 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| uniform 0.05 | 24+32 | 3.076 | 1.02 | 1.00 | 0.98 | 0.97 | **0.92** |
| uniform 0.10 | 12+16 | 2.999 | 1.05 | 1.05 | 1.04 | 1.05 | 1.00 |
| equalpop 15 | 15+15 | 2.984 | 1.06 | 1.05 | 1.07 | 1.06 | 1.11 |
| merge2 | 8+7 | 2.934 | 1.08 | 1.09 | 1.05 | 1.11 | 1.15 |
| equalpop 8 | 8+8 | 2.887 | 1.08 | 1.15 | 1.15 | 1.19 | 1.25 |
| merge4 | 4+4 | 2.781 | 1.16 | 1.17 | 1.15 | 1.30 | 1.33 |

## Findings

1. **The binning matters at the 5-30% level, and more bins help, with diminishing returns.** Merging the baseline into half as many bins weakens the limit by
   5-15%, and into a quarter by 15-33%. Doubling the resolution (uniform 0.05 GeV) leaves the limit near the minimum of the contour (dmsq about 1.5 eV^2)
   unchanged within 2% and **improves the high-dmsq limit by 3% to 13%** (dmsq 10-30 eV^2); the excluded area grows by about 1%.
2. **The baseline binning is already close to the best practical choice** at low dmsq: no variant beats it below 1.5 eV^2 beyond the 2% level.
3. **Equal-population ("adaptive") bins are worse than simply merging the baseline at the same bin count** (e.g. 9 bins: +11% to +28% against +7% to +11%;
   17 bins: up to +17%). They put a single wide bin over the whole high-energy tail (above about 1.0-1.2 GeV in the inclusive channel), which discards the
   high-energy shape information that carries the high-dmsq sensitivity. Equal statistical weight per bin is not equal oscillation information; the baseline
   keeps fine bins in the tail, which is where its advantage at high dmsq comes from.
4. **The cost of finer bins is large.** The scan time at the baseline is 121 s (inclusive) and 489 s (exclusive), and about 20-35 times more at 0.05 GeV
   (2,441 s and 16,842 s; both runs shared the machine with other jobs, so these are upper bounds). The process time is unchanged at 7-11 minutes.
   Uniform 0.025 GeV is impossible with the available statistics.

## Caveats

* Asimov sensitivity only; no fake-data realisations, and no Feldman-Cousins calibration of the delta-chi2 = 1.6 level (still deferred).
* Per-dmsq ratios above 3 eV^2 carry a few-per-cent scatter from the AMR contour resolution (the right-hand panels show it); the trends quoted are larger than that.
* Timings of runs that overlapped with other PROfit jobs are not comparable (`NOTE.txt` in those run directories); the 0.05 GeV scans overlapped, the
  others ran alone.
* The equal-population edges come from the Asimov population summed over the three detectors; a different weighting (for example per detector, or weighted by
  sensitivity) could move them. They are an illustration of a population-driven binning, not an optimised one.

## Suggested next steps

* A binning chosen for information, not population: uniform in log(E) or 1/E (constant oscillation phase resolution), or edges optimised against the delta-chi2 = 1.6
  limit directly.
* Nonuniform refinements of the baseline that keep its low-energy bins and add resolution only in the tail (1.0-2.0 GeV), where the 0.05 GeV gain appears.
* The same scans for the two-dimensional configurations.
