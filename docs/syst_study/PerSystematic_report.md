# Per-systematic investigation of the detector variations (SCE, WireMod, MattMod, DENT, Trigger)

Written 2026-10-10. Scripts: `plotting/syst_study/ps_*.py` (core analysis `ps_core.py`, cross-variation tests `ps_pairs.py`, fragility `ps_fragile.py`, trigger `ps_trigger.py`, duplicate fix `ps_fixdup.py`, tables `ps_summary.py`, figures `ps_plots.py`).
Numbers: `docs/syst_study/persyst/summary.json` (+ `core/<pair>.json`, `pairs_compare.json`, `trigger.json`, `fragile.json`, `fixdup.json`, `duplicates.json`). Figures: `docs/syst_study/persyst/*.png`.
cvwgt is not used anywhere. Charge scale / gain are not analysed here (used only as a comparison where stated).

## 0. Method, definitions and caveats (read first)

* Inputs: the paired event tables of `runs/syst_study/pairs/` (events common to nominal and varied sample, selected by the calculator selection on each side). `lost` = selected in the nominal only, `new` (gained) = selected in the varied only, `kept` = selected in both, `migration` = kept events that change energy bin. Ratio-1 = lost + new + migration per bin (checked: sum reproduces the ratio). Pair-table ratios reproduce the production histograms to <= 0.2 % in every bin (differences = duplicated keys, see Sec. 9; DENT 0.7 % in the last bin); the only real exception is ICARUS Run 4 WireModxThetaXW (Sec. 9, bug 2).
* Errors: paired Poisson bootstrap over events (200-300 replicas); this captures event sampling only, not the sensitivity of the reconstruction to a different random seed (no same-seed re-run exists, see Sec. 1).
* Species margins (approximate, my reconstruction of the argmax with 2x muon boost, other species such as electron/photon ignored): muon candidate `2 s_mu - max(s_pi, s_p)`, proton candidate `s_p - max(2 s_mu, s_pi)`. All selected events have margin > 0 on both, as expected, which supports the definition.
* Truth classes (from `true_cc, true_pdg, true_category`): signal-like numuCC (category 0/2), other numuCC, NC, nueCC, cosmic/other. For SBND WireModxThetaXW and WireModYZ the varied trees have no truth; I recovered it for 70 % of the `new` events from other variations' tables (same generated events), 30 % remain unknown, so class-resolved `new` rates for those two are not reliable and are not used.
* Channels: 1p / Mp defined by `reco_proton_multiplicity` on each side separately. Mp has no events in the first bin (0.3-0.4 GeV), so channel comparisons use 15 bins with >= 300 events per channel.
* "Floor": the lost fraction of events far from any species boundary (top five deciles of the margins).
* Production selection = tier 1 (event in the `selected` tree of one sample only); post cuts = tier 2 (fiducial, veto, crosser, end-x). From `s3_mechanism.json`: for every variation 95-99.6 % of lost/new events are tier 1 (DENT 89 %).

## 1. Cross-cutting findings that explain all variations

Table 1 (all numbers from `summary.json`; stat. errors on net are 0.03-0.12 %, bin errors 0.25-1.2 % except last bin 1.5-6 %).

| pair | N nominal | lost % | new % | net % | net/sqrt(L+G) | per-bin rms(r-1) % | chi2 vs unity (17 bins) | slope %/GeV | E_vis unchanged |
|---|---|---|---|---|---|---|---|---|---|
| SBND SCE0x | 346765 | 8.48 | 7.13 | -1.35 | -20 | 1.73 | 198 | +0.5 +- 0.5 | 0.00 |
| SBND SCE2x | 344192 | 6.94 | 8.47 | +1.53 | +23 | 1.60 | 222 | +0.1 +- 0.5 | 0.00 |
| SBND WireModxThetaXW | 282840 | 7.56 | 7.43 | -0.13 | -1.8 | 0.49 | **7** | +0.3 +- 0.5 | 0.03 |
| SBND WireModYZ (not used) | 323433 | 10.52 | 8.33 | -2.20 | -29 | 3.22 | 717 | **+7.7 +- 0.5** | 0.01 |
| SBND MattModYZ (used) | 335816 | 2.11 | 1.84 | -0.27 | -7.7 | 1.13 | 123 | +1.8 +- 0.2 | 0.92 |
| SBND DENT | 170423 | 13.00 | 12.12 | -0.87 | -7.2 | 1.05 | 44 | +2.4 +- 0.8 | 0.00 |
| ICARUS R2 SCE0x | 102238 | 6.51 | 5.92 | -0.58 | -5.3 | 1.09 | 40 | +0.9 +- 0.5 | 0.02 |
| ICARUS R2 SCE2x | 89070 | 5.31 | 5.89 | +0.59 | +5.2 | 1.00 | 39 | +0.6 +- 0.5 | 0.02 |
| ICARUS R2 WireModxThetaXW | 114800 | 6.64 | 7.13 | +0.49 | +4.5 | 0.94 | 29 | +1.3 +- 0.5 | 0.04 |
| ICARUS R2 WireModYZ | 113377 | 6.39 | 6.77 | +0.38 | +3.6 | 0.73 | 23 | -0.2 +- 0.5 | 0.03 |
| ICARUS R4 SCE0x | 162941 | 6.70 | 6.05 | -0.64 | -7.3 | 0.98 | 65 | +1.1 +- 0.4 | 0.02 |
| ICARUS R4 SCE2x | 164920 | 5.55 | 6.05 | +0.51 | +6.1 | 0.61 | 24 | +0.0 +- 0.4 | 0.02 |
| ICARUS R4 WireModxThetaXW (de-duplicated) | 147622 | 6.25 | 6.77 | +0.52 | +5.5 | 1.03 | 41 | +0.9 +- 0.4 | 0.05 |
| ICARUS R4 WireModYZ | 196008 | 6.59 | 6.46 | -0.13 | -1.6 | 0.40 | 11 | -0.0 +- 0.4 | 0.03 |

**F1. The net effect is 1-12 % of the churn.** Each variation re-selects 2-13 % of the events in each direction (lost + new 4-25 %), and the net is |net|/(lost+new) = 0.01-0.12. Chi2 of the applied ratio against unity (17 bins) is only 7-65 for 12 of 14 pairs, i.e. several variations (SBND WireModxThetaXW 7, ICARUS R4 WireModYZ 11, ICARUS R2 WireModYZ 23, ICARUS R4 SCE2x 24) are statistically compatible with no effect given the bootstrap errors. Only SBND SCE0x/2x (198/222), SBND WireModYZ (717), SBND MattModYZ (123) and DENT/ICARUS R4 SCE0x (44/65) show a significant systematic pattern.

**F2. The churn is dominated by a small population of "fragile" events, not by a smooth shift.** For nominal-selected events I counted in how many variations (sharing the SBND nominal: 9 variations; ICARUS: 5 per run) the event is dropped (`fragile.png`, `fragile.json`). The observed distribution is far wider than independent draws at the same mean loss rate: for SBND (9 variations) 9.1 % of events are dropped in >= 3 variations versus 2.1 % for independent binomial; for the 5 non-charge SBND variations 4.0 % versus 0.32 %; ICARUS R2/R4 3.3 %/3.5 % versus 0.18 %/0.20 %. 38-40 % of all loss incidents (61 % in the 9-variation count) come from this ~4 % population. These events are background-enriched (lost fractions by truth class: other-numuCC 19-30 % (DENT 41 %, MattModYZ 7 %), NC 14.5-19 % (DENT 30 %), signal 3.8-6 % (DENT 9.5 %, MattModYZ 1.5 %); signal share of lost events 0.61-0.68 versus 0.88-0.91 in the nominal) and sit near the species decision boundaries. Consequence: the per-bin ratio is the small imbalance between two large flows of mostly the same fragile events, so its noise is set by the fragile population, and its sign/shape is the residual of the physical response on top of it.

**F3. Where the flips happen (all variations, both detectors).**
* Standardized mean difference (lost vs kept, nominal features) is negative and large for the leading-proton margin (-0.55 to -0.83) and the muon margin (-0.40 to -1.05), moderately for proton calo-vs-range difference (-0.15 to -0.42), muon primary score (-0.33 to -0.56); `new` events show the same signs on the varied features (table in `summary.json`, heat map `overview.png`). The same pattern for every variation and both detectors means the flips are a property of the selection boundary, not of the individual physical effect.
* Lost fraction versus margin rises from a floor of 3-9 % (confident events; MattModYZ 1-1.5 %) to 17-36 % (proton margin) and 12.5-35 % (muon margin) in the lowest bin (`*_mech.png`; MattModYZ 6 %/10 %). The share of all lost events above the floor is 32-57 % for the leading-proton margin, 25-58 % for the muon margin (overlapping sets), and 9-21 % for the proton-KE bin edge near the 50 MeV threshold (the 50-60 MeV bin has 5.2-5.7 % of events and loses 14 % (ICARUS), 23 % (SBND SCE), 38 % (SBND DENT) versus floor 5-11 %). A muon-KE threshold (143 MeV) effect is weak (loss 10-18 % versus floor 6-12 %).
* A floor of 3-9 % lost (SBND 4-9 %, ICARUS 3-6 %; MattModYZ 1-1.5 %) exists even for events with margins near the maximum; it is 45-70 % of the total loss. I cannot attribute it with the available data: the failing side has no features (only selected events are stored). Candidates: sub-leading particles (a 25 MeV pion/photon/electron fragment appearing, second muon, proton-count fluctuations: more protons -> more loss, 8 % at 1p to 16 % at >= 5p) and clustering instabilities.
* Post cuts are minor except for DENT: SBND SCE 0.35/0.30 % (lost/new) of nominal, WireMod 0.25-0.45 %, ICARUS ~0.07-0.09 %, DENT 1.5/1.2 %.

**F4. Transform-type variations perturb every event; charge-type variations do not.** E_vis is a range-based quantity; 92 % of events keep exactly the same E_vis under MattModYZ (75 % under ChargeScale) but only 0-5 % do under SCE, DENT, WireModxThetaXW, WireModYZ (vertex jitter at the 3 mm voxel scale changes the reconstructed geometry of all events). The churn of the geometry-changing variations (6.3-13 %) is 3-6x the churn of MattModYZ (2.1 %), so a generic re-simulation noise of the production would be <~ 2 % (upper bound, MattMod itself is a perturbation). The remainder is sensitivity to input perturbation.

**What cannot be settled without new production:** whether the lost/new flows of 5-13 % are physics (a real, sign-definite detector response) or chaotic reaction of the network to any voxel-level perturbation. The odd/even and overlap tests below say the latter dominates the churn while the former is visible only in the net.

---

## 2. SBND SCE0x / SCE2x

Figures: `sbnd_sce_ratio.png`, `sbnd_sce_mech.png`, `sce_linear.png`, `sce_event_level.png`.
Not used in the SBND analysis (technote: DENT is used for the SBND E-field), computed for validation. Physical definition: reconstructed coordinates are displaced with the SBND SCE map applied backward (0x) / forward (2x) before clustering; nominal = no map (1x).

**What moves (kept events, `core/sbnd_run1__SCE*.json`).** Geometry only: median vertex shift 0.90 cm (84th percentile 1.6 cm), 1.1-1.2 cm near the cathode falling to ~0 at the anode (0.3 cm voxel quantization). Muon length: +0.82 % (0x) / -0.66 % (2x) near the cathode, ~ -0.1 % / +0.1 % near the anode, odd in sign. Median E_vis shift +0.17 % (0x) / -0.15 % (2x) (mean +0.35 % / -0.10 %, heavy tails); muon KE +0.23 % / -0.19 %; calo-vs-range difference -0.23 / +0.23 %-pt (muon); muon score +0.5 / -0.3 %-pt; direction median 0.47 deg. Heavy tails: 3D vertex shift std 23 cm (a few percent of events have a different candidate vertex). Proton lengths: median 0 (88 % changed, spread 15 %).

**Where events flip.** Tier 1 99.6 % (the net post-cut change is -0.05 %). SMD (lost vs kept) p-margin -0.63 / -0.69, mu-margin -0.42 / -0.58; excess above floor: p-margin 36-42 %, mu-margin 26-37 %, p_ke edge 14-17 %. Floor 5.5-6.4 % (0x), 4.1-4.4 % (2x). Position: lost fraction 11.0 % in the first decile of the distance to the cathode (5-27 cm) falling to 5.8 % at 160-190 cm (0x; deciles of the drift distance); 11.0/7.6 % at the y edges and 11.8/10.3 % at the z edges versus 7.3 % centrally, i.e. higher where SCE displacements are larger and near the TPC boundaries; weak track-angle dependence (10.1 % at theta_xz<7 deg to 7.3 % along drift).

**Net effect and its position dependence (`sbnd_sce_mech.png`).** Net per drift-distance bin for 0x: -4.6 % (first 20 cm from cathode), -3.8, -2.2, -1.8, -1.4, -0.6, 0.0, +0.2, +0.6, +1.5 % (180-200 cm); 2x is the mirror image (+3.9 % ... -1.5 %). Also edge enhancement in y (-2.9 % at y<-160 cm, -2.4 % at y>160 cm versus -0.85 % central) and z (-1.75 % at z<50 cm). Integrated -1.35 +- 0.07 % / +1.53 +- 0.07 %. This position pattern follows the geometric displacement and length change, a clear physical signature.

**Shape in E_vis.** Per-bin ratio-1: 0x -2.4 % at 0.3-0.4 GeV to about -1 % over 0.5-1.2 GeV (-3.8 +- 3 % in the last bin), 2x +0.8 % ... +2.5 %; the weighted slope is consistent with 0 (+0.5 +- 0.5 %/GeV). Decomposition (lost, new, migration) per bin is in `sbnd_sce_ratio.png`: lost and new are U-shaped (15 % at the lowest bin, 5-6 % at 0.8-1.0 GeV, 11-17 % at the top); the migration term is the odd, energy-scale-like part (0x: -2.2 % in bin 1 rising to +3.9 % in bin 16, mirrored for 2x) and the lost-new imbalance gives the offset. Variance shares of r-1 across bins are not clean (lost/new/migration +1.7/-1.6/+0.9 because lost and new are anticorrelated in E).

**Linearity / oddness (the key test).**
* Joint paired bootstrap (shared nominal, restricted to the common pool): odd part (r2x-r0x)/2 = +1.49 % mean, flat 1.2-1.8 % over the bins, chi2 = 695/17; even part (r2x+r0x)/2 = -0.01 % mean, chi2 = 31/17 (mildly high, likely correlated bin-to-bin structure). Integrated: odd 1.44 %, even 0.09 %. So the net ratio is odd to < 10 % (ICARUS: odd 0.58 %, even 0.00/-0.07 %).
* The event-level churn is **even**: 41 % of events lost in 0x are also lost in 2x (versus 6.9 % base rate, 5.9x); 47 % of 0x-new events are 2x-new (lower bound). Hamming distance of the selected sets: nominal-0x 15.7 %, nominal-2x 15.5 %, 0x-2x 17.5 % (a purely odd physical response would give 0x-2x = 2x nominal-0x; a purely random one gives equal distances). If the distances add linearly (d(n,x) = N + P, d(0x,2x) = N + 2P, N = sign-independent part, P = coherent part), the coherent part is P/(N+P) = d02/d_mean - 1 = 12 % of the nominal-to-varied selection difference (ICARUS 17-19 %); ~80-90 % is sign-independent. This is a model assumption (additivity of sets), not a measurement.
* Event-level shifts: central values are odd (vertex dx Spearman between 0x and 2x shifts -0.66, end-x -0.72, muon length -0.28) but the Pearson correlations are +0.5 because of common outliers, and PID score shifts are positively correlated (muon score +0.50, proton score +0.54, calo +0.39; i.e. even). Interpretation: the positions respond odd-linearly to the SCE map, the PID network responds to *any* displacement with a sign-independent random reshuffling of borderline events.
* Class level: signal-like numuCC net -1.0 % (0x), -0.3 % (2x) (even part -0.6 %); other-numuCC -5.5 % / +18.6 %; NC +2.9 % / +0.5 %. SCE2x's +1.53 % is almost entirely a +18.6 % growth of the other-numuCC contamination (15069 new versus 8956 lost background events, flat in all margins and E_vis, i.e. not a boundary effect; mechanism unverified: plausible is a pion/photon veto that fires less often when tracks are shortened, see Sec. 10 for the test). The ratio applied to all events therefore overstates the signal-specific response of 2x (+1.5 % applied, -0.3 % in signal) and understates it for 0x (-1.35 % versus -1.0 %).

**1p vs Mp.** Integrated all/1p/Mp: 0x -1.35/-1.39/-1.25 %, 2x +1.53/+1.46/+1.73 %; Mp-1p difference +0.13 +- 0.25 % / +0.27 +- 0.23 % (not significant). rms of (r_ch - r_all) over 15 bins 0.38 % (1p), 1.44 % (Mp; noise-subtracted 0.63 %) versus rms systematic 1.43 %; chi2 11 / 15 bins. The channel mismatch is statistically insignificant, below the systematic itself for 1p and comparable only through noise in Mp.

**Learned / uncertain / next.** SCE in SBND produces a real, odd, position-dependent (cathode-weighted) net change of +-1.4 % built on a 7-8 % even churn. The ratio applied to both channels is adequate. Uncertain: (i) whether a part of the net is an artifact of "1x = no transform" versus the transformed 0x/2x path (an identity-map transform sample, i.e. the nominal reconstructed through the same code path with the displacement set to zero, would measure the pure re-voxelization churn); (ii) the +18.6 % other-numuCC growth in 2x needs the failing-side information (pion/photon/electron multiplicities and KE near the 25 MeV thresholds for the nominal-selected and variation-failed events).

---

## 3. SBND WireModxThetaXW

Figures: `sbnd_wmxt_ratio.png`, `sbnd_wmxt_mech.png`. Physical definition (technote): hit width/area correction as function of (x, theta_xw) applied to the TPC signal; the signal-shape coverage of the analysis.

**What moves.** Calorimetry weakly: calo-vs-range -0.40 %-pt (muon), -0.26 (proton) (a charge-scale-like sensitivity predicts 3.5 %/GeV per %-pt, i.e. -1.4 %/GeV, observed +0.3 +- 0.5: not a charge-scale effect), muon score -0.16 %-pt, proton score ~0, E_vis median 0 (mean +0.10 %), E_vis unchanged in only 3 % of events (all others change by small amounts), median vertex shift 0.24 cm (84 %: 0.72 cm), median direction change 0.12 deg. The calo shift versus muon theta_xz is non-monotonic and changes sign between 10-deg bins (+0.65, -0.2, -1.3, -1.2, -1.0, +0.68, -0.22, -1.0, -0.1 %-pt), consistent with piecewise (x, theta_xw) maps not aligned with theta_xz.

**Where events flip.** Lost 7.56 %, new 7.43 % (as large as SCE0x), signal share of lost 0.61 (nominal 0.88), tier 2 only 0.25/0.45 %. SMD p-margin -0.72/-0.72, mu-margin -0.62/-0.58; excess above floor p-margin 44 %, mu-margin 37 %, p_ke 16 %; floor 4.3-4.8 %. Lost fraction falls with theta_xz from 9.7 % to 6.0 %.

**Shape/magnitude.** Net -0.13 +- 0.07 %; per-bin r-1 has rms 0.49 %, equal to the bootstrap errors (0.3-0.6 %), chi2 vs unity = 7 for 17 bins; no slope (+0.3 +- 0.5 %/GeV). This is a variation whose measurable signal is indistinguishable from event-sampling noise: the 7.5 % churn is balanced to 0.13 %. The spline from this ratio encodes bin-to-bin noise of ~0.5 %, not physics.

**1p vs Mp.** all/1p/Mp: -0.13/-0.05/-0.35 % (Mp-1p -0.30 +- 0.25 %); the rms channel difference (0.84 % / 1.48 %) is noise (noise-subtracted 0.36/0.39 %) and exceeds the (noise-level) systematic.

**Learned / uncertain / next.** In SBND this variation is statistically a null. Open: truth of `new` events is unavailable in the varied tree (30 % unrecoverable); a physical effect smaller than 0.3 % per bin is not resolvable with 2.8e5 selected events and 7.5 % churn; production of a second independent seed for the variation (or a larger sample) would show whether any real signal remains.

---

## 4. SBND WireModYZ (computed, not used) vs MattModYZ (used)

Figures: `sbnd_yz_ratio.png`, `sbnd_yz_compare.png`, `sbnd_yz_mech.png`.

| | WireModYZ | MattModYZ |
|---|---|---|
| net | -2.20 % | -0.27 % |
| lost / new | 10.5 / 8.3 % | 2.1 / 1.8 % |
| per-bin rms r-1 | 3.22 % | 1.13 % |
| ratio-1 low3 / mid / high3 bins | -3.7 / -1.0 / +5.2 % | -0.6 / +0.1 / +2.0 % |
| slope vs E_vis | +7.7 +- 0.5 %/GeV | +1.8 +- 0.2 %/GeV |
| events with identical E_vis | 1 % | 92 % |
| mean calo-vs-range shift (muon / proton) | +0.62 / +0.56 %-pt | +0.51 / +0.33 %-pt |
| mean muon-score shift | **-4.1 %-pt** | -1.1 %-pt |
| vertex median shift | 0.30 cm | 0.00 cm |

**Different physics, not a rescaled copy.** MattModYZ is a spatially dependent charge-response correction (data/MC maps) applied within the reconstruction: geometry is untouched (92 % identical E_vis, median vertex shift 0) and it acts purely as a small charge-scale shift (+0.5 %-pt in calo). Its energy trend obeys the same law as the SBND charge/gain variations (`sbnd_yz_compare.png`, bottom left): slope of ratio vs E_vis per %-pt of mean muon calo shift = 3.7 (GainHigh), 3.2 (ChargeScale), 3.9 (ChargeScale15), **3.6 (MattModYZ)** %/GeV; likewise score shift per %-pt: -2.5, -1.5, -1.8, -2.2. WireModYZ shifts the calorimetry by the same +0.6 %-pt but also modifies the signal shape (it changes the geometry of all events: 1 % identical E_vis), which drives the muon score ~3x further (-4.1 versus a charge-scale prediction of ~-1.4 %-pt), the churn 5x (10.5 versus 2.1 %) and the energy slope 3.5x above the charge-scale law (12.4 vs 3.6 %/GeV per %-pt). The two shifts of kept events are only weakly correlated event by event (Pearson 0.07-0.20; calo -0.10); 51 % of MattMod's lost events are also lost under WireModYZ (3369/6556) but only 10 % of WireModYZ's (3369/32537).

**Spatial structure.** In vertex (y,z) cells with > 300 events WireModYZ changes the count by up to -13.7 % (y in -200..-150, z < 50 cm), -9 % at z < 50, y~-120, and -8.3 % in a z~250-350 band at y~-50..0; MattModYZ varies between -1.7 % and +2.5 %. Cell-to-cell correlation of the two maps: 0.18; rms ratio 5.5. WireModYZ is concentrated in a few detector regions, which is why a spatially integrated net -2.2 % comes with a large churn.

**Energy shape (why both rise with E).** For WireModYZ the lost fraction falls with E (20 % at 0.3-0.4 GeV to 7 % above 0.8 GeV) while the new fraction rises above 1.0 GeV (7 % to 20 % in the last bin), net -3.7 % -> +7.9 %; migration is negligible (|mig| < 1 %). In MattModYZ lost 3.6 -> 1.5 %, new 3.1 % -> 5.4 % (net -0.5 % -> +3.7 %). The trend is the charge-scale response of the selection (high-E events are long muons whose track ends are driven by dE/dx PID); WireModYZ adds the extra score shift. The 1p and Mp channels have the same trend, Mp steeper (WireModYZ: Mp -11.6 %, -6.9 % ... +8.2 % versus 1p -3.7 % ... +4.4 %); channel difference Mp-1p integrated +0.53 +- 0.28 % (WireModYZ), +0.53 +- 0.11 % (MattModYZ: 1p -0.40 %, Mp +0.13 %: 4.8 sigma; here the 1p/Mp mismatch is as large as the combined effect, -0.27 %). rms of (r_ch - r_all): WireModYZ 0.89/2.37 % (noise-subtracted 0.50/1.96 %) versus systematic rms 2.59 %; MattModYZ 0.76/0.75 % (0.68/0.47 %) versus 0.71 %. Channel-difference chi2 35-36 (WireModYZ; ~3 sigma) and 25-26 (MattModYZ; ~2 sigma) for 15 bins.

**Double counting / consistency.** (i) The analysis allowlists only MattModYZ (WireModYZ is computed but not used; `xml/coverage/dp.xml` confirms), so there is no double counting of the two YZ variants. (ii) Partial overlap with the charge-type variations is plausible by construction: MattModYZ acts as a +0.5 %-pt charge shift, the same observable as ChargeScale (-1.8 %-pt) and Gain; 34 % of MattMod's lost events are also lost under ChargeScale (base 5.5 %); whether the MattMod map is independent of the data-driven ChargeScale calibration must be decided from how each is constrained, not from this study. (iii) Inconsistency across detectors: ICARUS uses WireModYZ (not MattMod), whose amplitude in ICARUS is a factor 4-8 smaller than SBND's (per-bin rms) (net +0.38 / -0.13 %, rms 0.73 / 0.40 %; mean calo shift +0.79 / -0.08 %-pt) while the SBND WireModYZ is the largest of the SBND geometry/shape-type effects (net -2.2 %, per-bin rms 3.2 %); the technote shows that the SBND WireModYZ overestimates the data/MC non-uniformity. (iv) WireModYZ and WireModxThetaXW lost sets overlap (38 % versus base 8 %), consistent with shared fragile events.

**Learned / uncertain / next.** MattModYZ is a controlled charge-scale-like perturbation whose energy trend (+3.7 % at 1.5-2.0 GeV, +1.8 +- 0.2 %/GeV) is the universal charge-scale response of the SBND selection and is a real feature, not noise (chi2 123/17). Not settled: whether the data/MC charge-map differences it encodes are already covered by ChargeScale (overlap). Next: a joint fit of ChargeScale + MattModYZ nuisances to a data-like pseudo-experiment, or check the correlation of the MattMod map with the ChargeScale shift per detector region.

---

## 5. SBND DENT

Figures: `sbnd_dent_ratio.png`, `sbnd_dent_mech.png`. Prior findings (`docs/dent_study/AUDIT.md`): end-point pile-up at the cathode is the robust DENT discriminator; crosser fraction depends on the generator model; the cuts (crosser veto + |x_end| >= 4 cm) act with opposite signs. Nominal is `sbnd_detvar_newcv` (different generation from the common CV), so net -0.87 % applies to the DENT/CV pair only.

**What is new here.**
* Churn in both TPC halves is the same: lost/new 13.5/12.9 % (vertex x<0) and 12.5/11.3 % (x>0); net -0.61 % versus -1.14 % (difference 0.5 %, ~2 sigma); median |dE_vis| 0.49 % on both sides and no event keeps E_vis unchanged (compare MattMod 92 %). The distortion described in the technote is localized to the east TPC field cage, yet the churn (13 %) is the largest of all variations and is not localized. Lost fraction versus distance to cathode: 17 % in the first 20 cm falling to 13 % and a flat 12-13 % beyond. So at least 12 % of the 13 % churn is not DENT-specific physics; it is the response of the selection to a pervasive geometry perturbation (the DENT transform is evidently applied to all events/TPCs or the CV and DENT samples differ otherwise). To be confirmed with production: does the DENT map act on the west TPC?
* Tier 2 (post cuts) is the largest of all variations: lost 2521 (1.48 % of nominal; veto 1291, end-x 731, fiducial 462, crosser 37), new 2073 (1.22 %; veto 1204, crosser 821, fiducial 27, end-x 21): the end-x cut removes the truncated tracks (731 lost), the crosser flag restores events (821 new), net post-cut -0.26 % (30 % of the -0.87 %); the rest, -0.61 %, is tier 1.
* Boundary dependence as elsewhere (p-margin SMD -0.62/-0.65, mu-margin -0.45/-0.50). Highest proton-threshold loss: the 50-60 MeV proton bin loses 38 % versus a floor of 11 %. Floor 8-12 %: the highest of all variations (reason unknown).
* Shape: -1.9 % (0.3-0.4 GeV) to ~0 above 0.9 GeV (+2.4 +- 0.8 %/GeV slope); E_vis median shift +0.11 % (mean +0.29 %), muon length +0.18 %. Channel: all/1p/Mp -0.87/-0.98/-0.56 % (Mp-1p +0.43 +- 0.40 %, not significant), chi2 of channel difference 6 for 14 bins.

**Learned / uncertain / next.** DENT shows the cleanest tier-2 mechanism (the end-x cut and crosser flag trade 2.5k lost and 2.1k new events), which supports the earlier audit, but its tier-1 churn is not localized. Missing: a test of the CV-versus-DENT sample difference in the west TPC (if the DENT field map is zero there, the observed 12.5 % churn is the processing noise floor for this pair, which is then the best available null for the other variations; ask production). The DENT result is a ratio of two samples whose events are not shared with the common CV, so its net cannot be combined with other SBND ratios event by event.

---

## 6. SBND TriggerEmulation

Figure: `trigger.png` (left column). Reproduced exactly (max |difference| = 0.0 for M1 and P1) from the nominal `sbnd_large` selected tree: an event passes if `0.647 x event_largest_flash_pe > 2000`, i.e. PE > 3091; M1/P1 shift the scale by -/+0.004, i.e. PE thresholds 3110/3072.

* Light yield: median flash PE per GeV of E_vis ~ 24.5k (rising to 25.9k in the first bin); PE median 9.7k at 0.3-0.4 GeV and 39k at 1.5-2 GeV; the 16th percentile at 0.3-0.4 GeV is 6.6k, twice the threshold. Only 0.18 % of selected events fail the nominal trigger (0.47 % in the first bin, < 0.1 % above 1 GeV) and < 0.01 % lie in the +-1 sigma window.
* Effect: |ratio - 1| <= 0.008 % in every bin; integrated +-0.001 %. Even a 15 % error on the PE scale (25 sigma) would change the selected count by only 0.02-0.04 %. Negligible by a wide margin, independent of E.

**Learned.** The SBND trigger systematic is irrelevant and can be dropped or fixed to zero; the splines from M1/P1 are exactly flat. Uncertain: nothing relevant (the quoted sigma 0.004 = 0.6 % cannot matter).

---

## 7. ICARUS Run 2 / Run 4 SCE0x / SCE2x

Figures: `icarus_sce_ratio.png`, `icarus_sce_mech.png`, `sce_linear.png`, `sce_event_level.png` (shared with Sec. 2).

**What moves.** Same observables as SBND, but with a smaller displacement: median vertex shift 0.34 cm (84 %: 0.6 cm), 0.3-0.4 cm in the first 100 cm of distance to cathode and ~0 beyond 120 cm (about one voxel; SBND 1.1-1.2 cm); muon length +0.53 % (0x) / -0.21 % (2x) near the cathode falling to ~0 at 100-150 cm; median E_vis +0.06 % / -0.05 %; calo shift -0.10 / +0.13 %-pt (muon); score shift < 0.1 %-pt. Run 2 and Run 4 are indistinguishable.

**Where events flip.** Lost 6.5/5.3 % (R2), 6.7/5.6 % (R4), new 5.9/5.9 %, 6.1/6.1 %; tier 2 0.07-0.09 %. SMD p-margin -0.66 to -0.74, mu-margin -0.40 to -0.49; the excess above floor is 45-52 % for the proton margin, 26-32 % for the muon margin, 13-16 % for the proton-KE threshold edge; floor 2.8-3.7 % (proton margin) and 3.7-5.1 % (muon margin). The 50-60 MeV proton bin loses 13-14.5 % versus floor 4.5-5.9 %. Class rates: signal lost 3.8-5.1 %, other-numuCC 19-22 %, NC 14.5-17 %; the 2x net comes from other-numuCC (+7.4 %/+8.3 %) with signal 0.0/-0.2 %, as in SBND.

**Net and shape.** Net -0.58 +- 0.12 % / +0.59 +- 0.11 % (R2), -0.65/+0.51 (R4); per-bin rms 0.6-1.1 % against bootstrap errors 0.5-1.2 % (chi2 vs unity: 40, 39, 65, 24 for 17 bins): the SCE0x are significant (R4 SCE0x 65), the SCE2x of Run 4 is compatible with noise. Slopes +0.9 +- 0.5 / +0.6 +- 0.5 (R2), +1.1 +- 0.4 / +0.0 +- 0.4 (R4) %/GeV: flat within errors. Net versus distance to cathode (0x, R2): -1.3, -1.2 ... +0.4 %, R4: -2.1 % in the first 15 cm to +0.8 % at 120 cm: same sign pattern as SBND at ~1/3 the amplitude.

**Oddness.** Joint bootstrap R2: odd +0.58 % (chi2 65/17), even 0.00 % (24/17); R4: odd +0.55 % (89/17), even -0.02 % (26/17). Same even churn: P(lost in 2x | lost in 0x) = 40 % (base 5.3-5.5 %), 35-40 % of 0x-new are 2x-new. Hamming: R2 14.3/12.1/15.7 %, R4 13.9/12.8/15.6 % (nominal-0x/nominal-2x/0x-2x): coherent part 19 % (R2) and 17 % (R4) of the nominal-to-varied selection difference, using the mean of the two nominal distances.

**1p vs Mp.** all/1p/Mp: R2 -0.58/-0.50/-0.82 %, +0.59/+0.64/+0.45 %; R4 -0.65/-0.64/-0.68 %, +0.50/+0.50/+0.51 %; Mp-1p within 1 sigma. Channel mismatch rms 0.3-0.6 % (1p) and 0.7-1.8 % (Mp, mostly noise) versus systematic 0.6-1.0 %; chi2 8-19 for 15 bins: not significant.

**Learned / uncertain / next.** SCE in ICARUS is a +-0.55-0.65 % odd net with the same churn structure as SBND; the odd net scales with the physical displacement (SBND/ICARUS odd-net ratio ~2.5, median vertex shift ratio ~2.6). Uncertain: the 5-6 % churn (non-odd part) as in SBND; no additional ICARUS-specific issue. The containment cut (Sec. 9) is not in the ratio calculation.

---

## 8. ICARUS WireModxThetaXW and WireModYZ (Run 2 and Run 4)

Figures: `icarus_wm_ratio.png`, `icarus_wm_mech.png`.

**WireModxThetaXW.** Physical: signal width/area as function of x and theta_xw. It acts on the calorimetry, strongly and angle dependent: the muon calo-vs-range difference shifts by +3.6 %-pt (R2), +1.95 %-pt (R4) on average, growing monotonically with theta_xz from +1.9 %-pt (tracks transverse to drift) to +5.8 %-pt (along drift), muon score -1.6 / -0.9 %-pt (to -3.6 %-pt along drift). E_vis median 0, mean +0.10 %, vertex 0.14 cm, direction 0.05 deg. In SBND the same variation shifts calorimetry by only -0.4 %-pt (erratic in angle). Net +0.49 +- 0.11 % (R2), +0.52 +- 0.10 % (R4, after de-duplication, see Sec. 9; the production histogram has +1.07 %); per-bin rms 0.9-1.0 % vs errors 0.5-0.8 %; chi2 vs unity 29 (R2), 41 (R4): weakly significant; slope +1.3 +- 0.5 / +0.9 +- 0.4 %/GeV. Lost/new 6.6/7.1 % (R2), 6.25/6.77 % (R4).
* **Channel dependence is significant here**: all/1p/Mp = +0.47/-0.07/+1.95 % (R2: Mp-1p = +2.02 +- 0.39 %, 5.2 sigma) and +0.52/+0.26/+1.20 % (R4: +0.94 +- 0.31 %, 3.0 sigma, de-duplicated). The applied (combined) ratio therefore gets Mp wrong by 1.5 % (R2) / 0.7 % (R4) and 1p by 0.5 % / 0.3 %, each comparable to or larger than the systematic itself (rms of channel difference 1.14/1.96 % (1p/Mp; noise-subtracted 0.98/1.43 %) versus rms systematic 0.99 %; chi2 47/45 for 15 bins).
* Sensitivity: the ICARUS selection responds much less to a calorimetric shift than SBND's: ratio slope per %-pt of muon calo shift 0.36 (R2) and 0.47 (R4) %/GeV versus 3.2-3.9 in SBND (a factor ~8-10); score shift per %-pt -0.44 (ICARUS WireModxThetaXW, both runs) versus -1.5 to -2.5 (SBND). That is why a +3.6 %-pt calorimetric shift produces only a 0.5 % net in ICARUS.

**WireModYZ.** Mean calo shift +0.79 %-pt (R2), -0.08 %-pt (R4), vertex 0.15 cm, score -0.2 / -0.07 %-pt; net +0.38 +- 0.10 % (R2) and -0.13 +- 0.08 % (R4); rms 0.73 / 0.40 %, chi2 vs unity 23 / 11 for 17 bins: compatible with noise (R4 certainly). No slope. Lost/new 6.4/6.8 % and 6.6/6.5 %. Channel: all/1p/Mp +0.37/+0.26/+0.68 % (R2), -0.13/-0.07/-0.29 % (R4), insignificant mismatch (chi2 7 / 14 for 15 bins). Same boundary pattern: p-margin SMD -0.78/-0.87 (R2), -0.83/-0.80 (R4), 55-57 % of the lost events above the margin floor (the highest of all variations). Calo shift is flat in theta_xz (+0.8 %-pt), consistent with a spatial (y,z) not angular map.

**ICARUS Run 2 versus Run 4.** Same churn, floors and boundaries; the physical shifts differ (WireModxThetaXW calo +3.6 versus +1.95 %-pt; WireModYZ +0.79 versus -0.08 %-pt) because the maps are run-specific, while the integrated nets are all within +-0.65 %.

**SBND versus ICARUS (same variation).**
* SCE: ICARUS displacement ~1/3, odd net 0.55-0.58 % vs 1.44-1.49 % (SBND), churn 5.3-6.7 % vs 6.9-8.5 %, floors lower (3-5 % vs 4-9 %), same boundary (proton margin) and class structure (the 2x net is carried by other-numuCC).
* WireModxThetaXW: much stronger calorimetric effect in ICARUS (+2..+6 %-pt) than in SBND (-0.4 %-pt), yet comparable small nets (SBND -0.13 % null, ICARUS +0.5 %) because of the ICARUS selection's lower calo sensitivity and the dominance of fragile churn.
* WireModYZ: SBND is a large (-2.2 %, slope +7.7 %/GeV) and different-in-kind perturbation (waveform-shape effect, 3.7x score shift), ICARUS is a null; the same variation name does not imply the same size.
* Trigger: SBND none, ICARUS +-2 % at 0.3-0.4 GeV (Sec. 9).

---

## 9. ICARUS TriggerEmulation (and SBND for comparison), and bugs/inconsistencies

Figure: `trigger.png`. Reproduced exactly for Run 2 and Run 4 (max |difference| = 0.0): pass if `scale x event_largest_flash_pe > threshold`, Run 2 0.615 +- 0.018, 6000 (PE > 9756; +-1 sigma window 9479-10050 PE = +-2.9 % in PE); Run 4 0.361 +- 0.013, 3000 (PE > 8310; window 8021-8621 = +-3.6 %).

| | SBND | ICARUS Run 2 | ICARUS Run 4 |
|---|---|---|---|
| nominal PE threshold | 3091 | 9756 | 8310 |
| median PE per GeV (>0.5 GeV) | 24.5k | 38.5k | 38.5k |
| selected events failing nominal trigger | 0.18 % | 2.07 % (20.6 % in 0.3-0.4 GeV, 13.2 %, 6.1 %, 2.2 %, 0.7 %) | 1.05 % (13.5 %, 6.9 %, 2.4 %, 0.7 %) |
| events inside the +-1 sigma window | < 0.01 % | 0.48 % (3.1 % in first bin) | 0.35 % (3.4 % in first bin) |
| M1 ratio-1, first 4 bins | 0 | -2.03, -1.55, -0.94, -0.50 % | -2.04, -1.50, -0.53, -0.25 % |
| P1 ratio-1, first 4 bins | 0 | +1.83, +1.57, +0.90, +0.41 % | +1.83, +1.06, +0.51, +0.17 % |
| integrated M1 / P1 | -0.001 / +0.001 % | -0.26 / +0.24 % | -0.20 / +0.16 % |

**Energy dependence.** Flash PE is proportional to E_vis (ICARUS 38.5k PE/GeV above 0.5 GeV; SBND 24.5k) with a ~25-30 % relative spread, so the PE distribution at 0.3-0.5 GeV (median 16-17k in ICARUS) overlaps the threshold (8-10k), while above 0.7 GeV it is more than 2x above it: the failure fraction and the fraction in the sigma window fall from 20.6 %/3.1 % (0.3-0.4 GeV) to < 0.1 % beyond 0.7 GeV. The ratio-1 equals minus/plus the fraction of the bin's events lying between the nominal and the shifted threshold and is close to linear (M1 ~ -P1); it is the cleanest systematic, deterministic subset of the same events (no re-simulation noise).

**Size.** ICARUS: +-2 % at E_vis < 0.4 GeV, +-1-1.5 % at 0.4-0.5 GeV, <= 0.5 % above 0.52 GeV, integrated +-0.2-0.25 %. This is larger than the bootstrap errors of the other ICARUS ratios (0.5-1.2 %) only in the first 2-3 bins; it is not negligible in the first bins but small elsewhere. The effect is linear in sigma (a sigma twice as large doubles it). Whether sigma = 0.018 (2.9 %) / 0.013 (3.6 %) is the right size cannot be tested here; it is determined by the data/MC PE comparison in the technote. SBND: negligible, see Sec. 6.

### Bugs and inconsistencies found

1. **ICARUS ratio calculator vs analysis selection**: `reco_user_containment_cut_filterEE` is not applied in the ratio calculation (in the pair tables 37.5 % of ICARUS Run 2 nominal-selected events have it = 0) although the analysis XMLs apply it for Run 2; flagged for the lead to quantify.
2. **Duplicated events in the ICARUS Run 4 WireModxThetaXW varied sample** (only sample affected): the `events` tree has 13,878 repeated (Run, Subrun, Evt) keys (0.53 %; runs 12965/12966/12968/12970) and the `selected` tree 1,107 extra copies (0.57 % versus 0.05-0.06 % in the other ICARUS trees); 89 % of the duplicate groups are identical rows. In the production ratio `WireModxThetaXW_icarus_Run4` every bin is inflated by +0.56 % on average (+0.84 to +1.32 % in the unique-key comparison); integrated production +1.07 % versus +0.52 % after de-duplication. The existing pair table (and `baseline.json`) excludes all copies on the varied side only, which biases that pair low (net -0.03 %, lost 6.77 % instead of 6.25 %). I used a de-duplicated table (`runs/syst_study/pairs_fixdup/`, script `ps_fixdup.py`). The production histogram should be regenerated from de-duplicated events (the effect, +0.56 % flat, is as large as the systematic itself). All SBND trees have ~0.46 % duplicated selected rows in nominal and variations alike (harmless: symmetric); no other sample has duplicates in `events`.
3. **Trigger requirement not applied to the detvar ratios** (the analysis weight includes `scale x PE > threshold`): the effect of applying it is <= 0.04 % integrated for all pairs; per bin (max over bins / rms over bins) SBND <= 0.20 % / 0.07 %, ICARUS R2 up to 0.58-1.13 % / 0.15-0.29 % (WireModYZ/xThetaXW), R4 <= 0.41 % / 0.12 %. For ICARUS Run 2 WireModYZ the rms shift (0.29 %) is 40 % of the systematic's rms (0.73 %).
4. **Truth not filled in the varied trees of SBND WireModxThetaXW and WireModYZ**; the `true_*` branches exist but are NaN (30 % of gained events cannot be classified even using other variations).
5. **Unclear calculator selection**: `detsys/calculator/configuration.h` (by_vis_energy for SBND) does not contain the |x_end| >= 4 cm cut that `common.selection_mask` (calculator_endx4) applies; the production histograms are reproduced with the latter, so `detsys/calculator_endx4` is presumably the production code; the older directory should not be used.
6. SBND SCE ratios (technote: not used) are present in `xml/coverage` allowlists only as commented lines (consistent).
7. Pair-table ratios differ from production by up to 0.2 % (0.7 % DENT last bin) from the exclusion of duplicate keys (symmetric duplicates drop out of the pair but not of the production count).

---

## 10. Cross-systematic comparison and open questions

| quantity | SBND SCE | SBND WireModxThetaXW | SBND WireModYZ | SBND MattModYZ | SBND DENT | ICARUS SCE | ICARUS WMxThetaXW | ICARUS WMYZ | Trigger |
|---|---|---|---|---|---|---|---|---|---|
| net (%) | -1.35 / +1.53 | -0.13 | -2.20 | -0.27 | -0.87 | -0.6 / +0.55 | +0.49 / +0.52 | +0.38 / -0.13 | 0 (SBND), +-0.2 (ICARUS) |
| lost + new (%) | 15.6 / 15.4 | 15.0 | 18.9 | 4.0 | 25.1 | 11.2-12.8 | 13-13.8 | 13 | -- |
| significance vs noise (chi2/17) | 198 / 222 | 7 | 717 | 123 | 44 | 24-65 | 29 / 41 | 23 / 11 | exact |
| energy trend | flat (odd offset) | none | -3.7 -> +7.9 % | -0.5 -> +3.7 % | -2 -> 0 % | flat | flat | flat | first 3 bins only |
| mechanism | odd geometry/length (cathode-weighted) + even churn | hit-shape; churn only | hit-shape + charge; PID score | pure charge scale (+0.5 %-pt) | tier-2 end-x/crosser + non-localized churn | as SBND, 1/3 amplitude | strong angle-dependent calo, low sensitivity | null | threshold crossing |
| 1p/Mp mismatch vs systematic | not significant | noise | <= 75 % of sys (chi2 36/15) | = sys, 4.8 sigma | n.s. | n.s. | R2 Mp +1.95 vs 1p -0.07 (5 sigma) | n.s. | none |

**What we learned (short).** (a) The detector ratios are dominated by sign-independent churn of a small, background-enriched, boundary-hugging population of fragile events; the physical response appears as a small residual: odd and position dependent for SCE, charge-scale-like for MattMod, null for WireMod xThetaXW (SBND) and ICARUS WireModYZ. (b) The energy trends are explained: WireModYZ/MattModYZ trends follow the SBND charge-scale law (3.6 %/GeV per %-pt of calo shift; WireModYZ 3.5x above it); SCE shape is a migration (energy scale, odd) plus flat offset; trigger trends are PE-versus-threshold geometry. (c) SBND and ICARUS differ because the selection in ICARUS is ~8-10x less sensitive to calorimetric shifts, because the SCE amplitude is 3x smaller, and because the ICARUS and SBND WireMod maps are different perturbations under one name.

**What remains uncertain, and the study that would settle it.**
1. Physical versus chaotic churn: an identity-transform reconstruction of the nominal (displacement set to zero through the SCE/DENT/WireMod code path) and a same-input different-seed re-simulation would give the noise floor directly. Without it, the churn-floor of 3-9 % per direction is unexplained.
2. Origin of the floor: dump the intermediate selection variables (pion/photon/electron multiplicities and KE near the 25 MeV thresholds, muon multiplicity, proton count at 50 MeV) for events lost in the varied sample (currently only selected events are stored), to attribute floor and boundary shares to specific thresholds; in particular the +18.6 % other-numuCC growth in SBND SCE2x.
3. Channel dependence: ICARUS WireModxThetaXW needs channel-specific ratios (Mp +1.95 % versus 1p -0.07 %, 5 sigma, Run 2) and SBND MattModYZ (1p -0.40 % versus Mp +0.13 %); production could supply 1p/Mp ratio histograms if the mismatch matters for the fit. The mismatch is statistically insignificant for the other variations.
4. Statistical content of the splines for null variations (SBND WireModxThetaXW chi2 7/17, ICARUS WireModYZ, ICARUS R4 SCE2x): the spline then encodes 0.3-1.2 % bin noise (see Q2); a smoothing or a flat-uncertainty replacement may be more honest.
5. DENT localization: check west-TPC behavior of the DENT map; if zero there, the west-TPC churn is a direct null for the pair.
6. Regenerate `WireModxThetaXW_icarus_Run4` from de-duplicated events; quantify the containment-cut mismatch (lead).

## Files

* Scripts: `plotting/syst_study/ps_common.py, ps_core.py, ps_pairs.py, ps_fragile.py, ps_trigger.py, ps_fixdup.py, ps_summary.py, ps_plots.py`
* Data: `docs/syst_study/persyst/summary.json` (all numbers), `core/*.json` (per pair), `pairs_compare.json`, `trigger.json`, `fragile.json`, `fixdup.json`, `duplicates.json`; `runs/syst_study/pairs_fixdup/`, `runs/syst_study/ps_truthmap_*.pkl` (truth lookup)
* Figures: `docs/syst_study/persyst/{overview, fragile, sbnd_sce_ratio, sbnd_sce_mech, sce_linear, sce_event_level, sbnd_wmxt_ratio, sbnd_wmxt_mech, sbnd_yz_ratio, sbnd_yz_compare, sbnd_yz_mech, trigger, sbnd_dent_ratio, sbnd_dent_mech, icarus_sce_ratio, icarus_sce_mech, icarus_wm_ratio, icarus_wm_mech}.png`
