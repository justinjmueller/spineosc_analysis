# Q3: Can the direct PROfit DetVar method reproduce the detector-systematic ratio histograms?

Short answer: **yes, to float precision (max deviation ~1e-7 per bin, 4e-4 once from a bin-edge event), but only with an intermediate-tree construction; the matched mode as it is used today cannot.** No cvwgt anywhere. Unmodified PROfit (`/Users/mueller/Projects/profit-asymmetric`, untouched) suffices. A one-line patch, tested in a copy (`/Users/mueller/Projects/profit-asymmetric-q3`, patch in `runs/syst_study/q3/q3_keep_zero_weight.patch`), makes the construction cleaner.

Figures: `docs/syst_study/q3/q3_fig*.png`. Machine-readable results: `docs/syst_study/q3/q3_summary.json`. Scripts: `plotting/syst_study/q3_*.py`. Runs: `runs/syst_study/q3/<case>/<variant>/`.

## 1. What PROfit's DetVar method does (code reading, src at c56f298)

1. **Config** (`src/PROconfig.cxx:1162-1460`): `<DetVarSection>` holds one `<cv>`, one or more `<variation name knobval>` and `<subchannel>`s; optional `cv_variation_matching_vars` (1215-1228). For every file a mini-config is built by `BuildDetVarConfig` (3339+); the DetVar subchannels inherit the `<branch>` definitions (selection formula as `additional_weight`, variables, bins) of the main MCFile (1400-1450). So the CV and variation trees are processed with the analysis' own formulas. `m_detvar_matching_vars` is copied into the mini-config only (3399).
2. **Reading** (`bin/PROfit_process.cxx:182-205`): each CV/variation tree goes through `PROcess_CAFAna` into its own `PROpeller` (per-event vectors `added_weights`, `variable_bin_indices`, ...). Cached in `<tag>_detvar_props.bin`.
3. **Weights / selection** (`src/PROcreate.cxx:1343-1385`): `mc_weight = product(branch weights, incl. additional_weight) * det_pot/mcpot` (POT per `<cv pot>`/`<variation pot>`). **`if(mc_weight == 0) return;` (line 1374): zero-weight events are never stored.** This is where "not selected" disappears. With matching vars set, key values are evaluated per TTree entry as float (1075-1095) and stored only if the event was stored (1092).
4. **Matching** (`PROfit_process.cxx:19-130`, `BuildDetVarMatchedSpecs`): key = vector of `static_cast<int>(round(float))` per matching var. A key set is built for CV and for every variation; `common_keys` = keys present in CV **and in all** variations. Only propeller events with a common key are histogrammed (`QuickFill` with the stored weight). Repeated keys are all kept (lists), if the key is common. Since stored events are by construction the selected ones, **matched mode = "selected in both" only**; the CV-only and variation-only events (selection migration) are dropped. This is the 3.4% / "half of the structure" loss in the DENT prototype.
5. **POT**: unmatched mode: the spectra carry the POT scaling (`det_pot/pot` per file). Matched mode: undone again (`PROfit_process.cxx:296-312`, spectra multiplied by `pot/det_pot`), i.e. raw weighted counts.
6. **Response** (`PROfit_process.cxx:326-420`, `PROsyst.cxx:1284-1440` `FillSpline`): bins where the matched CV is 0 are zeroed; `SystStruct` with `p_cv`, `p_multi_spec[knob]`; the knob=0 point is inserted as CV/CV=1; ratios `spec[knob]/p_cv` per bin of the chosen binning (`binning="varN"` of the allowlist) feed a cubic spline in knob (linear/"mirrored" for one variation; `mirror="false"` handled by the same `FillSpline`). **hist1d systematics go through the same `FillSpline`** (mode list at PROsyst.cxx:~150), so the interpolation, mirroring, and inflate behave identically to the production ratio histograms.
7. **Errors**: `PROspec::Error()` (sqrt sum w^2) is filled only for plotting. No statistical uncertainty of the response enters the spline, the covariance, or the chi2 (production ratio histograms also have zero errors; `calculator` sets bin errors 0). The MC-stat error comes from the main CV propeller only.
8. **Binning dependence**: the ratio is computed in the observable's own bins (any binning in the XML; `binning=` selects). Several DetVar sections with the same name are combined as CV-weighted averages (`dv_pieces`).

Why matched loses events, precisely: zero weight -> not stored (point 3) -> key intersection (point 4) -> only both-selected survive.

## 2. Tree constructions that work (tested)

Test harness: pair tables `runs/syst_study/pairs/` (union of selected unique-key events with `in_n`, `in_v`), written with uproot as real TTrees (note: `f["name"] = dict` writes an RNTuple that segfaults PROfit; use `mktree`). One-channel XML, one subchannel, bins = production bins, `process` then `plot` twice (without / with `--inject-systs-cv Q3_native 1.0`); response = injected / un-injected CV spectrum (`q3_native.py`). Equal POT on both sides.

| Variant | Trees | PROfit | Result (ChargeScale15, E_vis) |
|---|---|---|---|
| **A** pool-restricted, unmatched | selected rows in the common pool of CV and variation, equal POT | unmodified | = numpy pair ratio to 7e-8 |
| **Bm / Cz** matched | selected-only trees (Bm), or union trees with `selected` as 0/1 weight (Cz) | unmodified | max dev from production 0.116 (matched emulation agrees); matched counts: 325188 both, 21839 CV-only, 19598 var-only, as in my tables |
| **Ce** selected attribute, eps weight | every union event on both sides; `additional_weight="selected + 1e-9*(1-selected)"`; keys `Run,Subrun,Evt,Ei`; reco energy of unselected rows set to 0 (underflow) | unmodified | 7e-8; matching: 0 CV-only, 0 var-only |
| **Cp** selected attribute, plain 0/1 weight | as Cz | **patched** copy: `if(mc_weight==0 && inconfig.m_detvar_matching_vars.empty()) return;` | 7e-8 |
| **U** unmatched on the whole production trees | `events/<tree>/selected` of both files, cut as `additional_weight`, true POTs | unmodified | = numpy "unmatched full, POT" to 5e-4; differs from production by noise (section 3) |

Tested on SBND ChargeScale15, SBND DENT, ICARUS Run2 SCE2x, ICARUS Run4 WireMod xThetaXW (A, Bm, Ce, U; Cp on the first three). With **all repeated-key rows kept in the trees** the native result equals the production ratio histogram itself: 0.0 (ChargeScale15), 9e-8 (ICARUS Run2 SCE2x), 1.1e-7 (ICARUS Run4 WireMod), 3.9e-4 (DENT; one event on a float32 vs double bin edge) (`q3_exact.py`).

**Full event list** (the construction you described): all 6,274,248 pool events of ChargeScale15 on each side, `selected` from the selected tree: Ce and Cp both process in 36-38 s with ~3.7-3.9 GB RSS, plot 19 s each, same 7e-8 agreement. So the intermediate tree can be "events list + selected flag" directly; the union of selected events (370k rows, 2-4 s, <0.6 GB) is enough and much cheaper, because rows selected on neither side contribute nothing.

**Multi-universe** (SBND SCE0x = knob -1, SCE2x = knob +1, one 3-knot spline, matched + eps): native equals numpy on the common pool to <1e-7 at both knots. Matched mode requires keys in CV and all variations, so the pool is the intersection (6.155M vs 6.269M/6.223M pairwise); this moves the ratios by up to 0.39% (-1) and 0.15% (+1) relative to production (pairwise pools), consistent with statistical noise (production vs native: 0.58% / 0.18%). Intermediate knob values come from the same cubic `FillSpline` as for hist1d.

**Why eps works / what it costs**: weight 1e-9 is nonzero in float (POT scale 1e19/5e20 -> ~2e-11, still far from denormal), so the event is stored and keyed; its contribution to the bin sums is 1e-9 per event (<1e-3 events in total, if the energy is also set out of range it is exactly nil). Caveats: a nonzero-bin energy for unselected rows would add ~N_rows*1e-9; the main-MC branch must use the same formula or the unselected rows appear as 1e-9 events in the nominal spectrum (negligible but ugly); the patch (Cp) avoids both and is the clean option.

## 3. Numpy emulation of the three modes, 20 detector/variation pairs (E_vis bins)

Errors: paired Poisson bootstrap over union rows, B=300 (`q3_emulate.py`, `q3_fig2_summary_all.png`). "pair" is the ratio on unique-key rows (= selected-attribute = pool-restricted unmatched by algebra, verified natively). It differs from production only by the repeated-key rows.

| Mode | max abs dev (median / max over 20) | chi2/17 vs production (median / max) | Why |
|---|---|---|---|
| pair / selected-attribute / pool-restricted unmatched | 0.0015 / 0.013 | 0.12 / 62 | only repeated keys; the outlier is ICARUS Run4 WireModxThetaXW, see below |
| matched (both selected only) | 0.022 / 0.173 (GainHigh) | 79 / 3121 | drops 4-22% of rows (migration); error grows with migration fraction (`fig3`) |
| unmatched, full trees, POT-scaled | 0.0084 / 0.0245 | 8 / 48 | independent-ish extra events + POT/event-count mismatch |
| unmatched, full trees, raw counts (no POT) | up to 0.25 | up to 4e4 | POT ratio of samples is not 1 (ICARUS ~6%, DENT 5.5%) |

What "unmatched" means here: CV and variation are two processings of mostly the same events (pool = 92-99.5% of each sample), so even independent trees are strongly paired. The residual noise is from the non-pool events (0.5-8% of events) and from different POT/event ratios (mean bias of the POT-scaled ratio: |<dev>| <= 0.0075, typically 4e-4). If the samples were truly independent, the noise would be 1.9-5.2x (mean 2.9x) larger than the paired one (`fig6`): typical paired sigma per bin is 0.5%, independent 1.5%. So the "~1% statistical noise" of the prototype's unmatched mode is real but caused by a handful of unpaired events and by POT bookkeeping; it is not catastrophic here, yet it is avoidable.

**Discovered production artifact**: ICARUS Run4 WireModxThetaXW: the varied sample contains duplicated events (13,878 repeated keys in the events tree, 1,786 selected rows with repeated keys vs ~150 normally; 1,882 of the 2,090 duplicated selected rows are fully identical). The calculator counts them, inflating the production ratio by about +1.2% in every bin (production minus deduplicated ratio = -1.1..-1.3% in all 17 bins, 62 chi2). The direct method fed with the same rows reproduces the artifact; fed with unique keys it removes it. Worth deduplicating in the production tree.

## 4. Per-observable ratio files (v1.1.0v2) and the selection

The per-observable files `detvar/detsys_v1.1.0v2_<obs>.root` use the same cut as by_vis_energy (SBND crosser veto and |end_x|>=4), **not** the "fiducial && veto" cut of other calculator versions. With that cut: my "production-style" ratio (repeated keys kept) matches the file exactly (max dev 0.0) for muon_score and muon_endx_sbnd, and native PROfit (A, Ce) in those observable bins agrees with the repeated-key-free ratio (max dev 0.0008-0.0345, all from repeated-key rows, chi2 <= 0.3). 12 cases (`q3_obs_repro.py`, `fig4`). Matched mode fails badly in sparsely populated bins: muon score max dev 0.44 (ChargeScale15), 11.8 (GainHigh), 1.03 (WireModYZ); end-x DENT 0.088 with a spurious +7% at the cathode (`fig4`, third panel). Unmatched full-POT: up to 0.02 (SBND) but 0.06-0.13 in low-stat ICARUS bins.

## 5. Critical assessment

**Can the direct method reproduce the ratios?** Yes. Exact (<=1e-7) when the CV and variation trees carry the same events (pool) and the selected flag; this was shown through the real PROfit code, for SBND/ICARUS, three binnings, one- and two-sided knobs, and with 6.3M-event trees. Not with matched mode on selected-only trees (errors 1-17% E_vis, up to a factor 12 in sparse bins). Unmatched mode on the full production trees reproduces to 0.2-2.5% (median 0.8%) with correct POTs, 5-25% wrong without.

**Caveats found**
- Matching key goes through float, then `int(round())`. Run/Subrun/Evt/Ei are all < 2^24 (max Evt 591,635; max 1e6*E = 7.2e6), so exact. The prototype key `event_max_neutrino_energy` is cast to integer GeV, i.e. effectively not a key: ICARUS Run2 nominal has 400 (Run,Subrun,Evt) collisions between distinct events; with E rounded to GeV 230 stay ambiguous, the integer `Ei=round(1e6*E)` branch resolves all (`keys.json`). A pool_id branch also works. Add such an integer column to the intermediate tree.
- Duplicate keys: all rows are kept consistently if the key is common (matches calculator semantics, including its double counting).
- POT: matched mode cancels it (so equal POT is fine); unmatched mode needs the true per-file POT (read the `POT` histogram), otherwise 6-25% offsets.
- Stat uncertainty of the response: not propagated by either method (ratio histograms: errors set to 0; DetVar: Error() is plot-only). The paired bootstrap error is ~0.5% per bin typical (up to 6% in the last bins), comparable to the effect of many variations (DENT E_vis effect <= 2%, i.e. mostly noise-dominated). The native method could propagate it (patch to add an error band per spline bin) but doesn't today.
- Zero-CV-denominator bins are zeroed (ratio then 0 or 1 depending on path), whereas the calculator sets the ratio to 1; check empty bins (needed `--force` for empty muon-score bins).
- Memory/time: no problem (4 GB, <1 min for 6M rows/side). `std::map<vector<int>>` matching scales, but one DetVar section per variation and ~20 variations x several observables x ~20-40 s is minutes, not hours.

**Gains**
1. Any binning without re-running the calculator: the ratio is built per observable and per bin edge change directly from the trees; `binning=` per systematic.
2. Selection consistency: the response is evaluated with the analysis' own selection (subchannel truth cuts, proton multiplicity, flash cuts) on the varied sample. The calculator ratio is inclusive. For the exclusive channel reco_proton_multiplicity==1 the direct response differs from the inclusive production ratio by up to 17% (ChargeScale15, chi2 135/17), 8.5% (GainHigh), 3-4% others; `pm>=1` is identical because the selection already requires a proton (`q3_selection.py`, `fig5`).
3. Per-subchannel response, one consistent event list across all observables in a fit.
4. Same spline machinery as hist1d (mirroring, `force_0_cv`, inflate, knobs), so fit behaviour is unchanged.

**Losses / risks**
- The response is per DetVar subchannel and noisier (smaller samples, per-subchannel ratios); no stat error is propagated, so noise would be treated as systematic signal. Smoothing or a stat error would be wanted.
- Intersection pool across all variations of one knob (+-): small (0.2-0.4%) shifts.
- Trees must carry all branches used by the analysis selection for both CV and variation, and per-key CV/variation truth consistency is assumed.
- Matched-mode default is silently wrong for selection migration, nothing warns; the log line "num CV only/var only" is the only hint.
- Production reproduction includes its own artifacts (ICARUS Run4 WireMod duplicates).

## 6. Recommendation and steps

1. Keep ratio histograms as the validated reference; adopt native DetVar with the **selected-attribute** construction for observables where binning flexibility or exclusive selections matter (DENT, endx, exclusive proton channels).
2. Produce the intermediate trees from the existing machinery: for every variation, the union of keys selected in CV or variation within the pool; columns `Run,Subrun,Evt,Ei(int), selected, <all analysis branches>`; on the side where the event is not selected, a row with `selected=0` (reco branches may be dummy). Write with `mktree` (not dict assignment). Deduplicate repeated keys (or keep, matching calculator).
3. XML: `<DetVarSection cv_variation_matching_vars="Run,Subrun,Evt,Ei">` with equal POT, `additional_weight = <analysis cut> * selected` (or `selected + 1e-9*(1-selected)` with unmodified PROfit); set `reco` variables to an out-of-range value for unselected rows; use `--force` if empty bins.
4. Preferred PROfit change (in a branch of the user's checkout, tested in the copy): keep zero-weight events when `m_detvar_matching_vars` is non-empty (`src/PROcreate.cxx:1374`), see the patch file. Optional follow-ups: warn when matching drops CV-only/var-only events; propagate paired stat error of the response (not done).
5. Validate by reproducing the shipped ratio files (`q3_exact.py`, `q3_obs_repro.py`) before using native responses in fits.

## 7. Open points
- Not tested natively: more than two knots on a non-SBND detector, the multi-section combination, and full analysis XML (per-subchannel truth cuts) end to end; the selection-gain numbers (section 5) are numpy only, with no PROfit run.
- Flash/other cuts of the analysis XML were not evaluated (pair tables lack the needed truth-detector logic).
- `muon_score` empty bins required `--force` (PROfit aborts otherwise).
- The ICARUS Run4 WireMod duplicate cause is not traced to the production step.

## Files
Report: `docs/syst_study/Q3_direct_detvar.md`; summary `docs/syst_study/q3/q3_summary.json`; figures `docs/syst_study/q3/q3_fig1..6_*.png`. Scripts `plotting/syst_study/q3_{full,emulate,native,native_full,fullpool,multi,exact,obs_repro,observables,selection,keys,report,cmp,cmp2}.py`. Data `runs/syst_study/q3/` (`emulation.json`, `full_hists.json`, `obs_repro.json`, `selection.json`, `keys.json`, `observables.json`, `q3_keep_zero_weight.patch`, per-case run dirs). Patched copy: `/Users/mueller/Projects/profit-asymmetric-q3` (built).
Note: `observables.json` (E_vis ratio applied per event in other observables) was an exploratory comparison based on a wrong assumption (per-observable ratio files exist), and is not used in conclusions.
