# fdgen: blind fake-data weights from a YAML scenario file

Per MC event a scenario is a **ratio to the PROfit prediction**: the oscillation probability, evaluated the way PROfit
evaluates it, times any number of dial weights from the CVMC friend trees. Optionally the expected counts are Poisson
thrown. Generated outside PROfit; PROfit only reads the result as fake data.

```bash
PYTHONPATH=$(root-config --libdir) .venv/bin/python -m fdgen.generate fdgen/configs/closure_2d.yaml   # ~5 min
.venv/bin/python -m fdgen.data_xml output/fdgen/closure_2d osc_a -o xml/fakedata/data/fdgen_closure_2d_osc_a.xml
cd runs/fdgen_2d && ./fdgenjob.sh closure_2d_osc_a free          # global fit, MC caches shared (copies of the amr1 caches)
.venv/bin/python -m pytest fdgen/tests
```

## What follows PROfit, and where
| piece | PROfit | fdgen |
|---|---|---|
| probability | `SBL_2flav_numudis` (`src/PROmodels/PROmodelSine.cxx`): `1 - s sin²(1.266932679 Δm² L/E)`, `s` clamped to [0, 1], float32 | `osc.pmumu_profit` |
| where it is evaluated | at the **midpoint of the truth L/E bin** (`<bins unit="True L/E">`, 300 bins in [0, 3]); events at or above the last edge are dropped | same bins, same edge rule (`std::upper_bound`); out-of-grid events are reported and left at ratio 1 |
| which events oscillate | the `model_rule` of the XML branch the event falls in (1 = "Numu Dis", 0 = none) | the same XML, the same `additional_weight` expressions (evaluated with ROOT) |
| POT | weight x (detector POT / MC file POT) | the same factor, from the same XML |
| Poisson | `PROspec::PoissonVariation`: each reco bin ~ Poisson(mu_bin) | each event ~ Poisson(lambda_i), lambda_i = aw_i x pot_scale x ratio_i; a bin then holds Poisson(sum lambda_i) = Poisson(mu_bin) for any binning. numpy's stream, not mt19937, so the draws differ and the law is the same |

`dataset.poisson: false` (default) gives the Asimov dataset, so the first closure tests carry no fluctuations.

## Output
`<outdir>/<tree path>_mystery.root` per MC file: a copy of `selected` plus `<tree>_mysteryTree` (entry-aligned) with
`mysterywgtNN` and, with Poisson, `mysterycntNN`. The friend tree has to share the file with `selected` (PROfit ignores
`<friend filename=...>`), hence the copy. `manifest.yaml` is the shareable description, `KEY.yaml` the private truth
(`scenarios`, shuffled indices, expected yields).

## Worth knowing
* In `systematics_suite.xml` the NC branches carry `model_rule="1"`, so PROfit oscillates NC events with the numu
  disappearance probability. fdgen does the same because it follows the XML.
* Dials: `source: standard` (`selected_multisigmaTree`) or `rebruce` (`selected_fakedataTree`); the weight is the
  per-event spline/linear interpolation at the requested sigma over its value at 0; shifts outside the stored knots are refused.
  Dials are applied to the matched beam tree only (dirt, off-beam, non-matched stay at the CV). Dial composition is not
  yet exercised in a closure test.

## Cross-check against PROfit's own injection (Asimov, stat-only `plot`, `fdgen/crosscheck.py`)
`fdgen/configs/legacy_1d/crosscheck_dials.yaml` lists seven scenarios with the equivalent PROfit arguments; `runs/fakedata_data/fdgenxc.sh <sc> P|H`
runs PROfit's injection (P) and the by-hand dataset (H); `python -m fdgen.crosscheck` compares the total fake-data histograms.

| | reco (analysis) bins rms / max | true-E bins rms / max | total yield |
|---|---|---|---|
| oscillation only | 4e-4 / 0.3% | 1.2e-3 / 1.1% | 1e-5 |
| one or two standard dials (A-D) | 0.7-5% / 4-14% | 0.1-0.9% / 1-2.4% | <= 0.3% |
| oscillation + dials, four dials (E, F) | 3.4-3.7% / 10-12% | 0.5-0.7% / 2% | <= 0.4% |

PROfit's dial splines are binned in `var2`, which is the **true neutrino energy** (the third `<variable>` of each branch), one factor per
(detector, channel, subchannel, true-E bin), and reach the reco bins through the migration matrix. A by-hand per-event weight also depends on
the proton momentum, the interaction type and so on within a true-E bin, so the reco spectra differ by a few percent while the true-E spectra
agree to about 1%.

**Confirmed** with `project_true_e: true` (scenario option, `fdgen/configs/legacy_1d/crosscheck_proj.yaml`), which replaces every dial weight by PROfit's
representation of it, r_j = sum(aw w) / sum(aw) per branch and true-E bin. Against PROfit's injection the reco-bin difference then drops from
rms 0.7-5% (max 4-14%) to rms 0.05-0.11% (max 0.2-0.7%) for all six dial scenarios, and the true-E difference of the multi-dial scenarios
drops to the same ~1.2e-3 rms / 1% (one high-E bin) floor that the oscillation-only scenario already has. So the per-event dial weights are
exact and the differences are PROfit's model (dials as functions of true E only), not fdgen. The by-hand dataset therefore contains dial effects
PROfit cannot represent exactly, which is what a fake-data study wants; set `project_true_e` to get a dataset PROfit can reproduce.

## Downselected output and 2D fits
`dataset.output` in the YAML shrinks the slim files (about 40x; bit-identical PROfit histograms in the cross-check):
```yaml
output:
  selection: xml          # "xml" = events in any XML analysis branch, or a TFormula string (same syntax as the XML weights)
  branches: needed        # only the branches the XML expressions / selection reference + reco_proton_multiplicity + reco_pn (if present) + keep_extra
  keep_extra: []          # more branches; a missing one is an error
```
A selection that drops events which are in an analysis branch prints a warning with the lost weight (those events would be missing from the fake
data). The weight and count trees are filtered the same way, so they stay entry-aligned. `reco_pn` is kept by default, and the 2D analysis keeps it
anyway through `<variable2D>reco_visible_energy;reco_pn/1000</variable2D>`. For a 2D fit set `mc_xml: xml/contours/contours_base.xml`
(`fdgen/configs/closure_2d.yaml`) and run in `runs/fdgen_2d/` (its `fdmc.ini` points at that XML; caches are copies of `runs/amr_compare/amr1_*.bin`).

2D closure results (`closure_2d.yaml`, Asimov, full systematics, global fit): identity chi 9e-6, no signal; by-hand oscillation (2.0, 0.2) recovered
as 1.9999 / 0.20002 (PROfit's own `-i`: 2.0001 / 0.20002); MaCCRES +1 sigma by hand gives pull 0.45, chi 1.15 against PROfit's `--inject-systs`
pull 0.35, chi 0.34. The dial difference is larger than in 1D (0.30 vs 0.28), as expected if the per-event weight depends on pn, which PROfit's
true-E-only splines cannot represent; not yet confirmed with `project_true_e` in 2D.

## Truth-free export (blind hand-off)
`exports:` in the YAML writes, per detector, a data-like ROOT file of **Poisson-thrown events with reco variables only**, plus a data XML:
```yaml
exports:
  - name: pseudo_A            # files <name>_<stem>.root and <name>.xml in dataset.outdir
    scenario: O3              # which scenario's expectation is thrown; the id never appears in the files
    seed: 101
    keep: [reco_visible_energy, reco_pn, reco_proton_multiplicity]
    jitter_duplicates: 1e-6   # events thrown k>1 times get a relative 1e-6 jitter on the repeats (no identical rows); 0 = off
    detectors:
      SBND: {file: pseudo_A_sbnd.root, tree: events/sbnd_bnblight/selected, cuts: {events/sbnd_large/selected: "<reco cut>", events/sbnd_dirt/selected: "...", ...}}
```
- **Selection:** `cuts` (a TFormula string per MC tree, written next to the output file name) are the reco-only cuts of the analysis XML.
  `python -m fdgen.export --suggest-cuts xml/contours/contours_base.xml` prints them from the XML branch weights (truth terms and the channel's
  proton-multiplicity term dropped; `true_detector`/`Run`/flash terms of the flash cut are kept because the XML uses them). A thrown event that fails a
  cut is an error. The expectation of every event is still `aw x pot_scale x ratio` (so the data are what PROfit's MC prediction would give); the cuts
  only state and check the selection.
- **Truth-free:** only the `keep` branches are written. All MC components of a detector (beam, dirt, off-beam, non-matched) are merged and shuffled.
  No `true_*`, `cvwgt`, run/event ids, weights or scenario names are in the files; the XML uses the constant 0 for truth-based variables (the repo's
  real-data convention), the channel's proton-multiplicity term as the only weight, and the detector's target POT as the MCFile POT (factor 1).
- **Private:** the scenario, seed and expected yields go to `KEY.yaml`; do not share it. Set `dataset.write_slim: false` to skip the weight files.

Validation of the export (`fdgen/configs/export_test.yaml`: pseudo_A = O3, pseudo_B = null; 2D analysis; PROfit `plot` without `--scale-by-width`, so the histogram
contents are raw counts): pseudo_A against its own expectation (the by-hand Asimov O3 data) gives chi2/ndf = 98.6/90 (p = 0.25, pull mean -0.05, rms 1.05);
pseudo_B against the CV gives 86.2/90 (p = 0.59); the crossed comparisons are rejected (chi2 about 515 for 90 bins). The files hold only the three reco
branches, the XML contains no truth expression, weight or scenario name, and no two rows are identical. PROfit reads them through the blind XML like data.
Note: uproot-written nested directories crash ROOT (`TFile::Get` segfaults), so the export is written with ROOT (`export.write_tree`).

Oscillation points against PROfit's own `-i` injection (`fdgenosc.sh`, `python -m fdgen.osccheck`): all eight points O1-O8, 1D and 2D analyses: analysis-bin
rms 4e-4 (max 3e-3), total yield 1e-5 to 2e-5.

## Layout (Oct 2026 cleanup)
* `fdgen/configs/batch1.yaml` is the working copy for the first blind batch (2D analysis, 13 scenarios, 13 exports, `export_defaults` holding the cuts once);
  edit it directly. `dataset.keydir` keeps `KEY.yaml` and `manifest.yaml` (private) away from `dataset.outdir` (the shareable files).
* `fdgen/configs/closure_2d.yaml`, `export_test.yaml`: the 2D closure and export validation configs.
* `fdgen/configs/legacy_1d/`: configs on the 1D XML (`systematics_suite.xml`), kept because `crosscheck.py` and `validate_poisson.py` read two of them.
