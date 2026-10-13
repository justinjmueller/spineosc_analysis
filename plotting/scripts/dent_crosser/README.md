# dent_crosser scripts

Standalone analysis/plotting scripts from the DENT cathode-crosser investigation
(SBND Run 1, v1.1.0). Not wired into a pipeline -- each is a self-contained
script, run directly with `python3 <script>.py`. Paths to ntuples are
hardcoded at the top of each file; update them if the source files move.

- `cathode_crosser.py` -- ratio-histogram fractional uncertainty from DENT,
  with/without a cathode-crosser cut (the "fractional uncertainty with vs.
  without a cut" numbers: 2.04%/0.52% normalization, etc).
- `dent_migration.py` -- per-particle PID-score displacement between CV and
  DENT for common (truth-matched) events, used to demonstrate DENT migrates
  muon-score classification near the cathode.
- `dent_source.py` -- exploratory scan of DENT ratio histograms across
  variables, looking for which observables show strong shape effects.
- `calo_csda_crosser_split.py` -- DENT ratio split by cathode-crosser status,
  applied to muon_calo_csda.
- `muon_score_crosser.py` -- DENT ratio split by cathode-crosser status for
  muon_score; the script that found the real per-population DENT effect
  after the XML name-collision bug was fixed.
- `muon_score_dent_crosser_split.py` -- directional (east-to-west /
  west-to-east) crosser split for muon_score.
- `dent_crosser_ternary.py` -- ternary (3-class PID score) density plots of
  true muons, CV vs DENT, split by crosser/non-crosser (uses mutools.plotting.ternary).
- `dent_crosser_ternary_ratio.py` -- ternary ratio/displacement plot between
  CV and DENT (paired via common-event keys).
- `dent_crosser_true_particles.py` -- variant of the ternary ratio plot using
  v1.0.10 particle-level truth-matched ntuples instead of event-level selected
  trees.
- `dent_crosser_signal_ratio.py` -- variant of the ternary ratio plot using
  the `signal` tree instead of `selected`.
- `build_crosser_endx.py` -- generates
  `xml/detector_data_access/sbnd_run1/muon_endx_zoom_crosser{,_data}.xml`
  from `muon_endx_zoom_inclusive{,_data}.xml` by duplicating the channel into
  stay/xing populations. XML generator, not a plotting script, but kept here
  since it's part of the same investigation and the pattern is reusable for
  splitting other zoom XMLs the same way.
