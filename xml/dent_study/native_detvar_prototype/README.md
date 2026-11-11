# PROfit-native detector variations (prototype)

`<DetVarFiles>` lets PROfit build the detector-variation response itself from a CV tree and variation trees, in whatever binning the
channel uses, so no per-observable ratio histograms are needed. Both files are Tier-A-style (detvar CV as the MC) with a DENT section:

* `native_unmatched.xml` : `<DetVarSection filename=...>` with `<cv treename="events/sbnd_detvar_newcv/selected">` and
  `<variation name="DENT_native" treename="events/sbnd_detvar_newdent/selected" knobval="1">`; registered as `type="spline"`.
* `native_matched.xml`   : the same plus `cv_variation_matching_vars="Run,Subrun,Evt,event_max_neutrino_energy"`.

Result for `muon_abs_endx` (shape of the DENT-shifted prediction, vs the precomputed `DENT_sbnd_Run1` ratio):
unmatched agrees to within 1.1% per bin (statistical noise); matched is off by up to 3.4% and loses about half of the near-cathode
structure, because PROfit's matching keeps only events selected in BOTH samples (zero-weight events are never stored in the propeller).
Of 261,149 selected CV events and 242,605 selected DENT events only 208,183 are in both (CV-only 52,966; DENT-only 34,422), so the
selection migration that is a large part of DENT is discarded. Run: `PROfit --config <ini> process`, then
`plot --inject-systs-cv DENT_native 1.0`. The CV/variation `pot` attributes were set equal (1e19), which is harmless for area-normalised,
shape-only use but would matter for absolute rates.
