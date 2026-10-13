# schematics

Build PROfit XML configurations from blueprints, instead of copying a
1263-line file and editing four values in it.

A **template** is an existing XML with the values that vary between siblings
replaced by `{{placeholder}}` markers. A **blueprint** is a short TOML file
naming a template and supplying those values. `build.py` renders one into the
other.

```
schematics/
  build.py                          the generator
  templates/coverage.xml            102 lines: the study-specific structure
  partials/systematics_standard.xml 210 lines: the allowlists, shared by 16 XMLs today
  partials/detsys_histvarfiles.xml   62 lines: the HistVarSections, shared by 22
  partials/mcfiles_coverage.xml     892 lines: the MC files and their branches
  blueprints/coverage/*.toml         13 files, 4-5 values each
```

Templates compose two directives besides `{{placeholder}}`:

- `{{> partial}}` splices in a shared fragment, so the 208-line systematics
  block lives in one file rather than once per template.
- `{{#if SBND}} ... {{/if}}` keeps a region only when that detector is in the
  blueprint's `detectors` list. `ICARUS` is derived, for prose that belongs to
  the far detector generally rather than one of its runs.

Kept regions are emitted **verbatim**. That is deliberate: the XMLs carry
hand-tuned column alignment that no formatting rule reproduces, so the
generator never rewrites a line it keeps. Conditionals do not nest and there
is no `else`; anything needing more wants the raw escape hatch instead.

## The criterion: exact reproduction

The generator is only worth trusting if it reproduces what is already there,
so that is the test:

```
python3 schematics/build.py --check
```

renders every blueprint and compares it **byte-for-byte** with the committed
XML, printing a diff for any that disagree. It cannot quietly reformat,
reorder, or drop a comment: any such change is a failure.

It currently reports `13/13 reproduced exactly`.

Nothing under `xml/` is ever written. `--out DIR` renders into a directory of
your choosing; generating into the tree is a separate decision from proving
that the generator can.

## What one blueprint looks like

Thirteen coverage studies differ in four values. All of `dpt.xml` is:

```toml
template = "coverage"

[values]
study         = "dpt"
channel_label = "Reconstructed $\\delta p_{T}$ [MeV/c]"
fit_variable  = "reco_dpT"
fit_bins      = "<bins unit=\"dpT [MeV/c]\" min=\"0.0\" max=\"1000.0\" nbins=\"20\"/>"
```

Note `fit_variable`: `dp.xml` fits `reco_pn`, not `reco_dp`, so the branch
expression is a real parameter rather than something derivable from the study
name. That is the kind of detail that makes a naming convention a poor
substitute for a declaration.

A blueprint that names a value the template does not use, or omits one it
does, is an error rather than a silently half-rendered XML.

## What the check already found

`muon_primary_score.xml`, `proton_primary_score.xml` and `proton_score.xml`
each say

> (matching **muon_score_data.xml**'s `<detector>` block exactly)

in their header, where the other ten name their own data companion. It is a
stale reference left over from the file each was copied from -- harmless,
since it is a comment, but exactly the drift this directory exists to stop.

The committed XMLs are left alone, so those three blueprints carry an explicit
`detector_block_reference` override with a note. Correct the XMLs and the
override lines can be deleted; the check will confirm it.

## Not covered yet

The coverage family has three shapes, and only the largest is templated:

| shape | files | status |
|---|---|---|
| two-detector, 1263 lines | 13 | templated, exact |
| ICARUS-only, 932 lines | 6 | partials already render them; 12 lines of header prose to go |
| SBND-only, 631 lines | 6 | same, 14 lines |
| `zoom_icarus`, 915 lines | 1 | not yet |

The partials were built against all three shapes, so the detector
conditionals are already exercised: rendering them with
`detectors = ["ICARUSRun2", "ICARUSRun4"]` or `["SBND"]` reproduces the
systematics, HistVarFiles and MCFile blocks of the single-detector files
exactly. What remains for those twelve is a per-file header note explaining
the single-detector derivation, plus a couple of prose lines -- a uniform
12-line (ICARUS) or 14-line (SBND) delta. It needs one more idea than the
current directives have: a blueprint-declared flag, or a second template.

Also untouched: the `*_data.xml` companions, the `*_plot.ini` files, and the
`xml/contours/` family. The contours XMLs are the more valuable target in the
long run and the harder one -- see below.

## Where this should go next

The text substitution above is the cheap half. The expensive bugs in this
repo's XMLs are *index couplings* that no amount of copy-paste discipline
catches, and a generator is the natural place to eliminate them:

- the model needs true L/E at `var1` (`variable_index="1"`)
- the cross-section splines read `binning="var2"` (true energy)
- the `hist1d` detector systematics read `xvar=` the 1D reco-energy axis,
  which is `var3` in `contours_base.xml` and `var0` in `contours_recoe.xml`
- PROfit reads `plot="false"` flags by variable index off the **first**
  channel, so the same flag on a later channel does nothing

A blueprint that declared variables by role (`fit`, `le`, `truee`) and let the
generator assign the indices would make "switch the fit variable from 2D to
1D" a one-line change, rather than the 54-edit transform it currently is.

Worth enforcing there too, since the generator sees the whole configuration:

- every detector/channel/subchannel block is filled, which is what the
  `empty_bin_check` abort is about
- every `apply_to_subchannel` and norm pattern matches at least one fullname,
  **case-sensitively** -- the trap that would silently disable the CORSIKA
  norm if a detector were renamed
- every `hist1d` name resolves to a `HistVarSection` or a histogram that
  exists in the referenced detsys ROOT file
- no `--` inside an XML comment, which tinyxml2 rejects

One caveat worth knowing before generating into the tree: `PROconfig::CalcHash`
hashes fullnames, bin edges, weight names and branch variables -- not comments
or whitespace. A regeneration that only reformats will therefore **not**
invalidate anyone's multi-GB caches, while one that changes binning or weights
will. That is the desired behaviour, but it does mean the byte-for-byte check
here is stricter than PROfit itself would be.
