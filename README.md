# spineosc

PROfit XML configs for the SBND + ICARUS oscillation analysis.

## Repo layout

| Path | Purpose | Entry point / docs |
|---|---|---|
| `xml/` | PROfit XML and `.ini` configs, grouped by study (`contours`, `coverage`, `detector_data_access`, `dent_study`, `resolution`, `fakedata`, ...) | `xml/MANIFEST.toml`, `python scripts/xml_manifest.py --check` |
| `jobs/` | per-study job configs read by PROfit `--config` | see "Job configs" |
| `input/` | `sources.list.example` / `site.conf.example` templates and the machine file maps | see "First-time setup" |
| `scripts/` | deploy, XML generators and linter, notebook runner, study drivers (`scripts/dent_study/`) | `scripts/deploy.sh` |
| `plotting/studies/` | documented plotting scripts that replace the analysis notebooks | `plotting/studies/README.md` |
| `plotting/`, `detsys/plotting/` | legacy analysis notebooks, kept until the scripts are verified | `docs/RESIDUAL_TASKS.md` |
| `detsys/`, `detvar/` | detector-systematics calculator and the ratio-histogram ROOT files | `detsys/calculator*/` |
| `reweight/`, `reBruce/` | event reweighting (`reweight/Makefile`) and the `reBruce` submodule | submodule README |
| `spinefd/` | fake-data generation on the CV MC | `spinefd/README.md` |
| `docs/` | study reports, audits and `RESIDUAL_TASKS.md` (open work) | `docs/RESIDUAL_TASKS.md` |
| `runs/`, `figures/`, `work/` | generated run products, figures and symlinks (git-ignored) | |

Python environment: see `pyproject.toml` (Python >= 3.11; the study scripts also need the `mutools`
package, installed from `/Users/mueller/Projects/mutools` with `uv sync`).

Note: the PROfit wrappers `scripts/profit_run.sh` and `scripts/submit_pbs.sh` described below
are not present on this branch (they live in the commit history of another branch); the PROfit
build and `input/site.conf` settings are described in "Running PROfit".

The XML configs never hardcode an absolute, machine-specific path. Every
`filename="..."` attribute in `xml/*.xml` points at a relative
`work/<key>.root`, and `work/` is a directory of symlinks that
`scripts/deploy.sh` builds for you from `input/sources.list`. This is what
makes the same XML checkout runnable on your laptop, a collaborator's
machine, or a batch node like ALCF Aurora - only `input/sources.list`
differs per machine.

## First-time setup

1. Copy the template and point it at your local copies of the data files:

   ```
   cp input/sources.list.example input/sources.list
   ```

   Edit `input/sources.list`, replacing each `/path/to/...` with the real
   absolute path on your machine. See "File keys" below for what each line
   means.

2. Populate `work/`:

   ```
   scripts/deploy.sh
   ```

   This creates `work/<key>.root` -> `<path from sources.list>` for every
   line in the file, reports any source path that doesn't exist (and exits
   non-zero if so), and prunes any leftover symlink whose key is no longer
   in `sources.list`. Re-run it any time you edit `input/sources.list`.

3. Tell the scripts where your PROfit build and output area live:

   ```
   cp input/site.conf.example input/site.conf
   ```

   Edit it to point `profit_bin` at the directory containing the `PROfit`
   binary and `output_root` at a directory with room for the caches (a
   single tag's `.bin` files run to several GB, so keep this off the repo
   and off `$HOME`). The `pbs_*` settings only matter if you submit batch
   jobs.

`input/sources.list`, `input/site.conf` and `work/` are all gitignored -
they're machine-specific and fully reproducible from the checked-in
`.example` templates plus your local data, so there's nothing to commit.

## Running PROfit

Use `scripts/profit_run.sh`; it handles the working-directory problem for
you. The XML configs use paths relative to `work/`, but PROfit writes
multi-GB caches and plots that must not land in the repo, so the script
creates a run directory under `output_root` containing `work/` and `xml/`
symlinks back to the repo and runs from there. Relative paths resolve,
outputs stay out of the checkout.

```
scripts/profit_run.sh <jobs/NAME.ini> <subcommand> [extra PROfit args...]
```

**Build the cache.** This is the expensive step - it reads every MC file
the XML references and writes `<tag>_prop.bin` and `<tag>_syst.bin`:

```
scripts/profit_run.sh jobs/varset1_mc.ini process
```

For the full varset1 MC that takes about an hour, so it belongs in a batch
job rather than on a login node:

```
scripts/submit_pbs.sh jobs/varset1_mc.ini process 06:00:00 capacity
```

Mind the queue limits: `debug` allows 1-2 nodes and caps walltime at one
hour (and one job per user at a time), so anything longer needs `capacity`,
which allows up to seven days. The scheduler rejects an over-long request
outright rather than trimming it.

**Everything after that is cheap**, because it loads the cache instead of
rebuilding it:

```
scripts/profit_run.sh jobs/varset1_mc.ini plot
scripts/profit_run.sh jobs/varset1_mc.ini plot --with-splines
```

Add `--dry-run` to any invocation to run every check and print the command
without executing it - worth using on a login node, where a stray `process`
would start reading tens of GB:

```
scripts/profit_run.sh jobs/varset1_mc.ini process --dry-run
```

The run directory is keyed on the config's `tag`, not on the `.ini`
filename, because PROfit names its caches after the tag. Two job configs
sharing a tag therefore share one cache - which is how an MC-only and a
data/MC config can reuse the same hour of processing.

## Job configs (`jobs/*.ini`)

PROfit reads these natively via `--config`; there is no translation layer.
Keys are PROfit's own long option names without the leading dashes, and a
boolean is `key = true`:

```
xml = xml/varset1_mc.xml
tag = varset1_mc
nthread = 104
verbosity = 3
progress = true
```

Three things about this format will cost you time if you don't know them:

1. **Unrecognised keys are silently ignored.** A typo disables that setting
   with no warning at all. Check names against `PROfit --help`.

2. **Placement is significant and equally silent.** Everything listed under
   `Options:` in `PROfit --help` is *global* and belongs in the top-level
   block - including `max`, `output`, `scale-by-width`, `shapeonly` and
   `fit-variable`. A global key placed inside a `[section]` does nothing.

3. **A `[section]` activates that subcommand.** A config carrying both
   `[process]` and `[plot]` runs *both*, regardless of which one you name on
   the command line. Keep these files to global options only and pass
   subcommand options on the command line instead.

## Caches and when they go stale

PROfit embeds a hash of the XML in `<tag>_prop.bin` / `<tag>_syst.bin` and
refuses to load them if the XML has since changed. That check is doing real
work, so prefer bumping `tag` or re-running `process` over overriding it.

`--force` (or `force = true`) overrides the check. It is only safe when you
know the change cannot affect the cached content - a `plotname` or a
subchannel `color`, say. It is *not* safe for anything that enters the
per-event weights. In particular `<detector pot="...">` is folded into the
event weights at `process` time, so changing a POT genuinely requires a
re-`process`; forcing past it yields silently wrong normalization.

## File keys

`input/sources.list` uses the key pattern `<detector>_<runX>_<type>`:

- Detector: `sbnd` or `icarus`
- Run: `run1` for SBND (it only ever ran in one configuration), `run2` or
  `run4` for ICARUS
- Type: `cvmc` (central-value signal MC, "wsyst"), `dirt` (dirt MC),
  or `data` (beam-on data, "bnblight")

Plus a set of detector-independent `detsys_*` keys, one per ROOT file of
detector-systematics histograms referenced by `<HistVarFiles>` and the
`hist1d` allowlist entries in the XML configs:

- `detsys_main` - the histograms binned in reconstructed visible energy,
  used by the main analysis XMLs.
- One `detsys_<variable>` per coverage study under `xml/coverage/` -
  the same 27 systematics re-binned in that study's fit variable.
  Currently the PID scores (`muon_score`, `proton_score`,
  `muon_primary_score`, `proton_primary_score`), the range/momentum
  agreement variables (`muon_calo_csda`, `proton_calo_csda`,
  `muon_mcs_csda`), and the kinematic/TKI variables (`proton_ke`, `dp`,
  `dpt`, `dalphat`, `dphit`, `muon_azimuth`).

  The kinematic set is binned in the ntuples' native units - MeV, MeV/c
  and radians - rather than the GeV convention used by
  `xml/varset1_devsample.xml`. The detvar calculator hands the subchannel
  branch straight to RDataFrame as a column name, so a scaling expression
  like `reco_pn/1000` cannot be used there; keeping the XML channel
  binning in those same native units is what guarantees the fit variable
  and the systematic ratio axis agree.

Each is produced by `detsys/calculator` and then combined with e.g.

```
hadd detvar/detsys_v1.0.10_muon_score.root \
     detsys/profit_systematics/*_by_muon_score.root
```

That's 9 detector/run/type keys + 14 `detsys_*` keys = 23 lines in
`input/sources.list`.

**ICARUS Run 2 and Run 4 read out of the same underlying production
file** for a given type - `icarus_run2_dirt` and `icarus_run4_dirt`, for
instance, point at the identical dirt-MC file. This isn't a mistake:
PROfit needs each run to see its own distinct file path even when the
bytes are the same (it doesn't track run vs. content), so
`input/sources.list` aliases both run-specific keys to the one physical
file rather than relying on a manually-created "run4alias" symlink
sitting in the data directory.

## Updating after a data location changes

Edit `input/sources.list` and re-run `scripts/deploy.sh`. It's idempotent
and self-cleaning - safe to re-run any time, and it removes any
`work/*.root` symlink whose key you've since renamed or deleted so `work/`
never accumulates stale entries.

## Adding a new file key

If a new XML config introduces a `filename="work/some_new_key.root"`
reference:

1. Add `some_new_key=/path/to/the/file.root` to both
   `input/sources.list.example` (so collaborators know it's needed) and
   your own `input/sources.list`.
2. Run `scripts/deploy.sh`.

## Troubleshooting

**Segfault with no useful message.** Almost always a missing
`work/<key>.root` - PROfit crashes rather than reporting it. `profit_run.sh`
checks for this up front and points you at `scripts/deploy.sh`; you only see
the raw crash if you invoke PROfit directly.

**`Require data and MC to have same channels`.** The MC XML and the data XML
given as `data = ...` must declare identical `<detector pot="...">` values.
The `varset*_devsample*.xml` files carry the dev-sample exposure while the
`varset*_mc.xml` files carry full-exposure targets, so pairing them as-is
fails this check. Because POT is baked into the cache at `process` time, the
two exposures need separate configs and separate caches - not a runtime flag.
`--scale` does not help here: it scales subchannels in the loaded MC and
never touches the declared detector POT.

**`Invalid hex color format`.** Subchannel colours must be `#RRGGBB`;
eight-digit `#RRGGBBAA` is rejected. The failure surfaces at plot time, well
after processing has finished.

**`-m/--max` does nothing.** It appears in `PROfit --help` but is dead code
in the current build - parsed and never read. To truncate an event loop, use
the per-`<MCFile>` attributes `maxevents` or `partial_load_frac` instead.
