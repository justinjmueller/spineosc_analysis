# spineosc

PROfit XML configs for the SBND + ICARUS oscillation analysis.

## Repo layout

```
xml/                   PROfit XML configs (contours_base.xml, exclusive_selected_1d.xml, varset*_mc.xml, varset*_devsample*.xml)
input/
  sources.list.example  checked-in template listing every file key the XML configs expect
  sources.list           your machine's actual file paths (gitignored - you create this)
scripts/
  deploy.sh              populates work/ from input/sources.list
work/                   symlinks the XML configs read from (generated, gitignored)
detvar/detsys.root     detector systematics histograms (tracked in git)
```

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

3. Run PROfit **from the repo root** - the XML configs use paths relative
   to `work/`, so the working directory has to be the directory that
   contains `work/` (and `xml/`) for them to resolve:

   ```
   /path/to/PROfit -x xml/contours_base.xml -v 1 -w 4 --log log.txt plot
   ```

`input/sources.list` and `work/` are both gitignored - they're
machine-specific and fully reproducible from `input/sources.list.example`
plus your local data, so there's nothing to commit.

## File keys

`input/sources.list` uses the key pattern `<detector>_<runX>_<type>`:

- Detector: `sbnd` or `icarus`
- Run: `run1` for SBND (it only ever ran in one configuration), `run2` or
  `run4` for ICARUS
- Type: `cvmc` (central-value signal MC, "wsyst"), `dirt` (dirt MC),
  or `data` (beam-on data, "bnblight")

Plus one detector-independent key, `detsys_main`, for the ROOT file
holding the detector-systematics histograms referenced by
`<HistVarFiles>` and the `hist1d` allowlist entries in the XML configs.

That's 9 detector/run/type keys + `detsys_main` = 10 lines in
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
