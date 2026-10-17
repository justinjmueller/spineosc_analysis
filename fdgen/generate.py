#!/usr/bin/env python3
"""Generate fake-data weight (and optional Poisson count) branches from a YAML scenario file.

  PYTHONPATH=$(root-config --libdir) .venv/bin/python -m fdgen.generate fdgen/configs/batch1.yaml [--entry-stop N] [--dry-run]

For every MCFile of the PROfit MC XML, writes <outdir>/<tree path>_mystery.root with a copy of `selected` and a tree `<tree>_mysteryTree` (entry-aligned with
`selected`, in the same slim file as a copy of it): `mysterywgtNN` = the scenario's per-event ratio to the PROfit prediction (oscillation probability, as
PROfit evaluates it, times the dial weights), and with `poisson: true` also `mysterycntNN`, integer counts.

Poisson counts. PROfit throws each reco bin of the predicted spectrum as Poisson(mu_bin) (PROspec::PoissonVariation).
Here every MC event i gets k_i ~ Poisson(lambda_i), lambda_i = aw_i * pot_scale * ratio_i, the expected number of data
events it represents. A reco bin then holds a sum of independent Poissons, which is Poisson(sum lambda_i) = Poisson(mu_bin)
exactly, for any binning. The random stream is numpy's, not PROfit's mt19937, so the draws differ but the law is the same.
"""
import argparse
import os
import re
import sys

import numpy as np
import uproot
import yaml

from . import dials, events, export, osc, xmlinfo

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def out_name(treename):
    return treename.replace("/", "_") + "_mystery.root"


def normalize_exports(cfg):
    """Merge `export_defaults` into every entry of `exports`: detectors are merged per detector, other keys override; a detector `file` may use
    {name}. Each export needs `name`, `scenario` and `seed`."""
    import copy
    defaults = cfg.get("export_defaults", {}) or {}
    out = []
    for e in cfg.get("exports") or []:
        for need in ("name", "scenario", "seed"):
            if need not in e:
                sys.exit(f"export {e}: missing '{need}'")
        if not isinstance(e["scenario"], str):
            sys.exit(f"export {e['name']}: scenario must be a string id, got {e['scenario']!r} (quote ids like 'null' in YAML)")
        m = copy.deepcopy(defaults)
        for k, v in e.items():
            if k == "detectors":
                for det, spec in v.items():
                    m.setdefault("detectors", {}).setdefault(det, {}).update(copy.deepcopy(spec))
            else:
                m[k] = v
        if not m.get("detectors"):
            sys.exit(f"export {e['name']}: no detectors configured")
        for det, spec in m["detectors"].items():
            for need in ("file", "tree", "cuts"):
                if need not in spec:
                    sys.exit(f"export {e['name']} / {det}: missing '{need}'")
            spec["file"] = spec["file"].format(name=e["name"])
        out.append(m)
    names = [m["name"] for m in out]
    if len(set(names)) != len(names):
        sys.exit("export names must be unique")
    return out


def branch_names(path, treename):
    import ROOT
    f = ROOT.TFile.Open(path)
    names = {b.GetName() for b in f.Get(treename).GetListOfBranches()}
    f.Close()
    return names


def identifiers(texts):
    return {m for t in texts for m in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", t)}


def write_slim(src, treename, out, arrays, formula=None, keep=None):
    """Slim file: a copy of the MC file's `selected` tree and, next to it, `<tree>_mysteryTree` (entry-aligned).

    PROfit reads friend trees from the same file as the main tree (the <friend filename=...> attribute is parsed but
    ignored), so the weights must travel with a copy of the events. `formula` (a TTreeFormula string) keeps only the
    entries it selects; `keep` (a set of branch names) copies only those branches (None: all). `arrays` must already hold
    one row per kept entry.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    treedir, tree = treename.rsplit("/", 1)

    def directory(f):
        d = f
        for part in treedir.split("/"):        # one level at a time: nested mkdir/Snapshot paths land in the top directory
            d = d.GetDirectory(part) or d.mkdir(part)
        return d

    fin = ROOT.TFile.Open(src)
    st = fin.Get(treename)
    if keep is not None:
        present = {b.GetName() for b in st.GetListOfBranches()}
        missing = sorted(keep - present)
        if missing:
            raise ValueError(f"{treename}: branches to keep are not in the tree: {missing}")
        st.SetBranchStatus("*", 0)
        for name in keep:
            st.SetBranchStatus(name, 1)
    fo = ROOT.TFile.Open(out, "recreate")
    directory(fo).cd()
    if formula:
        copy = st.CopyTree(formula)
    else:
        copy = st.CloneTree(-1, "fast" if keep is None else "")
    n_copied = copy.GetEntries()
    copy.Write()
    fo.Close()
    fin.Close()
    n_rows = len(next(iter(arrays.values())))
    if n_copied != n_rows:
        raise RuntimeError(f"{treename}: copied {n_copied} entries but the weight arrays have {n_rows} rows")
    mname = f"{tree}_mysteryTree"
    tmp = out + ".tmp.root"
    if n_rows == 0:       # RDataFrame cannot snapshot empty arrays; build the empty tree directly
        from array import array
        ft = ROOT.TFile.Open(tmp, "recreate")
        t = ROOT.TTree(mname, "")
        bufs = {k: array("d" if v.dtype.kind == "f" else "i", [0]) for k, v in arrays.items()}
        for k, b in bufs.items():
            t.Branch(k, b, f"{k}/{'D' if b.typecode == 'd' else 'I'}")
        t.Write()
        ft.Close()
    else:
        ROOT.RDF.FromNumpy(arrays).Snapshot(mname, tmp, list(arrays))
    ft = ROOT.TFile.Open(tmp)
    fo = ROOT.TFile.Open(out, "update")
    directory(fo).cd()
    ft.Get(mname).CloneTree(-1, "fast").Write()
    fo.Close()
    ft.Close()
    os.remove(tmp)


def project_true_e(w, ev, te_edges):
    """Replace a per-event dial weight by PROfit's representation of it: one ratio per truth-energy bin, separately for every
    XML branch (detector, channel, subchannel), r_j = sum(aw w) / sum(aw) over the events of the bin; events outside the bins keep 1."""
    out = np.ones_like(w)
    idx = osc.le_bin_index(ev["true_e"], te_edges)             # same lower-edge-inclusive rule as PROfit's Binning::Bin
    sel = (ev["aw"] != 0) & (idx >= 0)
    key = ev["branch"][sel].astype(np.int64) * (len(te_edges) + 1) + idx[sel]
    num = np.bincount(key, weights=(ev["aw"] * w)[sel])
    den = np.bincount(key, weights=ev["aw"][sel])
    ratio = np.where(den != 0, num / np.where(den != 0, den, 1), 1.0)
    out[sel] = ratio[key]
    return out


def scenario_ratio(sc, ev, edges, f, entry_stop, rebruce_files, te_edges=None):
    """Per-event ratio to the PROfit prediction for one scenario and one MC file."""
    n = len(ev["aw"])
    r = np.ones(n)
    notes = []
    o = sc.get("oscillation")
    if o:
        p = osc.pmumu_profit(ev["le"], edges, o["dmsq"], o["sinsq2thmm"])
        use = ev["rule"] == 1
        bad = use & np.isnan(p)
        if bad.any():
            notes.append(f"{int(bad.sum())} rule-1 events outside the L/E grid (ratio set to 1)")
        use &= ~np.isnan(p)
        r[use] *= p[use]
    for e in sc.get("effects", []):
        treedir, tree = f["treename"].rsplit("/", 1)
        if tree != "selected":
            continue                  # only the matched beam-neutrino tree has GENIE weights
        path = f["_path"] if e["source"] == "standard" else rebruce_files.get(treedir)
        if path is None:
            continue
        want = f"{treedir}/selected{dials.TREES[e['source']][len('selected'):]}"
        if want not in uproot.open(path):
            continue                  # dirt / files without the friend tree stay at the CV
        w = dials.dial_weight(path, treedir, e["source"], e["dial"], e["sigma"], entry_stop)
        if sc.get("project_true_e"):
            w = project_true_e(w, ev, te_edges)
        r *= w
    return r, notes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config")
    ap.add_argument("--entry-stop", type=int, help="only the first N events of every file (testing)")
    ap.add_argument("--dry-run", action="store_true", help="evaluate and report expected yields, write nothing")
    ap.add_argument("--only", help="only MC files whose tree path contains this string")
    args = ap.parse_args(argv)
    cfg = yaml.safe_load(open(args.config))
    ds, scs = cfg["dataset"], cfg["scenarios"]
    base = os.path.join(REPO, ds.get("base_dir", "."))
    xml = os.path.join(REPO, ds["mc_xml"])
    info = xmlinfo.load(xml)
    edges = osc.le_bin_edges(*info["le_bins"])
    te_edges = osc.le_bin_edges(*info["te_bins"])      # truth-energy bins of the dial splines (binning=var2)
    outdir = os.path.join(REPO, ds["outdir"])
    keydir = os.path.join(REPO, ds.get("keydir", ds["outdir"]))        # KEY.yaml and manifest.yaml (private) can live elsewhere
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(keydir, exist_ok=True)
    poisson = bool(ds.get("poisson", False))
    out_cfg = ds.get("output", {}) or {}                       # selection / branches / keep_extra, see the README
    sel_cfg = out_cfg.get("selection")
    # reBruce dials live in the spinefd slim files, not in the CVMC file: tree directory -> file
    rebruce = {k: os.path.join(REPO, v) for k, v in ds.get("rebruce_files", {}).items()}
    ids = [s["id"] for s in scs]
    if len(set(ids)) != len(ids):
        sys.exit("scenario ids must be unique")
    order = np.random.default_rng(ds["seed"]).permutation(len(scs)) if ds.get("shuffle") else np.arange(len(scs))
    index = {sid: int(order[k]) for k, sid in enumerate(ids)}    # scenario id -> output index

    exports = normalize_exports(cfg)
    cut_cache, col_cache = {}, {}                          # (path, tree, cut) -> mask ; (path, tree, keep) -> columns, shared by all exports
    written_exports = set()
    export_rows = {}                                        # (export name, detector) -> list of per-file row dicts
    export_obs = {}                                         # (export name, detector) -> {"expected": x, "events": n}
    write_slim_files = bool(ds.get("write_slim", True))
    manifest = dict(name=ds["name"], mc_xml=ds["mc_xml"], poisson=poisson, n_scenarios=len(scs),
                    l_over_e_grid=list(info["le_bins"]), files=[])
    key = dict(name=ds["name"], seed=ds["seed"], scenarios=[dict(index=index[s["id"]], **s) for s in scs], expected={})
    for fi, f in enumerate(info["files"]):
        if args.only and args.only not in f["treename"]:
            continue
        f["_path"] = os.path.join(base, f["filename"])
        scale = xmlinfo.pot_scale(info, f)
        formula = None
        if sel_cfg:
            formula = " || ".join(f"(({b['weight']}) != 0)" for b in f["branches"]) if sel_cfg == "xml" else sel_cfg
        ev = events.evaluate(f["_path"], f["treename"], f["branches"], formula)
        n = ev["n"] if args.entry_stop is None else min(ev["n"], args.entry_stop)
        ev = {k: (v[:n] if isinstance(v, np.ndarray) else n) for k, v in ev.items()}
        if (ev["aw"] < 0).any() and poisson:
            sys.exit(f"{f['treename']}: negative event weights cannot be Poisson-thrown")
        mask = ev.get("sel")
        if mask is not None:
            lost = int(((ev["aw"] != 0) & ~mask).sum())
            if lost:
                print(f"  [{f['treename']}] WARNING: the selection drops {lost} events that are in an analysis branch "
                      f"(weight {float(np.abs(ev['aw'][(ev['aw'] != 0) & ~mask]).sum()):.4g})")
        keep = None
        if out_cfg.get("branches", "all") == "needed":
            # every expression PROfit evaluates for these branches: the weight and all <variable>/<variable2D> expressions (fitting variables)
            texts = [b["weight"] for b in f["branches"]] + [v for b in f["branches"] for v in b["variables"]] + [formula or ""]
            present = branch_names(f["_path"], f["treename"])
            # always kept when the tree has them: the proton multiplicity and reco_pn (the second fit variable of the 2D analysis),
            # so a 2D fit stays possible from a dataset generated with the 1D XML; keep_extra is mandatory (missing -> error)
            default = {"reco_proton_multiplicity", "reco_pn"} & present
            keep = (identifiers(texts) & present) | default | set(out_cfg.get("keep_extra", []))
        out = {}
        for sc in scs:
            k = index[sc["id"]]
            r, notes = scenario_ratio(sc, ev, edges, f, args.entry_stop, rebruce, te_edges)
            out[f"mysterywgt{k:02d}"] = r
            lam = ev["aw"] * scale * r
            key["expected"].setdefault(sc["id"], {})[f["treename"]] = float(lam.sum())
            if poisson:
                rng = np.random.default_rng(np.random.SeedSequence([ds["seed"], k, fi]))
                out[f"mysterycnt{k:02d}"] = rng.poisson(np.clip(lam, 0, None)).astype(np.int32)
            for note in notes:
                print(f"  [{f['treename']}] {sc['id']}: {note}")
        det = f["branches"][0]["detector"]
        for ex in ([] if args.dry_run else exports):          # a dry run writes and throws nothing
            if det not in ex["detectors"]:
                continue
            if args.entry_stop is not None:
                sys.exit("exports need the full files; do not use --entry-stop")
            sc_ids = {s_["id"]: s_ for s_ in scs}
            if ex["scenario"] not in sc_ids:
                sys.exit(f"export {ex['name']}: unknown scenario {ex['scenario']}")
            kx = index[ex["scenario"]]
            lam_x = ev["aw"] * scale * out[f"mysterywgt{kx:02d}"]
            cuts = ex["detectors"][det].get("cuts", {})
            if f["treename"] not in cuts:
                sys.exit(f"export {ex['name']}: no cut configured for {f['treename']} (python -m fdgen.export --suggest-cuts <xml> prints them)")
            ck = (f["_path"], f["treename"], cuts[f["treename"]])
            if ck not in cut_cache:
                cut_cache[ck] = events.evaluate_formula(*ck)
            passed = cut_cache[ck][:n]
            rng_x = np.random.default_rng(np.random.SeedSequence([ex["seed"], kx, fi, 4242]))
            cnt = rng_x.poisson(np.clip(lam_x, 0, None))
            bad = int(((cnt > 0) & ~passed).sum())
            if bad:
                sys.exit(f"export {ex['name']} [{f['treename']}]: {bad} thrown events fail the configured cuts (the cuts are stricter than the XML selection)")
            idx = np.repeat(np.nonzero(cnt > 0)[0], cnt[cnt > 0])
            keep = ex.get("keep", ["reco_visible_energy", "reco_pn", "reco_proton_multiplicity"])
            kk = (f["_path"], f["treename"], tuple(keep))
            if kk not in col_cache:
                col_cache[kk] = uproot.open(f["_path"])[f["treename"]].arrays(keep, library="np")
            arr = col_cache[kk]
            export_rows.setdefault((ex["name"], det), []).append({k: np.asarray(arr[k])[idx] for k in keep})
            o = export_obs.setdefault((ex["name"], det), dict(expected=0.0, events=0))
            o["expected"] += float(lam_x.sum())
            o["events"] += int(len(idx))
        path = os.path.join(outdir, out_name(f["treename"]))
        if not args.dry_run and write_slim_files:
            if args.entry_stop is not None:
                sys.exit("--entry-stop gives truncated weights that cannot be aligned with the copied tree; use --dry-run")
            rows = {k: (v[mask] if mask is not None else v) for k, v in out.items()}
            write_slim(f["_path"], f["treename"], path, rows, formula, keep)
        manifest["files"].append(dict(tree=f["treename"], file=out_name(f["treename"]), entries=int(n),
                                      mc_pot=f["mc_pot"], pot_scale=float(scale)))
        print(f"{f['treename']}: {n} events, scale {scale:.4f}, selected {int((ev['rule'] >= 0).sum())}, -> {path}")
    for ex in ([] if args.dry_run else exports):
        rng_p = np.random.default_rng(np.random.SeedSequence([ex["seed"], 99]))
        key.setdefault("exports", {})[ex["name"]] = dict(scenario=ex["scenario"], seed=ex["seed"], detectors={})
        for det, spec in ex["detectors"].items():
            parts = export_rows[(ex["name"], det)]
            cols = {k: np.concatenate([p_[k] for p_ in parts]) for k in parts[0]}
            perm = rng_p.permutation(len(next(iter(cols.values()))))     # merge the MC components and shuffle: no component order is visible
            cols = {k: v[perm] for k, v in cols.items()}
            jit = float(ex.get("jitter_duplicates", 1e-6))
            if jit > 0:        # an MC event thrown k > 1 times would otherwise appear as k identical rows, which real data never has
                cols = export.jitter_duplicates(cols, jit, rng_p)
            export.write_tree(os.path.join(outdir, spec["file"]), spec["tree"], cols)
            o = export_obs[(ex["name"], det)]
            key["exports"][ex["name"]]["detectors"][det] = dict(expected=o["expected"], events=o["events"])
            print(f"export {ex['name']} {det}: {o['events']} events (expected {o['expected']:.1f}) -> {spec['file']}")
        xml_text = open(xml).read()
        with open(os.path.join(outdir, ex.get("xml", ex["name"] + ".xml")), "w") as fx:
            fx.write(export.build_xml(xml_text, info, ex, info["det_pot"], outdir))
        written_exports.update([spec["file"] for spec in ex["detectors"].values()] + [ex.get("xml", ex["name"] + ".xml")])
    if exports and not args.dry_run:
        stale = sorted(set(os.listdir(outdir)) - written_exports)
        if stale:
            print(f"WARNING: {outdir} contains files this run did not write (left over from an earlier config?): {stale[:6]}"
                  f"{' ...' if len(stale) > 6 else ''}\n         Remove them before sharing the directory.")
    if args.dry_run:
        print(yaml.safe_dump(key["expected"]))
        return
    yaml.safe_dump(manifest, open(os.path.join(keydir, "manifest.yaml"), "w"), sort_keys=False)
    yaml.safe_dump(key, open(os.path.join(keydir, "KEY.yaml"), "w"), sort_keys=False)


if __name__ == "__main__":
    main()
