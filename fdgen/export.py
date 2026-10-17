"""Truth-free export: Poisson-thrown data events with reco variables only, in a data-like layout, plus the data XML to read them.

The export is built from the same per-event expectation as the weight datasets: lambda_i = aw_i * pot_scale * ratio_i (the XML branch
weight with the scenario's ratio and PROfit's POT factor), thrown as k_i ~ Poisson(lambda_i). Every event with k_i > 0 is written k_i times
with only the `keep` branches (reco-level), after the configured `cuts` (the selection, per MC tree) were required to hold. Components
(beam MC, dirt, off-beam, non-matched) of one detector are merged and shuffled, so nothing says which component an event came from.
"""
import re

import numpy as np

TRUTH = re.compile(r"\btrue_[A-Za-z0-9_]*|\bcvwgt\b")
# identifiers that look like truth but are file-level constants of the sample (detector id); allowed in a cut
ALLOWED = {"true_detector"}


def top_level_split(expr, sep="&&"):
    """Split on `sep` outside parentheses."""
    parts, depth, cur, i = [], 0, "", 0
    while i < len(expr):
        c = expr[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if depth == 0 and expr.startswith(sep, i):
            parts.append(cur.strip())
            cur, i = "", i + len(sep)
            continue
        cur += c
        i += 1
    parts.append(cur.strip())
    return [p for p in parts if p]


def unwrap(expr):
    expr = expr.strip()
    while expr.startswith("(") and expr.endswith(")"):
        depth = 0
        for i, c in enumerate(expr):
            depth += (c == "(") - (c == ")")
            if depth == 0 and i < len(expr) - 1:
                return expr
        expr = expr[1:-1].strip()
    return expr


def reco_cut(weight):
    """The reco-only selection of an XML branch weight: its top-level && terms without truth variables (true_detector is allowed), without
    the proton-multiplicity term (that defines the channel) and without the trailing `* cvwgt`."""
    w = weight.strip()
    w = re.sub(r"\*\s*cvwgt\s*$", "", w).strip()
    terms = top_level_split(unwrap(w))
    keep = []
    for t in terms:
        names = set(TRUTH.findall(t)) - ALLOWED
        if names or "reco_proton_multiplicity" in t:
            continue
        keep.append(t)
    return " && ".join(keep)


def channel_of(subchannel):
    return re.search(r"_(ccqelike|nproton|inclusive)_", subchannel).group(1)


def suggest_cuts(info):
    """tree path -> reco cut, taken from the XML; raises if the branches of one tree disagree."""
    out = {}
    for f in info["files"]:
        cuts = {reco_cut(b["weight"]) for b in f["branches"]}
        if len(cuts) != 1:
            raise ValueError(f"{f['treename']}: branches have different reco cuts: {cuts}")
        out[f["treename"]] = cuts.pop()
    return out


def channel_cut(info, f, channel):
    """The proton-multiplicity term of the channel in this file (e.g. '(reco_proton_multiplicity == 1)')."""
    for b in f["branches"]:
        if channel_of(b["subchannel"]) == channel:
            for t in top_level_split(unwrap(re.sub(r"\*\s*cvwgt\s*$", "", b["weight"].strip()).strip())):
                if "reco_proton_multiplicity" in t:
                    return t
    raise ValueError(f"{f['treename']}: no multiplicity term for channel {channel}")


def blind_variable(expr):
    """Truth-based variables become the constant 0, like the real-data XMLs do."""
    return "0" if TRUTH.search(expr) else expr


def write_tree(path, treename, cols):
    """Write numpy columns as a TTree at `treename` (e.g. events/sbnd_bnblight/selected) with ROOT itself: uproot's nested directories
    segfault when ROOT reads them, so the directories are made one level at a time and the tree is cloned into place."""
    import os
    import ROOT
    ROOT.gROOT.SetBatch(True)
    treedir, tree = treename.rsplit("/", 1)
    tmp = path + ".tmp.root"
    ROOT.RDF.FromNumpy({k: np.ascontiguousarray(v) for k, v in cols.items()}).Snapshot(tree, tmp, list(cols))
    ft = ROOT.TFile.Open(tmp)
    fo = ROOT.TFile.Open(path, "recreate")
    d = fo
    for part in treedir.split("/"):
        d = d.GetDirectory(part) or d.mkdir(part)
    d.cd()
    ft.Get(tree).CloneTree(-1, "fast").Write()
    fo.Close()
    ft.Close()
    os.remove(tmp)


def jitter_duplicates(cols, rel, rng):
    """Multiply the continuous columns of every repeated row (after its first copy) by 1 + rel * U(-1, 1). Integer-valued columns (e.g. the
    proton multiplicity) are untouched. rel = 1e-6 is far below any bin width, so the spectra do not change."""
    names = sorted(cols)
    stacked = np.stack([cols[k].astype(np.float64) for k in names], axis=1)
    _, first, inv = np.unique(stacked, axis=0, return_index=True, return_inverse=True)
    inv = np.ravel(inv)
    is_first = np.zeros(len(stacked), dtype=bool)
    is_first[first] = True
    dup = ~is_first
    out = {k: v.copy() for k, v in cols.items()}
    for k in names:
        v = cols[k]
        if v.dtype.kind == "f" and not np.all(v == np.round(v)):
            out[k] = np.where(dup, v * (1 + rel * rng.uniform(-1, 1, len(v))), v).astype(v.dtype)
    return out


def build_xml(xml_text, info, export, detector_pots, outdir):
    """The blind data XML for one export: same <mode>, <detector>, <model> and channels as the MC XML (every channel gets the single subchannel
    "data"), and one MCFile per detector reading the exported events. Weights are the channel's proton-multiplicity term alone (the cuts were
    already applied when the events were written); truth-based variables are the constant 0, as in the repo's real-data XMLs; the MCFile POT
    is the detector's target POT, so PROfit's POT factor is 1 and every event counts once."""
    import os
    s = xml_text
    out = ['<?xml version="1.0" ?>\n',
           f"<!--\nPseudo-data events, reco variables only (fdgen truth-free export '{export['name']}').\n"
           "The selection is already applied; no systematics, no real beam data.\n-->\n",
           re.search(r'<mode name="[^"]*"\s*/>', s).group(0) + "\n"]
    out += [m.group(0) + "\n" for m in re.finditer(r'<detector name="[^"]*"\s+pot="[^"]*"\s*/>', s)]
    out.append(re.search(r'<model tag=.*?</model>', s, re.S).group(0) + "\n")
    for m in re.finditer(r'<channel name=.*?</channel>', s, re.S):
        block = re.sub(r'\s*<subchannel [^>]*/>', "", m.group(0))
        out.append(block.replace("</channel>", '    <subchannel name="data" plotname="Data" color="#000000"/>\n\n</channel>') + "\n")
    for det, spec in export["detectors"].items():
        files = [f for f in info["files"] if f["branches"] and f["branches"][0]["detector"] == det]
        ref = next(f for f in files if f["treename"].endswith("/selected") and "dirt" not in f["treename"] and "offbeam" not in f["treename"])
        new = [f'<MCFile treename="{spec["tree"]}" filename="{os.path.join(os.path.abspath(outdir), spec["file"])}" scale = "1.0" pot="{detector_pots[det]:.6e}">']
        for channel in sorted({channel_of(b["subchannel"]) for b in ref["branches"]}):
            b = next(b for b in ref["branches"] if channel_of(b["subchannel"]) == channel)
            tags = []
            for v in b["variables"]:
                tags.append(f"        <variable{'2D' if ';' in v else ''}>{blind_variable(v)}</variable{'2D' if ';' in v else ''}>")
            new.append(f'    <branch\n        associated_subchannel = "nu_{det}_{channel}_data"\n        model_rule            = "0"\n'
                       f'        additional_weight     = "{channel_cut(info, ref, channel)}"\n    >\n' + "\n".join(tags) + "\n    </branch>")
        new.append("</MCFile>\n")
        out.append("\n".join(new))
    return "\n".join(out)


def main():
    import argparse
    import yaml
    from . import xmlinfo
    ap = argparse.ArgumentParser(description="print the per-tree reco cuts of an analysis XML as YAML (to paste into an export's `cuts`)")
    ap.add_argument("--suggest-cuts", metavar="XML", required=True)
    a = ap.parse_args()
    print(yaml.safe_dump(suggest_cuts(xmlinfo.load(a.suggest_cuts)), sort_keys=False, width=1000))


if __name__ == "__main__":
    main()
