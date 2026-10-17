#!/usr/bin/env python3
"""Per-detector, per-channel effect of every reBruce fake-data dial on the reco-energy spectrum, and the near/far double ratio.

For each dial (applied at sigma = +1 on the matched beam-neutrino MC; dirt/off-beam are unaffected and ignored here) and each of SBND,
ICARUS Run 2 and Run 4: f_bin = sum(aw * w) / sum(aw) - 1 in the analysis reco-energy bins of the channel (aw = the XML branch weight, so
cvwgt and the selection are included). The near/far double ratio is (1 + f_SBND) / (1 + f_ICARUS) bin by bin.

  PYTHONPATH=$(root-config --libdir) .venv/bin/python -m fdgen.nearfar
"""
import os
import re

import numpy as np
import uproot
import yaml

from . import dials, events, xmlinfo

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
XML = os.path.join(REPO, "xml/systematics/systematics_suite.xml")
BASE = os.path.join(REPO, "runs/fakedata_data")
SLIM = {"events/sbnd_large": "output/sbnd_run1_fakedata_slim.root", "events/icarus_run2": "output/icarus_run2_fakedata_slim.root",
        "events/icarus_run4": "output/icarus_run4_fakedata_slim.root"}
DET = {"events/sbnd_large": "SBND", "events/icarus_run2": "ICARUS Run 2", "events/icarus_run4": "ICARUS Run 4"}


def channel_edges():
    s = open(XML).read()
    out = {}
    for m in re.finditer(r'<channel name="(\w+)".*?</channel>', s, re.S):
        e = re.search(r'<bins unit="RecoE"\s+edges="([^"]+)"', m.group(0))
        out[m.group(1)] = np.array([float(x) for x in e.group(1).split()])
    return out


def main():
    info = xmlinfo.load(XML)
    edges = channel_edges()
    names = [k[:-6] for k in uproot.open(os.path.join(REPO, SLIM["events/sbnd_large"]))["events/sbnd_large/selected_fakedataTree"].keys() if k.endswith("_sigma")]
    res = {}                                     # (dial, det, channel) -> (cv bins, weighted bins)
    for treedir, slim in SLIM.items():
        f = next(x for x in info["files"] if x["treename"] == f"{treedir}/selected")
        path = os.path.join(BASE, f["filename"])
        ev = events.evaluate(path, f["treename"], f["branches"])
        reco = uproot.open(path)[f["treename"]]["reco_visible_energy"].array(library="np")
        chan = np.array([("ccqelike" if "ccqelike" in f["branches"][b]["subchannel"] else "nproton") if b >= 0 else "" for b in ev["branch"]])
        for name in names:
            w = dials.dial_weight(os.path.join(REPO, slim), treedir, "rebruce", name, 1.0)
            for ch, e in edges.items():
                m = (ev["aw"] != 0) & (chan == ch)
                cv = np.histogram(reco[m], e, weights=ev["aw"][m])[0]
                fk = np.histogram(reco[m], e, weights=(ev["aw"] * w)[m])[0]
                res[(name, DET[treedir], ch)] = (cv, fk)
        print("done", DET[treedir])
    np.save(os.path.join(REPO, "output/fdgen/nearfar/effects.npy"), {"res": res, "edges": edges, "names": names}, allow_pickle=True)


if __name__ == "__main__":
    main()
