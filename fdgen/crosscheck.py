#!/usr/bin/env python3
"""Compare the by-hand fdgen fake data with PROfit's own injection, bin by bin, from `plot` histogram files.

For every scenario of fdgen/configs/legacy_1d/crosscheck_dials.yaml, reads runs/fakedata_data/fdgen_xc{P,H}_<scenario>_xc{P,H}_<scenario>_PROplot_histograms.root
and compares the total fake-data histogram (prefix DATA, all subchannels) per detector, channel and variable.
Reports the largest |H/P - 1| per bin and the total-yield difference, and, relative to the injected effect (P/CV - 1), the
fraction of it that the by-hand data reproduces.
"""
import os

import numpy as np
import uproot
import yaml

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RUN = os.path.join(REPO, "runs", "fakedata_data")
TOTAL = 4294967295


def load(path):
    h = uproot.open(path)["hist1d"].arrays(library="np")
    h["prefix"] = np.asarray(h["prefix"]).astype(str)
    return h


def spec(h, prefix, det, ch, var):
    m = ((h["prefix"] == prefix) & (h["detector"] == det) & (h["channel"] == ch) & (h["variable"] == var) & (h["subchannel"] == TOTAL))
    _, o = np.unique(h["bin_index"][m], return_index=True)     # each histogram is written twice
    return h["bin_content"][m][o]


def main():
    cfg = yaml.safe_load(open(os.path.join(REPO, "fdgen/configs/legacy_1d/crosscheck_dials.yaml")))
    print(f"{'scenario':7s} {'injection':58s} {'max|H/P-1|':>10s} {'rms':>8s} {'yield H/P-1':>12s} {'max|P/CV-1|':>11s}  worst bin (det,ch,var,bin)")
    for s in cfg["scenarios"]:
        P = load(os.path.join(RUN, f"fdgen_xcP_{s['id']}_xcP_{s['id']}_PROplot_histograms.root"))
        H = load(os.path.join(RUN, f"fdgen_xcH_{s['id']}_xcH_{s['id']}_PROplot_histograms.root"))
        worst, rms, n, ypp, yhh, eff = 0.0, 0.0, 0, 0.0, 0.0, 0.0
        where = None
        for det in sorted(set(P["detector"].tolist())):
            for ch in sorted(set(P["channel"].tolist())):
                for var in sorted(set(P["variable"].tolist())):
                    p, h, cv = spec(P, "DATA", det, ch, var), spec(H, "DATA", det, ch, var), spec(P, "CV", det, ch, var)
                    if len(p) == 0 or len(p) != len(h):
                        continue
                    ok = p > 0
                    d = np.abs(h[ok] / p[ok] - 1)
                    if var == 0:                       # yields once per detector/channel (variables are different binnings of the same events)
                        ypp += p.sum(); yhh += h.sum()
                    if d.size and d.max() > worst:
                        worst, where = d.max(), (det, ch, var, int(np.argmax(np.where(ok, np.abs(h / np.where(ok, p, 1) - 1), 0))))
                    rms += (d ** 2).sum(); n += d.size
                    eff = max(eff, float(np.abs(p[ok] / np.where(cv[ok] > 0, cv[ok], np.nan) - 1).max(initial=0)))
        inj = " ".join(s["profit_inject"])
        print(f"{s['id']:7s} {inj:58s} {worst:10.2e} {np.sqrt(rms / max(n, 1)):8.1e} {yhh / ypp - 1:12.2e} {eff:11.3f}  {where}")


if __name__ == "__main__":
    main()
