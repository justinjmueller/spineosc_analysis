#!/usr/bin/env python3
"""Compare by-hand oscillation datasets with PROfit's own `-i` injection (stat-only `plot` histograms), bin by bin.

  .venv/bin/python -m fdgen.osccheck runs/fdgen_2d batch1_2d fdgen/configs/batch1.yaml   # prefix = the data-XML prefix; config = where the scenarios are
reads <rundir>/fdgen_osc{P,H}_<prefix>_<id>_fdgen_osc{P,H}_<prefix>_<id>_PROplot_histograms.root (made by fdgenosc.sh), for the scenarios of
fdgen/configs/<prefix>.yaml that have an oscillation block and no dial effects. Prints, for the analysis bins (variable 0; in the 2D analysis the plot is the reco-E projection) and the true-E
bins (variable 2; not plotted by contours_base.xml), the rms and maximum of |H/P - 1| over detectors, channels and bins, the total-yield difference and the largest
fractional deficit relative to the CV (the size of the injected effect).
"""
import os
import sys

import numpy as np
import yaml

from .crosscheck import load, spec

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def main():
    rundir, prefix = sys.argv[1], sys.argv[2]
    cfg_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(REPO, "fdgen/configs", f"{prefix}.yaml")
    cfg = yaml.safe_load(open(cfg_path))
    rows = [s for s in cfg["scenarios"] if s.get("oscillation") and not s.get("effects")]
    print(f"{'id':4s} {'dm2':>6s} {'s22':>6s} | {'analysis bins rms / max':>24s} | {'true-E bins rms / max':>22s} | {'yield H/P-1':>11s} | max|P/CV-1|")
    for s in rows:
        p = lambda r: os.path.join(REPO, rundir, f"fdgen_osc{r}_{prefix}_{s['id']}_fdgen_osc{r}_{prefix}_{s['id']}_PROplot_histograms.root")
        if not (os.path.exists(p("P")) and os.path.exists(p("H"))):
            print(f"{s['id']:4s} (missing)")
            continue
        P, H = load(p("P")), load(p("H"))
        out, yp, yh, eff = [], 0.0, 0.0, 0.0
        for var in (0, 2):
            d = []
            for det in (0, 1, 2):
                for ch in (0, 1):
                    a, b, cv = spec(P, "DATA", det, ch, var), spec(H, "DATA", det, ch, var), spec(P, "CV", det, ch, var)
                    ok = a > 0
                    d += list(np.abs(b[ok] / a[ok] - 1))
                    if var == 0:
                        yp += a.sum(); yh += b.sum()
                        eff = max(eff, float(np.abs(a[ok] / cv[ok] - 1).max()))
            d = np.array(d)
            out.append(f"{np.sqrt((d ** 2).mean()):.1e} / {d.max():.1e}" if d.size else "n/a (not plotted)")
        o = s["oscillation"]
        print(f"{s['id']:4s} {o['dmsq']:6.2f} {o['sinsq2thmm']:6.3f} | {out[0]:>24s} | {out[1]:>22s} | {yh / yp - 1:11.1e} | {eff:.3f}")


if __name__ == "__main__":
    main()
