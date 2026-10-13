#!/usr/bin/env python3
"""Runs the v1.1.0 data/MC coverage studies built by build_cov_v2.py (xml/coverage) with PROfit, N at a time.
Each study: PROfit --config <ini> -o v1 process  ->  runs/cov_v2/cov2_<stem>_v1_PROplot_histograms.root (area-norm, shape-only; see the ini).
Usage: python scripts/dent_study/run_cov_v2.py [-j 4] [stem ...]"""
import argparse, concurrent.futures as cf, pathlib, subprocess, sys, time
REPO = pathlib.Path(__file__).resolve().parents[2]; RUN = REPO / "runs" / "cov_v2"; PROFIT = REPO.parent / "profit-asymmetric" / "build" / "bin" / "PROfit"
XML = REPO / "xml" / "coverage"
def one(study):
    ini = XML / f"{study}_plot.ini"; tag = f"cov2_{study}"; out = RUN / f"{tag}_v1_PROplot_histograms.root"
    if out.exists() and out.stat().st_size > 1000: return study, "cached", 0
    t0 = time.time()
    with open(RUN / "logs" / f"{tag}.steps.log", "ab") as fh:
        rc = subprocess.run([str(PROFIT), "--config", str(ini), "-o", "v1", "-l", f"logs/{tag}_v1.log", "process"], cwd=RUN, stdout=fh, stderr=subprocess.STDOUT, timeout=3 * 3600).returncode
    return study, ("ok" if rc == 0 and out.exists() else f"FAILED rc={rc}"), time.time() - t0
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("-j", type=int, default=6); ap.add_argument("studies", nargs="*"); a = ap.parse_args()
    studies = a.studies or sorted(p.stem.replace("_plot", "") for p in XML.glob("*_plot.ini"))
    with cf.ThreadPoolExecutor(a.j) as ex:
        for s, st, dt in ex.map(one, studies): print(f"{time.strftime('%T')} {s}: {st} ({dt/60:.1f} min)", flush=True)
    print("ALLDONE", flush=True)
