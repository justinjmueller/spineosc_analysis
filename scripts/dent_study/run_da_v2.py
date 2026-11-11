#!/usr/bin/env python3
"""Runs the v1.1.0 detector-data-access studies built by build_da_v2.py (xml/detector_data_access/sbnd_run1) with PROfit, N at a time.
Each study: PROfit --config <ini> -o v1 process  ->  runs/da_v2/da2_sbnd_run1_<study>_v1_PROplot_histograms.root (area-norm, shape-only; see the ini).
Usage: python scripts/dent_study/run_da_v2.py [-j 6] [study ...]"""
import argparse, concurrent.futures as cf, pathlib, subprocess, sys, time
REPO = pathlib.Path(__file__).resolve().parents[2]; RUN = REPO / "runs" / "da_v2"; PROFIT = REPO.parent / "profit-asymmetric" / "build" / "bin" / "PROfit"
XML = REPO / "xml" / "detector_data_access" / "sbnd_run1"
def one(study):
    ini = XML / f"{study}_plot.ini"; tag = f"da2_sbnd_run1_{study}"; out = RUN / f"{tag}_v1_PROplot_histograms.root"
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
