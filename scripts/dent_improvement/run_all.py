#!/usr/bin/env python3
"""Runs every generated DENT-improvement study through the first study's driver, N at a time.

Each study goes through scripts/dent_study/run_study.py (process, no-shapeonly control, DENT
variants and covariance ablations, and fake-data pseudo-experiments). Results land in
runs/dent_study/results/<tag>.json. Studies whose result file already exists are skipped.

    python scripts/dent_improvement/run_all.py [-j 4] [--tier-c 50] [study ...]

Use `caffeinate -i -s -m -t <seconds>` first on a laptop: idle sleep stalls the jobs.
"""
import argparse
import concurrent.futures as cf
import re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
XML = REPO / "xml" / "dent_improvement" / "sbnd_run1"
LOGS = REPO / "runs" / "dent_study" / "logs"
DRIVER = REPO / "scripts" / "dent_study" / "run_study.py"


def run_one(ini: Path, tier_c: int, skip_variants: bool) -> tuple[str, str, float]:
    """Run one study; return its name, a status string, and the elapsed minutes."""
    name = ini.name.removesuffix("_plot.ini")
    tag = re.search(r"^tag\s*=\s*(\S+)", ini.read_text(), re.M).group(1)
    if (REPO / "runs" / "dent_study" / "results" / f"{tag}.json").exists():
        return name, "skipped (result exists)", 0.0
    start = time.time()
    with open(LOGS / f"{name}.driver.log", "ab") as log:
        cmd = [sys.executable, str(DRIVER), "--det", "sbnd_run1", "--study", name,
               "--ini", str(ini), "--tierC", str(tier_c)]
        if skip_variants:
            cmd.append("--skip-dent-variants")
        code = subprocess.run(cmd, cwd=REPO, stdout=log, stderr=subprocess.STDOUT).returncode
    return name, "ok" if code == 0 else f"FAILED rc={code}", (time.time() - start) / 60


def main() -> None:
    """Parse arguments and run the studies."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-j", type=int, default=4, help="studies run in parallel")
    parser.add_argument("--tier-c", type=int, default=50, help="pseudo-experiment seeds per study")
    parser.add_argument("--glob", default="*_plot.ini", help="which generated inis to consider")
    parser.add_argument("--skip-variants", action="store_true",
                        help="process and the no-shapeonly control only (Tier A)")
    parser.add_argument("studies", nargs="*", help="study names, e.g. muon_score_inclusive_s2")
    args = parser.parse_args()
    LOGS.mkdir(parents=True, exist_ok=True)
    inis = sorted(XML.glob(args.glob))
    if args.studies:
        inis = [p for p in inis if p.name.removesuffix("_plot.ini") in args.studies]
    with cf.ThreadPoolExecutor(args.j) as pool:
        futures = [pool.submit(run_one, ini, args.tier_c, args.skip_variants) for ini in inis]
        for future in cf.as_completed(futures):
            name, status, minutes = future.result()
            print(f"{time.strftime('%T')} {name}: {status} ({minutes:.1f} min)", flush=True)
    print("ALLDONE", flush=True)


if __name__ == "__main__":
    main()
