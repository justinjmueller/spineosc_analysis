#!/usr/bin/env python3
"""File-based job queue worker. A job is runs/dent_study/queue/<NNN>_<name>.job (one shell command line).
Claim = atomic rename to .running.<pid>; finish = rename to .done / .failed. Run several workers for parallelism.
Stops when queue/STOP exists, or after --idle-hours with nothing to do. Pauses (does not fail) if free disk < --min-free-gb."""
import argparse, os, pathlib, shutil, subprocess, sys, time

RUN = pathlib.Path(__file__).resolve().parents[2] / "runs" / "dent_study"
Q = RUN / "queue"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--idle-hours", type=float, default=10.0); ap.add_argument("--min-free-gb", type=float, default=25.0)
    a = ap.parse_args(); Q.mkdir(parents=True, exist_ok=True)
    pid = os.getpid(); idle_since = time.time(); (RUN / "logs").mkdir(exist_ok=True)
    wlog = open(RUN / "logs" / f"worker_{pid}.log", "a", buffering=1)
    print(f"[{time.strftime('%F %T')}] worker {pid} up", file=wlog)
    while not (Q / "STOP").exists():
        free = shutil.disk_usage(RUN).free / 1e9
        jobs = sorted(Q.glob("*.job")) if free >= a.min_free_gb else []
        if free < a.min_free_gb: print(f"[{time.strftime('%T')}] paused: only {free:.0f} GB free", file=wlog)
        claimed = None
        for j in jobs:
            tgt = j.with_suffix(f".running.{pid}")
            try: os.rename(j, tgt); claimed = tgt; break
            except FileNotFoundError: continue
        if claimed is None:
            if time.time() - idle_since > a.idle_hours * 3600: break
            time.sleep(20); continue
        idle_since = time.time(); cmd = claimed.read_text().strip()
        print(f"[{time.strftime('%T')}] START {claimed.name.split('.')[0]}: {cmd}", file=wlog)
        t0 = time.time(); rc = subprocess.run(cmd, shell=True, cwd=RUN).returncode
        final = claimed.with_name(claimed.name.split(".running")[0] + (".done" if rc == 0 else ".failed"))
        os.rename(claimed, final)
        print(f"[{time.strftime('%T')}] END   {final.name} rc={rc} {time.time()-t0:.0f}s", file=wlog)
    print(f"[{time.strftime('%F %T')}] worker {pid} exiting", file=wlog)

if __name__ == "__main__":
    main()
