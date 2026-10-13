#!/usr/bin/env python3
"""Moves jobs from runs/dent_study/queue_wait/ into queue/ once their requirements hold.
A waiting job's second line is `# REQUIRES dd:<observable>`: the merged per-observable ratio file detvar/detsys_dd_<obs>.root must
exist and contain every histogram the study XML references. Exits when queue_wait is empty."""
import pathlib, re, time, os
REPO = pathlib.Path(__file__).resolve().parents[2]
RUN = REPO / "runs" / "dent_study"
REQUIRED = ["DENT_sbnd_Run1", "ChargeScale15_sbnd_Run1", "WireModxThetaXW_sbnd_Run1", "MattModYZ_sbnd_Run1", "SCE0x_sbnd_Run1", "SCE2x_sbnd_Run1",
            "TriggerEmulationM1_sbnd_run1", "TriggerEmulationP1_sbnd_run1", "ToroidCalibration_generic_RunX", "ProtonMissLow_generic_RunX", "ProtonMissHigh_generic_RunX"]

def dd_ready(name):
    f = REPO / "detvar" / f"detsys_dd_{name}.root"
    if not f.exists() or time.time() - f.stat().st_mtime < 30:
        return False
    try:
        import uproot
        with uproot.open(f) as h:
            keys = {k.split(";")[0] for k in h.keys()}
        return all(r in keys for r in REQUIRED)
    except Exception:
        return False

if __name__ == "__main__":
    W, Q = RUN / "queue_wait", RUN / "queue"
    W.mkdir(exist_ok=True); Q.mkdir(exist_ok=True)
    while True:
        jobs = sorted(W.glob("*.job"))
        if not jobs: break
        for j in jobs:
            m = re.search(r"#\s*REQUIRES\s+dd:(\S+)", j.read_text())
            if m is None or dd_ready(m.group(1)):
                os.rename(j, Q / j.name)
        time.sleep(60)
