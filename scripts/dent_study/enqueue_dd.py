#!/usr/bin/env python3
"""Reprioritises the queue and adds the wide-net observables.  Priority classes (filename prefix, lexicographic):
   10 Tier A (detvar CV / DENT as MC) for the PRIMARY observables        - cheap, no ratio files needed
   20 Tier B + Tier C for the PRIMARY observables                         - wait for their detsys ratio file
   30 regenerated SBND data-access studies (existing XMLs)
   40 Tier A for the remaining observables      41 Tier B for the remaining observables (wait for ratio file)
   50 regenerated ICARUS studies"""
import json, os, pathlib, re
REPO = pathlib.Path(__file__).resolve().parents[2]; RUN = REPO / "runs" / "dent_study"
Q, W = RUN / "queue", RUN / "queue_wait"; Q.mkdir(exist_ok=True); W.mkdir(exist_ok=True)
PY = f"python3 {REPO}/scripts/dent_study/run_study.py"
DD = REPO / "xml/dent_study/sbnd_run1/dd"
obs = [o["name"] for o in json.loads((pathlib.Path(__file__).with_name("observables.json")).read_text())]
primary = json.loads((pathlib.Path(__file__).with_name("primary_selection.json")).read_text())["primary"]
TIERC_N = 50

# 1. re-prefix the pending baseline jobs
for j in sorted(Q.glob("*.job")):
    if re.match(r"^(10|20|30|40|41|50)_", j.name): continue
    pref = "50_" if "icarus_run2" in j.name else "30_"
    try: os.rename(j, Q / (pref + j.name))
    except FileNotFoundError: pass       # already claimed by a worker

def add(dirpath, name, cmd, requires=None):
    (dirpath / f"{name}.job").write_text(cmd + "\n" + (f"# REQUIRES dd:{requires}\n" if requires else ""))

for n in obs:
    isp = n in primary
    for which in ("cv", "dent"):
        ini = DD / f"tierA_{n}_{which}_plot.ini"
        add(Q, f"{'10' if isp else '40'}_A_{n}_{which}", f"{PY} --det sbnd_run1 --study tierA_{n}_{which} --ini {ini} --skip-dent-variants")
    add(W, f"{'20' if isp else '41'}_B_{n}", f"{PY} --det sbnd_run1 --sub dd --study {n}" + (f" --tierC {TIERC_N}" if isp else ""), requires=n)
print("queue:", len(list(Q.glob('*.job'))), "ready;", len(list(W.glob('*.job'))), "waiting on ratio files")
