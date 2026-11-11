#!/usr/bin/env python3
"""Builds the per-observable detector-variation ratio files for the DENT study.

Uses an ISOLATED copy of the detsys calculator (detsys/calculator_dd/) so the tracked configuration.h is untouched. Each
observable becomes one subchannel "dd_<name>" restricted to SBND, with the SAME selection the MC branches in the study XMLs
use (>=1 proton, cathode fiducialisation, high-y/z veto, flash-PE trigger emulation) so the ratio is computed on the
population it is applied to. Outputs: detsys/profit_systematics_dd/*.root, merged to detvar/detsys_dd_<name>.root.
Run with --stage compile|run|merge|all.
"""
import argparse, json, pathlib, re, shutil, subprocess, sys

REPO = pathlib.Path(__file__).resolve().parents[2]
SRC = REPO / "detsys" / "calculator"
DST = REPO / "detsys" / "calculator_dd"
OUTDIR = REPO / "detsys" / "profit_systematics_dd"
OBS = json.loads((pathlib.Path(__file__).with_name("observables.json")).read_text())
SEL = "reco_proton_multiplicity >= 1 && reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1 && 0.647 * event_largest_flash_pe > 2000"
EXE = pathlib.Path("/private/tmp/claude-550/-Users-mueller-Projects-spineosc/32857cee-ac1d-4b8e-b98a-74f527609260/scratchpad/run_detsys_dd")

# Detector systematics the study XMLs actually use. The template's Gain entries and WireModYZ/ChargeScale (non-15) are
# commented out of every allowlist, so their ratios are not computed.
UNUSED = ('"ChargeScale", ', '"WireModYZ", ', '"GainHigh", ', '"GainLow", ')

def gen_config(names=None, workdir=None):
    workdir = workdir or DST
    workdir.mkdir(parents=True, exist_ok=True)
    for f in ("classes.h", "main.cc"):
        shutil.copy(SRC / f, workdir / f)
    obs = [o for o in OBS if names is None or o["name"] in names]
    cfg = (SRC / "configuration.h").read_text()
    cfg = "\n".join(l for l in cfg.split("\n") if not any(f"detvar_syst({u}" in l for u in UNUSED))
    start = cfg.index("inline std::vector<SubchannelConfig> build_subchannels()")
    end = cfg.index("    /**\n     * @brief Builds the detector configuration for SBND.")
    entries = []
    for o in obs:
        edges = ", ".join(f"{e:.9g}" for e in o["edges"])
        entries.append(f'            {{"dd_{o["name"]}", "{SEL}", "{o["expr_cpp"]}", {{{edges}}}, {{"sbnd_run1"}}}}')
    body = ("inline std::vector<SubchannelConfig> build_subchannels()\n    {\n        return\n        {\n" +
            ",\n".join(entries) + "\n        };\n    }\n\n")
    (workdir / "configuration.h").write_text(cfg[:start] + body + cfg[end:])
    print("wrote", workdir / "configuration.h", f"({len(obs)} subchannels)")

def compile_(names=None, tag="", ):
    wd = DST / f"b{tag}" if tag else DST
    gen_config(names, wd)
    exe = pathlib.Path(str(EXE) + (f"_b{tag}" if tag else ""))
    r = subprocess.run(["clang++", "-std=c++17", "-O2", str(wd / "main.cc"), "-o", str(exe)] + subprocess.check_output(["root-config", "--cflags", "--libs"], text=True).split() + ["-lROOTDataFrame"],
                       capture_output=True, text=True)
    print("compile rc", r.returncode); print(r.stderr[-1500:])
    return r.returncode

def run(tag="", which_list=("sbnd",)):
    OUTDIR.mkdir(exist_ok=True)
    exe = pathlib.Path(str(EXE) + (f"_b{tag}" if tag else ""))
    for which in which_list:
        r = subprocess.run([str(exe), which, str(OUTDIR)], capture_output=True, text=True)
        (OUTDIR / f"run_{which}{tag}.log").write_text(r.stdout + r.stderr)
        print(which, "rc", r.returncode, "files now", len(list(OUTDIR.glob("*.root"))))

def merge():
    out = REPO / "detvar"
    for o in OBS:
        files = sorted(OUTDIR.glob(f"*_dd_{o['name']}.root"))
        if not files:
            print("NO FILES for", o["name"]); continue
        target = out / f"detsys_dd_{o['name']}.root"
        r = subprocess.run(["hadd", "-f", str(target), *map(str, files)], capture_output=True, text=True)
        print(f"{o['name']:20s} {len(files):3d} inputs -> {target.name} rc={r.returncode}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", default="all", choices=["compile", "run", "merge", "all", "generic"])
    ap.add_argument("--names", default=None, help="comma-separated observable subset (one parallel batch)"); ap.add_argument("--tag", default="")
    a = ap.parse_args(); names = a.names.split(",") if a.names else None
    if a.stage == "generic":
        # the flat (Toroid / ProtonMiss) ratios are observable-independent in value but binned per observable; needs every subchannel
        if compile_(None, "gen") != 0: sys.exit(1)
        run("gen", ("generic",)); sys.exit(0)
    if a.stage in ("compile", "all") and compile_(names, a.tag) != 0: sys.exit(1)
    if a.stage in ("run", "all"): run(a.tag)
    if a.stage in ("merge", "all"): merge()
