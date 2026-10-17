#!/usr/bin/env python3
"""Recomputes the detector-variation ratio files of the calorimetric observables with S2.

The first DENT study built detsys_dd_<observable>.root with the S0 selection (fiducial + high-y/z
veto + flash threshold). The S2 studies (crosser veto and a muon end point at least 4 cm from the
cathode) need the ratios computed on the same population. This reuses the first study's builder,
pointing it at isolated directories so that none of that study's files are touched. The S2 files
carry the suffix ``_s2``.

    python scripts/dent_improvement/build_ratios_s2.py --stage compile|run|merge|all
"""
import argparse
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRATCH = Path("/private/tmp/claude-550/-Users-mueller-Projects-spineosc/"
               "9dc16737-7847-471d-8c39-66d492c08963/scratchpad")
OBSERVABLES = ("muon_calo_csda", "proton_calo_csda")
S2_CUTS = " && reco_cathode_crosser == 0 && std::abs(reco_leading_muon_end_x) >= 4.0"


def load_builder():
    """Import the first study's builder and redirect it to isolated S2 directories."""
    spec = importlib.util.spec_from_file_location(
        "build_dd_detsys", REPO / "scripts" / "dent_study" / "build_dd_detsys.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.SEL = module.SEL + S2_CUTS
    module.OBS = [dict(o, name=o["name"] + "_s2") for o in module.OBS if o["name"] in OBSERVABLES]
    module.DST = REPO / "detsys" / "calculator_dd_s2"
    module.OUTDIR = REPO / "detsys" / "profit_systematics_dd_s2"
    module.EXE = SCRATCH / "run_detsys_dd_s2"
    return module


def main() -> int:
    """Run the requested stages."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--stage", default="all", choices=["compile", "run", "merge", "all"])
    stage = parser.parse_args().stage
    builder = load_builder()
    if stage in ("compile", "all") and builder.compile_(None, "") != 0:
        return 1
    if stage in ("run", "all"):
        builder.run("", ("sbnd", "generic"))
    if stage in ("merge", "all"):
        builder.merge()
    return 0


if __name__ == "__main__":
    sys.exit(main())
