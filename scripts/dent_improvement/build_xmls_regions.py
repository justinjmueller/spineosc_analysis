#!/usr/bin/env python3
"""Regional Tier A XMLs for the stage-3 contrast (docs/dent_improvement/PLAN.md).

For the muon calorimetric difference under S2, writes one Tier A study per region and sample: the
detvar CV and DENT samples, the older pair, the nominal sample and the nine other detector
variations (see ``build_xmls.TIER_A_SAMPLES``). A region is the half of the detector (sign of the
reco vertex x) times a drift-distance bin ``|x|``; it is applied as a cut on ``reco_vertex_x``, an
allowed variable, identically in the MC and data XMLs. The data XMLs are checked with the
data-access guard.

    python scripts/dent_improvement/build_xmls_regions.py
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_xmls as bx  # noqa: E402

STUDY = "muon_calo_csda"
DRIFT_EDGES = [5, 50, 100, 150, 205]
REGIONS = {}
for half, half_cut in (("xneg", "(reco_vertex_x < 0)"), ("xpos", "(reco_vertex_x >= 0)")):
    for lo, hi in zip(DRIFT_EDGES[:-1], DRIFT_EDGES[1:], strict=True):
        REGIONS[f"{half}_{lo}_{hi}"] = (
            f"{half_cut} && (abs(reco_vertex_x) >= {lo}) && (abs(reco_vertex_x) < {hi})")


def add_region(text: str, condition: str) -> str:
    """Append a region condition after the S2 cuts in every weight."""
    if text.count(bx.S2_CUTS) != text.count("additional_weight"):
        raise ValueError("every weight must contain the S2 cuts exactly once")
    return text.replace(bx.S2_CUTS, f"{bx.S2_CUTS} && {condition}")


def build(region: str, condition: str) -> None:
    """Write the data XML and every sample's Tier A XML and ini for one region."""
    source, stem, _, _ = bx.STUDIES[STUDY]
    mc, _ = bx.add_s2((source / f"{stem}.xml").read_text())
    data, _ = bx.add_s2((source / f"{stem}_data.xml").read_text())
    mc, data = add_region(mc, condition), add_region(data, condition)
    data_name = f"tierAR_{STUDY}_s2_{region}_data"
    (bx.DST / f"{data_name}.xml").write_text(data)
    for sample in bx.TIER_A_SAMPLES:
        name = f"tierAR_{STUDY}_s2_{region}_{sample}"
        bx.write(name, bx.tier_a(mc, sample), data_name, f"dimpAR_{STUDY}_s2_{region}_{sample}")


def main() -> int:
    """Build all regions and run the data-access guard on the data XMLs."""
    bx.DST.mkdir(parents=True, exist_ok=True)
    for region, condition in REGIONS.items():
        build(region, condition)
        print(f"{region}: {condition}")
    guard = bx.REPO / "scripts" / "dent_improvement" / "check_data_access.py"
    data_xmls = sorted(str(p) for p in bx.DST.glob("tierAR_*_data.xml"))
    return subprocess.run([sys.executable, str(guard), *data_xmls]).returncode


if __name__ == "__main__":
    sys.exit(main())
