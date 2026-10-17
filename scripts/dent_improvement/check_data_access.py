#!/usr/bin/env python3
"""Guard for the detector data access: fails if a data XML uses a branch outside the allowed set.

Allowed (technical note, "Detector data access"): the neutrino vertex position, the muon and
proton start and end point positions, their (E_calo - E_range)/E_range, their PID and primary
softmax scores, and the azimuthal angle. The polar angle and direction components, muon length,
and every energy, momentum or visible-energy quantity are forbidden. A handful of selection
flags used by the existing detector-access XMLs are also allowed.

This applies to the detector data access sample only. The separate coverage sample (xml/coverage/)
is not restricted this way, so do not run the guard on it.

    python scripts/dent_improvement/check_data_access.py xml/dent_improvement/*/*_data.xml
"""
import re
import sys
from pathlib import Path

PARTICLES = ("muon", "proton")
ALLOWED_EXACT = {
    "reco_vertex_x", "reco_vertex_y", "reco_vertex_z",
    "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser",
    "reco_proton_multiplicity",
    "reco_user_containment_cut_filterEE",
    "Run", "Subrun", "Evt", "event_max_neutrino_energy",
}
ALLOWED_PATTERNS = [
    *(rf"reco_leading_{p}_(start|end)_[xyz]" for p in PARTICLES),
    *(rf"reco_leading_{p}_calo_csda_diff" for p in PARTICLES),
    *(rf"reco_leading_{p}_(muon|pion|proton|mip|electron|photon|secondary|hadron)_softmax"
      for p in PARTICLES),
    *(rf"reco_leading_{p}_primary_softmax" for p in PARTICLES),
    *(rf"reco_leading_{p}_azimuthal_angle" for p in PARTICLES),
]
# Fragments that are always rejected, even if a pattern above were widened by mistake.
FORBIDDEN = ("polar_angle", "_dir_", "length", "len3d", "energy", "momentum", "_p_", "ke",
             "visible", "csda_range", "mcs")
TOKEN = re.compile(r"\b(?:reco|true)_[A-Za-z0-9_]+\b|\bevent_[A-Za-z0-9_]+\b")


def allowed(name: str) -> bool:
    """Return True if a branch name is on the allowed list and matches no forbidden fragment."""
    if name in ALLOWED_EXACT:
        return True
    if any(fragment in name for fragment in FORBIDDEN):
        return False
    return any(re.fullmatch(pattern, name) for pattern in ALLOWED_PATTERNS)


def violations(text: str) -> set[str]:
    """Branch-like tokens in an XML string that are not allowed."""
    return {t for t in TOKEN.findall(text) if not allowed(t)}


def main(argv: list[str]) -> int:
    """Check each XML given on the command line.

    Returns the number of files with violations.
    """
    bad = 0
    for path in map(Path, argv):
        found = violations(path.read_text())
        if found:
            bad += 1
            print(f"FORBIDDEN in {path}: {sorted(found)}")
    print(f"{len(argv) - bad}/{len(argv)} files clean")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main(sys.argv[1:]) else 0)
