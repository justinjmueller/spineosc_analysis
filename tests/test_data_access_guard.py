"""Tests for the detector-data-access branch whitelist."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "check_data_access",
    Path(__file__).resolve().parents[1] / "scripts" / "dent_improvement" / "check_data_access.py",
)
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)


def test_allowed_variables_pass():
    for name in ("reco_vertex_x", "reco_leading_muon_end_x", "reco_leading_proton_start_z",
                 "reco_leading_muon_muon_softmax", "reco_leading_proton_primary_softmax",
                 "reco_leading_muon_calo_csda_diff", "reco_leading_muon_azimuthal_angle",
                 "reco_cathode_crosser"):
        assert guard.allowed(name), name


def test_forbidden_variables_fail():
    for name in ("reco_visible_energy", "reco_leading_muon_length", "reco_leading_muon_polar_angle",
                 "reco_leading_muon_start_dir_x", "reco_leading_muon_azimuthal_angle_dent",
                 "reco_leading_muon_p", "reco_leading_muon_ke", "reco_leading_muon_mcs_ke",
                 "event_largest_flash_pe"):
        assert not guard.allowed(name), name


def test_expression_scan_finds_only_the_forbidden_names():
    text = ("(reco_veto_sbnd_highy_highz == 1) && abs(reco_leading_muon_end_x) >= 4 "
            "&& reco_visible_energy > 0.3")
    assert guard.violations(text) == {"reco_visible_energy"}
