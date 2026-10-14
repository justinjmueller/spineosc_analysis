"""Tests for the XML manifest classifier."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "xml_manifest", Path(__file__).resolve().parents[1] / "scripts" / "xml_manifest.py"
)
xml_manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(xml_manifest)

FAMILIES = [{"pattern": "xml/a/*"}, {"pattern": "xml/a/b/*"}, {"pattern": "xml/c/*"},
            {"pattern": "xml/c/d*"}, {"pattern": "xml/c/*x"}]


def test_most_specific_pattern_wins():
    family, tie = xml_manifest.match("xml/a/b/x.xml", FAMILIES)
    assert family["pattern"] == "xml/a/b/*" and not tie


def test_unmatched_path_returns_none():
    assert xml_manifest.match("xml/z/x.xml", FAMILIES) == (None, False)


def test_equal_specificity_is_reported_as_tie():
    family, tie = xml_manifest.match("xml/c/dx", FAMILIES)
    assert family is not None and tie


def test_current_manifest_is_consistent():
    result = xml_manifest.classify()
    assert xml_manifest.problems(result) == []
