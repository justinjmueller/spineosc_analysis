"""Tests for the figure-provenance sidecars (plotting/studies/provenance.py)."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plotting" / "studies"))
import provenance  # noqa: E402


def make(tmp_path):
    figure, data, script = tmp_path / "fig.pdf", tmp_path / "in.root", tmp_path / "make.py"
    figure.write_bytes(b"figure v1")
    data.write_bytes(b"input v1")
    script.write_text("# producer v1\n")
    sidecar = provenance.record(figure, [data], script, arguments={"detector": "sbnd"})
    return figure, data, script, sidecar


def test_fresh_figure_is_ok(tmp_path):
    _, _, _, sidecar = make(tmp_path)
    assert provenance.check(sidecar) == ("OK", [])
    rec = json.loads(sidecar.read_text())
    assert rec["arguments"] == {"detector": "sbnd"} and rec["inputs"][0]["sha256"]


def test_changed_input_script_or_figure_are_flagged(tmp_path):
    figure, data, script, sidecar = make(tmp_path)
    data.write_bytes(b"input v2")
    status, reasons = provenance.check(sidecar)
    assert status == "STALE" and "in.root" in reasons[0]
    data.write_bytes(b"input v1")
    script.write_text("# producer v2\n")
    assert provenance.check(sidecar)[0] == "STALE"
    script.write_text("# producer v1\n")
    figure.write_bytes(b"edited by hand")
    assert provenance.check(sidecar)[0] == "MODIFIED"


def test_touching_a_small_input_is_not_a_change(tmp_path):
    _, data, _, sidecar = make(tmp_path)
    data.write_bytes(b"input v1")  # new mtime, same content
    assert provenance.check(sidecar)[0] == "OK"


def test_missing_input_is_reported(tmp_path):
    _, data, _, sidecar = make(tmp_path)
    data.unlink()
    assert provenance.check(sidecar)[0] == "MISSING"
