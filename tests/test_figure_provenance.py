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


def make_profit_repo(tmp_path):
    """A miniature repository: ini -> XML -> work/ symlink -> data file, plus a PROfit output and its cache."""
    repo = tmp_path / "repo"
    (repo / "jobs").mkdir(parents=True)
    (repo / "xml").mkdir()
    (repo / "work").mkdir()
    (repo / "out").mkdir()
    (repo / "data").mkdir()
    (repo / "data" / "mc.root").write_bytes(b"mc v1")
    (repo / "work" / "sbnd_cvmc.root").symlink_to(repo / "data" / "mc.root")
    (repo / "xml" / "study.xml").write_text('<xml><file filename="work/sbnd_cvmc.root"/></xml>\n')
    (repo / "jobs" / "study.ini").write_text(
        "# Run:\n#   ../profit/build/bin/PROfit --config jobs/study.ini plot\nxml = xml/study.xml\ntag = mystudy\n\n[plot]\n")
    (repo / "out" / "mystudy_prop.bin").write_bytes(b"cache")
    output = repo / "out" / "mystudy_v1_PROplot_histograms.root"
    output.write_bytes(b"profit output")
    figure = repo / "fig.pdf"
    figure.write_bytes(b"figure")
    script = repo / "make.py"
    script.write_text("# producer\n")
    provenance._ini_index.cache_clear()
    sidecar = provenance.record(figure, [output], script, repo=repo)
    return repo, figure, sidecar


def test_profit_recipe_is_recorded(tmp_path):
    repo, _, sidecar = make_profit_repo(tmp_path)
    (entry,) = json.loads(sidecar.read_text())["profit"]
    assert entry["tag"] == "mystudy" and entry["xml"]["path"] == "xml/study.xml"
    assert [w["work_key"] for w in entry["work_files"]] == ["sbnd_cvmc.root"]
    assert entry["work_files"][0]["path"] == "data/mc.root"
    assert "[plot]" in entry["ini"]["content"] and entry["run_hints"][0].endswith("jobs/study.ini plot")
    assert "content_gz_b64" in entry["xml"]  # not in git here, so the XML is embedded
    assert provenance.check(sidecar, repo) == ("OK", [])


def test_changed_xml_or_work_file_makes_the_figure_stale(tmp_path):
    repo, _, sidecar = make_profit_repo(tmp_path)
    (repo / "data" / "mc.root").write_bytes(b"mc v2")
    status, reasons = provenance.check(sidecar, repo)
    assert status == "STALE" and "sbnd_cvmc.root" in " ".join(reasons)
    (repo / "data" / "mc.root").write_bytes(b"mc v1")
    (repo / "xml" / "study.xml").write_text("<xml/>\n")
    status, reasons = provenance.check(sidecar, repo)
    assert status == "STALE" and "xml/study.xml" in " ".join(reasons)


def test_restore_rewrites_the_recorded_configuration(tmp_path):
    repo, _, sidecar = make_profit_repo(tmp_path)
    original = (repo / "xml" / "study.xml").read_text()
    (repo / "xml" / "study.xml").write_text("<xml/>\n")
    written = provenance.restore(sidecar, tmp_path / "again", repo)
    assert {p.name for p in written} == {"study.ini", "study.xml"}
    assert (tmp_path / "again" / "xml" / "study.xml").read_text() == original
    steps = provenance.remake_steps(sidecar)
    assert any("PROfit --config jobs/study.ini plot" in s for s in steps) and steps[-1].startswith("python make.py")


def test_directory_inputs_are_tracked_by_count_size_and_mtime(tmp_path):
    scan = tmp_path / "scan"
    scan.mkdir()
    (scan / "a.root").write_bytes(b"a")
    figure, script = tmp_path / "fig.pdf", tmp_path / "make.py"
    figure.write_bytes(b"figure")
    script.write_text("# producer\n")
    sidecar = provenance.record(figure, [scan], script, repo=tmp_path)
    assert provenance.check(sidecar, tmp_path)[0] == "OK"
    (scan / "b.root").write_bytes(b"b")
    status, reasons = provenance.check(sidecar, tmp_path)
    assert status == "STALE" and "scan" in reasons[0]


def test_plain_and_phase2_inis_are_told_apart(tmp_path):
    inis = [tmp_path / "dp_cvrw_plot.ini", tmp_path / "dp_cvrw_phase2_plot.ini"]
    assert provenance._select_inis(inis, "cvrw") == [inis[0]]
    assert provenance._select_inis(inis, "cvrw_phase2") == [inis[1]]
    assert provenance._select_inis(inis[:1], "anything") == [inis[0]]


def test_prefix_matching_finds_the_tag_of_an_output_with_an_output_tag(tmp_path):
    repo, _, _ = make_profit_repo(tmp_path)
    provenance._ini_index.cache_clear()
    out = repo / "out" / "mystudy_cvrw_phase2_PROplot.root"
    assert [e["tag"] for e in provenance.profit_recipe(out, repo)] == ["mystudy"]
    assert provenance.profit_recipe(repo / "out" / "other_v1_PROplot.root", repo) == []


def add_manifest(repo, status="active", consumer="make.py", extra=""):
    (repo / "xml" / "MANIFEST.toml").write_text(
        f'[[family]]\npattern = "xml/study.xml"\nstatus = "{status}"\npurpose = "demo"\nconsumer = "{consumer}"\n'
        f'technote = "demo_section"\nregenerate = "python build.py"\n{extra}')
    provenance._manifest.cache_clear()


def test_manifest_family_is_recorded_and_checked(tmp_path):
    repo, figure, _ = make_profit_repo(tmp_path)
    add_manifest(repo)
    sidecar = provenance.record(figure, [repo / "out" / "mystudy_v1_PROplot_histograms.root"], repo / "make.py", repo=repo)
    (entry,) = json.loads(sidecar.read_text())["profit"]
    assert entry["xml"]["manifest"]["status"] == "active" and entry["xml"]["manifest"]["technote"] == "demo_section"
    assert provenance.check(sidecar, repo) == ("OK", [])
    assert any("python build.py" in step for step in provenance.remake_steps(sidecar))


def test_archived_family_makes_the_figure_stale_and_consumer_mismatch_is_a_note(tmp_path):
    repo, figure, _ = make_profit_repo(tmp_path)
    add_manifest(repo, consumer="other_script.py")
    sidecar = provenance.record(figure, [repo / "out" / "mystudy_v1_PROplot_histograms.root"], repo / "make.py", repo=repo)
    status, reasons = provenance.check(sidecar, repo)
    assert status == "OK" and any("manifest consumer" in r for r in reasons)
    add_manifest(repo, status="archived", extra='superseded_by = "xml/new.xml"\n')
    status, reasons = provenance.check(sidecar, repo)
    assert status == "STALE" and any("now archived" in r and "xml/new.xml" in r for r in reasons)
