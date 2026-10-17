import numpy as np

from fdgen import export


def test_reco_cut_drops_truth_multiplicity_and_cvwgt():
    w = ("(true_cc && (true_pdg == 14) && true_neutrino_id >= 0 && (reco_proton_multiplicity == 1) && (reco_fiducialize_cathode == 1) "
         "&& ((true_detector == 1 && (0.647 * event_largest_flash_pe > 2000)) || (true_detector == 2 && Run <= 10097))) * cvwgt")
    assert export.reco_cut(w) == ("(reco_fiducialize_cathode == 1) && ((true_detector == 1 && (0.647 * event_largest_flash_pe > 2000)) "
                                  "|| (true_detector == 2 && Run <= 10097))")
    # an XML weight without outer parentheses or cvwgt (dirt) and an off-beam branch whose only truth term is dropped
    assert export.reco_cut("true_cc && (reco_proton_multiplicity > 1) && (reco_veto == 1)") == "(reco_veto == 1)"
    assert export.reco_cut("(true_neutrino_id != true_neutrino_id) && (reco_proton_multiplicity == 1) && (reco_fiducialize_cathode == 1)") == "(reco_fiducialize_cathode == 1)"


def test_blind_variable_replaces_truth_with_zero():
    assert export.blind_variable("true_baseline / (1000 * true_neutrino_energy)") == "0"
    assert export.blind_variable("reco_visible_energy;reco_pn/1000") == "reco_visible_energy;reco_pn/1000"


def test_jitter_touches_only_repeated_continuous_rows():
    rng = np.random.default_rng(0)
    cols = {"e": np.array([0.5, 0.5, 0.7, 0.5]), "pn": np.array([100.3, 100.3, 50.1, 100.3]), "m": np.array([1., 1., 2., 1.])}
    out = export.jitter_duplicates(cols, 1e-6, rng)
    assert out["e"][0] == 0.5 and out["pn"][0] == 100.3 and out["e"][2] == 0.7        # first copies and unique rows are exact
    assert len(set(zip(out["e"], out["pn"]))) == 4                                    # no identical rows remain
    assert (out["e"][[1, 3]] != 0.5).all()                                           # repeats are jittered
    assert np.allclose(out["e"], cols["e"], rtol=2e-6)                                # and only by ~1e-6
    assert (out["m"] == cols["m"]).all()                                              # integer-valued column untouched


def test_write_tree_is_readable_by_root_in_nested_directories(tmp_path):
    import pytest
    ROOT = pytest.importorskip("ROOT")
    path = str(tmp_path / "x.root")
    export.write_tree(path, "events/sbnd_bnblight/selected", {"reco_visible_energy": np.linspace(0.3, 2.0, 50), "reco_proton_multiplicity": np.ones(50)})
    f = ROOT.TFile.Open(path)
    t = f.Get("events/sbnd_bnblight/selected")          # uproot's nested directories segfault here
    assert t.GetEntries() == 50 and {b.GetName() for b in t.GetListOfBranches()} == {"reco_visible_energy", "reco_proton_multiplicity"}


def test_normalize_exports_merges_defaults_and_formats_file_names():
    import pytest
    from fdgen.generate import normalize_exports
    cfg = {"export_defaults": {"keep": ["a"], "detectors": {"SBND": {"file": "{name}_s.root", "tree": "t", "cuts": {"x": "1"}}}},
           "exports": [{"name": "d1", "scenario": "O1", "seed": 1}, {"name": "d2", "scenario": "null", "seed": 2, "keep": ["b"]}]}
    out = normalize_exports(cfg)
    assert [e["detectors"]["SBND"]["file"] for e in out] == ["d1_s.root", "d2_s.root"]
    assert out[0]["keep"] == ["a"] and out[1]["keep"] == ["b"]
    with pytest.raises(SystemExit):          # an unquoted YAML `null` scenario id is None, not the string 'null'
        normalize_exports({**cfg, "exports": [{"name": "d3", "scenario": None, "seed": 3}]})
