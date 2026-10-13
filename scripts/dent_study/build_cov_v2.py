#!/usr/bin/env python3
"""Builds the data/MC coverage comparisons for the technote into xml/coverage/ from the current xml/coverage/ files (v1.1.0 selection).
Changes (and nothing else): MC blocks -> sbnd_large trees, v1.1.0 exposures, main beam MC from work/*_cvrw.root with `* cvwgt`, SBND veto + crosser + muon end-point
(>= 4 cm) cuts on every SBND branch (MC and data; ICARUS gets the veto and the Run 2 off-beam EE cut as in selected/*.xml); SBND WireMod YZ replaced by the MattMod YZ
non-uniformity variation; detector ratios from work/detsys_v2_<observable>.root for observables that include SBND; detector exposure = the v1.1.0 development-sample POT of dp_data.xml;
ini: area-norm and shapeonly both on. Run from the repository root:  python scripts/dent_study/build_cov_v2.py"""

# One-time migration (applied 2026-10-08): the v1.1.0 versions now ARE xml/coverage; the pre-migration originals were archived/removed.
# Re-running would transform already-migrated files, so this refuses unless --force is given (then restore the original tree from git history first).
import sys as _sys
if "--force" not in _sys.argv: raise SystemExit("build is a one-time migration that has already been applied; the XMLs in xml/coverage are the source of truth (see xml/MANIFEST.toml). Use --force only after restoring the originals.")
import pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import update_selected_signal_xmls as us
R = us.ROOT; SRC = R / "xml" / "coverage"; DST = R / "xml" / "coverage"
STEMS = ["muon_score", "proton_score", "muon_primary_score", "proton_primary_score", "proton_ke", "dp", "dpt", "dalphat", "dphit", "muon_azimuth",
         "vertexx_sbnd", "vertexy_sbnd", "vertexz_sbnd", "muon_endx_sbnd", "muon_endy_sbnd", "muon_endz_sbnd",
         "vertexx_icarus", "vertexy_icarus", "vertexz_icarus", "muon_endx_icarus", "muon_endy_icarus", "muon_endz_icarus", "zoom_icarus"]
DEV_POT = {"SBND": "4.514958e+18", "ICARUSRun2": "1.965837e+19", "ICARUSRun4": "1.887259e+19"}
CUTS = "(reco_veto_sbnd_highy_highz == 1) && (reco_cathode_crosser == 0) && (abs(reco_leading_muon_end_x) >= 4)"

def set_pots(t):
    for det, pot in DEV_POT.items(): t = re.sub(rf'(<detector name="{det}"\s+pot=")[^"]+(")', rf'\g<1>{pot}\2', t)
    return t

def data_cuts(t):
    n = [0]
    def blk(m):
        b = m.group(0); hdr = re.match(r"<MCFile[^>]*>", b).group(0)
        if "events/sbnd_bnblight" in hdr:
            b = re.sub(r'pot="[^"]+"', f'pot="{DEV_POT["SBND"]}"', b, count=1)
            def w(mm):
                x = mm.group(2).replace("(reco_veto_sbnd_highy_highz == 1) && ", "")
                assert "(reco_fiducialize_cathode == 1)" in x, x
                n[0] += 1; return f'{mm.group(1)}"' + x.replace("(reco_fiducialize_cathode == 1)", "(reco_fiducialize_cathode == 1) && " + CUTS, 1) + '"'
            b = re.sub(r'(additional_weight\s*=\s*)"([^"]*)"', w, b)
        elif "events/icarus_bnblight_run2" in hdr: b = re.sub(r'pot="[^"]+"', f'pot="{DEV_POT["ICARUSRun2"]}"', b, count=1)
        elif "events/icarus_bnblight_run4" in hdr: b = re.sub(r'pot="[^"]+"', f'pot="{DEV_POT["ICARUSRun4"]}"', b, count=1)
        return b
    return re.sub(r"<MCFile .*?</MCFile>", blk, t, flags=re.S), n[0]

def mattmod(t):
    def f(m):
        line = m.group(0)
        if 'name="WireModYZ_sbnd_Run1"' not in line: return line
        new = line.replace("WireModYZ_sbnd_Run1", "MattModYZ_sbnd_Run1").replace("WireMod yz SBND", "Non-uniformity SBND")
        if 'mirror=' not in new: new = new.replace(' plotname=', ' mirror="false" plotname=', 1)
        return "    <!-- SBND uses the MattMod non-uniformity variation in place of WireMod YZ (as in the analysis) -->\n" + new
    return re.sub(r"^.*<allowlist[^\n]*\n", f, t, flags=re.M)

def build(stem):
    ini_src = SRC / f"{stem}_plot.ini"
    mc = (SRC / f"{stem}.xml").read_text(); dt = (SRC / f"{stem}_data.xml").read_text()
    mc = us.fix_selected(mc).replace("events/sbnd/", "events/sbnd_large/")
    has_sbnd = "sbnd" in stem or stem in ("muon_score", "proton_score", "muon_primary_score", "proton_primary_score", "proton_ke", "dp", "dpt", "dalphat", "dphit", "muon_azimuth")
    if has_sbnd:
        mc = mattmod(mc); mc = re.sub(r"work/detsys_([A-Za-z0-9_]+)\.root", r"work/detsys_v2_\1.root", mc)
    mc = set_pots(mc); dt, nd = data_cuts(dt); dt = set_pots(dt)
    DST.mkdir(parents=True, exist_ok=True); (DST / f"{stem}.xml").write_text(mc); (DST / f"{stem}_data.xml").write_text(dt)
    ini = ini_src.read_text()
    ini = re.sub(r"^xml\s*=.*$", f"xml = {DST}/{stem}.xml", ini, flags=re.M); ini = re.sub(r"^data\s*=.*$", f"data = {DST}/{stem}_data.xml", ini, flags=re.M)
    ini = re.sub(r"^tag\s*=.*$", f"tag = cov2_{stem}", ini, flags=re.M); ini = re.sub(r"^#?\s*shapeonly\s*=.*$", "shapeonly = true", ini, flags=re.M)
    if "shapeonly" not in ini: ini = ini.replace("area-norm = true", "area-norm = true\nshapeonly = true")
    (DST / f"{stem}_plot.ini").write_text(ini)
    return nd

if __name__ == "__main__":
    for s in STEMS: print(s, "data branches with cuts:", build(s))
