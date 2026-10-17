#!/usr/bin/env python3
"""Brings xml/selected/{inclusive,exclusive}.xml and xml/signal/{inclusive,exclusive}.xml to v1.1.0.
Starts from the committed (HEAD) versions so that the systematics added to selected/inclusive.xml are dropped.

selected: v1.1.0 tree names and exposures, main beam MC read from work/*_cvrw.root with `* cvwgt`, reco_veto_sbnd_highy_highz in every
          branch (as in contours_base.xml), reco_cathode_crosser == 0 on SBND branches only.
signal:   v1.1.0 tree names and exposures, true_fiducialize_cathode, true_veto_sbnd_highy_highz in every branch and
          true_cathode_crosser == 0 on SBND branches only. No cvwgt (the signal tree is not reweighted).
Run from the repository root."""
import re, subprocess, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
POT = {"sbnd_large": "5.593437e+20", "icarus_run2": "1.040833e+21", "icarus_run4": "1.772140e+21", "sbnd_dirt": "3.021459e+20",
       "sbnd_offbeam": "3.983079e+20", "icarus_dirt_run2": "1.276535e+21", "icarus_dirt_run4": "1.774677e+21",
       "icarus_offbeam_run2": "1.979496e+20", "icarus_offbeam_run4": "3.422616e+19"}
MAIN = {"sbnd_large": "sbnd_run1", "icarus_run2": "icarus_run2", "icarus_run4": "icarus_run4"}

def head(path): return subprocess.check_output(["git", "show", f"HEAD:{path}"], cwd=ROOT, text=True)

def blocks(text):
    """split into (preamble, [(open_tag_match, body, close)...]) preserving everything"""
    return re.split(r"(?=<MCFile )", text)

def fix_selected(text):
    out = []
    for part in blocks(text):
        m = re.match(r'<MCFile treename="events/([^/]+)/(selected(?:_nonmatched)?)" filename="([^"]+)" scale = "1.0" pot="([^"]+)">', part)
        if not m: out.append(part); continue
        key, tree, fname, pot = m.groups(); key = "sbnd_large" if key == "sbnd" else key
        is_main = tree == "selected" and key in MAIN
        new_open = f'<MCFile treename="events/{key}/{tree}" filename="{"work/" + MAIN[key] + "_cvrw.root" if is_main else fname}" scale = "1.0" pot="{POT[key]}">'
        part = part.replace(m.group(0), new_open, 1)
        sbnd = key.startswith("sbnd")
        def fix(w):
            w0 = w.group(1)
            if "reco_veto_sbnd_highy_highz" not in w0:
                add = " && (reco_veto_sbnd_highy_highz == 1)" + (" && (reco_cathode_crosser == 0) && (abs(reco_leading_muon_end_x) >= 4)" if sbnd else "")
                w0 = w0.replace("(reco_fiducialize_cathode == 1)", "(reco_fiducialize_cathode == 1)" + add, 1)
            if key == "icarus_offbeam_run2" and "filterEE" not in w0:
                w0 += " && (reco_user_containment_cut_filterEE == 1)"
            return f'additional_weight     = "{("(" + w0 + ") * cvwgt") if is_main else w0}"'
        part = re.sub(r'additional_weight\s*=\s*"([^"]*)"', fix, part)
        out.append(part)
    return "".join(out)

def fix_signal(text):
    out = []
    for part in blocks(text):
        m = re.match(r'<MCFile treename="events/([^/]+)/signal" filename="([^"]+)" scale = "1.0" pot="([^"]+)">', part)
        if not m: out.append(part); continue
        key = m.group(1); key = "sbnd_large" if key == "sbnd" else key
        part = part.replace(m.group(0), f'<MCFile treename="events/{key}/signal" filename="{m.group(2)}" scale = "1.0" pot="{POT[key]}">', 1)
        extra = " && (true_fiducialize_cathode == 1) && (true_veto_sbnd_highy_highz == 1)" + (" && (true_cathode_crosser == 0) && (abs(true_leading_muon_end_x) >= 4)" if key.startswith("sbnd") else "")
        def fix(w):
            w0 = w.group(1)
            if "true_fiducialize_cathode" not in w0: w0 += extra
            return f'additional_weight     = "{w0}"'
        part = re.sub(r'additional_weight\s*=\s*"([^"]*)"', fix, part)
        out.append(part)
    return "".join(out)

NOTE_SEL = """<!--
v1.1.0 update. The beam MC is read from the CV-reweighted files (work/*_cvrw.root, written by reweight/) and its branches
multiply cvwgt into their weight; the reweighted files hold only the beam `selected` trees, so the dirt, off-beam and nonmatched blocks keep
their original files. Every branch applies the cathode fiducialization and the SBND high-y/high-z veto (reco_veto_sbnd_highy_highz);
SBND branches also remove cathode crossers (reco_cathode_crosser == 0) and muons whose end point is within 4 cm of the cathode
(abs(reco_leading_muon_end_x) >= 4). Neither is applied to ICARUS (two cryostats; the crosser flag is set for ~39% of its events). No systematics beyond mcstat are configured. Exposures are the v1.1.0 values from
plotting/studies/exposure.py.
-->
"""
NOTE_SIG = """<!--
v1.1.0 update. The signal definition uses truth quantities only: true_fiducialize_cathode and the SBND true_veto_sbnd_highy_highz
in every branch, true_cathode_crosser == 0 and abs(true_leading_muon_end_x) >= 4 on SBND branches only (not applied to ICARUS: two cryostats), and
true_user_containment_cut_filterEE on ICARUS Run 2. The signal tree is not CV-reweighted, so there is no cvwgt here.
-->
"""

def main():
    for name in ("inclusive", "exclusive"):
        t = fix_selected(head(f"xml/selected/{name}.xml"))
        i = t.index("<mode "); t = t[:i] + NOTE_SEL + "\n" + t[i:]
        (ROOT / f"xml/selected/{name}.xml").write_text(t)
        t = fix_signal(head(f"xml/signal/{name}.xml"))
        i = t.index("<mode "); t = t[:i] + NOTE_SIG + "\n" + t[i:]
        (ROOT / f"xml/signal/{name}.xml").write_text(t)
        print("wrote", name)

if __name__ == "__main__": main()
