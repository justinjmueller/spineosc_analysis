#!/usr/bin/env python3
"""Generates the wide-net DENT-study XMLs for the observables in observables.json.

Tier B  : production CV (cvrw + cvwgt) with the full systematic suite and the per-observable detector ratios
          (work/detsys_dd_<name>.root).                      -> xml/dent_study/sbnd_run1/dd/<name>{,_data}.xml
Tier A  : the detvar CV and DENT samples used DIRECTLY as the MC, MC-stat only (the DENT sample has no interaction-model
          weights, so no interaction-model systematics are possible there).
                                                              -> xml/dent_study/sbnd_run1/dd/tierA_<name>_{cv,dent}.xml
Both share the data XML, which carries the same selection as the MC (high-y/z veto included; no flash cut on data).
Template: the existing inclusive single-channel study xml/detector_data_access/sbnd_run1/muon_score_inclusive.xml.
"""
import json, pathlib, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import build_xmls as bx

REPO = bx.REPO
TEMPLATE = "muon_score_inclusive"
OBS = json.loads(pathlib.Path(__file__).with_name("observables.json").read_text())
OUT = REPO / "xml" / "dent_study" / "sbnd_run1" / "dd"

def substitute_observable(text, o):
    edges = " ".join(f"{e:g}" for e in o["edges"])
    text = re.sub(r'(<channel name="numu" unit=")[^"]*(")', lambda m: f'{m.group(1)}{o["label"].replace("$","")}{m.group(2)}', text, count=1)
    text = re.sub(r'<bins unit="Reco Muon Score"[^>]*/>', f'<bins unit="{o["name"]}" edges="{edges}"/>', text, count=1)
    def first_variable(m):
        br = m.group(0)
        return re.sub(r"<variable>.*?</variable>", f"<variable>{o['expr_xml']}</variable>", br, count=1, flags=re.S)
    text = re.sub(r"<branch\b.*?</branch>", first_variable, text, flags=re.S)
    text = text.replace(f"work/detsys_{TEMPLATE.replace('_inclusive','')}.root", f"work/detsys_dd_{o['name']}.root")
    # GainHigh/GainLow ratios are not computed (the Gain allowlist entry is commented out in the template), so drop the section
    text = re.sub(r"\s*<!-- Gain Variation -->\s*<HistVarSection name=\"Gain_sbnd_Run1\">.*?</HistVarSection>", "", text, flags=re.S)
    assert "Gain_sbnd_Run1\"> " not in text and 'HistVarSection name="Gain_sbnd_Run1"' not in text
    assert "Reco Muon Score" not in text and "detsys_muon_score" not in text
    return text

def tier_a(text, o, which):
    """detvar sample as the only MC; mcstat only."""
    tree = {"cv": "events/sbnd_detvar_newcv/selected", "dent": "events/sbnd_detvar_newdent/selected"}[which]
    blocks = re.findall(r"<MCFile .*?</MCFile>", text, flags=re.S)
    main = [b for b in blocks if 'treename="events/sbnd_large/selected"' in re.match(r"<MCFile[^>]*>", b).group(0)][0]
    main = re.sub(r"<MCFile[^>]*>", f'<MCFile treename="{tree}" filename="work/sbnd_run1_detvar.root" scale = "1.0" pot="1.0e+19">', main, count=1)
    main = re.sub(r"\s*<friend[^>]*/>", "", main)
    main = main.replace(' * cvwgt"', '"').replace('"(', '"(', 1)
    main = re.sub(r'incl_systematics\s*=\s*"true"', 'incl_systematics      = "false"', main)
    # drop all MC blocks, then append the single detvar block
    text = re.sub(r"<MCFile .*?</MCFile>\s*", "", text, flags=re.S)
    text = re.sub(r"<HistVarFiles>.*?</HistVarFiles>\s*", "", text, flags=re.S)
    text = re.sub(r"<variation_list>.*?</variation_list>",
                  '<variation_list>\n    <allowlist type="mcstat" plotname="MC Stats" tag="stat,MC Stats">mcstat</allowlist>\n</variation_list>', text, flags=re.S)
    return text.rstrip() + "\n\n" + main + "\n"

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = bx.SRC / "sbnd_run1"
    mc_t, _ = bx.convert_mc((src / f"{TEMPLATE}.xml").read_text(), "sbnd_run1")
    da_t, _ = bx.convert_data((src / f"{TEMPLATE}_data.xml").read_text(), "sbnd_run1")
    raw_mc = (src / f"{TEMPLATE}.xml").read_text()
    for o in OBS:
        n = o["name"]
        (OUT / f"{n}.xml").write_text(substitute_observable(mc_t, o))
        data = substitute_observable(da_t, o); (OUT / f"{n}_data.xml").write_text(data)
        (OUT / f"{n}_plot.ini").write_text(bx.INI.format(xml=OUT / f"{n}.xml", data=OUT / f"{n}_data.xml", tag=f"dent_sbnd_run1_dd_{n}"))
        for which in ("cv", "dent"):
            # Tier A starts from the ORIGINAL template (before cvwgt/cvrw conversion) so nothing from the production file leaks in
            txt = tier_a(substitute_observable(raw_mc, o), o, which)
            (OUT / f"tierA_{n}_{which}.xml").write_text(txt)
            (OUT / f"tierA_{n}_{which}_plot.ini").write_text(bx.INI.format(xml=OUT / f"tierA_{n}_{which}.xml", data=OUT / f"{n}_data.xml", tag=f"dentA_sbnd_run1_dd_{n}_{which}"))
    print(f"wrote {len(OBS)} observables x (B + A_cv + A_dent) to {OUT}")

if __name__ == "__main__":
    main()
