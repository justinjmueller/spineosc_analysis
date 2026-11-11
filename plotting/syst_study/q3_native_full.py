#!/usr/bin/env python3
"""Q3: unmatched mode on the WHOLE production selected trees (no intermediate tree, no pool restriction). Cut applied as an additional_weight formula.
usage: q3_native_full.py <det__variation> <potmode: equal|real>"""
import sys, json, subprocess, numpy as np, uproot, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from q3_native import *
def cutformula(det):
    f = "reco_fiducialize_cathode == 1 && reco_veto_sbnd_highy_highz == 1"
    return f + (" && reco_cathode_crosser == 0 && abs(reco_leading_muon_end_x) >= 4" if det == "sbnd_run1" else "")
def run_full(name, potmode):
    det = name.split("__")[0]; lab = name.split("__")[1]
    cfg = [c for c in configs()[det] if c[0] == lab][0]; _, _, vtree, (nfn, ntree), vfn = cfg
    pot = lambda fn, t: float(uproot.open(fn)[f"events/{t}/POT"].values()[0])
    pn, pv = (pot(nfn, ntree), pot(vfn, vtree)) if potmode == "real" else (1e19, 1e19)
    d = Q3 / name / f"U_{potmode}"; d.mkdir(parents=True, exist_ok=True); edges = " ".join(f"{x:g}" for x in BINS)
    xml = f'''<?xml version="1.0" ?>
<mode name="nu" />
<detector name="SBND" pot="1.0e+19" />
<model tag="numudis">
    <rule index="0" name="No Osc"/>
    <rule index="1" name="Numu Dis"/>
    <parameter name="L/E" variable_index="1"/>
</model>
<channel name="numu" unit="E_vis [GeV]">
    <bins unit="RecoE" edges="{edges}"/>
    <bins unit="True L/E" min="0.0" max="3.0" nbins="300" plot="false"/>
    <subchannel name="sel" plotname="Selected" color="#99CCFF"/>
</channel>
<DetVarFiles>
    <DetVarSection>
        <cv filename="{nfn}" treename="events/{ntree}/selected" pot="{pn:.9e}"/>
        <variation name="Q3_native" filename="{vfn}" treename="events/{vtree}/selected" pot="{pv:.9e}" knobval="1"/>
        <subchannel>nu_SBND_numu_sel</subchannel>
    </DetVarSection>
</DetVarFiles>
<variation_list>
    <allowlist type="mcstat" plotname="MC Stats" tag="stat,MC Stats">mcstat</allowlist>
    <allowlist type="spline" binning="var0" plotname="Q3 native" tag="DetVar">Q3_native</allowlist>
</variation_list>
<MCFile treename="events/{ntree}/selected" filename="{nfn}" scale="1.0" pot="{pn:.9e}">
    <branch associated_subchannel="nu_SBND_numu_sel" model_rule="0" additional_weight="{cutformula(det)}" incl_systematics="false">
        <variable>reco_visible_energy</variable>
        <variable>1.0</variable>
    </branch>
</MCFile>
'''
    (d / "q3.xml").write_text(xml); (d / "q3.ini").write_text(f"xml = {d}/q3.xml\ntag = q3\nseed = 1\nnthread = 4\nverbosity = 3\nlog = q3.log\nfile-verbosity = 2\nuse-fake-data = true\n")
    for cmd, pre in ((["process"], []), (["plot"], ["-o", "nom"]), (["plot"], ["-o", "inj", "--inject-systs-cv", "Q3_native", "1.0"])):
        with open(d / "steps.log", "ab") as fh: rc = subprocess.run([str(PROFIT), "--config", str(d / "q3.ini")] + pre + cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT, timeout=7200).returncode
        print(name, potmode, cmd[0], "rc", rc, flush=True)
        if rc: return
if __name__ == "__main__": run_full(sys.argv[1], sys.argv[2])
