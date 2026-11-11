#!/usr/bin/env python3
"""Q3: build intermediate trees (from the paired tables) and PROfit XML/ini for native DetVar tests; run PROfit; read back the response.
usage: q3_native.py <det__variation> <variant>   variants:
  A      pool-restricted selected trees, UNMATCHED mode (unmodified PROfit)
  Bm     matched mode on the same pool-restricted selected-only trees (selected-both only; the 'native matched' behaviour)
  Cz     union rows on both sides, `selected` as plain 0/1 weight, matched mode, UNMODIFIED PROfit (zero weights dropped -> expected failure)
  Ce     union rows on both sides, weight = selected + 1e-9*(1-selected), matched mode, UNMODIFIED PROfit
  Cp     like Cz (plain 0/1 weight) but with the PATCHED PROfit copy (profit-asymmetric-q3) that keeps zero-weight events when matching vars are set
"""
import sys, subprocess, pathlib, numpy as np, pandas as pd, uproot
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from common import *
Q3 = OUT / "q3"; PROFIT = pathlib.Path("/Users/mueller/Projects/profit-asymmetric/build/bin/PROfit"); PROFIT_PATCHED = pathlib.Path("/Users/mueller/Projects/profit-asymmetric-q3/build/bin/PROfit")
def make_trees(name):
    d = Q3 / name; d.mkdir(parents=True, exist_ok=True); fn = d / "trees.root"
    if fn.exists(): return fn
    P = pd.read_pickle(PAIRS / f"{name}.pkl"); inn = P.in_n.to_numpy(); inv = P.in_v.to_numpy()
    en = np.nan_to_num(P.n_reco_visible_energy.to_numpy()).astype(np.float32); ev = np.nan_to_num(P.v_reco_visible_energy.to_numpy()).astype(np.float32)
    k = {c: P[c].to_numpy().astype(np.int32) for c in ("Run", "Subrun", "Evt", "Ei")}
    def tree(m, e, sel, rows=None):
        r = np.ones(len(P), bool) if rows is None else rows
        out = {c: v[r] for c, v in k.items()}; out["reco_visible_energy"] = np.where(sel, e, 0).astype(np.float32)[r]; out["selected"] = sel.astype(np.int32)[r]; return out
    with uproot.recreate(fn) as f:
        for tn, br in (("cvA", tree(None, en, inn, inn)), ("varA", tree(None, ev, inv, inv)), ("cvU", tree(None, en, inn)), ("varU", tree(None, ev, inv))):   # A: selected-only; U: union rows with selected flag
            t = f.mktree(tn, {c: v.dtype for c, v in br.items()}); t.extend(br)           # real TTree (uproot's dict assignment would write an RNTuple)
    return fn
KEYS_ = "Run,Subrun,Evt,Ei"; NB_ = len(BINS) - 1
def make_xml(name, variant, outdir):
    fn = make_trees(name); tr = {"A": ("cvA", "varA"), "Bm": ("cvA", "varA")}.get(variant, ("cvU", "varU")); matched = variant != "A"
    wf = {"A": "selected", "Bm": "selected", "Cz": "selected", "Cp": "selected", "Ce": "selected + 1e-9*(1-selected)"}[variant]
    edges = " ".join(f"{x:g}" for x in BINS)
    mv = f' cv_variation_matching_vars="{KEYS_}"' if matched else ""
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
    <DetVarSection filename="{fn}"{mv}>
        <cv treename="{tr[0]}" pot="1.0e+19"/>
        <variation name="Q3_native" treename="{tr[1]}" pot="1.0e+19" knobval="1"/>
        <subchannel>nu_SBND_numu_sel</subchannel>
    </DetVarSection>
</DetVarFiles>
<variation_list>
    <allowlist type="mcstat" plotname="MC Stats" tag="stat,MC Stats">mcstat</allowlist>
    <allowlist type="spline" binning="var0" plotname="Q3 native" tag="DetVar">Q3_native</allowlist>
</variation_list>
<MCFile treename="{tr[0]}" filename="{fn}" scale="1.0" pot="1.0e+19">
    <branch associated_subchannel="nu_SBND_numu_sel" model_rule="0" additional_weight="{wf}" incl_systematics="false">
        <variable>reco_visible_energy</variable>
        <variable>1.0</variable>
    </branch>
</MCFile>
'''
    (outdir / "q3.xml").write_text(xml)
    (outdir / "q3.ini").write_text(f"xml = {outdir}/q3.xml\ntag = q3\nseed = 1\nnthread = 4\nverbosity = 3\nlog = q3.log\nfile-verbosity = 2\nuse-fake-data = true\n")
def run(name, variant):
    d = Q3 / name / variant; d.mkdir(parents=True, exist_ok=True); make_xml(name, variant, d)
    exe = PROFIT_PATCHED if variant == "Cp" else PROFIT
    steps = [(["process"], []), (["plot"], ["-o", "nom"]), (["plot"], ["-o", "inj", "--inject-systs-cv", "Q3_native", "1.0"])]
    for cmd, pre in steps:
        with open(d / "steps.log", "ab") as fh: rc = subprocess.run([str(exe), "--config", str(d / "q3.ini")] + pre + cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT, timeout=7200).returncode
        print(name, variant, cmd[0], pre[:2], "rc", rc, flush=True)
        if rc: break
def response(name, variant):
    """PROfit's response per bin: CV spectrum with Q3_native injected at +1 divided by the un-injected CV spectrum (spline at the knot = var/cv of the DetVar props)."""
    d = Q3 / name / variant; out = {}
    for t in ("nom", "inj"):
        a = uproot.open(d / f"q3_{t}_PROplot_histograms.root")["hist1d"].arrays(library="np"); m = (a["variable"] == 0) & (a["subchannel"] == 4294967295) & (a["prefix"] == "CV")
        out[t] = a["bin_content"][m].astype(float)[:NB_]   # the plot file lists the total CV twice; take the first set
    return np.where(out["nom"] > 0, out["inj"] / np.where(out["nom"] > 0, out["nom"], 1), 1.0), out["nom"], out["inj"]
if __name__ == "__main__": run(sys.argv[1], sys.argv[2])
