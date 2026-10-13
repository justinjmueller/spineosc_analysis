#!/usr/bin/env python3
"""Q3: multi-universe DetVar (SBND SCE0x = knob -1, SCE2x = knob +1) through matched mode with the selected-attribute (eps weight) construction on unmodified PROfit.
The matched key set is the intersection of ALL samples (CV, -1, +1), i.e. a common pool, unlike production where each variation has its own pairwise pool."""
import sys, subprocess, json, numpy as np, pandas as pd, uproot, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from q3_native import *
name = "sbnd_run1__SCE_pm"; d = Q3 / name; d.mkdir(parents=True, exist_ok=True); edges = " ".join(f"{x:g}" for x in BINS)
cf = {c[0]: c for c in configs()["sbnd_run1"]}; nfn, ntree = cf["SCE0x"][3]
ek = {"n": ekeys(nfn, ntree).drop_duplicates()}
for t, lab in (("0", "SCE0x"), ("2", "SCE2x")): ek[t] = ekeys(cf[lab][4], cf[lab][2]).drop_duplicates()
common = ek["n"].merge(ek["0"], on=KEYS + ["Ei"]).merge(ek["2"], on=KEYS + ["Ei"]); print("common pool", len(common), "pairwise pools", len(ek["n"].merge(ek["0"], on=KEYS + ["Ei"])), len(ek["n"].merge(ek["2"], on=KEYS + ["Ei"])), flush=True)
tabs = {t: pd.read_pickle(PAIRS / f"sbnd_run1__{lab}.pkl")[KEYS + ["Ei", "in_n", "in_v", "n_reco_visible_energy", "v_reco_visible_energy"]] for t, lab in (("0", "SCE0x"), ("2", "SCE2x"))}
S = {}
for t in tabs:
    m = tabs[t].merge(common, on=KEYS + ["Ei"], how="inner"); S[t] = m                                  # rows of this table inside the common pool
cvsel = pd.concat([S[t][S[t].in_n][KEYS + ["Ei", "n_reco_visible_energy"]] for t in S]).drop_duplicates(KEYS + ["Ei"])
rows = pd.concat([cvsel[KEYS + ["Ei"]]] + [S[t][S[t].in_v][KEYS + ["Ei"]] for t in S]).drop_duplicates(KEYS + ["Ei"]).reset_index(drop=True)
def side(tab, flagcol, ecol):
    t = tab[tab[flagcol]][KEYS + ["Ei", ecol]].rename(columns={ecol: "E"}); r = rows.merge(t, on=KEYS + ["Ei"], how="left"); sel = r.E.notna().to_numpy(); return sel, np.nan_to_num(r.E.to_numpy()).astype(np.float32)
cvs = pd.concat([S[t][S[t].in_n][KEYS + ["Ei", "n_reco_visible_energy"]] for t in S]).drop_duplicates(KEYS + ["Ei"]); r = rows.merge(cvs, on=KEYS + ["Ei"], how="left"); sel_n = r.n_reco_visible_energy.notna().to_numpy(); e_n = np.nan_to_num(r.n_reco_visible_energy.to_numpy()).astype(np.float32)
sel_0, e_0 = side(S["0"], "in_v", "v_reco_visible_energy"); sel_2, e_2 = side(S["2"], "in_v", "v_reco_visible_energy"); print("rows", len(rows), sel_n.sum(), sel_0.sum(), sel_2.sum(), flush=True)
fn = d / "trees.root"; k = {c: rows[c].to_numpy().astype(np.int32) for c in KEYS + ["Ei"]}
with uproot.recreate(fn) as f:
    for tn, sel, e in (("cvU", sel_n, e_n), ("m1U", sel_0, e_0), ("p1U", sel_2, e_2)):
        br = dict(k); br["reco_visible_energy"] = np.where(sel, e, 0).astype(np.float32); br["selected"] = sel.astype(np.int32); t = f.mktree(tn, {c: v.dtype for c, v in br.items()}); t.extend(br)
xml = f'''<?xml version="1.0" ?>
<mode name="nu" />
<detector name="SBND" pot="1.0e+19" />
<model tag="numudis"><rule index="0" name="No Osc"/><rule index="1" name="Numu Dis"/><parameter name="L/E" variable_index="1"/></model>
<channel name="numu" unit="E_vis [GeV]">
    <bins unit="RecoE" edges="{edges}"/>
    <bins unit="True L/E" min="0.0" max="3.0" nbins="300" plot="false"/>
    <subchannel name="sel" plotname="Selected" color="#99CCFF"/>
</channel>
<DetVarFiles>
    <DetVarSection filename="{fn}" cv_variation_matching_vars="Run,Subrun,Evt,Ei">
        <cv treename="cvU" pot="1.0e+19"/>
        <variation name="Q3_native" treename="m1U" pot="1.0e+19" knobval="-1"/>
        <variation name="Q3_native" treename="p1U" pot="1.0e+19" knobval="1"/>
        <subchannel>nu_SBND_numu_sel</subchannel>
    </DetVarSection>
</DetVarFiles>
<variation_list>
    <allowlist type="mcstat" plotname="MC Stats" tag="stat,MC Stats">mcstat</allowlist>
    <allowlist type="spline" binning="var0" plotname="Q3 native" tag="DetVar">Q3_native</allowlist>
</variation_list>
<MCFile treename="cvU" filename="{fn}" scale="1.0" pot="1.0e+19">
    <branch associated_subchannel="nu_SBND_numu_sel" model_rule="0" additional_weight="selected + 1e-9*(1-selected)" incl_systematics="false">
        <variable>reco_visible_energy</variable>
        <variable>1.0</variable>
    </branch>
</MCFile>
'''
(d / "q3.xml").write_text(xml); (d / "q3.ini").write_text(f"xml = {d}/q3.xml\ntag = q3\nseed = 1\nnthread = 4\nverbosity = 3\nlog = q3.log\nfile-verbosity = 2\nuse-fake-data = true\n")
steps = [(["process"], [])] + [(["plot"], ["-o", f"s{tag}"] + (["--inject-systs-cv", "Q3_native", str(v)] if v is not None else [])) for tag, v in (("nom", None), ("m1", -1.0), ("p1", 1.0), ("p05", 0.5), ("m05", -0.5))]
for cmd, pre in steps:
    with open(d / "steps.log", "ab") as fh: rc = subprocess.run([str(PROFIT), "--config", str(d / "q3.ini")] + pre + cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
    print(cmd[0], pre[:2], rc, flush=True)
def get(tag):
    a = uproot.open(d / f"q3_s{tag}_PROplot_histograms.root")["hist1d"].arrays(library="np"); m = (a["variable"] == 0) & (a["subchannel"] == 4294967295) & np.array([str(x) == "CV" for x in a["prefix"]]); return a["bin_content"][m].astype(float)[:len(BINS) - 1]
nom = get("nom"); res = {t: (get(t) / nom).tolist() for t in ("m1", "p1", "p05", "m05")}
# numpy common-pool ratios
def hist(sel, e): return np.histogram(e[sel], BINS)[0].astype(float)
Hn, H0, H2 = hist(sel_n, e_n), hist(sel_0, e_0), hist(sel_2, e_2); res["np_m1"] = (H0 / Hn).tolist(); res["np_p1"] = (H2 / Hn).tolist()
E_ = json.load(open(Q3 / "emulation.json")); res["prod_m1"] = E_["sbnd_run1__SCE0x"]["prod"]; res["prod_p1"] = E_["sbnd_run1__SCE2x"]["prod"]; res["pair_m1"] = E_["sbnd_run1__SCE0x"]["pair"]; res["pair_p1"] = E_["sbnd_run1__SCE2x"]["pair"]
res["n_common"] = len(common); res["n_rows"] = len(rows); json.dump(res, open(d / "result.json", "w"))
a = lambda k: np.array(res[k]); print("native vs numpy common pool: m1", np.abs(a("m1") - a("np_m1")).max(), "p1", np.abs(a("p1") - a("np_p1")).max())
print("native vs production: m1", np.abs(a("m1") - a("prod_m1")).max(), "p1", np.abs(a("p1") - a("prod_p1")).max()); print("common-pool numpy vs pairwise-pool pair: m1", np.abs(a("np_m1") - a("pair_m1")).max(), "p1", np.abs(a("np_p1") - a("pair_p1")).max())
print("p05", np.round(a("p05")[:5], 4), "m05", np.round(a("m05")[:5], 4), "midpoint of knots", np.round((1 + a("p1")[:5]) / 2, 4))
