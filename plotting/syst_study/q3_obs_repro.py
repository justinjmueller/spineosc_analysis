#!/usr/bin/env python3
"""Q3: reproduction of the PER-OBSERVABLE production ratio histograms (detvar/detsys_v1.1.0v2_<obs>.root; calculator cut = fiducial && veto, no SBND crosser cut) with
 (i) numpy emulation (pool-paired / matched-only / unmatched full trees), (ii) native PROfit in the observable's own binning (A = pool-restricted unmatched; Ce = selected-attribute matched, unmodified PROfit;
 Bm = matched on selected-only trees). cvwgt is not used.   usage: q3_obs_repro.py [obs ...]"""
import sys, json, subprocess, numpy as np, pandas as pd, uproot, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent)); from q3_native import *
OBS = {"muon_score": dict(col="reco_leading_muon_muon_softmax", fn="detvar/detsys_v1.1.0v2_muon_score.root", dets=("sbnd_run1", "icarus_run2")),
       "muon_endx_sbnd": dict(col="reco_leading_muon_end_x", fn="detvar/detsys_v1.1.0v2_muon_endx_sbnd.root", dets=("sbnd_run1",))}
CASES = {"muon_score": ["sbnd_run1__ChargeScale15", "sbnd_run1__DENT", "sbnd_run1__GainHigh", "sbnd_run1__SCE2x", "sbnd_run1__WireModYZ", "icarus_run2__SCE2x", "icarus_run2__ChargeScale", "icarus_run2__WireModYZ"],
         "muon_endx_sbnd": ["sbnd_run1__ChargeScale15", "sbnd_run1__DENT", "sbnd_run1__SCE2x", "sbnd_run1__GainHigh"]}
B = 200; rng = np.random.default_rng(11); cache = {}
def load(fn, tree, col, det):
    k = (str(fn), tree, col, det)
    if k not in cache:
        t = uproot.open(fn)[f"events/{tree}/selected"]; cols = sorted(set(KEYS + [EKEY, col, "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x"]) & set(t.keys())); d = t.arrays(cols, library="np")
        m = selection_mask(det, d)   # the v1.1.0v2 per-observable ratio files use the same cut as by_vis_energy (SBND: + crosser veto and |end_x| >= 4)
        df = pd.DataFrame({c: d[c][m] for c in KEYS}); df["Ei"] = np.round(d[EKEY][m] * 1e6).astype(np.int64); df["x"] = d[col][m].astype(np.float32)
        cache[k] = (df, float(uproot.open(fn)[f"events/{tree}/POT"].values()[0]), uproot.open(fn)[f"events/{tree}/events"].num_entries)
    return cache[k]
def ratio(a, b): return np.where(b > 0, a / np.where(b > 0, b, 1), 1.0)
res = json.load(open(Q3 / "obs_repro.json")) if (Q3 / "obs_repro.json").exists() else {}
for obs in (sys.argv[1:] or OBS):
    o = OBS[obs]; prodf = uproot.open(REPO / o["fn"])
    for name in CASES[obs]:
        det, lab = name.split("__"); cfg = [c for c in configs()[det] if c[0] == lab][0]; _, hname, vtree, (nfn, ntree), vfn = cfg
        key = f"{name}|{obs}"
        if hname not in prodf: print("no production hist", key); continue
        h = prodf[hname]; edges = h.axes[0].edges(); prod = h.values(); nb = len(prod)
        (dn, potn, nevn), (dv, potv, nevv) = load(nfn, ntree, o["col"], det), load(vfn, vtree, o["col"], det)
        pool = ekeys(nfn, ntree).drop_duplicates().merge(ekeys(vfn, vtree).drop_duplicates(), on=KEYS + ["Ei"])
        sn = dn.merge(pool, on=KEYS + ["Ei"]); sv = dv.merge(pool, on=KEYS + ["Ei"]); k4 = KEYS + ["Ei"]
        ndup_n = sn.duplicated(k4, keep=False).sum(); ndup_v = sv.duplicated(k4, keep=False).sum()
        un = sn[~sn.duplicated(k4, keep=False)]; uv = sv[~sv.duplicated(k4, keep=False)]
        P = un.rename(columns={"x": "xn"}).merge(uv.rename(columns={"x": "xv"}), on=k4, how="outer"); inn = P.xn.notna().to_numpy(); inv = P.xv.notna().to_numpy()
        xn = np.nan_to_num(P.xn.to_numpy()).astype(np.float32); xv = np.nan_to_num(P.xv.to_numpy()).astype(np.float32); both = inn & inv
        def bi(x, sel): i = np.digitize(x, edges) - 1; i[(i < 0) | (i >= nb)] = -1; i[x == edges[-1]] = nb - 1; return np.where(sel, i, -1)
        bn, bv = bi(xn, inn), bi(xv, inv)
        def H(b, w): m = b >= 0; return np.bincount(b[m], weights=w[m], minlength=nb)
        def meth(w): Hn, Hv = H(bn, w), H(bv, w); return ratio(Hv, Hn), ratio(H(np.where(both, bv, -1), w), H(np.where(both, bn, -1), w)), Hn
        pair, matched, Hn = meth(np.ones(len(P))); reps = np.array([meth(rng.poisson(1.0, len(P)).astype(float))[0] for _ in range(B)]); sig = reps.std(0, ddof=1)
        Fn = np.histogram(dn.x, edges)[0].astype(float); Fv = np.histogram(dv.x, edges)[0].astype(float); ufull = ratio(Fv / potv, Fn / potn)
        stat = lambda r: dict(max_abs=float(np.abs(r - prod).max()), rms=float(np.sqrt(((r - prod) ** 2).mean())), chi2=float((((r - prod) / np.maximum(sig, 1e-9)) ** 2).sum()), ndf=nb)
        # production-style ratio (repeated keys kept) for the exact calculator definition
        Hn_d = np.histogram(sn.x, edges)[0].astype(float); Hv_d = np.histogram(sv.x, edges)[0].astype(float); prodlike = ratio(Hv_d, Hn_d)
        # native PROfit
        d = Q3 / name / f"obs_{obs}"; d.mkdir(parents=True, exist_ok=True); fn_ = d / "trees.root"; kk = {c: P[c].to_numpy().astype(np.int32) for c in k4}
        with uproot.recreate(fn_) as f:
            for tn, sel, x, rows in (("cvA", inn, xn, inn), ("varA", inv, xv, inv), ("cvU", inn, xn, np.ones(len(P), bool)), ("varU", inv, xv, np.ones(len(P), bool))):
                br = {c: v[rows] for c, v in kk.items()}; br["obs"] = np.where(sel, x, edges[0] - 1e3).astype(np.float32)[rows]; br["selected"] = sel.astype(np.int32)[rows]
                t = f.mktree(tn, {c: v.dtype for c, v in br.items()}); t.extend(br)
        nat = {}
        for variant in ("A", "Bm", "Ce"):
            trn = ("cvA", "varA") if variant in ("A", "Bm") else ("cvU", "varU"); vd = d / variant; vd.mkdir(exist_ok=True)
            wf = "selected + 1e-9*(1-selected)" if variant == "Ce" else "selected"; mv = "" if variant == "A" else f' cv_variation_matching_vars="{KEYS_}"'
            ed = " ".join(f"{x:.9g}" for x in edges)
            xml = f'''<?xml version="1.0" ?>
<mode name="nu" />
<detector name="SBND" pot="1.0e+19" />
<model tag="numudis"><rule index="0" name="No Osc"/><rule index="1" name="Numu Dis"/><parameter name="L/E" variable_index="1"/></model>
<channel name="numu" unit="{obs}">
    <bins unit="obs" edges="{ed}"/>
    <bins unit="True L/E" min="0.0" max="3.0" nbins="300" plot="false"/>
    <subchannel name="sel" plotname="Selected" color="#99CCFF"/>
</channel>
<DetVarFiles>
    <DetVarSection filename="{fn_}"{mv}>
        <cv treename="{trn[0]}" pot="1.0e+19"/>
        <variation name="Q3_native" treename="{trn[1]}" pot="1.0e+19" knobval="1"/>
        <subchannel>nu_SBND_numu_sel</subchannel>
    </DetVarSection>
</DetVarFiles>
<variation_list>
    <allowlist type="mcstat" plotname="MC Stats" tag="stat,MC Stats">mcstat</allowlist>
    <allowlist type="spline" binning="var0" plotname="Q3 native" tag="DetVar">Q3_native</allowlist>
</variation_list>
<MCFile treename="{trn[0]}" filename="{fn_}" scale="1.0" pot="1.0e+19">
    <branch associated_subchannel="nu_SBND_numu_sel" model_rule="0" additional_weight="{wf}" incl_systematics="false">
        <variable>obs</variable>
        <variable>1.0</variable>
    </branch>
</MCFile>
'''
            (vd / "q3.xml").write_text(xml); (vd / "q3.ini").write_text(f"xml = {vd}/q3.xml\ntag = q3\nseed = 1\nnthread = 4\nverbosity = 3\nlog = q3.log\nfile-verbosity = 2\nuse-fake-data = true\n")
            (vd / "steps.log").unlink(missing_ok=True)
            for cmd, pre in ((["process"], []), (["plot"], ["-o", "nom"]), (["plot"], ["-o", "inj", "--inject-systs-cv", "Q3_native", "1.0"])):
                with open(vd / "steps.log", "ab") as fh: rc = subprocess.run([str(PROFIT), "--config", str(vd / "q3.ini"), "--force"] + pre + cmd, cwd=vd, stdout=fh, stderr=subprocess.STDOUT).returncode
                if rc: print("FAIL", key, variant, cmd, rc); break
            else:
                g = lambda t: (lambda a: a["bin_content"][(a["variable"] == 0) & (a["subchannel"] == 4294967295) & np.array([str(x) == "CV" for x in a["prefix"]])].astype(float)[:nb])(uproot.open(vd / f"q3_{t}_PROplot_histograms.root")["hist1d"].arrays(library="np"))
                nat[variant] = ratio(g("inj"), g("nom"))
        res[key] = dict(edges=edges.tolist(), prod=prod.tolist(), prodlike=prodlike.tolist(), pair=pair.tolist(), matched=matched.tolist(), unmatched_full_pot=ufull.tolist(), sig=sig.tolist(), Hn=Hn.tolist(),
                        native={k: v.tolist() for k, v in nat.items()}, stats=dict(prodlike=stat(prodlike), pair=stat(pair), matched=stat(matched), unmatched_full_pot=stat(ufull), **{"native_" + k: stat(v) for k, v in nat.items()}),
                        n=dict(both=int(both.sum()), only_n=int((inn & ~inv).sum()), only_v=int((~inn & inv).sum()), dup_rows_n=int(ndup_n), dup_rows_v=int(ndup_v)))
        print(key, {k: (round(v["max_abs"], 4), round(v["chi2"], 1)) for k, v in res[key]["stats"].items()}, flush=True)
        json.dump(res, open(Q3 / "obs_repro.json", "w"))
