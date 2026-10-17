"""Single source of truth for the DENT study report: an ordered list of ('md'|'code', text) cells.
build_report.py turns it into plotting/dent_study.ipynb and into a static HTML (executing the code cells)."""
import textwrap, pathlib
CELLS = []
def md(t): CELLS.append(("md", textwrap.dedent(t).strip("\n")))
def mdfile(name): CELLS.append(("md", (pathlib.Path(__file__).parent / "text" / f"{name}.md").read_text().strip("\n")))
def code(t): CELLS.append(("code", textwrap.dedent(t).strip("\n")))

mdfile("00_summary")
mdfile("01_data")

# ------------------------------------------------------------------------------------------------ setup
code('''
import sys, json, pathlib
HERE = next(p for p in (pathlib.Path.cwd() / d for d in ("dent_study", "plotting/dent_study", ".")) if (p / "dentstudy.py").exists())
sys.path.insert(0, str(HERE))
import numpy as np, pandas as pd
import dentstudy as ds
from mutools.plotting import use_style
use_style('rootlike')
REPO = ds.REPO
FIG = REPO / "docs" / "dent_study" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
res = ds.load_results(); tc = ds.load_tierC()
import re
LABEL = {n: o["label"].replace("$", "").replace("\\\\", "").replace("_{", "_").replace("}", "") for n, o in ds.OBS.items()}
def lab(s): return LABEL.get(s, s.replace("_", " "))
PRIM = [p for p in ds.PRIMARY if p in ds.OBS]
print(f"{res.tag.nunique()} study runs loaded; primary set = {len(PRIM)} observables; Tier C throws loaded for {tc.study.nunique() if len(tc) else 0} studies")
''')

mdfile("02_qc")
# ------------------------------------------------------------------------------------------------ QC
code('''
# Quality control per study: shape-only live, ndf = nbins - 1 (area-normalised), every step finished.
qc = res[res.variant == "v1"].groupby("tag").agg(det=("det", "first"), tier=("tier", "first"), live=("shapeonly_live", "first"),
        ndf_ok=("ndf", lambda s: True), n_rows=("chi2", "size")).reset_index()
bad_ndf = res[(res.variant == "v1") & (res.ndf != res.nbins - 1)]
print(f"studies with results: {len(qc)}  |  shape-only verified live (chi2 differs without the flag): {int(qc.live.fillna(False).sum())} of {len(qc)}")
print(f"rows where ndf != nbins-1 (would indicate no area normalisation): {len(bad_ndf)}")
summ = qc.groupby(["det", "tier"]).agg(studies=("tag", "size"), shapeonly_live=("live", lambda s: int(s.fillna(False).sum()))).reset_index()
ds.show_table(summ, "Studies by detector and tier")
# Audit: a DENT injection that PROfit silently ignored leaves chi2 identical to v1. List any such study.
w = res[res.variant.isin(["v1", "inj1p0"])].pivot_table(index=["tag", "channel"], columns="variant", values="chi2")
if "inj1p0" in w:
    inert = w[(w.v1 - w.inj1p0).abs() < 1e-9].reset_index().tag.unique().tolist()
    print(f"studies whose DENT injection left chi2 unchanged (invalid, must be re-run): {inert if inert else 'none'}")
''')

mdfile("03_mechanism")
code("""
for f, cap in (("mechanism_missing_length.png", "Track length missing from the flipped crossers: (a) distribution, (b) against the length beyond the cathode, (c) where the CV track ends, (d) probability of flipping."),
               ("mechanism_endpoint_cut.png", "(a) Distance of the muon endpoint from the cathode for the flipped events and for non-crossers. (b) DENT/CV ratio of selected events against reconstructed visible energy (MC only): without the crosser cut, with it, and with it plus |x_end| >= 5 cm. The extra cut leaves the ratio within 1.6 % of unity.")):
    ds.show_figure(FIG / f, cap)
""")

md("""
## 4. Tier A: the detvar samples against the data
""")

# ------------------------------------------------------------------------------------------------ Tier A
code('''
A = res[(res.tier == "A") & (res.variant == "v1")]
cv = A[A.which == "cv"].set_index("study"); dn = A[A.which == "dent"].set_index("study")
ta = pd.DataFrame({"chi2_cv": cv.chi2, "chi2_dent": dn.chi2, "ndf": cv.ndf}).dropna()
ta["dchi2"] = ta.chi2_cv - ta.chi2_dent; ta["label"] = [lab(s) for s in ta.index]; ta["primary"] = ta.index.isin(PRIM)
ta = ta.sort_values("dchi2", ascending=False)
tab = ta.reset_index()[["study", "primary", "ndf", "chi2_cv", "chi2_dent", "dchi2"]]
ds.show_table(tab.rename(columns={"study": "observable", "chi2_cv": "chi2 (detvar CV)", "chi2_dent": "chi2 (detvar DENT)", "dchi2": "chi2(CV) - chi2(DENT)"}), "Tier A: detvar samples used directly (MC statistics only)")
fig = ds.plot_tierA(ta[ta.primary], savepath=FIG / "tierA_primary.png", title="Tier A (detvar CV vs DENT)")
ds.show_figure(FIG / "tierA_primary.png", "Tier A, primary observables. Bars right of zero: the DENT sample agrees better with the data.")
''')

mdfile("04_tierA")

mdfile("05_tierB")
# ------------------------------------------------------------------------------------------------ Tier B
code('''
B = res[(res.tier == "B") & (res.channel == 0) & res.variant.isin(list(ds.VARIANTS_B))]
tb = pd.DataFrame(columns=["observable", "dchi2 (0 vs 1)"]); tcs = pd.DataFrame()
if B.empty: print("no Tier B results yet")
wide = B.pivot_table(index="study", columns="variant", values=["chi2", "p"]); ndf = B[B.variant == "v1"].set_index("study").ndf
tb = pd.DataFrame({"ndf": ndf})
for v in ds.VARIANTS_B:
    if ("chi2", v) in wide.columns: tb[f"chi2 {v}"] = wide[("chi2", v)]; tb[f"p {v}"] = wide[("p", v)]
if "chi2 v1" in tb and "chi2 inj1p0" in tb: tb["dchi2 (0 vs 1)"] = tb["chi2 v1"] - tb["chi2 inj1p0"]
tb = tb.reindex([p for p in PRIM if p in tb.index] + [s for s in tb.index if s not in PRIM]).reset_index().rename(columns={"study": "observable"})
tb.insert(1, "primary", tb.observable.isin(PRIM))
if len(tb): ds.show_table(tb, "Tier B: production CV (cvrw + cvwgt), full systematics, area-normalised and shape-only")
prim_b = B[B.study.isin(PRIM)]
if len(prim_b):
    ds.plot_pvalues(prim_b, [p for p in PRIM if p in set(prim_b.study)], [lab(p) for p in PRIM if p in set(prim_b.study)], savepath=FIG / "tierB_pvalues_primary.png")
    ds.show_figure(FIG / "tierB_pvalues_primary.png", "Tier B shape-only p-values for the primary observables, four prediction variants.")
''')

mdfile("06_ablation")
code(r"""
# Covariance ablation: how much of the DENT effect survives as other systematics are added. All configurations keep DENT in the
# covariance (so the injection is valid) and compare theta = 0 with theta = 1.0. Fixed configurations only; nothing is fitted.
CFG = {"stat + DENT": ("abl_statDENT_t0", "abl_statDENT_t1"), "stat + all detector": ("abl_det_t0", "abl_det_t1"),
       "all but other detector": ("abl_nodetother_t0", "abl_nodetother_t1"), "full (standard)": ("v1", "inj1p0")}
Bc = res[(res.tier == "B") & (res.channel == 0)]
rows = []
for name in [p for p in PRIM if p in set(Bc.study)] + [s for s in sorted(set(Bc.study)) if s not in PRIM]:
    g = Bc[Bc.study == name].set_index("variant").chi2
    r = {"observable": name, "primary": name in PRIM}
    for cfg, (a, b) in CFG.items():
        if a in g.index and b in g.index: r[f"{cfg}: chi2(0)"] = g[a]; r[f"{cfg}: chi2(1)"] = g[b]; r[f"{cfg}: dchi2"] = g[a] - g[b]
    rows.append(r)
abl = pd.DataFrame(rows)
if "stat + DENT: dchi2" in abl:
    ds.show_table(abl[["observable", "primary"] + [c for c in abl.columns if c.endswith("dchi2")]], "Covariance ablation: dchi2 = chi2(theta=0) - chi2(theta=1) (positive: the DENT shift helps), by covariance content")
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis
    ds.style(); pr = abl[abl.primary & abl["stat + DENT: dchi2"].notna()]
    if len(pr):
        fig, ax = plt.subplots(figsize=(11, 4.8)); x = np.arange(len(pr)); w = 0.2
        for i, (cfg, col) in enumerate(zip(CFG, ["C2", "C3", "C0", "C1"])):
            if f"{cfg}: dchi2" in pr: ax.bar(x + (i - 1.5) * w, pr[f"{cfg}: dchi2"], w, color=col, label=cfg)
        ax.axhline(0, color="0.3", lw=1); ax.set_yscale("symlog", linthresh=5); ax.set_xticks(x); ax.set_xticklabels([lab(o) for o in pr.observable], rotation=35, ha="right")
        ax.set_ylabel(r"$\chi^2(\theta{=}0)-\chi^2(\theta{=}1)$"); ax.legend(fontsize=9, ncol=2); mark_axis(ax, r"$\bf{SBND}$ Data Access", vadj=0.05)
        fig.savefig(FIG / "tierB_ablation.png", bbox_inches="tight", dpi=200); ds.show_figure(FIG / "tierB_ablation.png", "Tier B covariance ablation for the primary observables.")
else:
    print("ablation results not available yet")
""")

md("""
### Panels
PROfit-style panels for the ten primary observables (production CV, and the CV shifted by full DENT; area-normalised, shape-only).
""")

code('''
# PROfit-style panels (mutools, same configuration as plotting/studies/data_access.py) for the primary observables
import glob
made = []
for name in PRIM:
    o = ds.OBS[name]; e = o["edges"]
    for variant, tag_ in (("v1", "CV"), ("inj1p0", "CV + DENT")):
        f = ds.RUN / f"dent_sbnd_run1_dd_{name}_{variant}_PROplot_histograms.root"
        if not f.exists(): continue
        p = ds.render_panel(f, FIG / "panels", f"{name}_{variant}", lab(name), xlim=[e[0], e[-1]]); made.append((name, variant, p))
for name, variant, p in made:
    ds.show_figure(p, f"{lab(name)}: " + ("production CV" if variant == "v1" else "production CV shifted by 1.0 x DENT") + " (area-normalised, shape-only)")
print(len(made), "panels rendered")
''')

md("""
## 6. What the cathode and calorimetry results mean together
Tier A says the DENT sample matches the data's calorimetry and PID shapes better but its cathode geometry worse. Tier B says that, once the
other systematics are included, shifting the CV by DENT moves the chi-square by only a unit or two. These are consistent: the calorimetry/PID gain is covered by other detector
systematics, and the cathode-geometry loss is small in Tier B only because DENT is itself in the covariance. Section 7 asks the sharper question.
""")

mdfile("07_replica")

code('''
import replica as rp
rs = pd.read_csv(REPO / "docs" / "dent_study" / "replica_summary.csv")
rs = rs.set_index("study").loc[[p for p in PRIM if p in set(rs.study)]]
dmax = 0.0
for name in rs.index:
    v = rp.validate(rp.load(name)); dmax = max(dmax, max(abs(x) for x in v.values()))
print(f"replica versus PROfit chi-square, largest absolute difference over {len(rs)} primaries and 3 strengths: {dmax:.3f}")
t = rs.reset_index().copy(); t.insert(1, "label", [lab(s) for s in t.study])
t["nearer truth CV?"] = (abs(t.dchi2_fixed_data - t.mean_0) < abs(t.dchi2_fixed_data - t.mean_1))
ds.show_table(t[["label", "profit_dchi2", "dchi2_fixed_data", "mean_0", "sd_0", "p2_0", "mean_0.5", "sd_0.5", "p2_0.5", "mean_1", "sd_1", "p2_1", "dprime", "nearer truth CV?"]].rename(columns={
    "profit_dchi2": "PROfit change", "dchi2_fixed_data": "fixed-cov change (data)", "mean_0": "mean|CV", "sd_0": "sd|CV", "p2_0": "p2|CV", "mean_0.5": "mean|half", "sd_0.5": "sd|half", "p2_0.5": "p2|half", "mean_1": "mean|DENT", "sd_1": "sd|DENT", "p2_1": "p2|DENT", "dprime": "d'"}),
    "Fixed-covariance chi2(0) - chi2(1) in data, against 20,000 pseudo-experiments for truth = CV, half DENT, full DENT (two-sided p)")
print(f"data nearer the CV-truth mean than the DENT-truth mean: {int(t['nearer truth CV?'].sum())} of {len(t)}; smallest two-sided p over all truths: {t[['p2_0','p2_0.5','p2_1']].min().min():.3f}; d' range {t.dprime.min():.2f} to {t.dprime.max():.2f}")
import matplotlib.pyplot as plt
from mutools.plotting.helpers import mark_axis
ds.style(); fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(t) + 1.8)); y = np.arange(len(t))[::-1]
for off, col, key, nm in ((0.2, "C0", "0", "truth CV"), (0.0, "C2", "0.5", "truth half DENT"), (-0.2, "C1", "1", "truth full DENT")):
    ax.errorbar(t[f"mean_{key}"], y + off, xerr=t[f"sd_{key}"], fmt="o", color=col, ms=5, capsize=2, label=nm + " (mean +- 1 sd)")
ax.plot(t.dchi2_fixed_data, y, "k*", ms=11, label="data", zorder=5)
ax.axvline(0, color="0.4", lw=1); ax.set_yticks(y); ax.set_yticklabels(t.label)
ax.set_xlabel(r"fixed-covariance $\\chi^2(\\theta{=}0)-\\chi^2(\\theta{=}1)$" + "\\n[positive: DENT-shifted prediction closer to data]"); ax.legend(fontsize=8, loc="lower right")
mark_axis(ax, r"$\\bf{SBND}$ Data Access", vadj=0.05); fig.savefig(FIG / "replica_dchi2.png", bbox_inches="tight", dpi=200)
ds.show_figure(FIG / "replica_dchi2.png", "Data against the pseudo-experiment expectations for the three fixed truths (fixed covariance).")
''')

mdfile("08_replica_after")
mdfile("09_tierC")

# ------------------------------------------------------------------------------------------------ Tier C
code('''
rows = []
for name in PRIM:
    sub = tc[tc.study == name] if len(tc) else tc
    row_b = tb[tb.observable == name] if "observable" in tb else tb.iloc[0:0]
    if len(sub) == 0 or len(row_b) == 0: continue
    d_obs = float(row_b["dchi2 (0 vs 1)"].iloc[0]); r = ds.tierC_pvalues(sub, d_obs); r.update(study=name, label=lab(name), dchi2_data=d_obs); rows.append(r)
tcs = pd.DataFrame(rows)
if len(tcs):
    tcs["p_cv_holm"] = ds.holm(tcs.p_cv); tcs["p_dent_holm"] = ds.holm(tcs.p_dent)
    ds.show_table(tcs[["label", "n", "dchi2_data", "mean_cv", "sd_cv", "mean_dent", "sd_dent", "p_cv", "p_cv_holm", "p_dent", "p_dent_holm"]],
                  "Tier C: Delta chi2 = chi2(theta=0) - chi2(theta=1). p_CV = P(dchi2 >= data | truth CV); p_DENT = P(dchi2 <= data | truth DENT)")
    ds.plot_dchi2(tcs, savepath=FIG / "tierC_dchi2.png"); ds.show_figure(FIG / "tierC_dchi2.png", "Tier C: data Delta chi2 against the fake-data expectations.")
else:
    print("no Tier C throws available yet")
''')

code('''
import matplotlib.pyplot as plt
if len(tcs):
    n = len(tcs); ncol = 4; nrow = int(np.ceil(n / ncol)); fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow)); axes = np.atleast_1d(axes).ravel()
    for ax, (_, r) in zip(axes, tcs.iterrows()):
        ds.plot_tierC_hist(tc[tc.study == r.study], r.dchi2_data, r.label, ax=ax)
    for ax in axes[n:]: ax.axis("off")
    axes[0].legend(fontsize=8); fig.tight_layout(); fig.savefig(FIG / "tierC_hists.png", dpi=170, bbox_inches="tight")
    ds.show_figure(FIG / "tierC_hists.png", "Tier C: distributions of Delta chi2 for fake data from CV (blue) and from CV + DENT (orange); the black line is the data.")
''')

mdfile("10_baseline")
# ------------------------------------------------------------------------------------------------ baseline regeneration
code('''
base = res[(res.tier == "baseline") & res.variant.isin(["v1", "exDENT", "inj1p0"])]
rows = []
for (det, study, ch), g in base.groupby(["det", "study", "channel"]):
    d = g.drop_duplicates("variant", keep="last").set_index("variant"); r = dict(detector=det, study=study, channel=ch, ndf=int(d.loc["v1", "ndf"]), chi2=d.loc["v1", "chi2"], p=d.loc["v1", "p"])
    if "inj1p0" in d.index: r["chi2 CV+DENT"] = d.loc["inj1p0", "chi2"]; r["dchi2 (0 vs 1)"] = d.loc["v1", "chi2"] - d.loc["inj1p0", "chi2"]
    if "exDENT" in d.index: r["chi2 no DENT cov"] = d.loc["exDENT", "chi2"]
    rows.append(r)
bt = pd.DataFrame(rows)
ds.show_table(bt[bt.detector == "sbnd_run1"].sort_values("p"), "Regenerated SBND detector-access studies (shape-only), sorted by p-value")
ds.show_table(bt[bt.detector == "icarus_run2"].sort_values("p"), "Regenerated ICARUS Run 2 detector-access studies (shape-only), sorted by p-value")
''')

mdfile("11_limits")
mdfile("12_deviations")
