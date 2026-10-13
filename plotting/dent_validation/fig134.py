import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import dentval as dv, plotlib as pl
d = dv.load()
scan = pd.read_csv("scan_all_variables.csv")
NICE = {"reco_leading_muon_calo_csda_diff": "muon calo-vs-CSDA residual", "reco_leading_proton_calo_csda_diff": "proton calo-vs-CSDA residual",
        "reco_leading_muon_muon_softmax": "muon: muon score", "reco_leading_muon_proton_softmax": "muon: proton score", "reco_leading_muon_pion_softmax": "muon: pion score",
        "reco_leading_muon_primary_softmax": "muon: primary score", "reco_leading_proton_primary_softmax": "proton: primary score",
        "reco_leading_proton_secondary_softmax": "proton: secondary score", "reco_leading_proton_proton_softmax": "proton: proton score",
        "reco_leading_proton_mip_softmax": "proton: MIP score", "reco_leading_proton_muon_softmax": "proton: muon score", "reco_leading_proton_pion_softmax": "proton: pion score",
        "reco_leading_proton_hadron_softmax": "proton: hadron score", "abs_end_x": "|muon end x|", "muon_xing_geom": "muon vertex/end opposite sides",
        "reco_cathode_crosser": "cathode-crosser flag", "reco_neutrino2026_veto": "neutrino2026 veto (= not crosser)", "reco_leading_muon_end_x": "muon end x",
        "reco_vertex_x": "vertex x", "reco_leading_muon_start_x": "muon start x", "reco_leading_proton_start_x": "proton start x", "reco_leading_proton_end_x": "proton end x",
        "reco_leading_muon_start_dir_z": "muon direction z", "reco_leading_muon_start_dir_x": "muon direction x", "reco_vertex_y": "vertex y", "muon_len3d": "muon straight-line length"}
nice = lambda v: NICE.get(v, v.replace("reco_", "").replace("_", " "))

# ---------- Fig 1: ranked scan ----------
top = scan.drop_duplicates(subset=["sep", "chi_cv", "chi_dent"]).head(16).iloc[::-1]      # drop exact duplicate columns (flag / veto)
fig, ax = plt.subplots(1, 2, figsize=(12.5, 6.4), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1.25]})
y = np.arange(len(top)); lab = [f"{nice(v)}  (ndf {n})" for v, n in zip(top["var"], top["ndf"])]
ax[0].barh(y, top.sep / top.ndf, color="#7d7c78", height=0.62); ax[0].set_xscale("log")
ax[0].set_yticks(y); ax[0].set_yticklabels(lab, fontsize=8.8); ax[0].set_xlabel("CV vs DENT separation  χ²/ndf  (at data statistics)")
ax[0].set_title("(a) where DENT differs from CV", loc="left", fontsize=10.5)
h = 0.3
ax[1].barh(y + h / 2, top.chi_cv / top.ndf, height=h, color=pl.C_CV, label="data vs CV"); ax[1].barh(y - h / 2, top.chi_dent / top.ndf, height=h, color=pl.C_DENT, label="data vs DENT")
ax[1].set_xscale("log"); ax[1].set_yticks(y); ax[1].set_yticklabels([]); ax[1].set_xlabel("data vs MC  χ²/ndf  (shape, unit area)")
ax[1].set_title("(b) which one the data prefers (smaller is better)", loc="left", fontsize=10.5); ax[1].legend(loc="lower right")
for i, (c, e) in enumerate(zip(top.chi_cv / top.ndf, top.chi_dent / top.ndf)):
    ax[1].text(max(c, e) * 1.15, y[i], "DENT better" if e < c else "CV better", va="center", fontsize=8, color=pl.INK2)
fig.suptitle("DENT validation scan: every energy-free variable shared by data and the detvar MC, ranked by CV–DENT separation", x=0.005, ha="left", fontsize=11)
fig.savefig("fig1_scan_summary.png", dpi=150, facecolor=pl.SURF)

# ---------- helpers for grids ----------
def panel(fig, cell, col, logy=False, edges=None):
    a, r = pl.pair_axes(fig, cell)
    e = dv.edges_for(col, d["cv"]) if edges is None else edges
    S = pl.shapes(d, col, e); pl.draw(a, r, e, S, nice(col), "fraction / unit")
    if logy: a.set_yscale("log")
    m = dv.metrics({k: dv.hist(d[k], col, e) for k in d})
    a.set_title(nice(col), loc="left", fontsize=9.6)
    a.text(0.98, 0.95, f"χ²/ndf  CV {m['chi_cv']/m['ndf']:.0f}   DENT {m['chi_dent']/m['ndf']:.0f}", transform=a.transAxes, ha="right", va="top", fontsize=8, color=pl.INK2)
    rr = np.concatenate([S["dent"]["p"] / np.where(S["cv"]["p"] > 0, S["cv"]["p"], np.nan)]); lim = np.nanpercentile(np.abs(rr - 1), 95) * 1.6 + 0.02
    lim = min(max(lim, 0.06), 0.6); r.set_ylim(1 - lim, 1 + lim)
    return a, r

# ---------- Fig 3: calorimetry + PID features ----------
cols3 = [("reco_leading_muon_calo_csda_diff", False), ("reco_leading_proton_calo_csda_diff", False), ("reco_leading_muon_primary_softmax", True),
         ("reco_leading_muon_proton_softmax", True), ("reco_leading_muon_muon_softmax", True), ("reco_leading_muon_pion_softmax", True)]
fig = plt.figure(figsize=(14, 8.6)); gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.28, left=0.06, right=0.99, top=0.91, bottom=0.07)
for i, (c, lg) in enumerate(cols3):
    a, r = panel(fig, gs[i // 3, i % 3], c, lg)
    if i == 0: a.legend(loc="upper left", fontsize=8.5)
fig.suptitle("Calorimetry and PID features: data vs CV vs DENT   (markers = data ± stat; ratio panel is each curve over CV)", x=0.005, ha="left", fontsize=11)
fig.savefig("fig3_calorimetry_pid.png", dpi=150, facecolor=pl.SURF)

# ---------- Fig 4: next variables by separation ----------
shown = {c for c, _ in cols3} | {"abs_end_x", "muon_xing_geom", "reco_cathode_crosser", "reco_neutrino2026_veto", "reco_leading_muon_end_x", "reco_vertex_x"}
rest = [v for v in scan["var"] if v not in shown and not v.endswith("secondary_softmax") and "hadron" not in v][:12]
fig = plt.figure(figsize=(16, 11)); gs = fig.add_gridspec(3, 4, hspace=0.38, wspace=0.3, left=0.05, right=0.99, top=0.93, bottom=0.05)
for i, c in enumerate(rest):
    panel(fig, gs[i // 4, i % 4], c, logy=c.endswith("softmax"))
fig.suptitle("Next-ranked variables by CV–DENT separation   (markers = data ± stat; ratio panel is each curve over CV)", x=0.005, ha="left", fontsize=11)
fig.savefig("fig4_other_variables.png", dpi=140, facecolor=pl.SURF)
print("fig4 variables:", rest)
