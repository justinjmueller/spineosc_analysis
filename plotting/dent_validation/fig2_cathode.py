import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import dentval as dv, plotlib as pl
d = dv.load()
fig = plt.figure(figsize=(14, 8.6), constrained_layout=False)
gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.28, left=0.06, right=0.99, top=0.91, bottom=0.07)
# NOTE: the length panels were removed (muon length is not an allowed data-access variable); panels (a) and (b) are intentionally empty.
# (c) |end_x| zoom, non-crossers, fractions of all selected events per cm
a, r = pl.pair_axes(fig, gs[0, 2])
e = np.arange(0, 42, 2.0)
S = pl.shapes(d, "abs_end_x", e, sel=lambda x: x.reco_cathode_crosser.to_numpy() == 0, norm_to="all")
pl.draw(a, r, e, S, "|x$_{end}$| of the muon, non-crossers [cm]", "fraction of all events / cm"); a.legend(loc="upper right")
a.set_title("(c) distance of the endpoint from the cathode", loc="left", fontsize=10); r.set_ylim(0.5, 2.3)
# (d) signed end_x, xing (muon vertex/end opposite sides)
a, r = pl.pair_axes(fig, gs[1, 0]); e2 = np.arange(-60, 62, 4.0)
S = pl.shapes(d, "reco_leading_muon_end_x", e2, sel=lambda x: x.muon_xing_geom.to_numpy() == 1, norm_to="all")
pl.draw(a, r, e2, S, "muon endpoint x [cm]", "fraction of all events / cm"); a.set_title("(d) crossing muons (vertex and end on opposite sides)", loc="left", fontsize=10); r.set_ylim(0.6, 1.4)
# (e) signed end_x, stay
a, r = pl.pair_axes(fig, gs[1, 1])
S = pl.shapes(d, "reco_leading_muon_end_x", e2, sel=lambda x: x.muon_xing_geom.to_numpy() == 0, norm_to="all")
pl.draw(a, r, e2, S, "muon endpoint x [cm]", "fraction of all events / cm"); a.set_title("(e) non-crossing muons", loc="left", fontsize=10); r.set_ylim(0.6, 1.4)
# (f) vertex x zoom, all events
a, r = pl.pair_axes(fig, gs[1, 2])
S = pl.shapes(d, "reco_vertex_x", e2, norm_to="all")
pl.draw(a, r, e2, S, "reconstructed vertex x [cm]", "fraction of all events / cm"); a.set_title("(f) interaction vertex near the cathode", loc="left", fontsize=10); r.set_ylim(0.6, 1.4)
fig.suptitle("SBND cathode-region features: data vs CV vs DENT   (selection: ≥1 proton, cathode + high-y/z cuts, flash threshold; markers = data ± stat; ratio panel is each curve over CV)", x=0.005, ha="left", fontsize=11, color=pl.INK)
fig.savefig("fig2_cathode_features.png", dpi=150, facecolor=pl.SURF)
# numbers for the record
rows = []
def add(name, k, v, e): rows.append(dict(feature=name, sample=k, value=v, err=e))
for k in d:
    df = d[k]; n = len(df); cx = (df.reco_cathode_crosser == 1).mean(); pu = ((df.reco_cathode_crosser == 0) & (df.abs_end_x < 5)).mean()
    add("crosser fraction [%]", k, 100 * cx, 100 * np.sqrt(cx * (1 - cx) / n)); add("pile-up fraction (non-crosser, |x_end|<5cm) [%]", k, 100 * pu, 100 * np.sqrt(pu * (1 - pu) / n))
pd.DataFrame(rows).to_csv("cathode_feature_numbers.csv", index=False)
print(pd.DataFrame(rows).pivot(index="feature", columns="sample", values="value").round(3))
