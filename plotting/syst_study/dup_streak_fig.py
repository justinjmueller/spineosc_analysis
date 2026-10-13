import numpy as np, pandas as pd, uproot, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); D = pathlib.Path("/Users/mueller/data/spineosc/v1.1.0"); S = pathlib.Path("/Users/mueller/Projects/spineosc/docs/syst_study")
cases = [("ICARUS Run 4 dirt", "spineosc_numu_icarus_run4_v1.1.0_dirt_nosyst.root", "events/icarus_dirt_run4/events", "#0072B2"), ("ICARUS Run 4 WireModxThetaXW", "spineosc_numu_icarus_run4_v1.1.0_detvar_nosyst.root", "events/icarus_wiremod_xthetaxw_run4/events", "#D55E00")]
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2)); stats = {}
for lab, fn, tree, c in cases:
    df = uproot.open(D / fn)[tree].arrays(library="pd"); n = len(df); h = pd.util.hash_pandas_object(df, index=False).to_numpy(); first = {}; src = np.full(n, -1)
    for j, x in enumerate(h):
        if x in first: src[j] = first[x]
        else: first[x] = j
    rep = np.where(src >= 0)[0]; off = rep - src[rep]; cuts = np.r_[0, np.where((np.diff(rep) != 1) | (np.diff(off) != 0))[0] + 1, len(rep)]
    L = np.diff(cuts); j0 = rep[cuts[:-1]]; i0 = src[j0]; stats[lab] = (n, len(rep), len(L))
    axs[0].hist(L, bins=np.arange(0.5, 400, 5), histtype="step", lw=2, color=c, label=f"{lab}: {len(L):,} streaks")
    axs[1].plot(np.arange(20) * 5 + 2.5, np.histogram(rep, bins=20, range=(0, n))[0] / (n / 20) * 100, "-o", color=c, label=lab)
    axs[2].scatter(i0 / n * 100, j0 / n * 100, s=6, color=c, alpha=0.5, label=lab)
axs[0].set_yscale("log"); axs[0].set_xlabel("Streak length [consecutive repeated rows]"); axs[0].set_ylabel("Streaks"); axs[0].axvline(50, color="0.5", ls=":", lw=1); axs[0].text(52, 3, "50 = events per file", fontsize=10, color="0.35"); axs[0].legend(fontsize=9, loc="upper right")
axs[1].set_xlabel("Position in the tree [%]"); axs[1].set_ylabel("Repeated rows [% of rows in the slice]"); axs[1].legend(fontsize=9)
axs[2].plot([0, 100], [0, 100], color="0.6", lw=1); axs[2].set_xlabel("Position of the original block [% of tree]"); axs[2].set_ylabel("Position of the copy [% of tree]"); axs[2].legend(fontsize=9, loc="upper left")
mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[2], "events trees, storage order", alignment="right"); fig.tight_layout(); fig.savefig(S / "dup_streaks.png", dpi=95); print(stats)
