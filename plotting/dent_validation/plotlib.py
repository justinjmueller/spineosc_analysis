import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dentval as dv

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C_CV, C_DENT = "#2a78d6", "#eb6834"       # validated categorical slots 1 and 2; data is drawn in ink
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False})

def pair_axes(fig, gs_cell):
    """main + ratio axes in one gridspec cell"""
    sub = gs_cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
    a = fig.add_subplot(sub[0]); r = fig.add_subplot(sub[1], sharex=a)
    plt.setp(a.get_xticklabels(), visible=False)
    return a, r

def shapes(d, col, edges, sel=None, per_unit=True, norm_to="own"):
    """counts and unit-area densities for data/cv/dent. sel: function df->mask. norm_to 'own' = each sample's selected total;
    norm_to 'all' = divide by the sample's full (pre-sel) event count so fractions are of ALL selected events."""
    out = {}
    for k, df in d.items():
        full = len(df)
        x = df if sel is None else df[sel(df)]
        v = x[col].to_numpy()
        c = np.histogram(v, bins=edges)[0].astype(float)
        n = c.sum() if norm_to == "own" else full
        w = np.diff(edges) if per_unit else 1.0
        out[k] = dict(c=c, n=n, p=c / n / w, e=np.sqrt(np.maximum(c, 1)) / n / w)
    return out

def draw(a, r, edges, S, xlabel=None, ylabel="fraction / unit", label_unit=""):
    ctr = 0.5 * (edges[1:] + edges[:-1])
    for k, col in (("cv", C_CV), ("dent", C_DENT)):
        a.stairs(S[k]["p"], edges, color=col, lw=1.6, label={"cv": "MC CV", "dent": "MC DENT"}[k])
    a.errorbar(ctr, S["data"]["p"], yerr=S["data"]["e"], fmt="o", ms=3.6, color=INK, ecolor=INK, elinewidth=1, label="data", zorder=5)
    a.set_ylabel(ylabel)
    pc = np.where(S["cv"]["p"] > 0, S["cv"]["p"], np.nan)
    r.stairs(S["dent"]["p"] / pc, edges, color=C_DENT, lw=1.6)
    r.axhline(1, color=C_CV, lw=1.6)
    r.errorbar(ctr, S["data"]["p"] / pc, yerr=S["data"]["e"] / pc, fmt="o", ms=3.6, color=INK, ecolor=INK, elinewidth=1, zorder=5)
    r.set_ylabel("÷ CV")
    if xlabel: r.set_xlabel(xlabel)

def frac_vs(d, xcol, edges, num, den=None):
    """fraction num(df)/den(df) in bins of xcol; binomial errors. returns dict per sample."""
    out = {}
    for k, df in d.items():
        v = df[xcol].to_numpy(); m_den = np.ones(len(df), bool) if den is None else den(df)
        m_num = num(df) & m_den
        n = np.histogram(v[m_num], bins=edges)[0].astype(float); t = np.histogram(v[m_den], bins=edges)[0].astype(float)
        p = np.where(t > 0, n / t, np.nan); e = np.sqrt(np.where(t > 0, p * (1 - p) / np.maximum(t, 1), np.nan))
        out[k] = dict(p=p, e=e, n=n, t=t)
    return out

def draw_frac(a, r, edges, F, xlabel, ylabel, pct=True):
    ctr = 0.5 * (edges[1:] + edges[:-1]); s = 100 if pct else 1
    for k, col in (("cv", C_CV), ("dent", C_DENT)):
        a.plot(ctr, s * F[k]["p"], marker="s" if k == "cv" else "D", ms=4.5, lw=1.6, color=col, mec=SURF, mew=0.8, label={"cv": "MC CV", "dent": "MC DENT"}[k])
    a.errorbar(ctr, s * F["data"]["p"], yerr=s * F["data"]["e"], fmt="o", ms=4.2, color=INK, ecolor=INK, elinewidth=1, zorder=5, label="data")
    a.set_ylabel(ylabel)
    base = F["cv"]["p"]
    r.axhline(1, color=C_CV, lw=1.6); r.plot(ctr, F["dent"]["p"] / base, marker="D", ms=4.5, lw=1.6, color=C_DENT, mec=SURF, mew=0.8)
    r.errorbar(ctr, F["data"]["p"] / base, yerr=F["data"]["e"] / base, fmt="o", ms=4.2, color=INK, ecolor=INK, elinewidth=1, zorder=5)
    r.set_ylabel("÷ CV"); r.set_xlabel(xlabel)
