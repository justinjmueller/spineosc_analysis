"""
SBND data access request, phases 4 and 5: conditional-constraint overlays and
the global null fit. Companion of ``data_access.py`` (which dispatches here
for ``--phase 4`` and ``--phase 5``).

Phase 4: conditional constraint
-------------------------------
The observed E_nu distribution is used to condition each partner variable,
and the conditioned prediction is drawn over the data with its shrunken
uncertainty band. Inputs are the ``p4_enu_*_indep_joint`` runs, where the
11 hist1d detector systematics are split into independent per-observable
copies. The earlier shared-axis form gave the constrained half an
energy-binned detector shape, which inflated its chi-square by 11 (delta p)
and 16 (p_mu) units and produced a spurious tension. The 90 per-event
reweights are unaffected either way and carry most of the cross-block
correlation.

The constraint uses the whole partner block, both channels together, since
that is the information actually available. Each panel shows one channel,
and the chi-square quoted on it is restricted to the bins drawn, so the
number matches the picture. In a panel the red-outlined band is the
conditioned uncertainty, the grey hatch the unconditioned one and the
dashed line the conditioned prediction; the ratio panel is against the
unconditioned MC, so the conditioned prediction shows as a departure from
unity with its own band.

The joint covariance is symmetric, so both directions come from the SAME
runs (only the index arrays swap): ``flip=False`` is E_nu constraining the
partner (the request's Step 4); ``flip=True`` is the partner constraining
E_nu, which asks how much the interaction model pins down the spectrum the
oscillation fit uses. An internal check falls out of the reverse
direction: the unconditioned E_nu chi-square must agree across the three
joint runs (same channel, same systematics), up to the allcov_seed noise.

Caveat, sharp in the reverse direction: the statistical correlation between
the two blocks is neglected, yet both are filled by the same events.
Conditioning then treats a fluctuation already present in the constrained
block as independent information about it; the shift is over-applied and
the residual covariance over-shrunk, inflating the post-constraint
chi-square. The forward direction has the same defect but it does not bite
there. A reverse-direction p-value below 0.05 is not yet interpretable as a
tension until the 2D joint occupancy supplies the cross-block statistical
term.

Phase 5: global null fit
------------------------
Pull strengths of the fit to the null hypothesis (no oscillation), read
from ``p5_vis_energy_null_null_global_fit.txt``.

Outputs: ``constrained_<stem>_<enu|partner>_<singlep|multip>.<fmt>`` and
``pulls_p5_vis_energy_null.<fmt>`` in ``--output``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_access import CHANNEL_TOKENS, RUNS_DIR  # noqa: E402

import common  # noqa: E402

#: Selection value of the "total" subchannel in the histogram trees.
TOTAL_SUBCHANNEL = 4294967295

#: stem, number of channels in the joint config, and which are the
#: CONSTRAINED variable (E_nu always occupies the first two).
PHASE4_PAIRS: dict[str, dict] = {
    r"$\delta p$": dict(
        stem="p4_enu_dp_indep_joint",
        nch=4,
        partner=(2, 3),
        xlabel=r"Reconstructed $\delta p$ [MeV/c]",
    ),
    r"$p_{\mu}$": dict(
        stem="p4_enu_pmu_indep_joint",
        nch=4,
        partner=(2, 3),
        xlabel=r"Reconstructed $p_{\mu}$ [MeV/c]",
    ),
    r"$E_{ccqe}$": dict(
        stem="p4_enu_eccqe_indep_joint",
        nch=3,
        partner=(2,),
        xlabel=r"Reconstructed $E_{ccqe}$ [GeV]",
    ),
}

#: x label of the E_nu block when it is the constrained one.
ENU_XLABEL = "Reconstructed Visible Energy [GeV]"

#: Subchannel index -> (label, colour). The colour list is in XML
#: DECLARATION order (intime, cos, nc, nuecc, othercc, target), so target
#: takes the LAST colour; pairing the other way inverts the stack relative
#: to every other panel.
SUBCHANNEL_STYLE = {
    0: ("In-time Cosmics", "C0"),
    1: ("Cosmics", "C1"),
    2: (r"$\nu$ NC", "C2"),
    3: (r"$\nu_e$ CC", "C4"),
    4: (r"Other $\nu$ CC", "C3"),
    5: ("Target", "C5"),
}

#: The colour ``mark_axis`` uses (legend title).
ORANGE = "#d67a11"

#: Channel names within the constrained block, and in the direction table.
CHANNEL_NAME = {0: "Single-Proton", 1: "Multi-Proton"}
TABLE_CHANNEL = {0: "single-p", 1: "multi-p", 2: "single-p", 3: "multi-p"}

#: Global-fit pull file of phase 5 and the pull histogram range.
GLOBAL_FIT_FILE = RUNS_DIR / "p5_vis_energy_null_null_global_fit.txt"
PULL_RANGE = (-2.0, 2.0)
PULL_BINS = 20


def _bins(hist: pd.DataFrame, ch: int, prefix: str, sub: int):
    """
    Bin contents (times width) and edges of one SBND histogram.

    Parameters
    ----------
    hist : pandas.DataFrame
        The ``hist1d`` tree.
    ch : int
        Channel index.
    prefix : str
        ``"CV"`` or ``"DATA"``.
    sub : int
        Subchannel, `TOTAL_SUBCHANNEL` for the sum.

    Returns
    -------
    content, low, high : np.ndarray
        Content multiplied back by the bin width, and the bin edges.
    """
    b = hist[
        (hist["variable"] == 0)
        & (hist["channel"] == ch)
        & (hist["detector"] == 0)
        & (hist["subchannel"] == sub)
        & (hist["prefix"] == prefix)
    ]
    b = b.drop_duplicates("bin_index").sort_values("bin_index")
    width = (b["bin_high_edge"] - b["bin_low_edge"]).to_numpy()
    return (
        b["bin_content"].to_numpy() * width,
        b["bin_low_edge"].to_numpy(),
        b["bin_high_edge"].to_numpy(),
    )


def phase4_constraint(spec: dict, flip: bool = False, alpha: float = 1.0) -> dict:
    """
    Condition one block of a joint run on the other.

    Parameters
    ----------
    spec : dict
        Entry of `PHASE4_PAIRS`.
    flip : bool, optional
        False: E_nu conditions the partner; True: the reverse.
    alpha : float, optional
        Near-detector exposure scale; 1.0 is nominal.

    Returns
    -------
    dict
        ``hist``, ``N``, ``D``, ``C``, ``A`` (conditioning indices), ``B``
        (conditioned indices), ``nA``, ``mu_pre``, ``mu_post``,
        ``frac_pre``, ``frac_post`` and ``reduction``.
    """
    import uproot
    from mutools.profit import conditioned

    hist = pd.DataFrame(
        uproot.open(RUNS_DIR / f"{spec['stem']}_PROplot_histograms.root")["hist1d"].arrays(
            library="np"
        )
    )
    cv = [_bins(hist, c, "CV", TOTAL_SUBCHANNEL)[0] for c in range(spec["nch"])]
    da = [_bins(hist, c, "DATA", TOTAL_SUBCHANNEL)[0] for c in range(spec["nch"])]
    n_all, d_all = np.concatenate(cv), np.concatenate(da)
    cov = uproot.open(RUNS_DIR / f"{spec['stem']}_PROplot.root")[
        "Covariance/collapsed_total_frac_cov"
    ].values()
    n_a = len(cv[0]) + len(cv[1])  # E_nu occupies the first two channels
    enu, partner = np.arange(n_a), np.arange(n_a, len(n_all))
    nd, fd = (partner, enu) if flip else (enu, partner)  # nd conditions, fd is conditioned
    stat = np.diag(1.0 / np.maximum(n_all[nd], 1e-30))
    m = cov[np.ix_(nd, nd)] + stat / alpha
    shift = cov[np.ix_(fd, nd)] @ np.linalg.solve(m, (d_all[nd] - n_all[nd]) / n_all[nd])
    resid = cov[np.ix_(fd, fd)] - cov[np.ix_(fd, nd)] @ np.linalg.solve(m, cov[np.ix_(nd, fd)])
    return dict(
        hist=hist,
        N=n_all,
        D=d_all,
        C=cov,
        A=nd,
        B=fd,
        nA=n_a,
        mu_pre=n_all[fd],
        mu_post=n_all[fd] * (1 + shift),
        frac_pre=cov[np.ix_(fd, fd)],
        frac_post=resid,
        reduction=conditioned(cov, nd=nd, fd=fd, nd_stats=stat, alpha=alpha).reduction,
    )


def phase4_chi2(r: dict, mu: np.ndarray, frac: np.ndarray, idx: np.ndarray):
    """
    chi-square over a subset of the conditioned block.

    Restricted to the bins a panel actually draws, with a Poisson
    (data) statistical term added to the fractional covariance.

    Parameters
    ----------
    r : dict
        Output of `phase4_constraint`.
    mu : np.ndarray
        Prediction of the conditioned block.
    frac : np.ndarray
        Fractional covariance of the conditioned block.
    idx : np.ndarray
        Indices within the conditioned block.

    Returns
    -------
    value, ndf, pvalue : float, int, float
    """
    from scipy.stats import chi2 as chi2_dist

    d = r["D"][r["B"]][idx]
    cov = frac[np.ix_(idx, idx)] * np.outer(mu[idx], mu[idx]) + np.diag(np.maximum(d, 1e-30))
    res = d - mu[idx]
    value = float(res @ np.linalg.solve(cov, res))
    return value, len(idx), float(chi2_dist.sf(value, len(idx)))


def _pad(y: np.ndarray) -> np.ndarray:
    """
    Repeat the first and last entries, for ``step='mid'`` drawing on edge-padded x.

    Parameters
    ----------
    y : np.ndarray

    Returns
    -------
    np.ndarray
    """
    return np.r_[y[0], y, y[-1]]


def plot_constrained(r: dict, spec: dict, label: str, ch: int, drawn: tuple[int, ...],
                     offset: int, flip: bool):
    """
    One constrained-overlay panel (stack, bands, data, ratio) for a channel.

    Parameters
    ----------
    r : dict
        Output of `phase4_constraint`.
    spec : dict
        Entry of `PHASE4_PAIRS`.
    label : str
        Pair label (the constrainer when `flip`).
    ch : int
        Absolute channel index of the drawn block.
    drawn : tuple of int
        Channels of the constrained block.
    offset : int
        Index of this channel's first bin within the constrained block.
    flip : bool
        Direction, as in `phase4_constraint`.

    Returns
    -------
    fig : matplotlib.figure.Figure
    summary : tuple
        ``(label, channel, chi2, ndf, p, mean band pre, mean band post)``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from mutools.plotting.helpers import mark_axis
    from mutools.plotting.profit import add_error_band, construct_meta_handle, saver

    xlabel = ENU_XLABEL if flip else spec["xlabel"]
    by = label if flip else r"$E_{\nu}^{reco}$"
    content, lo, hi = _bins(r["hist"], ch, "CV", TOTAL_SUBCHANNEL)
    idx = np.arange(offset, offset + len(content))
    ctr, wid, edges = 0.5 * (lo + hi), hi - lo, np.r_[lo, hi[-1]]
    xs = np.r_[edges[0], ctr, edges[-1]]
    mu_pre, mu_post = r["mu_pre"][idx], r["mu_post"][idx]
    sig_pre = np.sqrt(np.diag(r["frac_pre"])[idx]) * mu_pre
    sig_post = np.sqrt(np.clip(np.diag(r["frac_post"])[idx], 0, None)) * mu_post
    dat = r["D"][r["B"]][idx]
    derr = np.sqrt(np.maximum(dat, 0))
    chi, ndf, pval = phase4_chi2(r, r["mu_post"], r["frac_post"], idx)

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    rx = fig.add_subplot(gs[1, 0], sharex=ax)

    # Stack in the declaration order mutools uses, with counter percentages.
    sub_vals, sub_lab, sub_col = [], [], []
    for s in (5, 4, 3, 2, 1, 0):
        v = _bins(r["hist"], ch, "CV", s)[0]
        if v.sum() <= 0:
            continue
        sub_vals.append(v)
        sub_lab.append(SUBCHANNEL_STYLE[s][0])
        sub_col.append(SUBCHANNEL_STYLE[s][1])
    total = sum(v.sum() for v in sub_vals)
    labels = [f"{n} ({v.sum() / total:.1%})" for n, v in zip(sub_lab, sub_vals, strict=True)]
    ax.hist(
        [ctr] * len(sub_vals),
        bins=edges,
        weights=sub_vals,
        stacked=True,
        color=sub_col,
        edgecolor="none",
        alpha=0.7,
        rasterized=saver.rasterized,
    )
    # Patches are built directly rather than via construct_proxy_stack, which
    # assigns colours in REVERSED subchannel order: correct for mutools' own
    # call site but it mispairs labels when used standalone.
    proxy = [Patch(facecolor=c, alpha=0.7, label=n) for n, c in zip(labels, sub_col, strict=True)]

    # Unconstrained prior band, then the constrained prediction over it. The
    # latter is a filled area with a red outline: the prior is already grey
    # hatch, and a second grey hatch on top is unreadable.
    prior_patch = add_error_band(
        ax, edges, mu_pre, yerr=[sig_pre, sig_pre], label="Unconstrained Uncertainty",
        color="0.55", hatch="\\\\\\",
    )
    ax.fill_between(xs, _pad(mu_post - sig_post), _pad(mu_post + sig_post), step="mid",
                    facecolor="red", alpha=0.12, zorder=3)
    for sgn in (-1, 1):
        ax.step(xs, _pad(mu_post + sgn * sig_post), where="mid", color="red", lw=1.1, zorder=3)
    post_patch = Patch(
        facecolor="red", alpha=0.12, edgecolor="red", label="Constrained Uncertainty"
    )
    ax.step(xs, _pad(mu_post), where="mid", color="k", ls="--", lw=1.5, zorder=4)
    post_line = Line2D([], [], color="k", ls="--", lw=1.5, label="Constrained MC")
    data_handle = ax.errorbar(ctr, dat, yerr=derr, xerr=0.5 * wid, fmt="ko", markersize=4,
                              linewidth=1.0, capsize=2, zorder=5, label="Data")
    chi_patch = Patch(facecolor="none", edgecolor="none", label=rf"$\chi^2$/ndf = {chi:.1f}/{ndf}")
    p_patch = Patch(facecolor="none", edgecolor="none", label=rf"$p$-value = {pval:.3f}")
    meta = construct_meta_handle(common.CODE_VERSION, common.SELECTION_VERSION)
    leg = ax.legend(
        handles=[*proxy, prior_patch, post_patch, post_line, data_handle, chi_patch, p_patch, meta],
        ncol=2,
        loc="upper right",
        fontsize=8,
        title=f"{CHANNEL_NAME.get(ch - drawn[0], '')}\nconstrained by {by}",
        title_fontsize=11,
    )
    leg.get_title().set_color(ORANGE)
    leg.get_title().set_fontweight("bold")
    ax.set_ylabel("Events / bin")
    ax.set_ylim(0, 1.75 * float(np.max(dat + derr)))
    ax.set_xlim(edges[0], edges[-1])
    ax.tick_params(labelbottom=False)
    mark_axis(ax, r"$\bf{SBND}$ Analysis-in-Progress", vadj=0.00)

    # Ratio: everything relative to the UNCONSTRAINED MC, so the constrained
    # prediction itself is visible as a departure from unity.
    add_error_band(rx, edges, np.ones_like(mu_pre), yerr=[sig_pre / mu_pre, sig_pre / mu_pre],
                   label="", color="0.55", hatch="\\\\\\")
    rlo, rhi = (mu_post - sig_post) / mu_pre, (mu_post + sig_post) / mu_pre
    rx.fill_between(xs, _pad(rlo), _pad(rhi), step="mid", facecolor="red", alpha=0.12, zorder=3)
    for y in (rlo, rhi):
        rx.step(xs, _pad(y), where="mid", color="red", lw=1.1, zorder=3)
    rx.step(xs, _pad(mu_post / mu_pre), where="mid", color="k", ls="--", lw=1.5, zorder=4)
    rx.errorbar(ctr, dat / mu_pre, yerr=derr / mu_pre, xerr=0.5 * wid, fmt="ko", markersize=4,
                linewidth=1.0, capsize=2, zorder=5)
    rx.axhline(1.0, color="red", linestyle="--", linewidth=1.0, zorder=1)
    rx.set_ylim(0.0, 2.0)
    rx.set_ylabel("Ratio to MC")
    rx.set_xlabel(xlabel)
    summary = (label, ch - drawn[0], chi, ndf, pval,
               float(np.mean(sig_pre / mu_pre)), float(np.mean(sig_post / mu_post)))
    return fig, summary


def direction_table(alpha: float = 1.0) -> list[tuple]:
    """
    Pre/post-constraint chi-square in both directions, from the same runs.

    Parameters
    ----------
    alpha : float, optional
        Near-detector exposure scale.

    Returns
    -------
    list of tuple
        ``(flip, pair, channel, pre, post, mean|shift|, band_pre, band_post)``
        where pre/post are ``(chi2, ndf, p)``.
    """
    rows = []
    for flip in (False, True):
        for lab, spec in PHASE4_PAIRS.items():
            if not (RUNS_DIR / f"{spec['stem']}_PROplot_histograms.root").exists():
                continue
            r = phase4_constraint(spec, flip, alpha)
            drawn, off = ((0, 1) if flip else spec["partner"]), 0
            for ch in drawn:
                n = len(_bins(r["hist"], ch, "CV", TOTAL_SUBCHANNEL)[0])
                idx = np.arange(off, off + n)
                off += n
                pre = phase4_chi2(r, r["mu_pre"], r["frac_pre"], idx)
                post = phase4_chi2(r, r["mu_post"], r["frac_post"], idx)
                shift = np.abs(r["mu_post"][idx] / r["mu_pre"][idx] - 1).mean()
                band_pre = float(np.mean(np.sqrt(np.clip(np.diag(r["frac_pre"])[idx], 0, None))))
                band_post = float(np.mean(np.sqrt(np.clip(np.diag(r["frac_post"])[idx], 0, None))))
                rows.append((flip, lab, TABLE_CHANNEL[ch], pre, post, shift, band_pre, band_post))
    return rows


def print_direction_table(rows: list[tuple]) -> None:
    """
    Print the both-directions table and the E_nu consistency check.

    Parameters
    ----------
    rows : list of tuple
        Output of `direction_table`.
    """
    for flip in (False, True):
        print("=" * 84)
        print("partner -> E_nu  (reverse)" if flip else "E_nu -> partner  (proposal Step 4)")
        print(f"{'pair':10s}{'channel':>10}{'chi2 pre':>13}{'p pre':>8}"
              f"{'chi2 post':>13}{'p post':>8}{'band pre':>10}{'band post':>10}{'|shift|':>9}")
        for f, lab, ch, pre, post, sh, bpre, bpost in rows:
            if f != flip:
                continue
            print(f"{lab:10s}{ch:>10}{pre[0]:9.2f}/{pre[1]:<3d}{pre[2]:8.3f}"
                  f"{post[0]:9.2f}/{post[1]:<3d}{post[2]:8.3f}{bpre:10.3f}{bpost:10.3f}{sh:9.3f}")
    # The reverse direction reuses one E_nu channel three times, so its
    # unconditioned chi2 must agree across pairs up to the allcov_seed noise.
    for ch in ("single-p", "multi-p"):
        v = [r[3][0] for r in rows if r[0] and r[2] == ch]
        if len(v) > 1:
            print(f"E_nu {ch:8s} unconditioned chi2 across pairs: "
                  f"{', '.join(f'{x:.2f}' for x in v)}  (spread {max(v) - min(v):.2f})")


def run_phase4(output: Path, fmt: str, flip: bool = False, alpha: float = 1.0) -> None:
    """
    Make the phase-4 panels and print both summary tables.

    Parameters
    ----------
    output : Path
        Figure directory.
    fmt : str
        Figure format.
    flip : bool, optional
        Direction of the panels (the table always reports both).
    alpha : float, optional
        Near-detector exposure scale; 1.0 is nominal.
    """
    direction = "partner -> E_nu" if flip else "E_nu -> partner"
    print(f"phase 4: {len(PHASE4_PAIRS)} pairs, alpha = {alpha}, direction {direction}")
    made = []
    for label, spec in PHASE4_PAIRS.items():
        stem = spec["stem"]
        if not (RUNS_DIR / f"{stem}_PROplot_histograms.root").exists():
            print(f"{label}: MISSING {stem}")
            continue
        r = phase4_constraint(spec, flip, alpha)
        drawn = (0, 1) if flip else spec["partner"]  # channels of the CONSTRAINED block
        offset = 0
        for ch in drawn:
            n = len(_bins(r["hist"], ch, "CV", TOTAL_SUBCHANNEL)[0])
            fig, summary = plot_constrained(r, spec, label, ch, drawn, offset, flip)
            offset += n
            token = CHANNEL_TOKENS.get(ch - drawn[0], f"ch{ch}")
            side = "enu" if flip else "partner"
            common.save_figure(fig, output / f"constrained_{stem}_{side}_{token}.{fmt}")
            made.append(summary)
    print(f"\n{len(made)} panels")
    print(f"{'pair':14s}{'channel':>14}{'chi2/ndf':>12}{'p':>8}{'band pre':>10}{'band post':>11}")
    for lab, ch, chi, ndf, pv, bpre, bpost in made:
        print(f"{lab:14s}{CHANNEL_NAME.get(ch, ch):>14}{chi:8.1f}/{ndf:<3d}{pv:8.3f}"
              f"{bpre:10.3f}{bpost:11.3f}")
    print_direction_table(direction_table(alpha))


def read_global_fit(path: Path = GLOBAL_FIT_FILE) -> pd.DataFrame:
    """
    Parse the ``name: value`` lines of a PROfit global-fit text file.

    Parameters
    ----------
    path : Path, optional
        The fit result; the first line (header) is skipped.

    Returns
    -------
    pandas.DataFrame
        Columns ``parameter`` and ``value``.
    """
    results: dict[str, list] = {"parameter": [], "value": []}
    for line in path.read_text().splitlines()[1:]:
        if ":" not in line:
            continue
        name, value = line.split(":", 1)
        results["parameter"].append(name.strip())
        results["value"].append(float(value))
    return pd.DataFrame(results)


def run_phase5(output: Path, fmt: str) -> None:
    """
    Histogram the pulls of the global null fit.

    Parameters
    ----------
    output : Path
        Figure directory.
    fmt : str
        Figure format.
    """
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis

    fits = read_global_fit()
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(fits["value"], bins=PULL_BINS, color="C0", edgecolor="k", alpha=0.7, range=PULL_RANGE)
    ax.set_xlim(*PULL_RANGE)
    ax.set_xlabel("Pull Strength")
    ax.set_ylabel("Count")
    mark_axis(ax, r"$\bf{SBND}$ Analysis-in-Progress")
    print(common.save_figure(fig, output / f"pulls_p5_vis_energy_null.{fmt}"))
