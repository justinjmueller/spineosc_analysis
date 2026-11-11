"""
Does the SBND DENT simulation improve data/MC agreement? Tables and figures.

A shape-only, no-fit study with the SBND release detector-access data.
This script replaces the notebook ``plotting/dent_study.ipynb`` (itself
generated from ``plotting/dent_study/report_cells.py``). It reads the
per-study result files that ``scripts/dent_study/run_study.py`` wrote to
``runs/dent_study`` (plus ``docs/dent_study/replica_summary.csv``),
prints the report tables and writes the figures. The analysis helpers
live in ``plotting/dent_study/`` (``dentstudy.py``, ``replica.py``) and
are imported, not copied.

The complete narrative, section by section, is kept in
``plotting/dent_study/text/*.md`` (00 summary, 01 data and rules, 02 QC,
03 mechanism, 04-06 Tier A and B, 07-08 fixed-covariance replica, 09
Tier C, 10 baseline, 11 limitations, 12 deviations) and the rendered
report is ``docs/dent_study/report.html``; the pre-registered plan is
``docs/dent_study/PREREGISTRATION.md``. The summary, limitations and
deviations are reproduced here so they travel with the code.

Question
--------
DENT is the detector variation that models the cathode-crossing (and
related) effects in SBND. After the cathode-crosser cut it is one of the
larger uncorrelated systematics. Does simulating it (the DENT detvar
sample, or the production CV shifted by the DENT ratio) describe the
data better than the CV does, and where does it match or fail to match?

Method
------
SBND release detector-access data (the file has no energy or momentum
column, and none is used) are compared with the SBND simulation, always
area-normalised and shape-only, in each observable's own binning, with
shape chi-square and p-value as the metrics. Nothing is fitted: DENT
enters only at fixed strengths 0, 0.5 and 1.0.

* Tier A uses the detvar CV and DENT samples directly (MC statistics
  only).
* Tier B uses the production CV (``cvrw`` plus ``cvwgt``) with the full
  systematic covariance.
* Tier C asks how well one observable could tell the two apart, using
  fake data.

Findings (as of 2026-10-07)
---------------------------
1. The data do not show that DENT applied to the central value improves
   the agreement, and they do not reject it. With the production
   covariance, shifting the CV by the full DENT effect changes the shape
   chi-square by less than 2 units for 9 of the 10 primary observables
   (5 improve, 5 worsen). The tenth, ``muon_calo_csda``, is dominated by
   a covariance re-centering artefact (see the replica stage). Across the
   111 SBND regenerated detector-access channels the largest changes are
   +5.8 and -3.1 units.
2. Tier A (MC statistics only) splits by physics. The DENT sample is
   closer to the data for the seven calorimetry and particle-ID
   observables (largest for ``muon_calo_csda``, 10,534 to 2,188 on 19
   degrees of freedom), and farther for the three cathode-geometry
   observables (``crosser_flag``, ``muon_abs_endx``, ``muon_endx_full``;
   for ``muon_abs_endx`` 51 to 338 on 19 degrees of freedom). The
   calorimetry/PID improvement disappears as the other systematics are
   added to the covariance (ablation stage), so it is not evidence that
   DENT is the cause of those data/MC differences.
3. The cathode features that make DENT large are only partly confirmed
   by the data. The fraction of events flagged as crossers is
   27.46 +- 0.18 % in data against 27.92 % (CV) and 25.36 % (DENT).
   Non-crossing muons ending within 5 cm of the cathode (the "pile-up")
   are 1.95 % of events in data against 2.07 % (CV) and 3.78 % (DENT):
   the data follow the CV, and the large DENT pile-up is not present.
   The deficit of crossing muons that end within 4 cm of the cathode is,
   however, in the direction DENT predicts (data 0.68 %, CV 1.14 %,
   DENT 0.55 %). These cathode comparisons were seen before the study
   plan was written, so they are hypothesis-prompted, not independent
   confirmations.
4. DENT as an uncertainty is not redundant. Removing DENT from the
   covariance makes the muon end-x observables fail (``muon_abs_endx``:
   p = 0.46 with DENT, 0.003 without; ``muon_endx_full``: 0.15 against
   2e-5). It absorbs a real data/MC difference.
5. No single observable has the power to decide. With pseudo-experiments
   through the fixed covariance, the separation between the "CV" and
   "full DENT" expectations is d' = 0.5 to 1.0 per observable, and all
   ten primaries have two-sided p >= 0.15 under all three truths (CV,
   half DENT, full DENT). The observables are correlated, so their
   p-values cannot be multiplied together.

Recommendation: keep DENT as an uncertainty, and do not shift the
central value by it on the strength of this study. The mechanism that
most plausibly makes DENT large at high energy is the near-cathode
endpoint pile-up, which the data lack. A near-cathode endpoint treatment
(an endpoint cut, or a dedicated near-cathode observable) is the natural
follow-up and is not part of this study.

What this study cannot say: it uses no energy or momentum quantity, so
it cannot test the effect of DENT on the reconstructed-energy
distributions that drive the analysis. The detvar samples are
neutrino-only and their CV is not the production nominal. ICARUS has no
DENT variation; its regenerated studies are reported only as a
baseline.

Why PROfit's own change in chi-square is not enough (replica stage)
-------------------------------------------------------------------
With ``--inject-systs-cv`` PROfit re-derives the systematic covariance
around the *shifted* prediction. For ``muon_calo_csda`` the first-bin
fractional DENT covariance falls from 0.84 to 0.51 and the
statistics-only chi-square from 14,606 to 3,556, so PROfit's chi-square
at strength 1.0 mixes the central-value shift with a change of the
covariance (PROfit reports -18.8; with the covariance held at its
strength-0 form it is +0.20). The numpy replica of PROfit's shape-only
chi-square, ``chi2 = (n-p)^T [diag(n) + F o (p p^T)]^-1 (n-p)``, with
``F`` the collapsed total fractional covariance and ``p`` the prediction
area-normalised to the data, reproduces PROfit at all three strengths.
Holding ``F`` at its strength-0 value for both predictions gives the
fixed-covariance change ``chi2(0) - chi2(1)`` (positive favours the
DENT-shifted prediction); 20,000 pseudo-experiments per truth
(strength 0, 0.5 or 1.0) calibrate it. The p-value is two-sided and
empirical, with a floor of 2/(N+1). No strength is estimated.

Baseline (regenerated detector-access studies)
----------------------------------------------
SBND (40 studies, 111 channel results) and ICARUS Run 2 (28 studies, 72
channel results): cvrw plus ``cvwgt`` on the main MC, the SBND data veto
added, area-normalised and shape-only. SBND: 44 of 111 channel results
have p < 0.05 and 26 have p < 1e-3; the worst are the primary-score
observables (muon primary score, p down to 1e-38). Adding DENT at full
strength changes any channel's chi-square by between -3.1 and +5.8 units
(median +0.19), so none of these failures is a DENT effect. ICARUS Run 2:
14 of 72 have p < 0.05 and 5 have p < 1e-3; ICARUS has no DENT
variation, so these are a baseline only. The p-values include the
production covariance but no look-elsewhere adjustment, and channels
within a study share events.

Limitations and what was not done
---------------------------------
* No energy. Energy and momentum quantities are blinded and embargoed
  for SBND data, so nothing here tests DENT on the energy
  distributions, which is where its uncertainty matters. A study of the
  high-energy region would need to be run by someone cleared to look.
* Power. Single observables have d' of about 0.5 to 1. The observables
  share events and cannot be combined into one test without a joint
  covariance, which was not built.
* Detvar CV is not the production nominal. Tier A chi-squares mix DENT
  with other data/MC differences. Tier B removes this by using the
  production CV and the DENT ratio.
* DENT enters as a one-sided ratio, taken from the detvar CV/DENT pair
  and applied to the production CV. Equivalence with using the DENT
  sample directly was checked for the native PROfit detvar mode
  (unmatched-event mode agrees); a Tier B run with the native detvar
  mode was not done.
* Cathode comparisons are hypothesis-prompted (seen before the plan).
* Trigger emulation: the flash-PE cut is applied to the MC only.
* Not examined: WireMod Y-Z (SBND uses MattMod Y-Z in the systematics),
  the near-cathode pile-up with an endpoint cut applied to the data
  comparison, and any ICARUS DENT-like variation.

Deviations from the plan (dated in ``docs/dent_study/PREREGISTRATION.md``)
--------------------------------------------------------------------------
* Removed three multiplicity observables (pion, photon, electron) that
  the selection forces to zero (24 to 21 observables). None primary.
* The first exploratory data dumps wrongly applied the MC flash cut to
  data. They were re-made without it; the headline cathode numbers moved
  by less than 0.1 sigma.
* After the first Tier C result (``crosser_flag``) a third truth (half
  DENT) and two-sided p-values were added: a third fixed point, not a
  fit.
* Two crossing-split studies register DENT under per-population names,
  so injections were silently ignored. The driver now detects this; both
  were re-run and the ``qc`` stage audits for any inert injection.
* Poisson-only Tier C was replaced by pseudo-experiments, then by the
  fixed-covariance replica, after finding that PROfit's injection
  re-derives the covariance.

Stages (``--only``)
-------------------
``qc``        quality control: shape-only live, ndf = nbins - 1, no inert injections
``tierA``     detvar CV against detvar DENT (table, ``tierA_primary`` figure)
``tierB``     production CV, four variants (table, ``tierB_pvalues_primary``)
``ablation``  covariance ablation (table, ``tierB_ablation``)
``panels``    PROfit-style data/MC panels of the primaries, CV and CV + DENT
``replica``   fixed-covariance replica and pseudo-experiments (``replica_dchi2``)
``tierC``     fake-data Tier C (``tierC_dchi2``, ``tierC_hists``)
``baseline``  regenerated detector-access studies (tables)

The cathode-mechanism figures of the report (``mechanism_*.png``) are
static images made by ``scripts/dent_study`` and are not regenerated.
The figures are written as PNG into ``--output`` (default
``figures/dent_study``, not the reviewed ``docs/dent_study/figures``).

Usage::

    python plotting/studies/dent_study.py --only tierA tierB
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

#: Helper modules (``dentstudy``, ``replica``) live here and import each
#: other as top-level modules, so the directory is put on ``sys.path``.
HELPER_DIR = common.REPO / "plotting" / "dent_study"

#: Replica summary written by ``plotting/dent_study/replica.py``.
REPLICA_SUMMARY = common.DOCS / "dent_study" / "replica_summary.csv"

#: Stage names, in report order.
STAGES = ("qc", "tierA", "tierB", "ablation", "panels", "replica", "tierC", "baseline")

#: Covariance-ablation configurations: label -> (theta = 0 variant, theta = 1 variant).
#: All keep DENT in the covariance so the injection is valid; nothing is fitted.
ABLATION_CONFIGS = {
    "stat + DENT": ("abl_statDENT_t0", "abl_statDENT_t1"),
    "stat + all detector": ("abl_det_t0", "abl_det_t1"),
    "all but other detector": ("abl_nodetother_t0", "abl_nodetother_t1"),
    "full (standard)": ("v1", "inj1p0"),
}

#: Bar colours of the ablation figure, one per configuration.
ABLATION_COLORS = ("C2", "C3", "C0", "C1")

#: Panels: (variant, caption stem) for the PROfit-style panels.
PANEL_VARIANTS = (("v1", "CV"), ("inj1p0", "CV + DENT"))

#: Watermark of the summary figures.
WATERMARK_DATA = r"$\bf{SBND}$ Data Access"


def load_helpers():
    """
    Import the study helper modules from ``plotting/dent_study``.

    The helpers are plain scripts (not a package) that refer to each
    other by bare module name, so their directory is put on
    ``sys.path`` here, at call time rather than import time.

    Returns
    -------
    dentstudy, replica : module
        The analysis helpers and the fixed-covariance replica.
    """
    if str(HELPER_DIR) not in sys.path:
        sys.path.insert(0, str(HELPER_DIR))
    import dentstudy
    import replica

    return dentstudy, replica


def show_table(df: pd.DataFrame, caption: str = "") -> None:
    """
    Print a table with a caption.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to print.
    caption : str, optional
        Caption printed above the table.
    """
    if caption:
        print(f"\n== {caption}")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


def make_labeller(ds):
    """
    Build the plain-text observable labeller used in tables and axes.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.

    Returns
    -------
    callable
        ``lab(name)`` -> display label with TeX markup stripped.
    """
    labels = {
        n: o["label"].replace("$", "").replace("\\", "").replace("_{", "_").replace("}", "")
        for n, o in ds.OBS.items()
    }

    def lab(s: str) -> str:
        return labels.get(s, s.replace("_", " "))

    return lab


def qc_report(ds, res: pd.DataFrame) -> None:
    """
    Quality control per study: shape-only live, ndf = nbins - 1, no inert injections.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    res : pandas.DataFrame
        Output of ``ds.load_results()``.
    """
    qc = (
        res[res.variant == "v1"]
        .groupby("tag")
        .agg(
            det=("det", "first"),
            tier=("tier", "first"),
            live=("shapeonly_live", "first"),
            ndf_ok=("ndf", lambda s: True),
            n_rows=("chi2", "size"),
        )
        .reset_index()
    )
    bad_ndf = res[(res.variant == "v1") & (res.ndf != res.nbins - 1)]
    print(
        f"studies with results: {len(qc)}  |  shape-only verified live (chi2 differs without "
        f"the flag): {int(qc.live.fillna(False).sum())} of {len(qc)}"
    )
    print(f"rows where ndf != nbins-1 (would indicate no area normalisation): {len(bad_ndf)}")
    summ = (
        qc.groupby(["det", "tier"])
        .agg(studies=("tag", "size"), shapeonly_live=("live", lambda s: int(s.fillna(False).sum())))
        .reset_index()
    )
    show_table(summ, "Studies by detector and tier")
    # Audit: a DENT injection that PROfit silently ignored leaves chi2 identical to v1.
    w = res[res.variant.isin(["v1", "inj1p0"])].pivot_table(
        index=["tag", "channel"], columns="variant", values="chi2"
    )
    if "inj1p0" in w:
        inert = w[(w.v1 - w.inj1p0).abs() < 1e-9].reset_index().tag.unique().tolist()
        print(
            "studies whose DENT injection left chi2 unchanged (invalid, must be re-run): "
            f"{inert if inert else 'none'}"
        )


def tier_a_table(res: pd.DataFrame, prim: list[str], lab) -> pd.DataFrame:
    """
    Tier A: detvar CV and DENT samples used directly (MC statistics only).

    Parameters
    ----------
    res : pandas.DataFrame
        Output of ``ds.load_results()``.
    prim : list of str
        Primary observable names.
    lab : callable
        Label function from :func:`make_labeller`.

    Returns
    -------
    pandas.DataFrame
        Indexed by study with ``chi2_cv``, ``chi2_dent``, ``ndf``,
        ``dchi2`` (CV minus DENT; positive means DENT is closer to the
        data), ``label`` and ``primary``, sorted by ``dchi2`` descending.
    """
    sel = res[(res.tier == "A") & (res.variant == "v1")]
    cv = sel[sel.which == "cv"].set_index("study")
    dn = sel[sel.which == "dent"].set_index("study")
    ta = pd.DataFrame({"chi2_cv": cv.chi2, "chi2_dent": dn.chi2, "ndf": cv.ndf}).dropna()
    ta["dchi2"] = ta.chi2_cv - ta.chi2_dent
    ta["label"] = [lab(s) for s in ta.index]
    ta["primary"] = ta.index.isin(prim)
    return ta.sort_values("dchi2", ascending=False)


def tier_b_table(ds, res: pd.DataFrame, prim: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Tier B: production CV, full systematics, area-normalised and shape-only.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    res : pandas.DataFrame
        Output of ``ds.load_results()``.
    prim : list of str
        Primary observable names.

    Returns
    -------
    tb : pandas.DataFrame
        One row per observable: ``ndf``, ``chi2``/``p`` per variant and
        ``dchi2 (0 vs 1)``; primaries first. Empty if there are no results.
    sel : pandas.DataFrame
        The underlying Tier B rows of ``res``.
    """
    sel = res[(res.tier == "B") & (res.channel == 0) & res.variant.isin(list(ds.VARIANTS_B))]
    if sel.empty:
        print("no Tier B results yet")
        return pd.DataFrame(columns=["observable", "dchi2 (0 vs 1)"]), sel
    wide = sel.pivot_table(index="study", columns="variant", values=["chi2", "p"])
    ndf = sel[sel.variant == "v1"].set_index("study").ndf
    tb = pd.DataFrame({"ndf": ndf})
    for v in ds.VARIANTS_B:
        if ("chi2", v) in wide.columns:
            tb[f"chi2 {v}"] = wide[("chi2", v)]
            tb[f"p {v}"] = wide[("p", v)]
    if "chi2 v1" in tb and "chi2 inj1p0" in tb:
        tb["dchi2 (0 vs 1)"] = tb["chi2 v1"] - tb["chi2 inj1p0"]
    order = [p for p in prim if p in tb.index] + [s for s in tb.index if s not in prim]
    tb = tb.reindex(order).reset_index().rename(columns={"study": "observable"})
    tb.insert(1, "primary", tb.observable.isin(prim))
    return tb, sel


def ablation_table(res: pd.DataFrame, prim: list[str]) -> pd.DataFrame:
    """
    Covariance ablation: how much of the DENT effect survives as systematics are added.

    Parameters
    ----------
    res : pandas.DataFrame
        Output of ``ds.load_results()``.
    prim : list of str
        Primary observable names.

    Returns
    -------
    pandas.DataFrame
        One row per observable with ``chi2(0)``, ``chi2(1)`` and
        ``dchi2`` columns per entry of :data:`ABLATION_CONFIGS`.
    """
    sel = res[(res.tier == "B") & (res.channel == 0)]
    rows = []
    for name in [p for p in prim if p in set(sel.study)] + [
        s for s in sorted(set(sel.study)) if s not in prim
    ]:
        g = sel[sel.study == name].set_index("variant").chi2
        r = {"observable": name, "primary": name in prim}
        for cfg, (a, b) in ABLATION_CONFIGS.items():
            if a in g.index and b in g.index:
                r[f"{cfg}: chi2(0)"] = g[a]
                r[f"{cfg}: chi2(1)"] = g[b]
                r[f"{cfg}: dchi2"] = g[a] - g[b]
        rows.append(r)
    return pd.DataFrame(rows)


def plot_ablation(ds, abl: pd.DataFrame, lab, path: Path):
    """
    Bar chart of dchi2 per covariance configuration for the primaries.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    abl : pandas.DataFrame
        Output of :func:`ablation_table`.
    lab : callable
        Label function from :func:`make_labeller`.
    path : Path
        Destination file.

    Returns
    -------
    Path or None
        The path written, or ``None`` if there is nothing to draw.
    """
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis

    ds.style()
    pr = abl[abl.primary & abl["stat + DENT: dchi2"].notna()]
    if not len(pr):
        return None
    fig, ax = plt.subplots(figsize=(11, 4.8))
    x = np.arange(len(pr))
    w = 0.2
    for i, (cfg, col) in enumerate(zip(ABLATION_CONFIGS, ABLATION_COLORS, strict=True)):
        if f"{cfg}: dchi2" in pr:
            ax.bar(x + (i - 1.5) * w, pr[f"{cfg}: dchi2"], w, color=col, label=cfg)
    ax.axhline(0, color="0.3", lw=1)
    ax.set_yscale("symlog", linthresh=5)
    ax.set_xticks(x)
    ax.set_xticklabels([lab(o) for o in pr.observable], rotation=35, ha="right")
    ax.set_ylabel(r"$\chi^2(\theta{=}0)-\chi^2(\theta{=}1)$")
    ax.legend(fontsize=9, ncol=2)
    mark_axis(ax, WATERMARK_DATA, vadj=0.05)
    return common.save_figure(fig, path)


def render_panels(ds, prim: list[str], lab, outdir: Path) -> list[tuple[str, str, Path]]:
    """
    PROfit-style data/MC panels (mutools) for the primary observables.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    prim : list of str
        Primary observable names.
    lab : callable
        Label function from :func:`make_labeller`.
    outdir : Path
        Output directory for the panels.

    Returns
    -------
    list of tuple
        ``(observable, variant, path)`` per panel rendered. Missing
        PROplot histogram files are skipped.
    """
    made = []
    for name in prim:
        edges = ds.OBS[name]["edges"]
        for variant, _ in PANEL_VARIANTS:
            f = ds.RUN / f"dent_sbnd_run1_dd_{name}_{variant}_PROplot_histograms.root"
            if not f.exists():
                continue
            p = ds.render_panel(
                f, outdir, f"{name}_{variant}", lab(name), xlim=[edges[0], edges[-1]]
            )
            made.append((name, variant, p))
    return made


def replica_table(rp, prim: list[str], lab) -> pd.DataFrame:
    """
    Fixed-covariance change in chi-square against pseudo-experiments.

    Parameters
    ----------
    rp : module
        The ``replica`` helper module.
    prim : list of str
        Primary observable names.
    lab : callable
        Label function from :func:`make_labeller`.

    Returns
    -------
    pandas.DataFrame
        Replica summary for the primaries with a ``label`` column and
        ``nearer truth CV?``; also prints the replica-versus-PROfit
        validation.
    """
    rs = pd.read_csv(REPLICA_SUMMARY)
    rs = rs.set_index("study").loc[[p for p in prim if p in set(rs.study)]]
    dmax = 0.0
    for name in rs.index:
        v = rp.validate(rp.load(name))
        dmax = max(dmax, max(abs(x) for x in v.values()))
    print(
        f"replica versus PROfit chi-square, largest absolute difference over {len(rs)} "
        f"primaries and 3 strengths: {dmax:.3f}"
    )
    t = rs.reset_index().copy()
    t.insert(1, "label", [lab(s) for s in t.study])
    t["nearer truth CV?"] = abs(t.dchi2_fixed_data - t.mean_0) < abs(t.dchi2_fixed_data - t.mean_1)
    return t


def print_replica(t: pd.DataFrame) -> None:
    """
    Print the replica table and its one-line summary.

    Parameters
    ----------
    t : pandas.DataFrame
        Output of :func:`replica_table`.
    """
    cols = [
        "label", "profit_dchi2", "dchi2_fixed_data", "mean_0", "sd_0", "p2_0", "mean_0.5",
        "sd_0.5", "p2_0.5", "mean_1", "sd_1", "p2_1", "dprime", "nearer truth CV?",
    ]  # fmt: skip
    names = {
        "profit_dchi2": "PROfit change",
        "dchi2_fixed_data": "fixed-cov change (data)",
        "mean_0": "mean|CV", "sd_0": "sd|CV", "p2_0": "p2|CV",
        "mean_0.5": "mean|half", "sd_0.5": "sd|half", "p2_0.5": "p2|half",
        "mean_1": "mean|DENT", "sd_1": "sd|DENT", "p2_1": "p2|DENT",
        "dprime": "d'",
    }  # fmt: skip
    show_table(
        t[cols].rename(columns=names),
        "Fixed-covariance chi2(0) - chi2(1) in data, against 20,000 pseudo-experiments for "
        "truth = CV, half DENT, full DENT (two-sided p)",
    )
    pmin = t[["p2_0", "p2_0.5", "p2_1"]].min().min()
    print(
        f"data nearer the CV-truth mean than the DENT-truth mean: "
        f"{int(t['nearer truth CV?'].sum())} of {len(t)}; smallest two-sided p over all "
        f"truths: {pmin:.3f}; d' range {t.dprime.min():.2f} to {t.dprime.max():.2f}"
    )


def plot_replica(ds, t: pd.DataFrame, path: Path):
    """
    Data fixed-covariance dchi2 against the three pseudo-experiment truths.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    t : pandas.DataFrame
        Output of :func:`replica_table`.
    path : Path
        Destination file.

    Returns
    -------
    Path
        The path written.
    """
    import matplotlib.pyplot as plt
    from mutools.plotting.helpers import mark_axis

    ds.style()
    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(t) + 1.8))
    y = np.arange(len(t))[::-1]
    truths = (
        (0.2, "C0", "0", "truth CV"),
        (0.0, "C2", "0.5", "truth half DENT"),
        (-0.2, "C1", "1", "truth full DENT"),
    )
    for off, col, key, nm in truths:
        ax.errorbar(
            t[f"mean_{key}"], y + off, xerr=t[f"sd_{key}"], fmt="o", color=col, ms=5, capsize=2,
            label=nm + " (mean +- 1 sd)",
        )  # fmt: skip
    ax.plot(t.dchi2_fixed_data, y, "k*", ms=11, label="data", zorder=5)
    ax.axvline(0, color="0.4", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(t.label)
    ax.set_xlabel(
        r"fixed-covariance $\chi^2(\theta{=}0)-\chi^2(\theta{=}1)$"
        "\n[positive: DENT-shifted prediction closer to data]"
    )
    ax.legend(fontsize=8, loc="lower right")
    mark_axis(ax, WATERMARK_DATA, vadj=0.05)
    return common.save_figure(fig, path)


def tier_c_table(ds, tc: pd.DataFrame, tb: pd.DataFrame, prim: list[str], lab) -> pd.DataFrame:
    """
    Tier C: fake-data p-values with Holm correction.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    tc : pandas.DataFrame
        Output of ``ds.load_tierC()``.
    tb : pandas.DataFrame
        Tier B table from :func:`tier_b_table` (supplies the data dchi2).
    prim : list of str
        Primary observable names.
    lab : callable
        Label function from :func:`make_labeller`.

    Returns
    -------
    pandas.DataFrame
        One row per primary with throws; empty if none are available.
    """
    rows = []
    for name in prim:
        sub = tc[tc.study == name] if len(tc) else tc
        row_b = tb[tb.observable == name] if "observable" in tb else tb.iloc[0:0]
        if len(sub) == 0 or len(row_b) == 0:
            continue
        d_obs = float(row_b["dchi2 (0 vs 1)"].iloc[0])
        r = ds.tierC_pvalues(sub, d_obs)
        r.update(study=name, label=lab(name), dchi2_data=d_obs)
        rows.append(r)
    tcs = pd.DataFrame(rows)
    if len(tcs):
        tcs["p_cv_holm"] = ds.holm(tcs.p_cv)
        tcs["p_dent_holm"] = ds.holm(tcs.p_dent)
    return tcs


def plot_tier_c_hists(ds, tc: pd.DataFrame, tcs: pd.DataFrame, path: Path, ncol: int = 4):
    """
    Grid of Tier C dchi2 histograms, one per primary observable.

    Parameters
    ----------
    ds : module
        The ``dentstudy`` helper module.
    tc : pandas.DataFrame
        Output of ``ds.load_tierC()``.
    tcs : pandas.DataFrame
        Output of :func:`tier_c_table`.
    path : Path
        Destination file.
    ncol : int, optional
        Number of panel columns.

    Returns
    -------
    Path
        The path written.
    """
    import matplotlib.pyplot as plt

    n = len(tcs)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (_, r) in zip(axes, tcs.iterrows(), strict=False):
        ds.plot_tierC_hist(tc[tc.study == r.study], r.dchi2_data, r.label, ax=ax)
    for ax in axes[n:]:
        ax.axis("off")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    return common.save_figure(fig, path, dpi=170)


def baseline_tables(res: pd.DataFrame) -> pd.DataFrame:
    """
    Regenerated detector-access studies (baseline), one row per channel.

    Parameters
    ----------
    res : pandas.DataFrame
        Output of ``ds.load_results()``.

    Returns
    -------
    pandas.DataFrame
        ``detector``, ``study``, ``channel``, ``ndf``, ``chi2``, ``p``
        and, where available, ``chi2 CV+DENT``, ``dchi2 (0 vs 1)`` and
        ``chi2 no DENT cov``.
    """
    base = res[(res.tier == "baseline") & res.variant.isin(["v1", "exDENT", "inj1p0"])]
    rows = []
    for (det, study, ch), g in base.groupby(["det", "study", "channel"]):
        d = g.drop_duplicates("variant", keep="last").set_index("variant")
        r = {
            "detector": det,
            "study": study,
            "channel": ch,
            "ndf": int(d.loc["v1", "ndf"]),
            "chi2": d.loc["v1", "chi2"],
            "p": d.loc["v1", "p"],
        }
        if "inj1p0" in d.index:
            r["chi2 CV+DENT"] = d.loc["inj1p0", "chi2"]
            r["dchi2 (0 vs 1)"] = d.loc["v1", "chi2"] - d.loc["inj1p0", "chi2"]
        if "exDENT" in d.index:
            r["chi2 no DENT cov"] = d.loc["exDENT", "chi2"]
        rows.append(r)
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; defaults to ``sys.argv[1:]``.
    """
    parser = common.base_parser(__doc__)
    parser.set_defaults(output=common.FIGURES / "dent_study", format="png")
    parser.add_argument(
        "--only",
        nargs="+",
        choices=STAGES,
        default=list(STAGES),
        help="stages to run (default: all)",
    )
    args = parser.parse_args(argv)
    common.set_provenance(__file__, [common.RUNS / "dent_study" / "results", REPLICA_SUMMARY],
                          {"stages": list(args.only or STAGES)})
    ext = args.format
    out = args.output
    run = set(args.only)

    ds, rp = load_helpers()
    ds.style()
    lab = make_labeller(ds)
    res = ds.load_results()
    tc = ds.load_tierC()
    prim = [p for p in ds.PRIMARY if p in ds.OBS]
    n_tc = tc.study.nunique() if len(tc) else 0
    print(
        f"{res.tag.nunique()} study runs loaded; primary set = {len(prim)} observables; "
        f"Tier C throws loaded for {n_tc} studies"
    )

    if "qc" in run:
        qc_report(ds, res)

    ta = tier_a_table(res, prim, lab)
    if "tierA" in run:
        tab = ta.reset_index()[["study", "primary", "ndf", "chi2_cv", "chi2_dent", "dchi2"]]
        show_table(
            tab.rename(
                columns={
                    "study": "observable",
                    "chi2_cv": "chi2 (detvar CV)",
                    "chi2_dent": "chi2 (detvar DENT)",
                    "dchi2": "chi2(CV) - chi2(DENT)",
                }
            ),
            "Tier A: detvar samples used directly (MC statistics only)",
        )
        path = out / f"tierA_primary.{ext}"
        out.mkdir(parents=True, exist_ok=True)
        ds.plot_tierA(ta[ta.primary], savepath=path, title="Tier A (detvar CV vs DENT)")
        print("wrote", path)

    tb, sel_b = tier_b_table(ds, res, prim)
    if "tierB" in run:
        if len(tb):
            show_table(
                tb,
                "Tier B: production CV (cvrw + cvwgt), full systematics, "
                "area-normalised and shape-only",
            )
        prim_b = sel_b[sel_b.study.isin(prim)]
        if len(prim_b):
            names = [p for p in prim if p in set(prim_b.study)]
            path = out / f"tierB_pvalues_primary.{ext}"
            out.mkdir(parents=True, exist_ok=True)
            ds.plot_pvalues(prim_b, names, [lab(p) for p in names], savepath=path)
            print("wrote", path)

    if "ablation" in run:
        abl = ablation_table(res, prim)
        if "stat + DENT: dchi2" in abl:
            show_table(
                abl[["observable", "primary"] + [c for c in abl.columns if c.endswith("dchi2")]],
                "Covariance ablation: dchi2 = chi2(theta=0) - chi2(theta=1) "
                "(positive: the DENT shift helps), by covariance content",
            )
            path = plot_ablation(ds, abl, lab, out / f"tierB_ablation.{ext}")
            if path:
                print("wrote", path)
        else:
            print("ablation results not available yet")

    if "panels" in run:
        made = render_panels(ds, prim, lab, out / "panels")
        print(len(made), "panels rendered")
        for _name, _variant, p in made:
            print("wrote", p)

    if "replica" in run:
        t = replica_table(rp, prim, lab)
        print_replica(t)
        print("wrote", plot_replica(ds, t, out / f"replica_dchi2.{ext}"))

    if "tierC" in run:
        tcs = tier_c_table(ds, tc, tb, prim, lab)
        if len(tcs):
            show_table(
                tcs[
                    [
                        "label", "n", "dchi2_data", "mean_cv", "sd_cv", "mean_dent", "sd_dent",
                        "p_cv", "p_cv_holm", "p_dent", "p_dent_holm",
                    ]  # fmt: skip
                ],
                "Tier C: Delta chi2 = chi2(theta=0) - chi2(theta=1). "
                "p_CV = P(dchi2 >= data | truth CV); p_DENT = P(dchi2 <= data | truth DENT)",
            )
            out.mkdir(parents=True, exist_ok=True)
            path = out / f"tierC_dchi2.{ext}"
            ds.plot_dchi2(tcs, savepath=path)
            print("wrote", path)
            print("wrote", plot_tier_c_hists(ds, tc, tcs, out / f"tierC_hists.{ext}"))
        else:
            print("no Tier C throws available yet")

    if "baseline" in run:
        bt = baseline_tables(res)
        for det, name in (("sbnd_run1", "SBND"), ("icarus_run2", "ICARUS Run 2")):
            show_table(
                bt[bt.detector == det].sort_values("p"),
                f"Regenerated {name} detector-access studies (shape-only), sorted by p-value",
            )


if __name__ == "__main__":
    main()
