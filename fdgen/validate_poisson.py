#!/usr/bin/env python3
"""Check that the per-event Poisson counts reproduce PROfit's per-bin Poisson(mu_bin) statistics.

PROfit throws every reco bin as Poisson(mu_bin). fdgen throws every MC event as Poisson(lambda_i); a bin then holds
sum k_i ~ Poisson(sum lambda_i). This script draws many datasets from the same event weights fdgen uses (ICARUS Run 2
beam MC with the by-hand oscillation applied) and compares, per reco-energy bin and for the total, the mean, the
variance and the full distribution (chi2 of the count histogram against Poisson(mu)) with the bin-level recipe.

  PYTHONPATH=$(root-config --libdir) .venv/bin/python -m fdgen.validate_poisson [--throws 2000]
"""
import argparse
import os

import numpy as np
import uproot
import yaml
from scipy import stats

from . import events, osc, xmlinfo

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--throws", type=int, default=2000)
    ap.add_argument("--tree", default="events/icarus_run2/selected")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(os.path.join(REPO, "fdgen/configs/legacy_1d/poisson_osc.yaml")))
    info = xmlinfo.load(os.path.join(REPO, cfg["dataset"]["mc_xml"]))
    f = next(x for x in info["files"] if x["treename"] == a.tree)
    path = os.path.join(REPO, cfg["dataset"]["base_dir"], f["filename"])
    ev = events.evaluate(path, f["treename"], f["branches"])
    edges = osc.le_bin_edges(*info["le_bins"])
    p = osc.pmumu_profit(ev["le"], edges, 2.0, 0.2)
    ratio = np.where((ev["rule"] == 1) & ~np.isnan(p), p, 1.0)
    lam = ev["aw"] * xmlinfo.pot_scale(info, f) * ratio
    reco = uproot.open(path)[a.tree]["reco_visible_energy"].array(library="np")
    sel = lam > 0
    # 8 reco bins of equal expected content
    cum = np.cumsum(lam[sel][np.argsort(reco[sel])])
    cuts = np.interp(np.linspace(0, cum[-1], 9)[1:-1], cum, np.sort(reco[sel]))
    ib = np.digitize(reco, cuts)
    bins = [sel & (ib == k) for k in range(8)]
    mu = np.array([lam[b].sum() for b in bins])
    rng = np.random.default_rng(12345)
    counts = np.zeros((a.throws, 8))
    for t in range(a.throws):
        k = rng.poisson(np.clip(lam, 0, None))
        for j, b in enumerate(bins):
            counts[t, j] = k[b].sum()
    print(f"{a.tree}: {sel.sum()} events with lambda > 0, {a.throws} throws; sum lambda = {lam.sum():.1f}")
    print("bin   mu       mean     (mean-mu)/sem   var/mu   chi2/ndf of count histogram vs Poisson(mu)   p")
    ok = True
    for j in range(8):
        c = counts[:, j]
        z = (c.mean() - mu[j]) / (np.sqrt(mu[j] / a.throws))
        vr = c.var(ddof=1) / mu[j]
        lo, hi = max(0, int(mu[j] - 4 * np.sqrt(mu[j]))), int(mu[j] + 4 * np.sqrt(mu[j]))
        edges_h = np.arange(lo, hi + 2) - 0.5
        obs, _ = np.histogram(c, edges_h)
        exp = a.throws * np.diff(stats.poisson.cdf(edges_h[1:-0 or None] if False else edges_h, mu[j]))
        exp *= obs.sum() / exp.sum()
        m = exp > 5
        chi2 = ((obs[m] - exp[m]) ** 2 / exp[m]).sum()
        pv = stats.chi2.sf(chi2, m.sum() - 1)
        ok &= abs(z) < 4 and 0.9 < vr < 1.1 and pv > 1e-3
        print(f"{j}  {mu[j]:8.1f} {c.mean():8.1f}   {z:+6.2f}        {vr:6.3f}   {chi2:6.1f}/{m.sum() - 1:<3d}   {pv:.3f}")
    tot = counts.sum(1)
    print(f"total  mu={mu.sum():.1f} mean={tot.mean():.1f} var/mu={tot.var(ddof=1) / mu.sum():.3f}")
    print("PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
