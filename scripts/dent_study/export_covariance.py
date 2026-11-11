#!/usr/bin/env python3
"""Exports, for one dd observable, what is needed to REPRODUCE PROfit's shape-only chi2 outside PROfit and to hold the covariance fixed:
data counts n, the area-normalised prediction p at DENT strength theta = 0, 0.5, 1.0 (fixed injections), and the collapsed total fractional
covariance F at each theta, plus PROfit's own chi2 for each (the check). Validated: chi2 = (n-p)^T [diag(n) + F o (p p^T)]^-1 (n-p).
Note --inject-systs-cv re-derives the systematic covariance around the SHIFTED prediction, so PROfit's chi2(theta=1) mixes the shift with a
covariance change; holding F at its theta=0 form separates them. Output: runs/dent_study/replica/<obs>.npz"""
import os, pathlib, re, subprocess, sys, glob
import numpy as np, uproot
REPO = pathlib.Path(__file__).resolve().parents[2]; RUN = REPO / "runs" / "dent_study"; PROFIT = REPO.parent / "profit-asymmetric/build/bin/PROfit"
obs = sys.argv[1]; tag = f"dent_cov_{obs}"
ini_src = (REPO / f"xml/dent_study/sbnd_run1/dd/{obs}_plot.ini").read_text()
(RUN / "inis").mkdir(exist_ok=True); (RUN / "replica").mkdir(exist_ok=True)
ini = RUN / "inis" / f"{tag}.ini"; ini.write_text(re.sub(r"^log = .*$", f"log = {tag}.log", re.sub(r"^tag = .*$", f"tag = {tag}", ini_src, flags=re.M), flags=re.M))
names = sorted(set(re.findall(r'<allowlist[^>]*type="hist1d"[^>]*\sname="(DENT[^"]*)"', re.sub(r"<!--.*?-->", "", (REPO / f"xml/dent_study/sbnd_run1/dd/{obs}.xml").read_text(), flags=re.S))))
def run(label, extra, sub):
    subprocess.run([str(PROFIT), "--config", str(ini), "-o", label, *extra, *sub], cwd=RUN, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, timeout=3600)
run("x", [], ["process"])
out = {}
for theta, label, extra in ((0.0, "t0", []), (0.5, "t05", ["--inject-systs-cv", *[x for n in names for x in (n, "0.5")]]), (1.0, "t1", ["--inject-systs-cv", *[x for n in names for x in (n, "1.0")]])):
    run(label, extra, ["plot", "--with-covar"])
    h = uproot.open(RUN / f"{tag}_{label}_PROplot_histograms.root"); hh = h["hist1d"].arrays(library="np")
    def spec(prefix):
        m = (hh["prefix"] == prefix) & (hh["variable"] == 0) & (hh["channel"] == 0) & (hh["detector"] == 0) & (hh["subchannel"] == 4294967295)
        _, first = np.unique(hh["bin_index"][m], return_index=True)
        return hh["bin_low_edge"][m][first], hh["bin_high_edge"][m][first], hh["bin_content"][m][first], hh["bin_error"][m][first]
    lo, hi, cv, _ = spec("CV"); _, _, dat, dat_e = spec("DATA"); w = hi - lo
    n = (dat / dat_e) ** 2; p = cv * w / np.sum(cv * w) * n.sum()
    F = uproot.open(RUN / f"{tag}_{label}_PROplot.root")["Covariance/collapsed_total_frac_cov"].values()
    c = h["chi2"].arrays(library="np"); out.update({f"n": n, f"lo": lo, f"hi": hi, f"p_{theta:g}": p, f"F_{theta:g}": F, f"chi_profit_{theta:g}": float(c["chi2"][0])})
np.savez(RUN / "replica" / f"{obs}.npz", **out)
for f in glob.glob(str(RUN / f"{tag}_*")): os.remove(f)
for f in glob.glob(str(RUN / "logs" / f"{tag}_*")): os.remove(f)
print("saved", obs)
