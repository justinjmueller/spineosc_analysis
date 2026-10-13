#!/usr/bin/env python3
"""Runs one PROfit study end to end and harvests its chi2 numbers.

Steps (each logged; a failing step is recorded and the remaining independent steps still run):
  1. process   area-norm + shapeonly (from the generated ini)                      -> label v1
  2. noshape   same ini with the shapeonly line removed; chi2 must DIFFER from v1   -> proof the flag is live
  3. SBND only, DENT variants (fixed injections, never a fit):
       exDENT  DENT removed from the covariance, central value = CV
       inj0p5 / inj1p0   central value shifted by 0.5 / 1.0 x DENT, DENT kept in the covariance
Results are written to runs/dent_study/results/<tag>.json. Caches (.bin) and PDFs are deleted at the end unless --keep-bins.
"""
import argparse, json, os, pathlib, re, subprocess, sys, time, glob

REPO = pathlib.Path(__file__).resolve().parents[2]
RUN = REPO / "runs" / "dent_study"
PROFIT = REPO.parent / "profit-asymmetric" / "build" / "bin" / "PROfit"
DENT_NAME = "DENT_sbnd_Run1"

def dent_names(xml_path):
    """Names of the active (non-commented) hist1d DENT allowlist entries. Most studies have DENT_sbnd_Run1; the crossing-split
    studies have one per population (_stay, _xing, ...). Injecting/excluding the wrong name is silently ignored by PROfit."""
    txt = re.sub(r"<!--.*?-->", "", pathlib.Path(xml_path).read_text(), flags=re.S)
    return sorted(set(re.findall(r'<allowlist[^>]*type="hist1d"[^>]*\sname="(DENT[^"]*)"', txt)))

# The other detector-response entries that sit alongside DENT in every study's allowlist (hist1d names, active ones are read from the XML)
def detector_names(xml_path):
    txt = re.sub(r"<!--.*?-->", "", pathlib.Path(xml_path).read_text(), flags=re.S)
    return sorted(set(re.findall(r'<allowlist[^>]*type="hist1d"[^>]*\sname="([^"]*)"', txt)))

def sh(cmd, log, timeout):
    t0 = time.time()
    with open(log, "ab") as fh:
        fh.write(("\n$ " + " ".join(map(str, cmd)) + "\n").encode())
        try:
            r = subprocess.run(cmd, cwd=RUN, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = -9
    return rc, time.time() - t0

def read_chi2(path):
    import uproot
    with uproot.open(path) as f:
        t = f["chi2"].arrays(library="np")
    return [dict(detector=int(t["detector"][i]), channel=int(t["channel"][i]), variable=int(t["variable"][i]),
                 chi2=float(t["chi2"][i]), nbins=int(t["nbins"][i]), ndf=int(t["ndf"][i])) for i in range(len(t["chi2"]))
            if str(t["prefix"][i]) == "CV"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True); ap.add_argument("--study", required=True); ap.add_argument("--sub", default="")
    ap.add_argument("--keep-bins", action="store_true"); ap.add_argument("--skip-dent-variants", action="store_true")
    ap.add_argument("--ini", default=None, help="explicit ini (default derived from det/sub/study)")
    ap.add_argument("--tag-prefix", default="")
    ap.add_argument("--wait-for", default=None, help="comma-separated files that must exist (and be stable) before starting")
    ap.add_argument("--tierC", type=int, default=0, help="number of PSEUDO-EXPERIMENTS (Poisson + systematic nuisance throws) per (truth, test) combination; 0 = skip")
    ap.add_argument("--tierC-poisson", type=int, default=0, help="number of Poisson-only fake-data throws (no nuisance throws; too narrow to calibrate p-values); 0 = skip")
    a = ap.parse_args()
    if a.wait_for:                                   # ratio files are produced by a parallel job; poll up to 8 h
        deadline = time.time() + 8 * 3600
        while time.time() < deadline:
            fs = [pathlib.Path(x) for x in a.wait_for.split(",")]
            if all(f.exists() and f.stat().st_size > 1000 and time.time() - f.stat().st_mtime > 60 for f in fs): break
            time.sleep(60)
    ini = pathlib.Path(a.ini) if a.ini else REPO / "xml" / ("detector_data_access" if a.det.startswith("icarus") else "dent_study") / a.det / a.sub / f"{a.study}_plot.ini"
    tag = re.search(r"^tag\s*=\s*(\S+)", ini.read_text(), re.M).group(1)
    (RUN / "results").mkdir(exist_ok=True); (RUN / "logs").mkdir(exist_ok=True); (RUN / "inis").mkdir(exist_ok=True)
    slog = RUN / "logs" / f"{tag}.steps.log"
    res = dict(tag=tag, det=a.det, study=a.study, sub=a.sub, ini=str(ini), started=time.strftime("%F %T"), steps={}, chi2={})

    def hist_path(label): return RUN / f"{tag}_{label}_PROplot_histograms.root"
    def step(label, extra, timeout, subcmd="plot", ini_path=ini):
        rc, dt = sh([PROFIT, "--config", ini_path, "-o", label, "-l", f"logs/{tag}_{label}.log"] + extra + [subcmd], slog, timeout)
        ok = rc == 0 and hist_path(label).exists()
        res["steps"][label] = dict(rc=rc, seconds=round(dt, 1), ok=ok, cmd=" ".join(map(str, extra + [subcmd])))
        if ok:
            try: res["chi2"][label] = read_chi2(hist_path(label))
            except Exception as e: res["steps"][label]["harvest_error"] = repr(e); res["steps"][label]["ok"] = False
        return ok

    # 1. process (+ the [plot] section) with the generated ini
    ok = step("v1", [], 150 * 60, subcmd="process")
    if not ok:
        ok = step("v1", [], 150 * 60, subcmd="process")          # one retry (stalls under memory pressure have been seen)
    if ok:
        # 2. same ini without shapeonly: chi2 must change, ndf may not (ndf cannot prove shape-only)
        txt = re.sub(r"^shapeonly\s*=.*\n", "", ini.read_text(), flags=re.M)
        ini_ns = RUN / "inis" / f"{tag}_noshape.ini"; ini_ns.write_text(txt)
        step("noshape", [], 30 * 60, ini_path=ini_ns)
        names = dent_names(res["ini"].replace("_plot.ini", ".xml")); res["dent_names"] = names
        # 3. DENT variants (SBND only), fixed injections
        if a.det.startswith("sbnd") and not a.skip_dent_variants and names:
            pairs = lambda v: [x for n in names for x in (n, v)]
            step("exDENT", ["--exclude-systs", *names], 30 * 60)
            step("inj0p5", ["--inject-systs-cv", *pairs("0.5")], 30 * 60)
            step("inj1p0", ["--inject-systs-cv", *pairs("1.0")], 30 * 60)
            # Covariance ablation (fixed configurations; DENT stays in the covariance so that the injection is valid):
            #   stat+DENT, stat+all detector entries, everything except the OTHER detector entries. Each at strength 0 and 1.0.
            alln = detector_names(res["ini"].replace("_plot.ini", ".xml")); others = [n for n in alln if n not in names]
            res["ablation_names"] = dict(all_hist1d=alln, others=others)
            for cfg, cflags in (("abl_statDENT", ["--syst-list", "stat", *names]), ("abl_det", ["--syst-list", "stat", *alln]), ("abl_nodetother", ["--exclude-systs", *others])):
                if cfg == "abl_nodetother" and not others: continue
                step(cfg + "_t0", cflags, 30 * 60)
                step(cfg + "_t1", cflags + ["--inject-systs-cv", *pairs("1.0")], 30 * 60)
            # PROfit ignores unknown names but still exits 0: read the logs and demand that the variants really differ from v1
            for lab_, bad in (("exDENT", "matches no systematic"), ("inj0p5", "Unrecognized spline"), ("inj1p0", "Unrecognized spline")):
                lp = RUN / "logs" / f"{tag}_{lab_}.log"
                if lp.exists() and bad in lp.read_text(errors="ignore"):
                    res["steps"][lab_]["ok"] = False; res["steps"][lab_]["ignored"] = bad; res["chi2"].pop(lab_, None)
            if "v1" in res["chi2"] and "inj1p0" in res["chi2"] and all(abs(x["chi2"] - y["chi2"]) < 1e-9 for x, y in zip(res["chi2"]["v1"], res["chi2"]["inj1p0"])):
                res["steps"]["inj1p0"]["ok"] = False; res["steps"]["inj1p0"]["ignored"] = "chi2 identical to v1"; res["chi2"].pop("inj1p0", None)
    # 4. Tier C: fake-data throws. truth CV, DENT x0.5 or DENT x1.0 (fake data), prediction at test strength 0 or 1.0 (fixed injections;
    #    nothing is fitted). The same seed gives the same throw, so chi2(test 0) - chi2(test 1) is a paired difference.
    #    "pe": --pseudo-experiment (Poisson + Gaussian nuisance throws via the total covariance): the calibrated version.
    #    "po": --poisson-throw only: the nuisance parameters stay at their central values, so the spread is much too narrow.
    for key, nthrows, mflag, mlab in (("tierC_pe", a.tierC, "--pseudo-experiment", "PE"), ("tierC", a.tierC_poisson, "--poisson-throw", "C")):
        if not (ok and nthrows and a.det.startswith("sbnd") and names): continue
        res[key] = []; nfail = 0
        for seed in range(1, nthrows + 1):
            if nfail >= 4 and not res[key]:
                res["steps"][key + "_aborted"] = "first %d fake-data runs failed" % nfail; break
            for truth, tflags in (("cv", []), ("half", ["--inject-systs", *[x for n in names for x in (n, "0.5")]]), ("dent", ["--inject-systs", *[x for n in names for x in (n, "1.0")]])):
                for test, pflags in ((0.0, []), (1.0, ["--inject-systs-cv", *[x for n in names for x in (n, "1.0")]])):
                    lab = f"{mlab}_{truth}_{test:g}_{seed}"
                    rc, dt = sh([PROFIT, "--config", ini, "-o", lab, "-l", f"logs/{tag}_{lab}.log", "-s", str(seed), "--use-fake-data", mflag] + tflags + pflags + ["plot"], slog, 20 * 60)
                    hp = hist_path(lab)
                    if not (rc == 0 and hp.exists()): nfail += 1
                    else:
                        try:
                            c = read_chi2(hp)
                            res[key].append(dict(seed=seed, truth=truth, test=test, chi2=[x["chi2"] for x in c], ndf=[x["ndf"] for x in c]))
                        except Exception as e: res["steps"].setdefault(key + "_errors", []).append(repr(e))
                    for q in glob.glob(str(RUN / f"{tag}_{lab}_*")) + glob.glob(str(RUN / f"logs/{tag}_{lab}.log")): os.remove(q)
    res["finished"] = time.strftime("%F %T")
    # verification that shapeonly is live: chi2 differs between v1 and noshape in at least one row
    if "v1" in res["chi2"] and "noshape" in res["chi2"]:
        d = [abs(x["chi2"] - y["chi2"]) for x, y in zip(res["chi2"]["v1"], res["chi2"]["noshape"])]
        res["shapeonly_live"] = bool(max(d) > 1e-6)
    old = RUN / "results" / f"{tag}.json"
    if old.exists():
        try:
            prev = json.loads(old.read_text())
            for k in ("tierC", "tierC_pe"):
                if k not in res and prev.get(k): res[k] = prev[k]; res.setdefault("carried_over", []).append(k)
        except Exception: pass
    old.write_text(json.dumps(res, indent=1))
    if not a.keep_bins:
        for p in glob.glob(str(RUN / f"{tag}_*.bin")) + glob.glob(str(RUN / f"{tag}_*.pdf")): os.remove(p)
    sys.exit(0 if res["steps"].get("v1", {}).get("ok") else 1)

if __name__ == "__main__":
    main()
