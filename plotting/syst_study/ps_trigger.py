#!/usr/bin/env python3
"""Trigger emulation: reproduce TriggerEmulationM1/P1 from the nominal `selected` tree (event passes if scale*event_largest_flash_pe > threshold), explain the energy dependence, and
quantify the interplay with the other (non-trigger) detvar ratios, which are computed without the trigger requirement.  Output docs/syst_study/persyst/trigger.json"""
from ps_common import *
from common import sel_tree, SBND_NOS, icarus_nos
import uproot
CFG = {"sbnd_run1": dict(scale=0.647, sigma=0.004, thr=2000, fn=SBND_NOS, tree="sbnd_large", hist="sbnd_run1"),
       "icarus_run2": dict(scale=0.615, sigma=0.018, thr=6000, fn=icarus_nos("run2"), tree="icarus_run2", hist="icarus_run2"),
       "icarus_run4": dict(scale=0.361, sigma=0.013, thr=3000, fn=icarus_nos("run4"), tree="icarus_run4", hist="icarus_run4")}
f = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root"); out = {}
def hist_ratio(E, pe, scale, thr, mask=None):
    ok = np.isfinite(pe) & (scale * pe > thr) if mask is None else mask & np.isfinite(pe) & (scale * pe > thr)
    return np.histogram(E[ok], BINS)[0].astype(float)
for det, c in CFG.items():
    d, m = sel_tree(c["fn"], c["tree"], det); E = d["reco_visible_energy"].to_numpy(float); pe = d["event_largest_flash_pe"].to_numpy(float); m = np.asarray(m)
    N = hist_ratio(E, pe, c["scale"], c["thr"], m); res = dict(cfg={k: v for k, v in c.items() if k in ("scale", "sigma", "thr")}, n_rows=int(len(d)), n_sel=int(m.sum()))
    for tag, s in (("M1", c["scale"] - c["sigma"]), ("P1", c["scale"] + c["sigma"])):
        V = hist_ratio(E, pe, s, c["thr"], m); r = V / N; prod = f[f"TriggerEmulation{tag}_{c['hist']}"].values(); res[tag] = dict(ratio=r.tolist(), prod=prod.tolist(), maxdiff=float(np.max(np.abs(r - prod))))
    res["n_nom_bin"] = N.tolist()
    # PE-threshold geometry
    pe_thr = c["thr"] / c["scale"]; res["pe_threshold"] = pe_thr; res["pe_window_pm1sigma"] = [c["thr"] / (c["scale"] + c["sigma"]), c["thr"] / (c["scale"] - c["sigma"])]
    sel = m & np.isfinite(pe); res["frac_fail_nominal_trigger"] = float((c["scale"] * pe[sel] <= c["thr"]).mean()); res["frac_in_1sigma_window"] = float(((pe[sel] > res["pe_window_pm1sigma"][0]) & (pe[sel] <= res["pe_window_pm1sigma"][1])).mean())
    # per-bin: fraction failing the nominal trigger (events that the detvar ratios include but the analysis would drop), PE quantiles and PE per GeV
    b = binidx(E); rows = []
    for i in range(NBIN):
        k = sel & (b == i); p = pe[k]
        rows.append(dict(n=int(k.sum()), fail=float((c["scale"] * p <= c["thr"]).mean()), pe_q=np.percentile(p, [5, 16, 50, 84, 95]).tolist(), pe_per_gev=float(np.median(p / E[k])),
                         in_win=float(((p > res["pe_window_pm1sigma"][0]) & (p <= res["pe_window_pm1sigma"][1])).mean())))
    res["bins"] = rows
    # PE distribution vs visible energy (2D, coarse) for the figure
    ok = sel & (b >= 0); res["pe_hist2d"] = np.histogram2d(E[ok], np.log10(np.clip(pe[ok], 1, None)), bins=[BINS, np.linspace(2.5, 5.5, 61)])[0].tolist()
    # linear-response check: ratio-1 for +-1 sigma should be opposite and proportional to sigma
    rM = np.array(res["M1"]["ratio"]) - 1; rP = np.array(res["P1"]["ratio"]) - 1
    res["asym"] = dict(M1_int=float(hist_ratio(E, pe, c["scale"] - c["sigma"], c["thr"], m).sum() / N.sum() - 1), P1_int=float(hist_ratio(E, pe, c["scale"] + c["sigma"], c["thr"], m).sum() / N.sum() - 1))
    # sensitivity: what sigma would give |r-1|=1% in the lowest bin
    out[det] = res; print(det, "maxdiff M1/P1", res["M1"]["maxdiff"], res["P1"]["maxdiff"], "fail", round(res["frac_fail_nominal_trigger"], 4), "r-1 M1[0:4]", np.round(rM[:4] * 100, 2), "P1", np.round(rP[:4] * 100, 2), "int", res["asym"])
# ---- effect of the missing trigger requirement on the OTHER detvar ratios (pair tables; trigger applied with nominal scale on each side's own flash PE)
inter = {}
for name in PAIRS:
    det = name.split("__")[0]; c = CFG[det]; df = load(name); bn, bv = prep(df)
    pen = col(df, "n", "event_largest_flash_pe"); pev = col(df, "v", "event_largest_flash_pe")
    tn = np.isfinite(pen) & (c["scale"] * pen > c["thr"]); tv = np.isfinite(pev) & (c["scale"] * pev > c["thr"])
    bn_t = np.where(tn, bn, -1); bv_t = np.where(tv | (df.in_n.to_numpy() & ~df.in_v.to_numpy() & False), bv, -1)
    n0, v0 = counts(bn, bv); n1, v1 = counts(bn_t, bv_t); r0 = v0 / n0; r1 = v1 / n1
    # same trigger decision for both sides from the NOMINAL flash PE (the analysis applies it to the nominal MC event): events gained (no nominal PE) take the varied PE
    pe_any = np.where(np.isfinite(pen), pen, pev); t_any = np.isfinite(pe_any) & (c["scale"] * pe_any > c["thr"])
    bn_a = np.where(t_any, bn, -1); bv_a = np.where(t_any, bv, -1); n2, v2 = counts(bn_a, bv_a); r2 = v2 / n2
    pe_same = float(np.nanmean(pen[df.in_n & df.in_v] == pev[df.in_n & df.in_v]))
    inter[name] = dict(frac_fail_nom=float((~tn[df.in_n.to_numpy()]).mean()), r_all=r0.tolist(), r_trig_own=r1.tolist(), r_trig_nomPE=r2.tolist(), max_dr_own=float(np.max(np.abs(r1 - r0))), max_dr_nomPE=float(np.max(np.abs(r2 - r0))),
                       rms_dr_nomPE=float(np.sqrt(np.mean((r2 - r0) ** 2))), rms_sys=float(np.sqrt(np.mean((r0 - 1) ** 2))), int_all=float(v0.sum() / n0.sum() - 1), int_trig=float(v2.sum() / n2.sum() - 1), pe_identical_frac_kept=pe_same)
    print(name, "fail nom", round(inter[name]["frac_fail_nom"], 4), "int all/trig", round(100 * inter[name]["int_all"], 2), round(100 * inter[name]["int_trig"], 2), "max|dr|", round(inter[name]["max_dr_nomPE"], 4), "rms sys", round(inter[name]["rms_sys"], 4))
out["interplay"] = inter
json.dump(out, open(PSFIG / "trigger.json", "w"))
