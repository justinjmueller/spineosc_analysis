#!/usr/bin/env python3
"""How much does each detector systematic matter for the oscillation sensitivity? Uses the PROfit chi2 surfaces (runs/legacy_root/contours_*_surf.root, produced BEFORE the v1.1.0v2 detector
ratios and cuts, so they describe the previous systematic configuration). 90% CL contour at delta chi2 = 1.6.  Metric: sin^2(2theta) at which chi2 crosses 1.6, per dm^2 row; and the excluded
area in (log sin^2 2theta, log dm^2).  'no X' = systematic X removed; improvement = (limit_full - limit_noX) / limit_full."""
import json, numpy as np, uproot, pathlib
P = pathlib.Path("/Users/mueller/Projects/spineosc/runs/legacy_root"); OUT = pathlib.Path("/Users/mueller/Projects/spineosc/docs/syst_study"); OUT.mkdir(exist_ok=True)
NAMES = {"full": "contours_full_v1", "noSBNDChargeScale": "contours_full_nochargescale", "noICARUSChargeScale": "contours_full_noicchargescale", "noSBNDDENT": "contours_full_nodent",
         "noSBNDWireModYZ": "contours_full_nowmyz", "noSBNDWireModXTh": "contours_full_nowmxtxw", "noICARUSWireModYZ": "contours_full_noicwmyz", "noICARUSWireModXTh": "contours_full_noicwmxtxw",
         "noICARUSSCE": "contours_full_noicsce", "noSBNDTrigger": "contours_full_notrigemu", "noICARUSTrigger": "contours_full_noictrigemu", "blind": "contours_blind_v1", "sbndunblind": "contours_sbndunblind_v1"}
def load(n):
    h = uproot.open(P / f"{n}_surf.root")["surf"]; return h.values(), h.axes[0].centers(), h.axes[1].centers()
res = {}; limits = {}
for k, n in NAMES.items():
    try: v, x, y = load(n)
    except Exception as e: continue
    lim = np.full(len(y), np.nan)
    # which axis index is sin^2(2theta)? the TH2D is (x=sin^2 2theta, y=dm^2) per the notebook: values shape (nx, ny)
    for j in range(len(y)):
        col = v[:, j]; idx = np.where(col > 1.6)[0]
        if len(idx) == 0: continue
        i = idx[0]
        if i == 0: lim[j] = x[0]; continue
        # log-linear interpolation in chi2 between i-1 and i
        c0, c1 = col[i - 1], col[i]; t = (1.6 - c0) / (c1 - c0) if c1 != c0 else 0.5; lim[j] = 10 ** (np.log10(x[i - 1]) + t * (np.log10(x[i]) - np.log10(x[i - 1])))
    area = float(np.sum(v > 1.6) / v.size); limits[k] = lim
    res[k] = dict(area_frac=area, limit=[None if not np.isfinite(l) else float(l) for l in lim])
dm2 = load("contours_full_v1")[2]
sel = (dm2 > 0.5) & (dm2 < 20)
out = {"dm2": dm2.tolist(), "systems": {}}
full = limits["full"]
for k, lim in limits.items():
    ok = sel & np.isfinite(lim) & np.isfinite(full)
    impr = float(np.nanmean((full[ok] - lim[ok]) / full[ok])) if k != "full" else 0.0
    out["systems"][k] = dict(area_frac=res[k]["area_frac"], mean_limit_improvement_vs_full_0p5_to_20eV2=impr, limit=res[k]["limit"])
    print(f"{k:22s} excluded-area fraction {res[k]['area_frac']:.4f}   mean limit improvement vs full (0.5<dm2<20): {100*impr:+.2f}%")
json.dump(out, open(OUT / "sensitivity_surfaces.json", "w"))
E=0
