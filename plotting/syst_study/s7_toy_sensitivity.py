#!/usr/bin/env python3
"""S7: ILLUSTRATIVE Asimov sensitivity toy for the near/far detector design, to see how much each detector systematic matters for the 2-flavour numu-disappearance limit and why
removing the SBND or the ICARUS charge-scale nuisance can give a similar improvement.  NOT a replacement for PROfit.
Model: reco-visible-energy spectra (17 analysis bins) of selected nominal events (v1.1.0 nosyst, calculator selection, no cvwgt) scaled to the 'full' scenario exposures (SBND 1e20, ICARUS Run2 2e20, Run4 3e20 POT);
numu CC events oscillate with P = 1 - s2 sin^2(1.27 dm2 L/E_true) using the stored true baseline and energy; Asimov data = no oscillation. Nuisances: all active detector systematics per detector/run
(ratio histograms of detvar/detsys-v1.1.0v2.root; one-sided hist1d with linear mirrored response, SCE and trigger two-sided with the +-1 histograms), a flux/xsec normalisation shared by all detectors (10%),
and an uncorrelated 3% normalisation per detector/run. chi2 profiled over nuisances; contour level 1.6."""
import json, sys, os, numpy as np, uproot, pathlib
from scipy.optimize import least_squares
from common import *
POT = {"sbnd_run1": (5.593437e20, 1e20), "icarus_run2": (1.040833e21, 2e20), "icarus_run4": (1.772140e21, 3e20)}
NOS = {"sbnd_run1": (SBND_NOS, "sbnd_large"), "icarus_run2": (icarus_nos("run2"), "icarus_run2"), "icarus_run4": (icarus_nos("run4"), "icarus_run4")}
COLS = ["reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", "true_cc", "true_pdg", "true_baseline", "true_neutrino_energy"]
DM2 = np.logspace(-1, 2, 30); S2 = np.logspace(-2, 0, 30)
H = uproot.open(REPO / "detvar" / "detsys-v1.1.0v2.root")
def hv(n): return H[n].values()
SYST = {"sbnd_run1": [("ChargeScale", "one", "ChargeScale15_sbnd_Run1"), ("WireModxThetaXW", "one", "WireModxThetaXW_sbnd_Run1"), ("MattModYZ", "one", "MattModYZ_sbnd_Run1"), ("DENT", "one", "DENT_sbnd_Run1"),
                      ("Trigger", "two", ("TriggerEmulationM1_sbnd_run1", "TriggerEmulationP1_sbnd_run1"))]}
for run, R in (("run2", "Run2"), ("run4", "Run4")):
    SYST[f"icarus_{run}"] = [("ChargeScale", "one", f"ChargeScale_icarus_{R}"), ("SCE", "two", (f"SCE0x_icarus_{R}", f"SCE2x_icarus_{R}")), ("WireModxThetaXW", "one", f"WireModxThetaXW_icarus_{R}"),
                             ("WireModYZ", "one", f"WireModYZ_icarus_{R}"), ("Trigger", "two", (f"TriggerEmulationM1_icarus_{run}", f"TriggerEmulationP1_icarus_{run}"))]
FLUX_SIG = float(os.environ.get('FLUX_SIG', 0.10)); NORM_SIG = float(os.environ.get('NORM_SIG', 0.03)); TAG = os.environ.get('TAG', ''); ONLY = os.environ.get('ONLY', '')
DET = list(POT); spec = {}
for det in DET:
    fn, tree = NOS[det]; d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="np"); m = selection_mask(det, d); E = d["reco_visible_energy"][m]
    sig = (d["true_cc"][m] > 0) & (np.abs(d["true_pdg"][m]) == 14); LE = d["true_baseline"][m] / (1000.0 * d["true_neutrino_energy"][m]); w = POT[det][1] / POT[det][0]
    b = np.digitize(E, BINS) - 1; ok = (b >= 0) & (b < len(BINS) - 1); bb = b[ok]; nb = len(BINS) - 1
    tot = np.bincount(bb, minlength=nb) * w
    S = np.array([np.bincount(bb[sig[ok]], weights=np.sin(1.27 * dm * LE[ok][sig[ok]]) ** 2, minlength=nb) * w for dm in DM2])    # (ndm2, nbins): nonosc-signal events lost per unit s2
    spec[det] = (tot, S); print(det, "expected events in range", int(tot.sum()), flush=True)
def response(kind, h, a):
    if kind == "one": return 1 + a * (h - 1)
    hm, hp = h; return 1 + 0.5 * a * (hp - hm) + 0.5 * a * a * (hp + hm - 2)      # quadratic through (-1,hm),(0,1),(1,hp)
def build(excl):
    names = []; cfg = []
    for det in DET:
        for lab, kind, hh in SYST[det]:
            if (det, lab) in excl or (det.split("_")[0], lab) in excl: continue
            h = hv(hh) if kind == "one" else (hv(hh[0]), hv(hh[1])); cfg.append((det, kind, h)); names.append(f"{det}:{lab}")
    return names, cfg
def chi2_grid(excl):
    names, cfg = build(excl); ns = len(cfg); nflux = 1; nnorm = len(DET); npar = ns + nflux + nnorm; chi = np.zeros((len(S2), len(DM2))); x0 = np.zeros(npar)
    for j, dm in enumerate(DM2):
        for i, s2 in enumerate(S2):
            def resid(p):
                out = []; fl = 1 + FLUX_SIG * p[ns]
                for k, det in enumerate(DET):
                    tot, S = spec[det]; mu = tot - s2 * S[j]; f = np.ones_like(tot)
                    for q, (dd, kind, h) in enumerate(cfg):
                        if dd == det: f = f * response(kind, h, p[q])
                    mu = mu * f * fl * (1 + NORM_SIG * p[ns + 1 + k]); out.append((mu - tot) / np.sqrt(tot))
                return np.concatenate(out + [p])
            sol = least_squares(resid, x0, method="lm"); chi[i, j] = np.sum(sol.fun ** 2)
    return names, chi
def limit(chi, level=1.6):
    lim = np.full(len(DM2), np.nan)
    for j in range(len(DM2)):
        idx = np.where(chi[:, j] > level)[0]
        if len(idx):
            i = idx[0]
            if i == 0: lim[j] = S2[0]
            else: t = (level - chi[i - 1, j]) / (chi[i, j] - chi[i - 1, j]); lim[j] = 10 ** (np.log10(S2[i - 1]) + t * (np.log10(S2[i]) - np.log10(S2[i - 1])))
    return lim
SCEN = {"all": set(), "no systematics": {(d, l) for d in DET for l, _, _ in SYST[d]}, "noSBNDCS": {("sbnd_run1", "ChargeScale")}, "noICARUSCS": {("icarus_run2", "ChargeScale"), ("icarus_run4", "ChargeScale")},
        "noSBNDMattMod": {("sbnd_run1", "MattModYZ")}, "noSBNDWireModXTh": {("sbnd_run1", "WireModxThetaXW")}, "noSBNDDENT": {("sbnd_run1", "DENT")}, "noSBNDTrigger": {("sbnd_run1", "Trigger")},
        "noICARUSSCE": {("icarus_run2", "SCE"), ("icarus_run4", "SCE")}, "noICARUSWireModXTh": {("icarus_run2", "WireModxThetaXW"), ("icarus_run4", "WireModxThetaXW")},
        "noICARUSWireModYZ": {("icarus_run2", "WireModYZ"), ("icarus_run4", "WireModYZ")}, "noICARUSTrigger": {("icarus_run2", "Trigger"), ("icarus_run4", "Trigger")}}
res = {}
for k, ex in SCEN.items():
    if ONLY and k not in ONLY.split(','): continue
    names, chi = chi2_grid(ex); lim = limit(chi); res[k] = dict(chi=chi.tolist(), limit=lim.tolist()); print(k, "n_nuisance_syst", len(names), "mean log10 limit (1<dm2<10):", np.nanmean(np.log10(lim[(DM2 > 1) & (DM2 < 10)])).round(3), flush=True)
res["dm2"] = DM2.tolist(); res["s2"] = S2.tolist(); json.dump(res, open(FIG / f"s7_toy_sensitivity{TAG}.json", "w"))
base = np.array(res["all"]["limit"]); sel = (DM2 > 0.5) & (DM2 < 20)
for k in [k for k in SCEN if k in res]:
    l = np.array(res[k]["limit"]); print(f"{k:22s} mean limit improvement vs all (0.5<dm2<20): {100*np.nanmean((base[sel]-l[sel])/base[sel]):+.2f}%   median chi2 ratio in band: {np.median(np.array(res[k]['chi'])[(np.array(res['all']['chi'])>0.5)&(np.array(res['all']['chi'])<20)]/np.array(res['all']['chi'])[(np.array(res['all']['chi'])>0.5)&(np.array(res['all']['chi'])<20)]):.3f}")
