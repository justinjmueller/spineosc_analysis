#!/usr/bin/env python3
"""Per-pair core analysis -> docs/syst_study/persyst/core/<pair>.json.
For every variation: ratio (+paired bootstrap), lost/gained/migration flows, channel (1p / Mp) ratios, truth composition, event-level shifts of kept events,
feature-conditioned lost/gained fractions, standardized mean differences, and position/direction dependence."""
import sys, time
from ps_common import *
CORE = PSFIG / "core"; CORE.mkdir(exist_ok=True)
only = sys.argv[1:]
COND = ["mu_score", "mu_margin", "p_score", "p_margin", "p_mip", "mu_len", "p_len", "mu_ke", "p_ke", "E_vis", "mu_calo_csda", "p_calo_csda", "opening", "mu_primary", "p_primary",
        "drift", "vy", "vz", "mu_dirx", "mu_diry", "mu_dirz", "theta_xz", "n_p"]
SHIFT = ["E_vis", "mu_ke", "mu_len", "p_len", "p_ke", "mu_score", "p_score", "mu_calo_csda", "p_calo_csda", "mu_primary", "p_primary", "vx", "vy", "vz", "mu_endx", "mu_endy", "mu_endz", "opening"]
def L(a): return [None if not np.isfinite(x) else float(x) for x in np.asarray(a, float).ravel()]
def summarize_shift(d):
    d = d[np.isfinite(d)]
    if len(d) == 0: return None
    return dict(n=int(len(d)), mean=float(d.mean()), sem=float(d.std() / np.sqrt(len(d))), std=float(d.std()), median=float(np.median(d)), frac_nonzero=float((np.abs(d) > 1e-6).mean()),
                p16=float(np.percentile(d, 16)), p84=float(np.percentile(d, 84)))
for name in PAIRS:
    if only and not any(o in name for o in only): continue
    t0 = time.time(); det, var = name.split("__"); df = load(name); inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); both = inn & inv; lost = inn & ~inv; gain = ~inn & inv
    R = {}
    bn, bv = prep(df); n0, v0 = counts(bn, bv); r = v0 / n0; prodr = prod_hist(prod_name(name)); R["prod_match_maxdiff"] = float(np.max(np.abs(prodr - r)))
    R["integrated"] = dict(n_nom=int(inn.sum()), n_var=int(inv.sum()), lost=int(lost.sum()), gained=int(gain.sum()), both=int(both.sum()), net=float(inv.sum() / inn.sum() - 1),
                           lost_frac=float(lost.sum() / inn.sum()), gained_frac=float(gain.sum() / inn.sum()),
                           net_sigma_churn=float((inv.sum() - inn.sum()) / np.sqrt(lost.sum() + gain.sum())))
    fl = flows(df); R["flows"] = {k: L(v) for k, v in fl.items()}
    bc = boot_channels(df, nb=200); R["ratio"] = L(bc["ratio"][0]); R["err"] = L(bc["err"][0]); R["r1p"] = L(bc["ratio"][1]); R["e1p"] = L(bc["err"][1]); R["rMp"] = L(bc["ratio"][2]); R["eMp"] = L(bc["err"][2])
    R["n1p"] = L(bc["nom"][1]); R["nMp"] = L(bc["nom"][2]); R["nall"] = L(bc["nom"][0]); R["vall"] = L(bc["var"][0])
    # channel mismatch: differences per bin with bootstrap error of the difference
    d1 = bc["ratio"][1] - bc["ratio"][0]; dM = bc["ratio"][2] - bc["ratio"][0]; e1 = np.nanstd(bc["boots"][:, 1] - bc["boots"][:, 0], axis=0); eM = np.nanstd(bc["boots"][:, 2] - bc["boots"][:, 0], axis=0)
    ok = (bc["nom"][1] >= 300) & (bc["nom"][2] >= 300); sysrms = float(np.sqrt(np.mean((bc["ratio"][0][ok] - 1) ** 2)))
    def rms(d, e): return float(np.sqrt(np.mean(d[ok] ** 2))), float(np.sqrt(max(0.0, np.mean(d[ok] ** 2) - np.mean(e[ok] ** 2))))
    R["channel"] = dict(d1p=L(d1), dMp=L(dM), e_d1p=L(e1), e_dMp=L(eM), rms_sys=sysrms, rms_d1p=rms(d1, e1)[0], rms_dMp=rms(dM, eM)[0], rms_d1p_sub=rms(d1, e1)[1], rms_dMp_sub=rms(dM, eM)[1], nbin_ok=int(ok.sum()),
                        chi2_d1p=float(np.nansum(((d1 / e1)[ok]) ** 2)), chi2_dMp=float(np.nansum(((dM / eM)[ok]) ** 2)), nbin=NBIN,
                        int_1p=float(bc["var"][1].sum() / bc["nom"][1].sum() - 1), int_Mp=float(bc["var"][2].sum() / bc["nom"][2].sum() - 1), int_all=float(bc["var"][0].sum() / bc["nom"][0].sum() - 1))
    R["channel"].update(boot_int_channels(df))
    mn = col(df, "n", "reco_proton_multiplicity"); mv = col(df, "v", "reco_proton_multiplicity")
    mig = np.zeros((2, 2), int)
    for a in (0, 1):
        for b in (0, 1): mig[a, b] = int((both & ((mn > 1) == a) & ((mv > 1) == b)).sum())
    R["mult_migration"] = mig.tolist()
    # truth composition (nominal-side truth where available, else varied side)
    T = attach_truth(df, det); tt = truth_class(T["true_cc"], T["true_pdg"], T["true_category"]); R["truth_missing_gained"] = float((tt[gain] < 0).mean()) if gain.any() else 0.0
    R["truth_missing_lost"] = float((tt[lost] < 0).mean()); comp = {}
    for lab, m in (("nominal", inn), ("lost", lost), ("gained", gain), ("kept", both), ("varied", inv)):
        k = tt[m]; comp[lab] = [float((k == c).mean()) for c in range(5)] + [float((k < 0).mean())]
    R["composition"] = comp
    # class-resolved lost / gained rates (fraction of the nominal events of that class that are lost; of the varied-selected events of that class that are new), unknown-truth excluded
    rates = {}
    for c in range(5):
        mn_ = inn & (tt == c); mv_ = inv & (tt == c)
        rates[CLASSES[c]] = dict(n_nom=int(mn_.sum()), lost=float((lost & (tt == c)).sum() / max(1, mn_.sum())), gained_over_nom=float((gain & (tt == c)).sum() / max(1, mn_.sum())),
                                 net=float((mv_.sum() - mn_.sum()) / max(1, mn_.sum())))
    R["class_rates"] = rates
    sig = (tt == 0); known = tt >= 0; unk = tt < 0
    def cnt(mask_n, mask_v):
        return np.bincount(bn[(bn >= 0) & mask_n], minlength=NBIN).astype(float), np.bincount(bv[(bv >= 0) & mask_v], minlength=NBIN).astype(float)
    ns, vs = cnt(sig & inn, sig & inv); nb_, vb = cnt(~sig & known & inn, ~sig & known & inv); nu_, vu = cnt(unk & inn, unk & inv)
    # bounds: unknown-truth events (mostly gained ones) counted as signal (hi for r_sig, lo for r_bkg) or as background
    R["sig"] = dict(r_sig=L(vs / np.maximum(ns, 1)), r_bkg=L(vb / np.maximum(nb_, 1)), r_sig_hi=L((vs + vu) / np.maximum(ns + nu_ * 0, 1)), r_bkg_hi=L((vb + vu) / np.maximum(nb_, 1)),
                    bkg_frac=L(nb_ / np.maximum(ns + nb_, 1)), n_sig=L(ns), n_bkg=L(nb_), unknown_var=L(vu))
    cls_lost = [L(np.bincount(bn[(bn >= 0) & lost & (tt == c)], minlength=NBIN) / n0) for c in range(5)]
    cls_gain = [L(np.bincount(bv[(bv >= 0) & gain & (tt == c)], minlength=NBIN) / n0) for c in range(5)]
    R["class_flows"] = dict(lost=cls_lost, gained=cls_gain)
    # features
    Fn = feature_table(df, "n"); Fv = feature_table(df, "v")
    for F, s in ((Fn, "n"), (Fv, "v")):
        F["drift"] = drift_dist(det, F["vx"])
    # kept-event shifts
    sh = {}
    for k in SHIFT:
        d = (Fv[k] - Fn[k])[both]; s = summarize_shift(d); fr = None
        if k in ("E_vis", "mu_ke", "mu_len", "p_len", "p_ke"):
            den = Fn[k][both]; ok = den > 0; fr = summarize_shift((d[ok] / den[ok]))
        sh[k] = dict(abs=s, rel=fr)
    dvx, dvy, dvz = (Fv["vx"] - Fn["vx"])[both], (Fv["vy"] - Fn["vy"])[both], (Fv["vz"] - Fn["vz"])[both]; sh["vertex_3d"] = summarize_shift(np.sqrt(dvx ** 2 + dvy ** 2 + dvz ** 2))
    de = np.sqrt(((Fv["mu_endx"] - Fn["mu_endx"]) ** 2 + (Fv["mu_endy"] - Fn["mu_endy"]) ** 2 + (Fv["mu_endz"] - Fn["mu_endz"]) ** 2)[both]); sh["mu_end_3d"] = summarize_shift(de)
    dn = np.stack([col(df, "n", f"reco_leading_muon_start_dir_{a}") for a in "xyz"], 1)[both]; dvv = np.stack([col(df, "v", f"reco_leading_muon_start_dir_{a}") for a in "xyz"], 1)[both]
    ang = np.degrees(np.arccos(np.clip((dn * dvv).sum(1) / (np.linalg.norm(dn, axis=1) * np.linalg.norm(dvv, axis=1)), -1, 1))); sh["mu_dir_angle_deg"] = summarize_shift(ang)
    R["shifts"] = sh
    R["E_unchanged_frac"] = float((Fv["E_vis"][both] == Fn["E_vis"][both]).mean())
    # conditioned lost / gained fractions and SMD
    cond = {}; smds = {}
    for k in COND:
        x = Fn[k]; xv = Fv[k]; xk = x[inn & np.isfinite(x)]
        if len(xk) < 100: continue
        if k == "p_ke": edges = np.array([50., 60, 70, 80, 100, 130, 170, 250, 400, 1000, 3000])
        elif k == "mu_ke": edges = np.array([143., 170, 200, 250, 300, 400, 600, 1000, 1500, 3000])
        elif k == "n_p": edges = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 20])
        else:
            edges = np.unique(np.nanpercentile(xk, np.linspace(0, 100, 11))); edges[0] -= 1e-6; edges[-1] += 1e-6
        if len(edges) < 3: continue
        lf, lfe, den = cond_frac(x, lost, inn, edges); gf, gfe, denv = cond_frac(xv, gain, inv, edges)
        cond[k] = dict(edges=L(edges), lost=L(lf), lost_err=L(lfe), n_nom=L(den), gained=L(gf), gained_err=L(gfe), n_var=L(denv))
        smds[k] = dict(lost_vs_kept=smd(x[lost], x[both]), gained_vs_kept=smd(xv[gain], xv[both]))
    # signal-only lost fraction (nominal features) and true kinematics
    sigc = {}
    for k in ["p_margin", "mu_margin", "n_p", "E_vis", "drift", "mu_len"]:
        x = Fn[k]; xk = x[inn & sig & np.isfinite(x)]
        if len(xk) < 100: continue
        edges = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 20]) if k == "n_p" else np.unique(np.nanpercentile(xk, np.linspace(0, 100, 9))); 
        if k != "n_p": edges[0] -= 1e-6; edges[-1] += 1e-6
        lf, lfe, den = cond_frac(x, lost, inn & sig, edges); gf, gfe, denv = cond_frac(Fv[k], gain, inv & sig, edges)
        sigc[k] = dict(edges=L(edges), lost=L(lf), lost_err=L(lfe), n=L(den), gained=L(gf), gained_err=L(gfe), n_var=L(denv))
    for k, key, edges in (("t_p_ke", "true_leading_proton_ke", np.array([0, 30, 50, 60, 80, 100, 150, 250, 400, 1000, 5000.])), ("t_mu_ke", "true_leading_muon_ke", np.array([0, 100, 143, 170, 200, 300, 500, 800, 1500, 5000.])),
                          ("t_Enu", "true_neutrino_energy", np.array([0, 0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0, 3.0, 20.]))):
        x = T[key]
        if not np.isfinite(x[inn]).any(): continue
        lf, lfe, den = cond_frac(x, lost, inn, edges); gf, gfe, denv = cond_frac(x, gain, inv, edges)
        sigc[k] = dict(edges=L(edges), lost=L(lf), lost_err=L(lfe), n=L(den), gained=L(gf), gained_err=L(gfe), n_var=L(denv))
    R["cond_extra"] = sigc
    # background-only (non-signal, known truth) lost/gained fractions and SMDs
    bkg = known & ~sig; sigcb = {}; smdb = {}
    for k in ["p_margin", "mu_margin", "n_p", "E_vis", "mu_len", "p_len", "mu_calo_csda", "p_calo_csda", "mu_primary", "opening", "drift"]:
        x = Fn[k]; xk = x[inn & bkg & np.isfinite(x)]
        if len(xk) < 100: continue
        edges = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 20]) if k == "n_p" else np.unique(np.nanpercentile(xk, np.linspace(0, 100, 7)))
        if k != "n_p": edges[0] -= 1e-6; edges[-1] += 1e-6
        lf, lfe, den = cond_frac(x, lost, inn & bkg, edges); gf, gfe, denv = cond_frac(Fv[k], gain, inv & bkg, edges)
        sigcb[k] = dict(edges=L(edges), lost=L(lf), lost_err=L(lfe), n=L(den), gained=L(gf), gained_err=L(gfe), n_var=L(denv))
        smdb[k] = dict(lost_vs_kept=smd(x[lost & bkg], x[both & bkg]), gained_vs_kept=smd(Fv[k][gain & bkg], Fv[k][both & bkg]))
        if k in ("n_p",): pass
    R["cond_bkg"] = sigcb; R["smd_bkg"] = smdb; R["n_bkg_lost"] = int((lost & bkg).sum()); R["n_bkg_gained"] = int((gain & bkg).sum())
    R["cond"] = cond; R["smd"] = smds
    # threshold-boundary fractions: lost/gained/kept with mu_ke< 0.2 GeV or p_ke<0.08 GeV
    thr = {}
    for lab, m, F in (("lost", lost, Fn), ("kept_n", both, Fn), ("gained", gain, Fv), ("kept_v", both, Fv)):
        thr[lab] = dict(p_ke_lt_80MeV=float(np.nanmean(F["p_ke"][m] < 80)), mu_ke_lt_200MeV=float(np.nanmean(F["mu_ke"][m] < 200)), mu_margin_lt_0=float(np.nanmean(F["mu_margin"][m] < 0)),
                        p_margin_lt_0=float(np.nanmean(F["p_margin"][m] < 0)), mu_score_lt_0p4=float(np.nanmean(F["mu_score"][m] < 0.4)), p_score_lt_0p4=float(np.nanmean(F["p_score"][m] < 0.4)))
    R["thr"] = thr
    # position dependence of the shifts of kept events (SCE-like); per position bin mean fractional E shift, mean mu_len shift, vertex displacement, lost/gained fractions
    pos = {}
    for k, edges in (("drift", np.linspace(0, 150 if not det.startswith("sbnd") else 200, 11)), ("vy", np.linspace(-200, 200, 11) if det.startswith("sbnd") else np.linspace(-180, 140, 11)),
                     ("vz", np.linspace(0, 500, 11) if det.startswith("sbnd") else np.linspace(-900, 900, 13)), ("theta_xz", np.linspace(0, 90, 10))):
        xk = Fn[k]; b = np.digitize(xk, edges) - 1; ok = both & (b >= 0) & (b < len(edges) - 1) & np.isfinite(xk); rows = []
        for i in range(len(edges) - 1):
            m = ok & (b == i); mm = m & (Fn["E_vis"] > 0)
            dE = ((Fv["E_vis"] - Fn["E_vis"]) / Fn["E_vis"])[mm]; dL = ((Fv["mu_len"] - Fn["mu_len"]) / Fn["mu_len"])[m & (Fn["mu_len"] > 0)]; dv = np.sqrt(dvx ** 2 + dvy ** 2 + dvz ** 2)
            vv = np.sqrt((Fv["vx"] - Fn["vx"]) ** 2 + (Fv["vy"] - Fn["vy"]) ** 2 + (Fv["vz"] - Fn["vz"]) ** 2)[m]
            rows.append(dict(n=int(m.sum()), dE=float(np.nanmean(dE)) if len(dE) else None, dE_sem=float(np.nanstd(dE) / np.sqrt(max(1, len(dE)))) if len(dE) else None,
                             dL=float(np.nanmean(dL)) if len(dL) else None, dL_sem=float(np.nanstd(dL) / np.sqrt(max(1, len(dL)))) if len(dL) else None,
                             dcalo=float(np.nanmean((Fv["mu_calo_csda"] - Fn["mu_calo_csda"])[m])) if m.any() else None, dscore=float(np.nanmean((Fv["mu_score"] - Fn["mu_score"])[m])) if m.any() else None,
                             dvtx=float(np.nanmedian(vv)) if len(vv) else None, dvx_signed=float(np.nanmean((Fv["vx"] - Fn["vx"])[m])) if m.any() else None))
        bvv = np.digitize(Fv[k], edges) - 1; nn_ = np.array([int((inn & (b == i)).sum()) for i in range(len(edges) - 1)], float); nv_ = np.array([int((inv & (bvv == i)).sum()) for i in range(len(edges) - 1)], float)
        pos[k] = dict(edges=L(edges), rows=rows, r_pos=L(np.where(nn_ > 50, nv_ / np.maximum(nn_, 1) - 1, np.nan)), r_pos_err=L(np.sqrt(nv_ + nn_ - 0) / np.maximum(nn_, 1) * 0 + (np.sqrt(np.maximum(nn_ * 0.2, 1))) * 0), n_nom=L(nn_), n_var=L(nv_))
    R["pos"] = pos
    # side (x<0 / x>0 = SBND east/west TPC or ICARUS east/west cryostat): churn and shift
    side = {}; vxn_, vxv_ = Fn["vx"], Fv["vx"]
    for lab, fn_, fv_ in (("x<0", vxn_ < 0, vxv_ < 0), ("x>0", vxn_ > 0, vxv_ > 0)):
        N = int((inn & fn_).sum()); dEr = ((Fv["E_vis"] - Fn["E_vis"]) / Fn["E_vis"])[both & fn_]
        side[lab] = dict(N=N, lost=float((lost & fn_).sum() / N), gained=float((gain & fv_).sum() / N), net=float((inv & fv_).sum() / N - 1), dE_mean=float(np.nanmean(dEr)), dE_absmed=float(np.nanmedian(np.abs(dEr))),
                         E_unchanged=float(np.mean(dEr == 0)), dmu_len_mean=float(np.nanmean((Fv["mu_len"] - Fn["mu_len"])[both & fn_])))
    R["side"] = side
    # per-bin mean fractional E shift of events kept on both sides (by nominal bin) and the migration it predicts (flux through the bin edges)
    bk = bn[both]; dEk = ((Fv["E_vis"] - Fn["E_vis"]) / Fn["E_vis"])[both]; mu_, med_, sd_ = [], [], []
    for i in range(NBIN):
        x = dEk[bk == i]; x = x[np.isfinite(x)]; mu_.append(float(x.mean()) if len(x) else np.nan); med_.append(float(np.median(x)) if len(x) else np.nan); sd_.append(float(x.std()) if len(x) > 1 else np.nan)
    R["Eshift"] = dict(mean=L(mu_), median=L(med_), std=L(sd_))
    (CORE / f"{name}.json").write_text(json.dumps(R)); print(name, f"{time.time()-t0:.0f}s", "prodmatch", R["prod_match_maxdiff"], "net", round(R["integrated"]["net"] * 100, 2), flush=True)
    del df, Fn, Fv
