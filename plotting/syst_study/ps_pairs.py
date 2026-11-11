#!/usr/bin/env python3
"""Cross-variation comparisons: (1) SCE 0x vs 2x odd/even decomposition and event-level linearity; (2) SBND WireModYZ vs MattModYZ (and vs ChargeScale / WireModxThetaXW) event-level overlap and shift correlation.
Output: docs/syst_study/persyst/pairs_compare.json"""
from ps_common import *
from scipy import stats
KEYK = KEYS + ["Ei"]
def shifts(df):
    """event-keyed shifts for events selected on both sides"""
    both = (df.in_n & df.in_v).to_numpy(); d = df.loc[both, KEYK].copy()
    def rel(c, den_min=0):
        a, b = col(df, "n", c)[both], col(df, "v", c)[both]; return np.where(a > den_min, (b - a) / np.where(a > den_min, a, 1), np.nan)
    d["dE"] = rel("reco_visible_energy"); d["dmu_len"] = rel("reco_leading_muon_length"); d["dp_len"] = rel("reco_leading_proton_length")
    for nm, c in (("mu_score", "reco_leading_muon_muon_softmax"), ("p_score", "reco_leading_proton_proton_softmax"), ("mu_calo", "reco_leading_muon_calo_csda_diff"), ("p_calo", "reco_leading_proton_calo_csda_diff"),
                  ("mu_primary", "reco_leading_muon_primary_softmax"), ("p_primary", "reco_leading_proton_primary_softmax"), ("vx", "reco_vertex_x"), ("vy", "reco_vertex_y"), ("vz", "reco_vertex_z"), ("mu_endx", "reco_leading_muon_end_x")):
        d["d" + nm] = col(df, "v", c)[both] - col(df, "n", c)[both]
    return d
def state(df):
    """per event: nominal-selected flag, varied-selected flag"""
    return df[KEYK + ["in_n", "in_v"]].copy()
def corr_table(a, b, names, la, lb):
    m = a.merge(b, on=KEYK, suffixes=("_a", "_b")); out = {"n_common_kept": int(len(m))}
    for n in names:
        x = m[n + "_a"].to_numpy(float); y = m[n + "_b"].to_numpy(float); ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 50: continue
        r = stats.pearsonr(x[ok], y[ok])[0]; rs = stats.spearmanr(x[ok], y[ok])[0]
        # robust slope: regression of b on a on the central 98%
        lo, hi = np.percentile(x[ok], [1, 99]); c = ok & (x >= lo) & (x <= hi); sl = np.polyfit(x[c], y[c], 1)[0] if c.sum() > 50 else np.nan
        out[n] = dict(pearson=float(r), spearman=float(rs), slope_b_on_a=float(sl), mean_a=float(np.nanmean(x[ok])), mean_b=float(np.nanmean(y[ok])), std_a=float(np.nanstd(x[ok])), std_b=float(np.nanstd(y[ok])))
    return out
def contingency(sa, sb):
    """among events selected in nominal and present in both pools: lost-in-A vs lost-in-B ; and for gained"""
    m = sa.merge(sb, on=KEYK, suffixes=("_a", "_b"), how="inner"); res = {}
    nom = m.in_n_a.to_numpy() & m.in_n_b.to_numpy()
    la = nom & ~m.in_v_a.to_numpy(); lb = nom & ~m.in_v_b.to_numpy()
    res["nom_common"] = int(nom.sum()); res["lost_a"] = int(la.sum()); res["lost_b"] = int(lb.sum()); res["lost_both"] = int((la & lb).sum())
    exp = la.sum() * lb.sum() / max(1, nom.sum()); res["lost_both_expected_indep"] = float(exp); res["lost_overlap_ratio"] = float((la & lb).sum() / max(1e-9, exp))
    # fraction of the lost-in-A that is also lost in B (vs the base loss rate in B)
    res["P_lostB_given_lostA"] = float((la & lb).sum() / max(1, la.sum())); res["P_lostB_base"] = float(lb.sum() / max(1, nom.sum()))
    mo = sa.merge(sb, on=KEYK, suffixes=("_a", "_b"), how="outer")
    ga = (~mo.in_n_a.fillna(True).astype(bool) & mo.in_v_a.fillna(False).astype(bool)).to_numpy(); gb = (~mo.in_n_b.fillna(True).astype(bool) & mo.in_v_b.fillna(False).astype(bool)).to_numpy()
    res["gained_a"] = int(ga.sum()); res["gained_b"] = int(gb.sum()); res["gained_both"] = int((ga & gb).sum())
    res["P_gainedB_given_gainedA_lowerbound"] = float((ga & gb).sum() / max(1, ga.sum())); res["P_gainedA_given_gainedB_lowerbound"] = float((ga & gb).sum() / max(1, gb.sum()))
    return res
out = {}
# ---- (1) SCE odd/even
for det in ("sbnd_run1", "icarus_run2", "icarus_run4"):
    c0 = json.loads((PSFIG / "core" / f"{det}__SCE0x.json").read_text()); c2 = json.loads((PSFIG / "core" / f"{det}__SCE2x.json").read_text())
    r0 = np.array(c0["ratio"]) - 1; r2 = np.array(c2["ratio"]) - 1; e0 = np.array(c0["err"]); e2 = np.array(c2["err"])
    odd = 0.5 * (r2 - r0); even = 0.5 * (r2 + r0); e = 0.5 * np.sqrt(e0 ** 2 + e2 ** 2)   # independent errors assumed (common nominal makes them positively correlated -> this is conservative for odd and optimistic for even)
    w = 1 / np.array([e0, e2]).max(0) ** 2; sl = float(np.sum(w * r2 * (-r0)) / np.sum(w * r0 ** 2))   # fit r2 = s * (-r0)
    res = dict(odd=odd.tolist(), even=even.tolist(), err=e.tolist(), chi2_even=float(np.sum((even / e) ** 2)), chi2_odd=float(np.sum((odd / e) ** 2)), nbin=len(e), slope_r2_vs_minus_r0=sl,
               int0=c0["integrated"]["net"], int2=c2["integrated"]["net"], int_odd=0.5 * (c2["integrated"]["net"] - c0["integrated"]["net"]), int_even=0.5 * (c2["integrated"]["net"] + c0["integrated"]["net"]),
               lost0=c0["integrated"]["lost_frac"], gain0=c0["integrated"]["gained_frac"], lost2=c2["integrated"]["lost_frac"], gain2=c2["integrated"]["gained_frac"])
    # class level
    res["class_net"] = {k: dict(x0=c0["class_rates"][k]["net"], x2=c2["class_rates"][k]["net"], n=c0["class_rates"][k]["n_nom"]) for k in c0["class_rates"]}
    d0 = load(f"{det}__SCE0x"); d2 = load(f"{det}__SCE2x"); s0, s2 = shifts(d0), shifts(d2)
    res["event_corr"] = corr_table(s0, s2, ["dE", "dmu_len", "dp_len", "dmu_score", "dp_score", "dmu_calo", "dp_calo", "dvx", "dvy", "dvz", "dmu_endx", "dmu_primary", "dp_primary"], "0x", "2x")
    res["states"] = contingency(state(d0), state(d2))
    # odd/even of the lost/gained rates as function of nominal features (position): fraction lost vs drift distance for 0x and 2x
    res["lost_vs_drift"] = dict(edges=c0["cond"]["drift"]["edges"], lost0=c0["cond"]["drift"]["lost"], lost2=c2["cond"]["drift"]["lost"], gain0=c0["cond"]["drift"]["gained"], gain2=c2["cond"]["drift"]["gained"])
    out[f"SCE_{det}"] = res; print("SCE", det, "slope", round(sl, 2), "chi2 odd/even", round(res["chi2_odd"]), round(res["chi2_even"]), "int", res["int0"], res["int2"], "corr dE", res["event_corr"].get("dE", {}).get("pearson"), "dvx", res["event_corr"].get("dvx", {}).get("pearson"), flush=True)
    del d0, d2
# ---- (2) SBND YZ comparison
dW = load("sbnd_run1__WireModYZ"); dM = load("sbnd_run1__MattModYZ"); dX = load("sbnd_run1__WireModxThetaXW"); dC = load("sbnd_run1__ChargeScale"); dC15 = load("sbnd_run1__ChargeScale15"); dG = load("sbnd_run1__GainHigh")
SH = {k: shifts(d) for k, d in (("WireModYZ", dW), ("MattModYZ", dM), ("WireModxThetaXW", dX), ("ChargeScale", dC), ("ChargeScale15", dC15), ("GainHigh", dG))}
ST = {k: state(d) for k, d in (("WireModYZ", dW), ("MattModYZ", dM), ("WireModxThetaXW", dX), ("ChargeScale", dC))}
yz = dict(); names = ["dE", "dmu_len", "dp_len", "dmu_score", "dp_score", "dmu_calo", "dp_calo", "dmu_primary", "dp_primary", "dvx", "dvy", "dvz"]
for a, b in (("WireModYZ", "MattModYZ"), ("WireModYZ", "WireModxThetaXW"), ("WireModYZ", "ChargeScale"), ("MattModYZ", "ChargeScale")):
    yz[f"{a}_vs_{b}"] = dict(corr=corr_table(SH[a], SH[b], names, a, b), states=contingency(ST[a], ST[b]))
# mean/rms shifts per variation (kept events) for charge-like observables
summ = {}
for k, s in SH.items():
    summ[k] = {n: dict(mean=float(np.nanmean(s[n])), std=float(np.nanstd(s[n])), sem=float(np.nanstd(s[n]) / np.sqrt(s[n].notna().sum()))) for n in ["dE", "dmu_len", "dp_len", "dmu_calo", "dp_calo", "dmu_score", "dp_score", "dmu_primary", "dp_primary"]}
    summ[k]["frac_E_unchanged"] = float((s["dE"] == 0).mean()); summ[k]["n"] = int(len(s))
yz["shift_summary"] = summ
# spatial structure of the YZ effects (kept events): mean dE and (lost-gained) fraction vs vertex (z,y) 2D coarse grid, WireModYZ vs MattModYZ, and E_vis-weighted charge proxy (calo) shift
def grid(df, edges_y, edges_z):
    inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); both = inn & inv
    vy, vz = col(df, "n", "reco_vertex_y"), col(df, "n", "reco_vertex_z"); vyv, vzv = col(df, "v", "reco_vertex_y"), col(df, "v", "reco_vertex_z")
    by = np.digitize(vy, edges_y) - 1; bz = np.digitize(vz, edges_z) - 1; byv = np.digitize(vyv, edges_y) - 1; bzv = np.digitize(vzv, edges_z) - 1
    nb = (len(edges_y) - 1) * (len(edges_z) - 1); idx = np.where((by >= 0) & (by < len(edges_y) - 1) & (bz >= 0) & (bz < len(edges_z) - 1), by * (len(edges_z) - 1) + bz, -1)
    idxv = np.where((byv >= 0) & (byv < len(edges_y) - 1) & (bzv >= 0) & (bzv < len(edges_z) - 1), byv * (len(edges_z) - 1) + bzv, -1)
    nn = np.bincount(idx[inn & (idx >= 0)], minlength=nb).astype(float); nv = np.bincount(idxv[inv & (idxv >= 0)], minlength=nb).astype(float)
    return (nv / np.maximum(nn, 1) - 1).reshape(len(edges_y) - 1, -1), nn.reshape(len(edges_y) - 1, -1)
ey = np.linspace(-200, 200, 9); ez = np.linspace(0, 500, 11)
gW, nW = grid(dW, ey, ez); gM, nM = grid(dM, ey, ez)
mask = (nW > 300) & (nM > 300)
yz["grid"] = dict(edges_y=ey.tolist(), edges_z=ez.tolist(), WireModYZ=np.where(nW > 300, gW, np.nan).tolist(), MattModYZ=np.where(nM > 300, gM, np.nan).tolist(), n=nW.tolist(),
                  corr_cells=float(np.corrcoef(gW[mask], gM[mask])[0, 1]), ratio_rms=float(np.sqrt(np.mean(gW[mask] ** 2)) / np.sqrt(np.mean(gM[mask] ** 2))))
out["YZ"] = yz
json.dump(out, open(PSFIG / "pairs_compare.json", "w"))
print("WireModYZ vs MattModYZ shift corr:", {k: round(v["pearson"], 2) for k, v in yz["WireModYZ_vs_MattModYZ"]["corr"].items() if isinstance(v, dict)})
print(yz["WireModYZ_vs_MattModYZ"]["states"])
print({k: (round(v["dmu_calo"]["mean"], 4), round(v["dp_calo"]["mean"], 4), round(v["dE"]["mean"], 4)) for k, v in summ.items()})

# ---- (3) three-way Hamming distances among nominal, 0x, 2x selections (events present in both pair tables; absent in a table => not selected there)
ham = {}
for det in ("sbnd_run1", "icarus_run2", "icarus_run4"):
    a = state(load(f"{det}__SCE0x")); b = state(load(f"{det}__SCE2x")); m = a.merge(b, on=KEYK, how="outer", suffixes=("_0", "_2"))
    n = ((m.in_n_0.fillna(False).astype(bool)) | (m.in_n_2.fillna(False).astype(bool))).to_numpy(); s0 = m.in_v_0.fillna(False).astype(bool).to_numpy(); s2 = m.in_v_2.fillna(False).astype(bool).to_numpy()
    # restrict to events whose key is present in both tables OR selected in nominal-only-of-one... (nominal-selected events are in both pools except pool differences)
    both_t = (m.in_n_0.notna() & m.in_n_2.notna()).to_numpy(); sel = both_t | n
    # pure comparison on events in both tables
    k = ~(n & ~both_t); N = int(n[k].sum())   # drop nominal-selected events absent from one table (pool mismatch); not-selected-anywhere events are implicit zeros
    ham[det] = dict(n_nominal=N, d_n0=int((n[k] ^ s0[k]).sum()), d_n2=int((n[k] ^ s2[k]).sum()), d_02=int((s0[k] ^ s2[k]).sum()), n_0=int(s0[k].sum()), n_2=int(s2[k].sum()),
                    d_n0_over_N=float((n[k] ^ s0[k]).sum() / N), d_n2_over_N=float((n[k] ^ s2[k]).sum() / N), d_02_over_N=float((s0[k] ^ s2[k]).sum() / N))
    print("hamming", det, ham[det])
out["hamming"] = ham; json.dump(out, open(PSFIG / "pairs_compare.json", "w"))

# ---- (4) joint paired bootstrap of the SCE 0x/2x odd and even parts (shared nominal events -> correlated errors), restricted to the pool common to nominal, 0x and 2x
from common import configs, ekeys
jb = {}
for det in ("sbnd_run1", "icarus_run2", "icarus_run4"):
    cf = {c[0]: c for c in configs()[det]}; lab0, lab2 = cf["SCE0x"], cf["SCE2x"]
    ek = ekeys(*lab0[3]).drop_duplicates(); e0 = ekeys(lab0[4], lab0[2]).drop_duplicates(); e2 = ekeys(lab2[4], lab2[2]).drop_duplicates()
    pool = ek.merge(e0, on=KEYK).merge(e2, on=KEYK)[KEYK]; pool["inpool"] = True
    a = load(f"{det}__SCE0x"); b = load(f"{det}__SCE2x")
    fa = a[KEYK].copy(); fa["bn"] = np.where(a.in_n, binidx(col(a, "n", "reco_visible_energy")), -1); fa["b0"] = np.where(a.in_v, binidx(col(a, "v", "reco_visible_energy")), -1)
    fb = b[KEYK].copy(); fb["bn_b"] = np.where(b.in_n, binidx(col(b, "n", "reco_visible_energy")), -1); fb["b2"] = np.where(b.in_v, binidx(col(b, "v", "reco_visible_energy")), -1)
    m = fa.merge(fb, on=KEYK, how="outer").merge(pool, on=KEYK, how="inner")
    bnn = np.where(m.bn.notna(), m.bn, m.bn_b).astype(int); b0 = m.b0.fillna(-1).astype(int).to_numpy(); b2 = m.b2.fillna(-1).astype(int).to_numpy()
    rng = np.random.default_rng(11); NB = 300; O = np.empty((NB, NBIN)); Ev = np.empty((NB, NBIN)); R0 = np.empty((NB, NBIN)); R2 = np.empty((NB, NBIN))
    def rr(w):
        n = np.bincount(bnn[bnn >= 0], weights=w[bnn >= 0], minlength=NBIN); v0 = np.bincount(b0[b0 >= 0], weights=w[b0 >= 0], minlength=NBIN); v2 = np.bincount(b2[b2 >= 0], weights=w[b2 >= 0], minlength=NBIN); return v0 / n - 1, v2 / n - 1
    for i in range(NB):
        w = rng.poisson(1.0, len(bnn)).astype(float); x0, x2 = rr(w); R0[i], R2[i] = x0, x2; O[i] = 0.5 * (x2 - x0); Ev[i] = 0.5 * (x2 + x0)
    x0, x2 = rr(np.ones(len(bnn))); odd = 0.5 * (x2 - x0); even = 0.5 * (x2 + x0)
    jb[det] = dict(odd=odd.tolist(), even=even.tolist(), e_odd=O.std(0).tolist(), e_even=Ev.std(0).tolist(), chi2_odd=float(np.sum((odd / O.std(0)) ** 2)), chi2_even=float(np.sum((even / Ev.std(0)) ** 2)), nbin=NBIN,
                   e0=R0.std(0).tolist(), e2=R2.std(0).tolist(), corr_r0_r2=float(np.mean([np.corrcoef(R0[:, j], R2[:, j])[0, 1] for j in range(NBIN)])), n_common_pool_rows=int(len(m)),
                   r0_common=x0.tolist(), r2_common=x2.tolist(), int0=float(np.sum(np.bincount(b0[b0 >= 0], minlength=NBIN)) / np.sum(np.bincount(bnn[bnn >= 0], minlength=NBIN)) - 1),
                   int2=float(np.sum(np.bincount(b2[b2 >= 0], minlength=NBIN)) / np.sum(np.bincount(bnn[bnn >= 0], minlength=NBIN)) - 1))
    print("joint SCE", det, "chi2 odd/even", round(jb[det]["chi2_odd"]), round(jb[det]["chi2_even"]), "avg corr(r0,r2) boot", round(jb[det]["corr_r0_r2"], 2), "int0/int2", round(jb[det]["int0"], 4), round(jb[det]["int2"], 4), flush=True)
out["joint_boot"] = jb; json.dump(out, open(PSFIG / "pairs_compare.json", "w"))
