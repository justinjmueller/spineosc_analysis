#!/usr/bin/env python3
"""S9: figures for the charge-scale / energy-trend sections of the report."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ana import *
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); C = {"sbnd": "#D55E00", "ic2": "#0072B2", "ic4": "#56B4E9"}
E6 = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0]); ctr6 = 0.5 * (E6[1:] + E6[:-1])
# ---------- Fig A: composition of the nominal selected sample versus energy
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
for lab, name, c in (("SBND Run 1", "sbnd_run1__ChargeScale15", C["sbnd"]), ("ICARUS Run 2", "icarus_run2__ChargeScale", C["ic2"]), ("ICARUS Run 4", "icarus_run4__ChargeScale", C["ic4"])):
    df = load(name); m = df.in_n.to_numpy(); E = df.n_reco_visible_energy.to_numpy(float)[m]; mult = df.n_reco_proton_multiplicity.to_numpy(float)[m]; muke = df.n_reco_leading_muon_ke.to_numpy(float)[m] / 1000
    fMp = []; mm = []; fmu = []
    for a, b in zip(E6[:-1], E6[1:]):
        s = (E >= a) & (E < b); fMp.append((mult[s] > 1).mean()); mm.append(mult[s].mean()); fmu.append((muke[s] / E[s]).mean())
    axs[0].plot(ctr6, np.array(fMp) * 100, "-o", color=c, label=lab); axs[1].plot(ctr6, mm, "-o", color=c, label=lab); axs[2].plot(ctr6, np.array(fmu) * 100, "-o", color=c, label=lab)
for ax, yl in zip(axs, ("Multi-proton events [%]", "Mean reco proton multiplicity", "Muon KE / $E_{vis}$ [%]")): ax.set_xlabel("Reco visible energy [GeV]"); ax.set_ylabel(yl)
axs[0].legend(fontsize=11, loc="lower right"); mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[2], "Nominal selected sample", alignment="right"); fig.tight_layout(); fig.savefig(FIG / "s9_composition_vs_energy.png", dpi=100); plt.close(fig)
# ---------- Fig B: hazard of loss/gain versus proton multiplicity
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4)); mb = [1, 2, 3, 4, 8]; x = np.arange(4)
for ax, (grp, title) in zip(axs, (("lost", "Lost (selected in nominal only)"), ("gained", "Gained (selected in the variation only)"))):
    for k, (name, lab, c, ls) in enumerate((("sbnd_run1__ChargeScale15", "SBND ChargeScale15", C["sbnd"], "-"), ("sbnd_run1__ChargeScale", "SBND ChargeScale", C["sbnd"], "--"), ("icarus_run2__ChargeScale", "ICARUS Run2 ChargeScale", C["ic2"], "-"),
                                          ("icarus_run4__ChargeScale", "ICARUS Run4 ChargeScale", C["ic4"], "-"), ("sbnd_run1__MattModYZ", "SBND MattModYZ", "#009E73", "-"), ("sbnd_run1__DENT", "SBND DENT", "#CC79A7", "-"))):
        df = load(name); inn = df.in_n.to_numpy(); inv = df.in_v.to_numpy(); side, base, ev = ("n", inn, inn & ~inv) if grp == "lost" else ("v", inv, ~inn & inv)
        mult = df[f"{side}_reco_proton_multiplicity"].to_numpy(float); p = [ev[base & (mult >= a) & (mult < b)].mean() * 100 for a, b in zip(mb[:-1], mb[1:])]; ax.plot(x + 1, p, ls, marker="o", color=c, label=lab)
    ax.set_xticks(x + 1); ax.set_xticklabels(["1", "2", "3", "4-7"]); ax.set_xlabel("Reco proton multiplicity"); ax.set_ylabel("Fraction of selected events [%]"); ax.text(0.04, 0.93, title, transform=ax.transAxes, fontweight="bold", va="top")
axs[0].legend(fontsize=10, loc="center left"); mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[1], "Selection flips vs multiplicity", alignment="right"); fig.tight_layout(); fig.savefig(FIG / "s9_hazard_vs_multiplicity.png", dpi=100); plt.close(fig)
# ---------- Fig C: calorimetric response change versus KE
def calo_prof(name, part, ked):
    df = load(name); b = (df.in_n & df.in_v).to_numpy(); L = df[f"n_reco_leading_{part}_length"].to_numpy(float); Lv = df[f"v_reco_leading_{part}_length"].to_numpy(float); same = b & (np.abs(L - Lv) < 1e-6)
    ke = df[f"n_reco_leading_{part}_ke"].to_numpy(float); d = df[f"v_reco_leading_{part}_calo_csda_diff"].to_numpy(float) - df[f"n_reco_leading_{part}_calo_csda_diff"].to_numpy(float); o = []
    for lo, hi in zip(ked[:-1], ked[1:]):
        m = same & (ke >= lo) & (ke < hi) & np.isfinite(d); q = np.percentile(d[m], [16, 50, 84]); o.append(q)
    return np.array(o)
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
for ax, (part, ked) in zip(axs, (("muon", np.array([143, 200, 300, 400, 550, 750, 1000, 2000])), ("proton", np.array([50, 70, 100, 140, 200, 300, 450, 2000])))):
    cc = 0.5 * (ked[1:] + ked[:-1])
    for name, lab, c, ls in (("sbnd_run1__ChargeScale", "SBND ChargeScale", C["sbnd"], "--"), ("sbnd_run1__ChargeScale15", "SBND ChargeScale15", C["sbnd"], "-"), ("icarus_run2__ChargeScale", "ICARUS Run2 ChargeScale", C["ic2"], "-"), ("icarus_run4__ChargeScale", "ICARUS Run4 ChargeScale", C["ic4"], "-"), ("sbnd_run1__GainLow", "SBND Gain $-$5%", "0.4", ":")):
        o = calo_prof(name, part, ked) * 100; ax.plot(cc, o[:, 1], ls, marker="o", ms=4, color=c, label=lab)
        if name in ("sbnd_run1__ChargeScale15", "icarus_run2__ChargeScale"): ax.fill_between(cc, o[:, 0], o[:, 2], color=c, alpha=0.15, lw=0)
    ax.axhline(0, color="0.5", lw=0.8); ax.set_xscale("log"); ax.set_xlabel(f"Leading {part} KE (range) [MeV]"); ax.set_ylabel("Change in calorimetric/range KE [%]"); ax.text(0.04, 0.93, part.capitalize(), transform=ax.transAxes, fontweight="bold", va="top")
axs[0].legend(fontsize=10, loc="lower right"); mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[1], "Events with unchanged track length", alignment="right"); fig.tight_layout(); fig.savefig(FIG / "s9_calo_response.png", dpi=100); plt.close(fig)
print("figs done")
