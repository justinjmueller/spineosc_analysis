#!/usr/bin/env python3
"""Per-systematic deep dive (SBND, final selection: fiducial, veto, flash, crosser == 0, |x_end| >= 4 cm).
For each active detector variation, on the events common to the nominal and varied samples: split the change in the selected count in each reco-energy bin into
 (lost) selected in the nominal only, (gained) selected in the variation only, (migration) selected in both but moved between energy bins;
then show where in the PID scores / track length the lost and gained events sit.  The nominal of DENT is the detvar CV, the others use sbnd_large (nosyst).
Output: docs/syst_deepdive/*.  Run: python plotting/dent_study/syst_deepdive.py"""
import json, pathlib, sys
import numpy as np, pandas as pd, uproot
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).parent)); import cathode_scan as cs
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike")
OUT = cs.REPO / "docs" / "syst_deepdive"; OUT.mkdir(parents=True, exist_ok=True); K = cs.K
F = ["reco_leading_muon_muon_softmax", "reco_leading_muon_pion_softmax", "reco_leading_muon_proton_softmax", "reco_leading_muon_length", "reco_leading_muon_ke", "reco_leading_proton_ke", "reco_leading_proton_proton_softmax"]
COLS = K + ["reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz", "reco_cathode_crosser", "reco_leading_muon_end_x", "reco_proton_multiplicity", "event_largest_flash_pe", "true_cc", "true_pdg", "true_category"] + F
VARS = [("Charge scale", "sbnd_detvar_smear15", (cs.NOS, "sbnd_large")), ("WireMod X-$\\theta_{xW}$", "sbnd_detvar_wiremod_xthetaxw", (cs.NOS, "sbnd_large")),
        ("Non-uniformity", "sbnd_detvar_mattmod_yz", (cs.NOS, "sbnd_large")), ("DENT", "sbnd_detvar_newdent", (cs.DET, "sbnd_detvar_newcv"))]
EDG = {"1p": np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0]), "Mp": np.array([0.4, 0.6, 0.8, 1.0, 1.3, 2.0])}
FEAT = {"reco_leading_muon_muon_softmax": ("Muon score of the muon", np.linspace(0, 1, 11)), "reco_leading_muon_pion_softmax": ("Pion score of the muon", np.array([0, .02, .05, .1, .2, .3, .5, 1.0])),
        "reco_leading_muon_length": ("Muon length [cm]", np.array([0, 25, 50, 75, 100, 150, 200, 300, 600]))}
def load(fn, tree):
    d = uproot.open(fn)[f"events/{tree}/selected"].arrays(COLS, library="pd").drop_duplicates(K, keep=False).set_index(K)
    d = d[(d.reco_fiducialize_cathode == 1) & (d.reco_veto_sbnd_highy_highz == 1) & (0.647 * d.event_largest_flash_pe > 2000) & (d.reco_cathode_crosser == 0) & (d.reco_leading_muon_end_x.abs() >= 4) & (d.reco_proton_multiplicity >= 1)]
    d = d.copy(); d["sig"] = (d.true_cc > 0) & (d.true_pdg.abs() == 14) & d.true_category.isin([0, 2]); return d
def binidx(E, edges): return np.digitize(E, edges) - 1
if __name__ == "__main__":
    res = {}; fig, axs = plt.subplots(4, 4, figsize=(19, 15))
    for ri, (label, tree, (nfn, ntree)) in enumerate(VARS):
        v = load(cs.DET, tree); n = load(nfn, ntree); com = cs.keys(cs.DET, tree).intersection(cs.keys(nfn, ntree)); v = v[v.index.isin(com)]; n = n[n.index.isin(com)]
        both = n.index.intersection(v.index); lost = n.index.difference(v.index); gain = v.index.difference(n.index)
        res[label] = dict(n_nom=len(n), n_var=len(v), lost=len(lost), gained=len(gain), both=len(both), ratio=len(v) / len(n) - 1,
                          lost_frac=len(lost) / len(n), gained_frac=len(gain) / len(n),
                          lost_signal_share=float(n.loc[lost].sig.mean()), gained_signal_share=float(v.loc[gain].sig.mean()), nominal_signal_share=float(n.sig.mean()), both_signal_share=float(n.loc[both].sig.mean()))
        res[label]["bins"] = {}
        for ch, pm in (("1p", lambda d: d.reco_proton_multiplicity == 1), ("Mp", lambda d: d.reco_proton_multiplicity > 1)):
            e = EDG[ch]; nb = len(e) - 1; nn = n[pm(n)]; vv = v[pm(v)]
            Nn = np.bincount(np.clip(binidx(nn.reco_visible_energy.to_numpy(), e), -1, nb)[(binidx(nn.reco_visible_energy.to_numpy(), e) >= 0) & (binidx(nn.reco_visible_energy.to_numpy(), e) < nb)], minlength=nb)
            def bc(d):
                b = binidx(d.reco_visible_energy.to_numpy(), e); ok = (b >= 0) & (b < nb); return np.bincount(b[ok], minlength=nb)
            lo = bc(nn.loc[nn.index.difference(v.index)]); ga = bc(vv.loc[vv.index.difference(n.index)])
            bn = nn.loc[nn.index.intersection(vv.index)]; bv = vv.loc[bn.index]; mig = bc(bv) - bc(bn)
            res[label]["bins"][ch] = dict(edges=e.tolist(), nominal=Nn.tolist(), lost=(-lo / Nn).tolist(), gained=(ga / Nn).tolist(), migration=(mig / Nn).tolist(), total=((bc(vv) / Nn) - 1).tolist())
        # feature dependence of the losses and gains
        for ci, key in enumerate(FEAT):
            nm, ed = FEAT[key]; ax = axs[ri, 1 + ci]; ctr = 0.5 * (ed[1:] + ed[:-1]); wd = np.diff(ed)
            bn_ = np.digitize(n[key], ed) - 1; tot = np.array([(bn_ == b).sum() for b in range(len(ed) - 1)])
            lost_ = np.array([((bn_ == b) & n.index.isin(lost)).sum() for b in range(len(ed) - 1)])
            bv_ = np.digitize(v[key], ed) - 1; totv = np.array([(bv_ == b).sum() for b in range(len(ed) - 1)]); gain_ = np.array([((bv_ == b) & v.index.isin(gain)).sum() for b in range(len(ed) - 1)])
            fl = lost_ / np.maximum(tot, 1); fg = gain_ / np.maximum(totv, 1)
            ax.errorbar(ctr, fl * 100, xerr=wd / 2, yerr=100 * np.sqrt(fl * (1 - fl) / np.maximum(tot, 1)), fmt="s", color="#D55E00", label="Lost (nominal sample)", ms=5)
            ax.errorbar(ctr, fg * 100, xerr=wd / 2, yerr=100 * np.sqrt(fg * (1 - fg) / np.maximum(totv, 1)), fmt="o", color="#0072B2", label="Gained (varied sample)", ms=5)
            ax.set_xlabel(nm); ax.set_ylabel("Fraction of events [%]"); ax.set_yscale("symlog", linthresh=1); ax.set_ylim(0, None)
            res[label].setdefault("features", {})[key] = dict(edges=ed.tolist(), lost_frac=fl.tolist(), gained_frac=fg.tolist(), n_nom=tot.tolist())
        # decomposition panel (1p and Mp, bins in sequence)
        ax = axs[ri, 0]; x0 = 0; ticks = []; lbl = []
        for ch in ("1p", "Mp"):
            b = res[label]["bins"][ch]; k = len(b["total"]); x = np.arange(k) + x0
            ax.bar(x, np.array(b["lost"]) * 100, color="#D55E00", label="Lost" if ch == "1p" else None); ax.bar(x, np.array(b["gained"]) * 100, color="#0072B2", label="Gained" if ch == "1p" else None)
            ax.bar(x, np.array(b["migration"]) * 100, width=0.35, color="#009E73", label="Energy migration" if ch == "1p" else None); ax.plot(x, np.array(b["total"]) * 100, "k_", ms=14, mew=2.5, label="Net ratio $-$ 1" if ch == "1p" else None)
            ticks += list(x); lbl += [f"{0.5*(e0+e1):.2f}" for e0, e1 in zip(b["edges"][:-1], b["edges"][1:])]; x0 += k + 1
        ax.set_xticks(ticks); ax.set_xticklabels(lbl, rotation=60, fontsize=9); ax.axhline(0, color="0.4", lw=1); ax.set_ylabel("Fraction of nominal [%]"); ax.set_xlabel("Reco visible energy [GeV]  (1p | Mp)")
        ax.text(0.03, 0.95, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")
        r = res[label]; print(label, f"N={r['n_nom']} ratio-1={100*r['ratio']:+.2f}% lost={100*r['lost_frac']:.2f}% gained={100*r['gained_frac']:.2f}% | signal share: nominal {r['nominal_signal_share']:.3f} lost {r['lost_signal_share']:.3f} gained {r['gained_signal_share']:.3f}", flush=True)
    axs[0, 0].legend(fontsize=9, loc="lower right"); axs[0, 1].legend(fontsize=9, loc="upper center")
    mark_axis(axs[0, 0], r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(axs[0, 3], "SBND Run 1", alignment="right"); fig.tight_layout()
    fig.savefig(OUT / "deepdive_sbnd.png", dpi=100); fig.savefig(OUT / "deepdive_sbnd.pdf"); (OUT / "deepdive.json").write_text(json.dumps(res))
