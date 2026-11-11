#!/usr/bin/env python3
"""Event 'fragility': across all SBND variations that share the nominal (sbnd_large), how many variations drop a given nominal-selected event? Compared with an independent-binomial expectation at the
same mean loss rate.  Also per truth class.  Output docs/syst_study/persyst/fragile.json (+ figure fragile.png)."""
from ps_common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from mutools.plotting import use_style; use_style("rootlike"); plt.rcParams.update({"axes.labelsize": 10, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8, "axes.titlesize": 9})
KEYK = KEYS + ["Ei"]; out = {}
def run(det, vars_, tag_):
    base = None; cols = []
    for v in vars_:
        d = load(f"{det}__{v}", v2=False); n = d.loc[d.in_n, KEYK + ["n_true_cc", "n_true_pdg", "n_true_category"]].copy(); n[f"lost_{v}"] = ~d.loc[d.in_n, "in_v"].to_numpy()
        base = n if base is None else base.merge(n[KEYK + [f"lost_{v}"]], on=KEYK, how="inner"); del d
    L = np.column_stack([base[f"lost_{v}"].to_numpy() for v in vars_]); K = L.shape[1]; cnt = L.sum(1); p = L.mean()
    obs = np.bincount(cnt, minlength=K + 1) / len(cnt); exp = stats.binom.pmf(np.arange(K + 1), K, p)
    tc = truth_class(base["n_true_cc"].to_numpy(float), base["n_true_pdg"].to_numpy(float), base["n_true_category"].to_numpy(float))
    res = dict(K=K, n_events=int(len(cnt)), mean_loss_rate=float(p), observed=obs.tolist(), binomial=exp.tolist(), frac_lost_ge3=float((cnt >= 3).mean()), binom_ge3=float(exp[3:].sum()),
               share_of_loss_incidents_from_events_lost_ge3=float(cnt[cnt >= 3].sum() / cnt.sum()), never_lost=float((cnt == 0).mean()), binom_never=float(exp[0]),
               per_class={CLASSES[c]: dict(n=int((tc == c).sum()), mean_lost=float(cnt[tc == c].mean() / K) if (tc == c).any() else None, frac_ge3=float((cnt[tc == c] >= 3).mean()) if (tc == c).any() else None) for c in range(5)})
    # by-class split by number of lost variations for signal only
    sig = tc == 0; res["signal_only"] = dict(observed=(np.bincount(cnt[sig], minlength=K + 1) / sig.sum()).tolist(), binomial=stats.binom.pmf(np.arange(K + 1), K, L[sig].mean()).tolist(), mean=float(L[sig].mean()))
    print(det, tag_, "K", K, "N", len(cnt), "mean loss", round(p, 4), "P(>=3) obs/binom", round(res["frac_lost_ge3"], 4), round(res["binom_ge3"], 5), "never lost obs/binom", round(res["never_lost"], 3), round(exp[0], 3))
    return res
sb = ["SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ", "MattModYZ", "ChargeScale", "ChargeScale15", "GainHigh", "GainLow"]
out["sbnd_all9"] = run("sbnd_run1", sb, "9 variations (4 charge/gain-type included only as a fragility probe)")
out["sbnd_geo5"] = run("sbnd_run1", ["SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ", "MattModYZ"], "5 non-charge variations")
for run_ in ("icarus_run2", "icarus_run4"):
    out[run_] = run(run_, ["SCE0x", "SCE2x", "WireModxThetaXW", "WireModYZ", "ChargeScale"], "5 variations")
json.dump(out, open(PSFIG / "fragile.json", "w"))
fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
for ax, k, t in zip(axs, ("sbnd_all9", "icarus_run2", "icarus_run4"), ("SBND (9 variations)", "ICARUS Run 2 (5)", "ICARUS Run 4 (5)")):
    r = out[k]; x = np.arange(r["K"] + 1); ax.bar(x - 0.2, r["observed"], 0.4, color="#0072B2", label="observed"); ax.bar(x + 0.2, r["binomial"], 0.4, color="#D55E00", label="independent binomial, same mean")
    ax.set_yscale("log"); ax.set_xlabel("number of variations that drop the (nominal-selected) event"); ax.set_ylabel("fraction of events"); ax.set_title(t); ax.legend()
from mutools.plotting.helpers import mark_axis; mark_axis(axs[0], r"$\bf{SBN}$ Analysis-in-Progress")
fig.tight_layout(); fig.savefig(PSFIG / "fragile.png", dpi=110)
