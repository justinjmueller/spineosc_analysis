#!/usr/bin/env python3
"""Tables and figure for the duplicate study (reads docs/syst_study/duplicates_all.json, duplicates_subruns.json, duplicates_ratio_impact.json)."""
import json, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mutools.plotting import use_style; from mutools.plotting.helpers import mark_axis
use_style("rootlike"); S = pathlib.Path("/Users/mueller/Projects/spineosc/docs/syst_study"); r = json.loads((S / "duplicates_all.json").read_text())
def short(k): f, t = k.split("::"); return f.replace("spineosc_numu_", "").replace("_v1.1.0", "").replace(".root", ""), t.replace("events/", "")
rows = []
for k, v in r.items():
    f, t = short(k)
    if v.get("index_only") or t.endswith(("/selected", "/signal")) or "wsyst" in f: continue
    if "dirt_wsyst" in f: continue
    rows.append((f, t, v))
lines = ["| file | tree | rows | K1 excess | K2 excess | K1 degeneracy broken by energy | full-row identical excess | K2-duplicate rows identical |", "|---|---|---|---|---|---|---|---|"]
for f, t, v in rows:
    nz = v["K1"]["excess"]; lines.append(f"| {f} | {t[:-7]} | {v['rows']:,} | {nz:,} ({100*nz/v['rows']:.3f}%) | {v['K2']['excess']:,} | " + (f"{100*v['energy_breaks_frac']:.1f}%" if v.get('energy_breaks_frac') is not None else "n/a") + f" | {v['full_row_excess']:,} | {v['K2_dup_rows_all_identical']:,} of {v['K2_dup_rows']:,} |")
(S / "duplicates_events_table.md").write_text("\n".join(lines))
sl = ["| file | tree | rows | K2 excess (distinct interactions in one event + identical copies) | identical-row excess | share of K2 excess that is identical |", "|---|---|---|---|---|---|"]
for k, v in r.items():
    f, t = short(k)
    if v.get("index_only") or not t.endswith("/selected") or "wsyst" in f.replace("wsyst_cvrw", "x") and "cvrw" not in f: continue
    sl.append(f"| {f} | {t[:-9]} | {v['rows']:,} | {v['K2']['excess']:,} ({100*v['K2']['excess']/v['rows']:.3f}%) | {v['full_row_excess']:,} | {100*v['full_row_excess']/max(1,v['K2']['excess']):.1f}% |")
(S / "duplicates_selected_table.md").write_text("\n".join(sl))
# figure: events trees, K1 / K2 / full-row excess as % of rows
sel = [(f, t, v) for f, t, v in rows if ("nosyst" in f) and not t.startswith(("icarus_detvar_sce", "icarus_wiremod_yz", "sbnd_detvar_sce", "sbnd_detvar_m5", "sbnd_detvar_p5", "sbnd_detvar_wiremod_yz", "sbnd_detvar_mattmod", "icarus_detvar_charge")) and "offbeam" not in t]
labels = [f"{t[:-7]}" for f, t, v in sel]
fig, ax = plt.subplots(figsize=(13, 6)); x = np.arange(len(sel)); w = 0.27
for i, (key, lab, c) in enumerate((("K1", "Run, Subrun, Evt", "#999999"), ("K2", "+ max neutrino energy", "#0072B2"), ("full", "identical in every branch", "#D55E00"))):
    vals = [100 * (v["K1"]["excess"] if key == "K1" else v["K2"]["excess"] if key == "K2" else v["full_row_excess"]) / v["rows"] for f, t, v in sel]; ax.bar(x + (i - 1) * w, np.maximum(vals, 1e-5), w, color=c, label=lab)
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=9); ax.set_ylabel("Excess rows [% of tree rows]"); ax.legend(fontsize=10, loc="upper left"); mark_axis(ax, r"$\bf{SBN}$ Analysis-in-Progress"); mark_axis(ax, "events trees (nosyst files)", alignment="right"); fig.tight_layout(); fig.savefig(S / "dup_events_excess.png", dpi=95)
print("tables written;", len(rows), "event trees;", len(sl) - 2, "selected trees")
