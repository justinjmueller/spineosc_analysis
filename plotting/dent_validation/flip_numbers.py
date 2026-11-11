"""Numbers quoted in the report's mechanism section, from the event-matched CV/DENT dump (MC only).
Requires flip.csv (made by dump_matched_events.C, selection: >=1 proton, cathode + high-y/z cuts, flash threshold on MC)."""
import numpy as np, pandas as pd, json
d = pd.read_csv("/private/tmp/claude-550/-Users-mueller-Projects-spineosc/32857cee-ac1d-4b8e-b98a-74f527609260/scratchpad/flip.csv", dtype={"cls": str})
cr = d[d.cls.isin(["10", "11"])].copy(); cr["flip"] = cr.cls == "10"; cr["dL"] = cr.Lcv - cr.Ldt; cr["post"] = cr.excv.abs()
cr["mu_cv"] = (cr.vxcv * cr.excv) < 0; cr["mu_dt"] = (cr.vxdt * cr.exdt) < 0
fl = cr[cr.flip]; st = cr[~cr.flip]; nc = d[d.cls == "00"]
out = dict(
    n_matched=int(len(d)), n_crossers_cv=int(len(cr)), n_flipped=int(len(fl)), flip_fraction=float(len(fl) / len(cr)),
    reverse_flips=int((d.cls == "01").sum()), true_crosser_fraction_flipped=float((fl.tcx != 0).mean()),
    median_dL_flipped=float(fl.dL.median()), mean_dL_flipped=float(fl.dL.mean()),
    frac_dL_gt50=float((fl.dL > 50).mean()), frac_dL_lt2=float((fl.dL.abs() < 2).mean()),
    tail_share_of_missing_length=float(fl[fl.dL > 50].dL.sum() / fl.dL.sum()),
    muon_driven=float((fl.mu_cv & ~fl.mu_dt).mean()), non_muon_flag=float((~fl.mu_cv).mean()),
    endpoint_within5_flipped_dent=float((fl.exdt.abs() < 5).mean()), endpoint_within5_flipped_cv=float((fl.excv.abs() < 5).mean()),
    endpoint_within5_noncrossers=float((nc.exdt.abs() < 5).mean()), median_end_x_flipped_dent=float(fl.exdt.abs().median()))
print(json.dumps(out, indent=1)); json.dump(out, open("flip_numbers.json", "w"), indent=1)
