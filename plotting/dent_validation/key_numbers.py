"""Headline cathode-feature numbers (data vs CV vs DENT), regenerated with the corrected data selection (no flash cut on data)."""
import sys; sys.path.insert(0, '.')
import numpy as np, pandas as pd, dentval as dv
d = dv.load(); N = {k: len(v) for k, v in d.items()}
rows = []
def row(name, mask):
    c = {k: int(mask(d[k]).sum()) for k in d}; f = {k: c[k] / N[k] for k in d}; err = np.sqrt(c["data"]) / N["data"]
    rows.append(dict(feature=name, data=100 * f["data"], cv=100 * f["cv"], dent=100 * f["dent"], data_over_cv=f["data"] / f["cv"], dent_over_cv=f["dent"] / f["cv"],
                     err_over_cv=err / f["cv"], sig_vs_cv=(f["data"] - f["cv"]) / err, sig_vs_dent=(f["data"] - f["dent"]) / err, n_data=c["data"]))
x = lambda lo, hi, xing: (lambda df: (df.muon_xing_geom == xing) & (df.reco_leading_muon_end_x >= lo) & (df.reco_leading_muon_end_x < hi))
row("crossing muons, |end_x| < 4 cm", lambda df: (df.muon_xing_geom == 1) & (df.abs_end_x < 4))
row("non-crossing muons, |end_x| < 4 cm", lambda df: (df.muon_xing_geom == 0) & (df.abs_end_x < 4))
row("all muons, |end_x| < 4 cm", lambda df: df.abs_end_x < 4)
row("crosser flag, all", lambda df: df.reco_cathode_crosser == 1)
row("non-crosser flag, |end_x| < 5 cm (pile-up)", lambda df: (df.reco_cathode_crosser == 0) & (df.abs_end_x < 5))
t = pd.DataFrame(rows); pd.set_option("display.width", 220)
print(f"N data={N['data']} cv={N['cv']} dent={N['dent']}"); print(t.round(3).to_string(index=False)); t.to_csv("key_numbers.csv", index=False)
