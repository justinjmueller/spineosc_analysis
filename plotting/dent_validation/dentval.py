"""Shared helpers for the DENT data/MC validation suite (SBND detector-access data vs CV and DENT detvar MC)."""
import numpy as np, pandas as pd

SCRATCH = "/private/tmp/claude-550/-Users-mueller-Projects-spineosc/32857cee-ac1d-4b8e-b98a-74f527609260/scratchpad/dv_new"
# energy-like / bookkeeping columns never scanned or plotted
SKIP = {"Run", "Subrun", "Evt", "event_largest_flash_pe", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz"}
DISCRETE = {"reco_cathode_crosser", "reco_neutrino2026_veto", "reco_electron_multiplicity", "reco_muon_multiplicity",
            "reco_photon_multiplicity", "reco_pion_multiplicity", "reco_proton_multiplicity"}

def load(path=SCRATCH):
    d = {k: pd.read_csv(f"{path}/{k}.csv") for k in ("data", "cv", "dent")}
    for df in d.values():
        df["muon_len3d"] = np.sqrt((df.reco_leading_muon_end_x - df.reco_leading_muon_start_x) ** 2 +
                                   (df.reco_leading_muon_end_y - df.reco_leading_muon_start_y) ** 2 +
                                   (df.reco_leading_muon_end_z - df.reco_leading_muon_start_z) ** 2)
        df["abs_end_x"] = df.reco_leading_muon_end_x.abs()
        df["muon_xing_geom"] = ((df.reco_vertex_x * df.reco_leading_muon_end_x) < 0).astype(float)
    return d

def edges_for(col, ref):
    v = ref[col].dropna().to_numpy()
    if col in DISCRETE:
        return np.arange(np.floor(v.min()) - 0.5, min(np.ceil(v.max()), 6) + 1.5)
    lo, hi = np.percentile(v, [0.5, 99.5])
    if col.endswith("softmax"):
        lo, hi = 0.0, 1.0
    if hi <= lo:
        hi = lo + 1
    return np.linspace(lo, hi, 31)

def hist(df, col, edges):
    v = df[col].dropna().to_numpy()
    return np.histogram(np.clip(v, edges[0] + 1e-9, edges[-1] - 1e-9), bins=edges)[0].astype(float)

def metrics(cnt):
    """cnt: dict with 'data','cv','dent' count arrays. All compared as shapes (unit area), with proper
    Poisson-type variances; the CV-vs-DENT separation is evaluated at DATA statistics."""
    nd, nc, ne = cnt["data"].sum(), cnt["cv"].sum(), cnt["dent"].sum()
    pd_, pc, pe = cnt["data"] / nd, cnt["cv"] / nc, cnt["dent"] / ne
    def chi(a, na, b, nb):
        var = a / na + b / nb          # approx: p/N with p~(a+b)/2 would also do; keep the larger/safer sum
        var = np.where(var > 0, var, np.inf)
        return float(np.sum((a * 1.0 - b) ** 2 / var / 1.0) / 1.0) * 1.0
    # all as p-differences: chi2 = sum (pa-pb)^2 / (pa/Na + pb/Nb), i.e. Poisson error on each normalised bin
    def chi2(pa, na, pb, nb):
        var = pa / na + pb / nb
        m = var > 0
        return float(np.sum((pa[m] - pb[m]) ** 2 / var[m])), int(m.sum()) - 1
    s_cd, ndf = chi2(pc, nd, pe, nd)       # CV vs DENT with both scaled to the data statistics
    x_c, _ = chi2(pd_, nd, pc, nc)         # data vs CV
    x_e, _ = chi2(pd_, nd, pe, ne)         # data vs DENT
    return dict(sep=s_cd, chi_cv=x_c, chi_dent=x_e, ndf=ndf)
