#!/usr/bin/env python3
"""Defines the wide-net DENT-study observables and writes observables.json.

Bin edges come from the detvar CV MC ONLY (equal-population quantiles), never from data. Each observable is energy-free and
exists in the SBND detector-access data tree. `expr_xml` is the TTreeFormula spelling used in the PROfit XML, `expr_cpp` the
C++ spelling used by the detsys calculator (they differ only for abs).
"""
import json, numpy as np, pandas as pd, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "plotting" / "dent_validation"))
import dentval as dv

L3D = "sqrt(pow(reco_leading_muon_end_x-reco_leading_muon_start_x,2)+pow(reco_leading_muon_end_y-reco_leading_muon_start_y,2)+pow(reco_leading_muon_end_z-reco_leading_muon_start_z,2))"
# name, column or expression, label, kind, (nbins), notes
SPECS = [
    # --- cathode-region / DENT-mechanism observables (primary candidates, chosen from MC-only DENT sensitivity)
    ("muon_calo_csda",     "reco_leading_muon_calo_csda_diff",  r"Muon calo-vs-CSDA residual",              "cont", 20),
    ("muon_len3d",         L3D,                                  r"Muon straight-line length [cm]",          "cont", 20),
    ("muon_abs_endx",      "abs(reco_leading_muon_end_x)",       r"Muon $|x_{end}|$ [cm]",                   "cont", 20),
    ("crosser_flag",       "reco_cathode_crosser",               r"Cathode-crosser flag",                    "disc", 2),
    ("muon_endx_full",     "reco_leading_muon_end_x",            r"Muon end $x$ [cm]",                       "cont", 24),
    ("vertexx_full",       "reco_vertex_x",                      r"Vertex $x$ [cm]",                         "cont", 24),
    ("proton_endx_full",   "reco_leading_proton_end_x",          r"Leading-proton end $x$ [cm]",             "cont", 24),
    ("muon_startx_full",   "reco_leading_muon_start_x",          r"Muon start $x$ [cm]",                     "cont", 24),
    # --- PID / calorimetry
    ("proton_calo_csda",   "reco_leading_proton_calo_csda_diff", r"Proton calo-vs-CSDA residual",            "cont", 20),
    ("muon_proton_score",  "reco_leading_muon_proton_softmax",   r"Muon track: proton score",                "soft", 16),
    ("muon_pion_score",    "reco_leading_muon_pion_softmax",     r"Muon track: pion score",                  "soft", 16),
    ("proton_mip_score",   "reco_leading_proton_mip_softmax",    r"Proton track: MIP score",                 "soft", 16),
    ("proton_pion_score",  "reco_leading_proton_pion_softmax",   r"Proton track: pion score",                "soft", 16),
    ("proton_muon_score",  "reco_leading_proton_muon_softmax",   r"Proton track: muon score",                "soft", 16),
    # --- other geometry (controls and diagnostics)
    ("muon_endy",          "reco_leading_muon_end_y",            r"Muon end $y$ [cm]",                       "cont", 20),
    ("vertexy",            "reco_vertex_y",                      r"Vertex $y$ [cm]",                         "cont", 20),
    ("proton_endy",        "reco_leading_proton_end_y",          r"Leading-proton end $y$ [cm]",             "cont", 20),
    ("proton_endz",        "reco_leading_proton_end_z",          r"Leading-proton end $z$ [cm]",             "cont", 20),
    ("muon_dir_x",         "reco_leading_muon_start_dir_x",      r"Muon direction $x$",                      "cont", 20),
    ("muon_dir_y",         "reco_leading_muon_start_dir_y",      r"Muon direction $y$",                      "cont", 20),
    ("muon_dir_z",         "reco_leading_muon_start_dir_z",      r"Muon direction $z$",                      "cont", 20),
]

def cpp(expr):
    return expr.replace("abs(", "std::abs(")

def sig3(x):
    return float(f"{x:.3g}")

def main():
    cv = dv.load()["cv"]               # MC only; the data frame is loaded by dentval but never touched here
    cols = {c: cv[c].to_numpy(dtype=float) for c in cv.columns}
    obs = []
    for name, expr, label, kind, nb in SPECS:
        v = np.asarray(eval(expr, {"sqrt": np.sqrt, "pow": np.power, "abs": np.abs}, cols), dtype=float)
        v = v[np.isfinite(v)]
        if kind == "disc":
            hi = {"reco_cathode_crosser": 1, "reco_pion_multiplicity": 3, "reco_photon_multiplicity": 3, "reco_electron_multiplicity": 2}[expr]
            edges = [i - 0.5 for i in range(0, hi + 2)]
        else:
            lo, hi = np.percentile(v, [0.5, 99.5]) if kind == "cont" else (0.0, 1.0)
            if kind == "soft":
                q = np.quantile(v, np.linspace(0.02, 0.98, nb - 1))
                edges = sorted(set([0.0] + [sig3(x) for x in q] + [1.0]))
            else:
                q = np.quantile(np.clip(v, lo, hi), np.linspace(0, 1, nb + 1))
                edges = sorted(set(sig3(x) for x in q))
        edges = [float(e) for e in edges]
        n_mc = np.histogram(np.clip(v, edges[0], edges[-1]), bins=edges)[0]
        obs.append(dict(name=name, expr_xml=expr, expr_cpp=cpp(expr), label=label, kind=kind, edges=edges,
                        nbins=len(edges) - 1, mc_min_bin=int(n_mc.min()), mc_events=int(len(v))))
    out = pathlib.Path(__file__).with_name("observables.json")
    out.write_text(json.dumps(obs, indent=1))
    for o in obs:
        print(f"{o['name']:20s} {o['nbins']:3d} bins  [{o['edges'][0]:g}, {o['edges'][-1]:g}]  min MC bin {o['mc_min_bin']:6d}")
    print("wrote", out)

if __name__ == "__main__":
    main()
