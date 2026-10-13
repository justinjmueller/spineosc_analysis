#!/usr/bin/env python3
"""GUMPLE-style nuisance pull matrix for the 2D fake-data campaign (rows: parameters with the largest pull in any dataset,
columns: fake datasets, one panel per analysis; nuisance-only fit and oscillation-floating fit side by side).

Pulls are the optimiser's best-fit values in prior sigma (see fd2d.read_fit). Datasets are the data-XML ("ds") route: the injected dial is
never a nuisance parameter. Writes plotting/fakedata_2d_figures/pull_matrix_<cfg>.pdf and tables/pull_matrix_<cfg>_<fit>.csv.

  /Users/mueller/Projects/mutools/.venv/bin/python plotting/scripts/fd_pull_matrix.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import fd2d as f2

OUT = Path("/Users/mueller/Projects/spineosc/plotting/fakedata_2d_figures")
(OUT / "tables").mkdir(exist_ok=True, parents=True)
# (stem, short label)
DATASETS = [("ctrl", "CV (closure)"), ("cpi_subthreshold_minus100", "sub-thr. CPI -100%"), ("cpi_subthreshold_plus100", "sub-thr. CPI +100%"),
            ("mec_bdt", "asym. MEC (BDT)"), ("jaesung_lowq2_pi_enhancement_prefsi", "low-Q2 pi, pre-FSI"),
            ("jaesung_lowq2_pi_enhancement_postfsi", "low-Q2 pi, post-FSI"), ("pi_fsi_ha2025", "pi FSI hA2025"),
            ("t2k_nc1pi_p_costh", "T2K NC1pi"), ("ub_cc1p0pi_dpt_dat_loW", "uB CC1p0pi loW"), ("ub_cc1p0pi_dpt_dat_midW", "uB CC1p0pi midW"),
            ("ub_cc1p0pi_dpt_dat_hiW", "uB CC1p0pi hiW"), ("ub_cc2p0pi_delta_pt", "uB CC2p0pi"), ("qe_zexp_mva_to_lqcd", "LQCD"),
            ("minerva_3dqelike_bnb", "MINERvA 3D QE-like")]
TOP = 24


def matrix(cfg, variant):
    cols, vals = [], {}
    for stem, lab in DATASETS:
        if f2.available("ds", cfg, variant, stem):
            fit = f2.read_fit("ds", cfg, variant, stem)
            cols.append(lab); vals[lab] = {n: fit["values"][n] for n in f2.nuisance_names(fit)}
    df = pd.DataFrame(vals)
    return df.loc[df.abs().max(axis=1).sort_values(ascending=False).index[:TOP]] if len(df) else df


if __name__ == "__main__":
    for cfg, name in f2.CONFIGS.items():
        fig, axes = plt.subplots(1, 2, figsize=(15, 8.5), sharey=True)
        for ax, variant, ttl in zip(axes, ("fix", "free"), ("oscillation fixed (null)", "oscillation floating")):
            m = matrix(cfg, variant)
            if m.empty:
                continue
            m.to_csv(OUT / "tables" / f"pull_matrix_{cfg}_{variant}.csv")
            im = ax.imshow(m.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
            for i in range(m.shape[0]):
                for j in range(m.shape[1]):
                    ax.text(j, i, f"{m.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="k" if abs(m.values[i, j]) < 0.7 else "w")
            ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=60, ha="right", fontsize=8)
            ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=8); ax.set_title(f"{name}: {ttl}", fontsize=11)
        fig.colorbar(im, ax=axes, fraction=0.02, label=r"post-fit nuisance value [prior $\sigma$]")
        fig.savefig(OUT / f"pull_matrix_{cfg}.pdf", bbox_inches="tight"); plt.close(fig)
        print("wrote pull_matrix", cfg)
