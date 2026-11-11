"""Oscillation weight, as PROfit computes it (SBL_2flav_numudis, src/PROmodels/PROmodelSine.cxx).

PROfit does not evaluate the probability at each event's L/E. Events are binned into the model's truth L/E grid
(`<bins unit="True L/E" min max nbins>`, lower edge inclusive, values at or above the last edge dropped) and the
probability is evaluated at the bin midpoint:

    P = 1 - s * sin(1.266932679 * dmsq * le_mid)^2,   s = clamp(sinsq2thmm, 0, 1),   all in float32.
"""
import numpy as np

K = np.float32(1.266932679)


def le_bin_edges(lo, hi, nbins):
    return np.linspace(lo, hi, nbins + 1).astype(np.float32)


def le_bin_index(le, edges):
    """PROfit's Binning::Bin: std::upper_bound, so a lower edge belongs to its bin; -1 outside [edges[0], edges[-1])."""
    le = np.asarray(le, dtype=np.float32)
    idx = np.searchsorted(edges, le, side="right") - 1
    idx[(le < edges[0]) | (le >= edges[-1]) | ~np.isfinite(le)] = -1
    return idx


def pmumu_profit(le, edges, dmsq, sinsq2thmm):
    """Per-event survival probability. Events outside the grid get NaN (PROfit drops them from the prediction)."""
    idx = le_bin_index(le, edges)
    mid = ((edges[1:] + edges[:-1]) / np.float32(2))[np.clip(idx, 0, None)]
    s = np.float32(min(max(sinsq2thmm, 0.0), 1.0))
    sinterm = np.sin(K * np.float32(dmsq) * mid)
    p = np.float32(1) - s * sinterm * sinterm
    p = p.astype(np.float64)
    p[idx < 0] = np.nan
    return p
