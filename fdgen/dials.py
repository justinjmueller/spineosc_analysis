"""Dial weights from the CVMC friend trees, interpolated the way PROfit builds its splines.

Both the standard trees (selected_multisigmaTree) and the reBruce tree (selected_fakedataTree) store, per event, a vector
of knot weights `<name>` and the knot positions in sigma `<name>_sigma`. The standard trees list their knots unsorted
(-1, +1, -2, +2, -3, +3) and have no sigma = 0 entry. PROfit sorts the knots, adds the point (0, 1) when it is missing
(the ratio to the central value) and interpolates with the CAFAna cubic (PROsyst.cxx, "This is cubic interpolation"):
  * interior segments: cubic Hermite with slopes (y[i+1] - y[i-1]) / 2 at both ends (knot-index units),
  * first segment: quadratic through y0, y1 with slope (y2 - y0) / 2 at y1,
  * last segment: y[n-2] + d t + (y[n-1] - y[n-2] - d) t^2 with d = (y[n-1] - y[n-3]) / 2,
  * two points only: linear.
fdgen applies this per event; PROfit applies it to per-bin sums of weights. The two agree exactly at the knots, so
closure tests should use shifts on a knot (integer sigma); between knots they differ at second order, as they do for
products of several dials. Shifts outside the knot range are refused.
"""
import numpy as np
import uproot

TREES = {"standard": "selected_multisigmaTree", "rebruce": "selected_fakedataTree"}


def knots(arr):
    """(n_events, n_knots) array for the modal knot count, and a boolean mask of the events that have it.

    Events without a GENIE interaction record (non-neutrino, or too few weights) carry an empty or short vector; they get
    ratio 1 in dial_weight.
    """
    import awkward as ak
    n = ak.to_numpy(ak.num(arr, axis=1))
    if len(n) == 0:
        return np.ones((0, 0)), np.zeros(0, dtype=bool)
    modal = int(np.bincount(n).argmax())
    ok = n == modal
    out = np.ones((len(n), modal))
    if ok.any():
        out[ok] = ak.to_numpy(arr[ok])
    return out, ok


def cafana_spline(x, s, y):
    """Value at sigma = x of the CAFAna cubic through knots s (sorted, shape (K,)) with per-event values y (n, K)."""
    n = y.shape[1]
    if n < 2:
        raise ValueError("need at least two knots")
    if x < s[0] - 1e-9 or x > s[-1] + 1e-9:
        raise ValueError(f"sigma={x} outside the knot range [{s[0]}, {s[-1]}]")
    seg = min(max(int(np.searchsorted(s, x, side="right")) - 1, 0), n - 2)    # last knot <= x, the last segment closes at s[-1]
    t = (x - s[seg]) / (s[seg + 1] - s[seg])
    if n == 2:
        return y[:, 0] + (y[:, 1] - y[:, 0]) * t
    if seg == 0:
        y1, y2, d = y[:, 0], y[:, 1], (y[:, 2] - y[:, 0]) / 2
        return y1 + (-2 * y1 + 2 * y2 - d) * t + (y1 - y2 + d) * t * t
    if seg == n - 2:
        y5, y6, d = y[:, n - 2], y[:, n - 1], (y[:, n - 1] - y[:, n - 3]) / 2
        return y5 + d * t + (y6 - y5 - d) * t * t
    y0, y1, y2, y3 = y[:, seg - 1], y[:, seg], y[:, seg + 1], y[:, seg + 2]
    d1, d2 = (y2 - y0) / 2, (y3 - y1) / 2
    return y1 + d1 * t + (-3 * y1 + 3 * y2 - 2 * d1 - d2) * t ** 2 + (2 * y1 - 2 * y2 + d1 + d2) * t ** 3


def dial_weight(path, treedir, source, name, sigma, entry_stop=None):
    t = uproot.open(path)[f"{treedir}/selected{TREES[source][len('selected'):]}"]
    w, ok = knots(t[name].array(entry_stop=entry_stop))
    s, ok_s = knots(t[name + "_sigma"].array(entry_stop=entry_stop))
    ok &= ok_s
    if len(ok) == 0:
        return np.ones(0)
    if not np.all(s[ok] == s[ok][0]):
        raise ValueError(f"{name}: knot positions vary between events")
    order = np.argsort(s[ok][0])
    s, w = s[ok][0][order], w[:, order]
    if not np.any(s == 0.0):                       # PROfit's added knob = 0 point: the ratio is 1 at the central value
        at = int(np.searchsorted(s, 0.0))
        s = np.insert(s, at, 0.0)
        w = np.insert(w, at, 1.0, axis=1)
    r = np.ones(len(ok))
    r[ok] = cafana_spline(float(sigma), s, w[ok])
    r[~np.isfinite(r)] = 1.0
    if not ok.all():
        print(f"  [{name}] {int((~ok).sum())} events without the full knot vector keep ratio 1")
    return r
