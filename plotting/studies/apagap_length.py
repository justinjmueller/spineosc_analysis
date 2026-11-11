"""
How long is the SBND APA gap? Direct probes in the detector data access, area-normalised, with no sense of the data's scale.

Companion of ``apagap.py`` (same data, selection, allowed variables and no-scale rules: only ratios, lengths, fractions and
likelihood differences are printed, plotted or stored; the data are never counted out loud). Two probes of the width of the
dead region in z:

1. **Equivalent length** (no model). With ``R(z)`` the data/simulation ratio, the length ``integral (1 - R) dz`` over a
   window is the stretch of z that the missing events would fill if the detector were uniformly efficient. Reported for the
   deficit side and for the pile-up just above it, with their sum (the net loss, i.e. what is neither displaced nor
   compensated).
2. **Box scan with the simulation's own reconstruction.** The simulation's true vertex is smooth across the gap, and its
   reconstruction resolves it to a robust 0.35 cm. A no-gap pseudo-reconstruction is built from the true vertex plus the
   residuals measured away from the gap. Events whose *true* vertex lies in a box ``[a, b)`` are then lost (fraction ``loss``),
   displaced to the +z edge of the box (fraction ``snap``) or reconstructed normally; the reconstructed-vertex histogram is
   compared with the data in 0.5 cm bins (Poisson likelihood, the prediction area-normalised to the data, so a loss only
   changes the shape). ``L = b - a`` is the gap length; its uncertainty comes from the profile of ``-2 ln L`` over
   everything else. Variants: loss only, snap only, and both.

Outputs (``--output``, default ``docs/apa_gap``)::

    length_edges.<fmt>    data/MC in 0.5 cm bins for the vertex, the proton end and the muon end, with the best-fit gap
    length_scan.<fmt>     -2 ln L over the box edges (best model) and the profile in the gap length for the three variants
    length_summary.json   best-fit edges, length and uncertainty, fractions, likelihood differences, equivalent lengths

Usage::

    /Users/mueller/Projects/mutools/.venv/bin/python plotting/studies/apagap_length.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

import apagap as ag  # noqa: E402
import common  # noqa: E402

#: Fit range in vertex z [cm] and the bin width of the compared histograms.
FIT_RANGE = (238.0, 270.0)
BIN = 0.5

#: Resolution of the precomputed (true z, reco z) grid [cm].
GRID = 0.1

#: Edges scanned [cm]: lower edge a, upper edge b (b > a), in steps of ``STEP``.
A_RANGE, B_RANGE, STEP = (249.0, 255.0), (251.0, 258.0), 0.25

#: Number of resolution draws per simulated event when building the pseudo-reconstruction.
N_DRAWS = 12

#: Poisson toys used for the length uncertainty and the closure test.
N_TOYS = 40


def residual_pool(mc: pd.DataFrame, gap_centre: float = 253.0, margin: float = 8.0) -> np.ndarray:
    """Reco - true vertex z of simulated events away from the gap: the resolution without the gap."""
    away = (mc.true_vertex_z.sub(gap_centre).abs() > margin) & mc.true_vertex_z.between(226, 282)
    return (mc.reco_vertex_z - mc.true_vertex_z)[away].to_numpy()


def pseudo_reco(mc: pd.DataFrame, pool: np.ndarray, seed: int = 7):
    """No-gap reconstruction: true vertex z plus a resolution draw from the pool; returns (true, reco, weight)."""
    rng = np.random.default_rng(seed)
    true = np.repeat(mc.true_vertex_z.to_numpy(), N_DRAWS)
    reco = true + rng.choice(pool, size=true.size)
    weight = np.repeat(mc.w.to_numpy(), N_DRAWS) / N_DRAWS
    return true, reco, weight


class GapModel:
    """Predicted reconstructed-vertex histogram for a dead box in true z, from the no-gap pseudo-reconstruction."""

    def __init__(self, mc: pd.DataFrame, pool: np.ndarray, fit_range=FIT_RANGE):
        true, reco, weight = pseudo_reco(mc, pool)
        self.t_edges = np.arange(226.0, 282.0 + GRID / 2, GRID)
        n_fine = int(round((fit_range[1] - fit_range[0]) / GRID))
        self.r_edges = fit_range[0] + GRID * np.arange(n_fine + 1)
        self.fit_edges = np.arange(fit_range[0], fit_range[1] + BIN / 2, BIN)
        self.table, _, _ = np.histogram2d(true, reco, [self.t_edges, self.r_edges], weights=weight)
        self.per = int(round(BIN / GRID))
        self.pool = pool
        self._kernels: dict[float, np.ndarray] = {}

    def rebin(self, fine: np.ndarray) -> np.ndarray:
        """Fine (0.1 cm) reconstructed bins to the compared (0.5 cm) bins."""
        return fine.reshape(-1, self.per).sum(axis=1)

    def box(self, a: float, b: float) -> np.ndarray:
        """Fine reconstructed histogram of the events whose true z lies in [a, b)."""
        i, j = int(round((a - self.t_edges[0]) / GRID)), int(round((b - self.t_edges[0]) / GRID))
        return self.table[i:j].sum(axis=0)

    def kernel(self, edge: float) -> np.ndarray:
        """Fine reconstructed histogram (unit area) of an event displaced to ``edge``: the edge plus a resolution draw."""
        key = round(float(edge), 4)
        if key not in self._kernels:
            hist, _ = np.histogram(edge + self.pool, self.r_edges)
            self._kernels[key] = hist / hist.sum() if hist.sum() else hist.astype(float)
        return self._kernels[key]

    def predict(self, a: float, b: float, loss: float, snap: float) -> np.ndarray:
        """Prediction in the compared bins for a gap [a, b) with the given loss and snap fractions."""
        total = self.table.sum(axis=0)
        inside = self.box(a, b)
        weight_in = self.table[int(round((a - self.t_edges[0]) / GRID)):int(round((b - self.t_edges[0]) / GRID))].sum()
        fine = total - inside + (1 - loss - snap) * inside + snap * weight_in * self.kernel(b)
        return self.rebin(fine)


def deviance(data: np.ndarray, prediction: np.ndarray) -> float:
    """Poisson -2 ln L (up to a constant), the prediction area-normalised to the data."""
    mu = prediction * data.sum() / prediction.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(data > 0, data * np.log(data / mu), 0.0)
    return float(2 * np.sum(mu - data + term))


def bin_deviance(data: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    """Per-bin Poisson deviance (the prediction area-normalised to the data)."""
    mu = prediction * data.sum() / prediction.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(data > 0, data * np.log(data / mu), 0.0)
    return 2 * (mu - data + term)


def fractions(variant: str):
    """The (loss, snap) grid of a model variant."""
    values = np.round(np.arange(0.0, 1.0001, 0.1), 2)
    if variant == "loss":
        return [(x, 0.0) for x in values]
    if variant == "snap":
        return [(0.0, x) for x in values]
    return [(x, y) for x in values for y in values if x + y <= 1.0001]


def scan(model: GapModel, data: np.ndarray, variant: str) -> dict:
    """-2 ln L over the box edges, profiled over the fractions; returns the map and the best fractions."""
    a_grid = np.arange(A_RANGE[0], A_RANGE[1] + STEP / 2, STEP)
    b_grid = np.arange(B_RANGE[0], B_RANGE[1] + STEP / 2, STEP)
    grid = np.full((len(a_grid), len(b_grid)), np.inf)
    best = {}
    for i, a in enumerate(a_grid):
        for j, b in enumerate(b_grid):
            if b - a < 2 * STEP:
                continue
            candidates = [(deviance(data, model.predict(a, b, ls, sn)), ls, sn) for ls, sn in fractions(variant)]
            value, ls, sn = min(candidates)
            grid[i, j] = value
            best[(i, j)] = (ls, sn)
    return {"a": a_grid, "b": b_grid, "dev": grid, "fractions": best}


def profile_length(result: dict) -> tuple[np.ndarray, np.ndarray]:
    """-2 ln L minimised over everything but the length L = b - a."""
    lengths = np.round(result["b"][None, :] - result["a"][:, None], 4)
    values = sorted(set(lengths[np.isfinite(result["dev"])]))
    curve = [result["dev"][(lengths == length) & np.isfinite(result["dev"])].min() for length in values]
    return np.array(values), np.array(curve)


def interval(lengths: np.ndarray, curve: np.ndarray, delta: float = 1.0):
    """Best length and the range where the profile stays within ``delta`` of its minimum (the edges of the scan if it does)."""
    best = int(np.argmin(curve))
    ok = curve <= curve[best] + delta
    lo, hi = best, best
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(curve) - 1 and ok[hi + 1]:
        hi += 1
    return float(lengths[best]), float(lengths[lo]), float(lengths[hi])


def best_of(result: dict) -> dict:
    """Best edges, length, profile interval and fractions of one scan."""
    finite = np.isfinite(result["dev"])
    i, j = np.unravel_index(np.argmin(np.where(finite, result["dev"], np.inf)), result["dev"].shape)
    lengths, curve = profile_length(result)
    best, low, high = interval(lengths, curve)
    loss, snap = result["fractions"][(i, j)]
    return {"a": float(result["a"][i]), "b": float(result["b"][j]), "length": float(result["b"][j] - result["a"][i]),
            "length_best_profile": best, "length_1sigma": [low, high], "loss": loss, "snap": snap,
            "minus2lnL": float(result["dev"][i, j])}


def toys(model: GapModel, prediction: np.ndarray, n_total: float, n_toys: int, seed: int = 11) -> dict:
    """Refit the snap model to Poisson toys of a prediction: the spread of the fitted length and the closure of its mean."""
    rng = np.random.default_rng(seed)
    lengths, lows = [], []
    for _ in range(n_toys):
        counts = rng.poisson(prediction * n_total / prediction.sum())
        fit = best_of(scan(model, counts, "snap"))
        lengths.append(fit["length"])
        lows.append(fit["a"])
    return {"length_mean": float(np.mean(lengths)), "length_std": float(np.std(lengths)),
            "lower_edge_mean": float(np.mean(lows)), "lower_edge_std": float(np.std(lows)), "n_toys": n_toys}


def equivalent_lengths(data: pd.DataFrame, mc: pd.DataFrame) -> dict:
    """integral (1 - R) dz over the deficit window and the pile-up above it, R = data/MC in 0.5 cm bins; no counts returned."""
    edges = np.arange(244.0, 262.01, 0.5)
    nd, nm, _ = ag.binned(data.reco_vertex_z, mc.reco_vertex_z, mc.w, edges)
    ratio = nd / nm
    sigma = np.sqrt(nd) / nm
    centres = (edges[:-1] + edges[1:]) / 2
    out = {}
    for name, (lo, hi) in {"deficit": (251.5, 254.5), "pile_up": (254.5, 256.5)}.items():
        sel = (centres > lo) & (centres < hi)
        sign = 1.0 if name == "deficit" else -1.0
        out[name] = {"window": [lo, hi], "length_cm": float(sign * np.sum((1 - ratio[sel]) * BIN)),
                     "error_cm": float(np.sqrt(np.sum((sigma[sel] * BIN) ** 2)))}
    net = out["deficit"]["length_cm"] - out["pile_up"]["length_cm"]
    out["net_loss_cm"] = {"length_cm": float(net),
                          "error_cm": float(np.hypot(out["deficit"]["error_cm"], out["pile_up"]["error_cm"]))}
    return out


def plot_edges(data, mc, best: dict, path: Path):
    """Data/MC for the vertex, proton end and muon end in 0.5 cm bins; the best-fit gap shaded."""
    import matplotlib.pyplot as plt

    edges = np.arange(244.0, 262.01, 0.5)
    centres = (edges[:-1] + edges[1:]) / 2
    names = {"reco_vertex_z": "vertex", "reco_leading_proton_end_z": "proton end", "reco_leading_muon_end_z": "muon end"}
    fig, axes = plt.subplots(len(names), 1, figsize=(8, 8), sharex=True)
    for ax, (column, label) in zip(axes, names.items(), strict=True):
        nd, nm, _ = ag.binned(data[column], mc[column], mc.w, edges)
        ax.errorbar(centres, nd / nm, yerr=np.sqrt(nd) / nm, fmt="o", color=ag.DATA_COLOR, ms=3.5)
        ax.axhline(1.0, color="0.5", lw=0.8)
        ax.axvspan(best["a"], best["b"], color="C3", alpha=0.15, label="best-fit gap (true z)")
        ax.set_ylim(0.0, 2.2)
        ax.set_ylabel(f"{label}\nData / MC")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[-1].set_xlabel("z [cm]")
    common.save_figure(fig, path, arguments={"kind": "length_edges"})


def plot_scan(results: dict, best_variant: str, path: Path):
    """-2 ln L map over (a, b) for the best variant and the length profiles of all variants."""
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.8))
    result = results[best_variant]
    delta = result["dev"] - np.nanmin(result["dev"][np.isfinite(result["dev"])])
    image = left.pcolormesh(result["b"], result["a"], np.where(np.isfinite(delta), np.minimum(delta, 25), np.nan),
                            cmap="viridis_r", shading="nearest")
    left.contour(result["b"], result["a"], np.where(np.isfinite(delta), delta, 99), levels=[1, 4], colors="w")
    left.set_xlabel("upper edge b [cm]")
    left.set_ylabel("lower edge a [cm]")
    left.set_title(f"{best_variant}: $\\Delta(-2\\ln L)$ (white: 1 and 4)", fontsize=10)
    fig.colorbar(image, ax=left)
    floor = min(np.min(profile_length(r)[1]) for r in results.values())
    for name, res in results.items():
        lengths, curve = profile_length(res)
        right.plot(lengths, curve - floor, label=name)
    right.axhline(1.0, color="0.5", lw=0.8)
    right.set_ylim(-0.5, 12)
    right.set_xlabel("gap length L = b - a [cm]")
    right.set_ylabel("$-2\\ln L$ (profile) relative to the best model")
    right.legend(fontsize=8)
    common.save_figure(fig, path, arguments={"kind": "length_scan"})


def main(argv: list[str] | None = None) -> None:
    """Command-line entry point."""
    parser = common.base_parser(__doc__)
    parser.set_defaults(output=common.DOCS / "apa_gap")
    args = parser.parse_args(argv)
    common.configure_style(fmt=args.format)
    inputs = [common.DATA_DIR / ag.DATA_FILE, common.DATA_DIR / ag.MC_FILE]
    common.set_provenance(__file__, inputs, {"fit_range": list(FIT_RANGE), "bin": BIN},
                          restricted=[common.DATA_DIR / ag.DATA_FILE])

    data, mc = ag.load_data(), ag.load_mc(with_truth=True)
    mc = ag.normalise(data, mc)
    pool = residual_pool(mc)
    model = GapModel(mc, pool)
    counts, _ = np.histogram(data.reco_vertex_z, model.fit_edges)

    null = deviance(counts, model.predict(253.0, 253.5, 0.0, 0.0))
    results = {variant: scan(model, counts, variant) for variant in ("loss", "snap", "loss+snap")}
    summary = {"resolution_robust_sigma_cm": float(1.4826 * np.median(np.abs(pool[np.abs(pool) < 5] - np.median(pool[np.abs(pool) < 5])))),
               "equivalent_length": equivalent_lengths(data, mc), "variants": {}}
    floor = min(np.min(r["dev"][np.isfinite(r["dev"])]) for r in results.values())
    for name, res in results.items():
        finite = np.isfinite(res["dev"])
        i, j = np.unravel_index(np.argmin(np.where(finite, res["dev"], np.inf)), res["dev"].shape)
        lengths, curve = profile_length(res)
        best, low, high = interval(lengths, curve)
        loss, snap = res["fractions"][(i, j)]
        summary["variants"][name] = {
            "a": float(res["a"][i]), "b": float(res["b"][j]), "length": float(res["b"][j] - res["a"][i]),
            "length_best_profile": best, "length_1sigma": [low, high], "loss": loss, "snap": snap,
            "delta_minus2lnL_vs_no_gap": float(null - res["dev"][i, j]), "delta_vs_best_variant": float(res["dev"][i, j] - floor),
        }
    best_name = min(summary["variants"], key=lambda k: summary["variants"][k]["delta_vs_best_variant"])
    summary["best_variant"] = best_name
    best = summary["variants"][best_name]
    ndf = len(counts) - 1 - 3  # area-normalised: -1; free parameters of the snap model: a, b, snap
    summary["goodness_of_fit"] = {"minus2lnL": best["minus2lnL"] if "minus2lnL" in best else None, "ndf": ndf}
    summary["goodness_of_fit"]["minus2lnL"] = float(deviance(counts, model.predict(best["a"], best["b"], best["loss"], best["snap"])))
    best_prediction = model.predict(best["a"], best["b"], best["loss"], best["snap"])
    summary["toys_from_best_fit"] = toys(model, best_prediction, counts.sum(), N_TOYS)
    # where the misfit sits: per-bin deviance inside the gap region and in the rest of the range (expected ~1 per bin)
    centres = (model.fit_edges[:-1] + model.fit_edges[1:]) / 2
    inside = (centres >= 248.0) & (centres < 259.0)
    scatter = {}
    for label, prediction in (("no_gap", model.predict(253.0, 253.5, 0.0, 0.0)), ("best_model", best_prediction)):
        per_bin = bin_deviance(counts, prediction)
        scatter[label] = {"inside_248_259_per_bin": float(per_bin[inside].mean()), "outside_per_bin": float(per_bin[~inside].mean())}
    summary["deviance_per_bin"] = scatter
    factor = max(1.0, scatter["best_model"]["outside_per_bin"])
    summary["length_uncertainty_cm"] = {"toys": summary["toys_from_best_fit"]["length_std"], "scatter_factor": factor,
                                        "inflated": summary["toys_from_best_fit"]["length_std"] * float(np.sqrt(factor))}
    # closure: the same scan applied to the simulation's own reconstruction (what its built-in gap looks like)
    mc_hist, _ = np.histogram(mc.reco_vertex_z, model.fit_edges, weights=mc.w)
    asimov = best_of(scan(model, mc_hist * counts.sum() / mc_hist.sum(), "snap"))
    summary["simulation_reconstruction_as_data"] = {"asimov": asimov, "toys": toys(model, mc_hist, counts.sum(), N_TOYS)}
    summary["robustness"] = {}
    for label, kwargs in {"fit_range_240_268": {"fit_range": (240.0, 268.0)}, "fit_range_236_272": {"fit_range": (236.0, 272.0)},
                          "resolution_margin_5": {"margin": 5.0}, "resolution_margin_12": {"margin": 12.0}}.items():
        alt_pool = residual_pool(mc, margin=kwargs.get("margin", 8.0))
        alt_model = GapModel(mc, alt_pool, kwargs.get("fit_range", FIT_RANGE))
        alt_counts, _ = np.histogram(data.reco_vertex_z, alt_model.fit_edges)
        fit = best_of(scan(alt_model, alt_counts, "snap"))
        summary["robustness"][label] = {k: fit[k] for k in ("a", "b", "length", "length_1sigma", "snap")}
    plot_edges(data, mc, summary["variants"][best_name], args.output / f"length_edges.{args.format}")
    plot_scan(results, best_name, args.output / f"length_scan.{args.format}")
    (args.output / "length_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main(sys.argv[1:])
