#!/usr/bin/env python3
"""
How much of each detector-model uncertainty is stochastic (finite MC statistics)?

The detector systematics enter the fit as per-bin ratios of a varied sample to a nominal sample.
Each ratio carries the true effect of the variation plus statistical noise from the finite number of
simulated events, and the noise is currently treated as systematic. This study measures the noise.

Method
------
The two samples of each pair are matched event by event on the key (Run, Subrun, Evt,
event_max_neutrino_energy). Every event of the union gets one Poisson(1) weight per bootstrap
replica, the *same* weight in both samples, so the pairing, and therefore the correlation between
the numerator and the denominator, is preserved. Events that appear in only one sample are
resampled independently like any other event. In each replica the ratio ``varied / nominal`` is
recomputed in every bin; the spread over replicas is the statistical uncertainty ``sigma`` of the
ratio. Per bin and variation the output gives the observed shift ``r - 1``, ``sigma``, the
significance ``z = (r - 1) / sigma``, and a noise-subtracted size
``sqrt(max((r - 1)^2 - sigma^2, 0))``. Per variation it gives the fraction of the summed squared
shift that is noise, ``sum(sigma^2) / sum((r - 1)^2)``.

Nothing here smooths or modifies the ratios; it only quantifies how much of each is noise.

Usage
-----
    python plotting/studies/detvar_stat_bootstrap.py [--n-boot 500] [--output docs/detvar_stat]
    python plotting/studies/detvar_stat_bootstrap.py --null-test     # calibration

Outputs: ``bootstrap.csv`` (per bin), ``summary.csv`` (per variation), and a figure.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, REPO, configure_style, save_figure  # noqa: E402

DATA = DATA_DIR
NOMINAL_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_nosyst.root"
DETVAR_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"
KEY = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
END_CUT = 4.0
EDGES = np.array([0.3, 0.4, 0.46, 0.52, 0.58, 0.64, 0.70, 0.76, 0.82, 0.88, 0.94, 1.0, 1.06, 1.13,
                  1.2, 1.3, 1.5, 2.0])
SEED = 20261014

# systematic -> (nominal file, nominal tree, varied tree). DENT has its own paired nominal sample;
# the others are compared with the production nominal sample, as in the detsys calculator.
PRODUCTION = (NOMINAL_FILE, "sbnd_large")
DETVAR_CV = (DETVAR_FILE, "sbnd_detvar_newcv")
SYSTEMATICS = {
    "DENT": (DETVAR_CV, "sbnd_detvar_newdent"),
    "ChargeScale": (PRODUCTION, "sbnd_detvar_charge_scale"),
    "ChargeScale15": (PRODUCTION, "sbnd_detvar_smear15"),
    "SCE0x": (PRODUCTION, "sbnd_detvar_sce0x"),
    "SCE2x": (PRODUCTION, "sbnd_detvar_sce2x"),
    "WireModxThetaXW": (PRODUCTION, "sbnd_detvar_wiremod_xthetaxw"),
    "WireModYZ": (PRODUCTION, "sbnd_detvar_wiremod_yz"),
    "MattModYZ": (PRODUCTION, "sbnd_detvar_mattmod_yz"),
    "GainHigh": (PRODUCTION, "sbnd_detvar_p5_gain"),
    "GainLow": (PRODUCTION, "sbnd_detvar_m5_gain"),
}
# The DENT pair is made of two samples with different POT; the analysis uses only the events common
# to both (as the detsys calculator does), so the bootstrap must too, or the POT ratio becomes an
# offset of about -5% in every bin.
COMMON_ONLY = {"DENT"}
COLUMNS = [*KEY, "reco_visible_energy", "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
           "reco_cathode_crosser", "reco_leading_muon_end_x"]


def load_energy(path: Path, tree: str) -> pd.Series:
    """Reco visible energy of the events selected under S2, indexed by event key."""
    frame = uproot.open(path)[f"events/{tree}/selected"].arrays(COLUMNS, library="pd")
    frame = frame.drop_duplicates(KEY, keep=False).set_index(KEY)
    keep = ((frame.reco_fiducialize_cathode == 1) & (frame.reco_veto_sbnd_highy_highz == 1)
            & (frame.reco_cathode_crosser == 0) & (frame.reco_leading_muon_end_x.abs() >= END_CUT))
    return frame.loc[keep, "reco_visible_energy"]


def common_keys(path: Path, nominal_tree: str, varied_tree: str) -> pd.MultiIndex:
    """Keys of the events present in the ``events`` trees of both samples."""
    def keys(tree: str) -> pd.MultiIndex:
        frame = uproot.open(path)[f"events/{tree}/events"].arrays(KEY, library="pd")
        return pd.MultiIndex.from_frame(frame.drop_duplicates(KEY))

    return keys(nominal_tree).intersection(keys(varied_tree))


def bin_index(energy: np.ndarray) -> np.ndarray:
    """Bin number of each energy, -1 outside the binning (and for NaN)."""
    index = np.digitize(energy, EDGES) - 1
    index[(index < 0) | (index >= len(EDGES) - 1) | np.isnan(energy)] = -1
    return index


def bootstrap_ratio(nominal: pd.Series, varied: pd.Series, n_boot: int,
                    rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bin counts and the replica ratios of a paired sample.

    Returns the nominal and varied counts (nominal weights) and the ratio in every replica,
    with the replica axis first.
    """
    union = nominal.index.union(varied.index)
    nom_bin = bin_index(nominal.reindex(union).to_numpy())
    var_bin = bin_index(varied.reindex(union).to_numpy())
    weights = rng.poisson(1.0, size=(n_boot + 1, len(union))).astype(np.float32)
    weights[0] = 1.0                                           # row 0 is the observed sample
    n_bins = len(EDGES) - 1
    nom_counts = np.stack([weights[:, nom_bin == b].sum(axis=1) for b in range(n_bins)], axis=1)
    var_counts = np.stack([weights[:, var_bin == b].sum(axis=1) for b in range(n_bins)], axis=1)
    return nom_counts[0], var_counts[0], var_counts / np.maximum(nom_counts, 1)


def run(n_boot: int, output: Path) -> None:
    """Bootstrap every systematic and write the tables and the figure."""
    rng = np.random.default_rng(SEED)
    cache: dict[tuple[Path, str], pd.Series] = {}

    def energy(spec: tuple[Path, str]) -> pd.Series:
        if spec not in cache:
            cache[spec] = load_energy(*spec)
        return cache[spec]

    rows = []
    for name, (nominal_spec, varied_tree) in SYSTEMATICS.items():
        nominal = energy(nominal_spec)
        varied = energy((DETVAR_FILE, varied_tree))
        if name in COMMON_ONLY:
            common = common_keys(DETVAR_FILE, nominal_spec[1], varied_tree)
            nominal, varied = nominal[nominal.index.isin(common)], varied[varied.index.isin(common)]
        n_nom, n_var, ratios = bootstrap_ratio(nominal, varied, n_boot, rng)
        shift = ratios[0] - 1
        sigma = ratios[1:].std(axis=0)
        for b in range(len(EDGES) - 1):
            rows.append({"systematic": name, "bin": b, "lo": EDGES[b], "hi": EDGES[b + 1],
                         "n_nominal": int(n_nom[b]), "n_varied": int(n_var[b]),
                         "delta": shift[b], "sigma_stat": sigma[b],
                         "z": shift[b] / sigma[b] if sigma[b] > 0 else np.nan,
                         "shift_noise_subtracted": np.sqrt(max(shift[b] ** 2 - sigma[b] ** 2, 0))})
        print(f"{name}: {len(nominal)} nominal, {len(varied)} varied events selected")
    table = pd.DataFrame(rows)
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "bootstrap.csv", index=False)

    summary = table.groupby("systematic", sort=False).apply(summarize, include_groups=False)
    summary.to_csv(output / "summary.csv")
    print("\n" + summary.round(3).to_string())
    plot(table, output)


def summarize(group: pd.DataFrame) -> pd.Series:
    """Per-systematic noise summary over the energy bins."""
    total = float((group.delta ** 2).sum())
    noise = float((group.sigma_stat ** 2).sum())
    return pd.Series({
        "rms_shift": float(np.sqrt((group.delta ** 2).mean())),
        "rms_sigma_stat": float(np.sqrt((group.sigma_stat ** 2).mean())),
        "noise_fraction": noise / total if total > 0 else np.nan,
        "bins_with_abs_z_lt_2": int((group.z.abs() < 2).sum()),
        "bins": len(group),
    })


def plot(table: pd.DataFrame, output: Path) -> None:
    """Shift and statistical uncertainty per bin for every systematic."""
    import matplotlib.pyplot as plt

    configure_style()
    names = list(dict.fromkeys(table.systematic))
    fig, axes = plt.subplots(2, 5, figsize=(20, 7), sharex=True)
    for ax, name in zip(axes.ravel(), names, strict=False):
        sub = table[table.systematic == name]
        centers = 0.5 * (sub.lo + sub.hi)
        ax.errorbar(centers, 100 * sub.delta, yerr=100 * sub.sigma_stat, fmt="o", ms=3,
                    capsize=2)
        ax.axhline(0, color="0.5", lw=1)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("Reconstructed visible energy [GeV]")
    axes[0, 0].set_ylabel("Varied / nominal - 1 [%]")
    axes[1, 0].set_ylabel("Varied / nominal - 1 [%]")
    save_figure(fig, output / "bootstrap.pdf")


def null_test(n_repeats: int = 300, n_boot: int = 150, min_energy: float = 1.3) -> None:
    """Calibrate the bootstrap: split one sample at random into two halves (true ratio 1).

    If the replica spread is the right statistical uncertainty, ``z = (ratio - 1) / sigma`` has unit
    rms and about 4.6% of values exceed 2. Only the highest-energy bins are used, where the counts
    are smallest and the bootstrap is least trustworthy.
    """
    rng = np.random.default_rng(SEED)
    nominal = load_energy(*PRODUCTION)
    nominal = nominal[nominal >= min_energy]
    z_by_bin: dict[int, list[float]] = {}
    for _ in range(n_repeats):
        in_a = rng.random(len(nominal)) < 0.5
        _, _, ratios = bootstrap_ratio(nominal[in_a], nominal[~in_a], n_boot, rng)
        z = (ratios[0] - 1) / ratios[1:].std(axis=0)
        for b in np.nonzero(EDGES[:-1] >= min_energy)[0]:
            z_by_bin.setdefault(int(b), []).append(float(z[b]))
    for b, values in z_by_bin.items():
        z = np.array(values)
        print(f"{EDGES[b]:.2f}-{EDGES[b + 1]:.2f} GeV: z rms {np.sqrt((z ** 2).mean()):.2f}, "
              f"fraction |z| > 2: {(np.abs(z) > 2).mean():.3f} (ideal: 1.00, 0.046)")


def main() -> None:
    """Parse arguments and run."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--n-boot", type=int, default=500)
    parser.add_argument("--output", type=Path, default=REPO / "docs" / "detvar_stat")
    parser.add_argument("--null-test", action="store_true",
                        help="calibrate the bootstrap on random half-splits instead of running it")
    args = parser.parse_args()
    if args.null_test:
        null_test()
    else:
        run(args.n_boot, args.output)


if __name__ == "__main__":
    main()
