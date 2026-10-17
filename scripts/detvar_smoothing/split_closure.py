#!/usr/bin/env python3
"""Closure test 2 with a three-way event split: how well does a ratio represent its variation?

Events are assigned to one of three parts by a deterministic hash of (Run, Subrun, Evt), so the
parts are independent and the same expression works in PROfit weights and in pandas:

    ratio part   the varied / nominal ratio of the systematic under test is made from these events
    MC part      the nominal (unweighted) sample, with the full systematic suite, is the prediction
    data part    the variation sample itself is the pseudo-data

A global fit with every systematic floating (absolute normalization) then reports the fitted value
of the injected parameter (1 if the ratio represents the variation at one sigma) and the chi-square.
The same fit with the pseudo-data taken from the *reference* sample of the pair (no variation) is
the null. PROfit bakes the ratio histograms into its processed binaries and its cache hash ignores
them, so every ratio file gets its own tag. See docs/detvar_smoothing/PLAN.md.

    python scripts/detvar_smoothing/split_closure.py muon_score_inclusive_s2 DENT --perm 0
"""
import argparse
import functools
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import closure  # noqa: E402

RUN = closure.RUN
DATA = Path("/Users/mueller/data/spineosc/v1.1.0")
NOMINAL_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_nosyst.root"
DETVAR_FILE = DATA / "spineosc_numu_sbnd_v1.1.0_detvar_nosyst.root"
KEY = ["Run", "Subrun", "Evt", "event_max_neutrino_energy"]
DETECTOR_POT = 1.0e20
SPLIT = "((Evt + 7*Subrun + 13*Run) % 3)"
FLASH = "(0.647 * event_largest_flash_pe > 2000)"

# systematic -> (reference file, reference tree, varied tree)
PAIRS = {
    "DENT": (DETVAR_FILE, "sbnd_detvar_newcv", "sbnd_detvar_newdent"),
    "ChargeScale15": (NOMINAL_FILE, "sbnd_large", "sbnd_detvar_smear15"),
    "WireModxThetaXW": (NOMINAL_FILE, "sbnd_large", "sbnd_detvar_wiremod_xthetaxw"),
    "MattModYZ": (NOMINAL_FILE, "sbnd_large", "sbnd_detvar_mattmod_yz"),
}
# study -> (observable branch, ratio file of the systematics that are not under test, ratio
# file name the study XML refers to)
STUDIES = {
    "muon_score_inclusive_s2": ("reco_leading_muon_muon_softmax",
                                "work/detsys_v2_muon_score.root"),
    "reco_energy_inclusive_s2": ("reco_visible_energy", "work/detsys_main.root"),
}


def part_of(frame: pd.DataFrame) -> np.ndarray:
    """Part (0, 1, 2) of each event, the same hash as the PROfit weight expression."""
    return ((frame.Evt + 7 * frame.Subrun + 13 * frame.Run) % 3).to_numpy()


@functools.cache
def pot(path: Path, tree: str) -> float:
    """Protons on target of one sample."""
    return float(uproot.open(path)[f"events/{tree}/POT"].values()[0])


@functools.cache
def load_sample(path: Path, tree: str, variable: str) -> pd.DataFrame:
    """Events passing the calculator selection of the S2 ratio files, with their part."""
    columns = [*KEY, variable, "reco_fiducialize_cathode", "reco_veto_sbnd_highy_highz",
               "reco_cathode_crosser", "reco_leading_muon_end_x"]
    frame = uproot.open(path)[f"events/{tree}/selected"].arrays(columns, library="pd")
    frame = frame.drop_duplicates(KEY, keep=False)
    keep = ((frame.reco_fiducialize_cathode == 1) & (frame.reco_veto_sbnd_highy_highz == 1)
            & (frame.reco_cathode_crosser == 0) & (frame.reco_leading_muon_end_x.abs() >= 4))
    frame = frame[keep].copy()
    frame["part"] = part_of(frame)
    return frame


def ratio_from_part(syst: str, variable: str, edges: np.ndarray,
                    part: int) -> tuple[np.ndarray, np.ndarray]:
    """Varied / nominal ratio per bin and its sigma from one part (1 and 1 where undefined)."""
    ref_file, ref_tree, var_tree = PAIRS[syst]
    reference = load_sample(ref_file, ref_tree, variable)
    varied = load_sample(DETVAR_FILE, var_tree, variable)
    if syst == "DENT":                         # the analysis uses only the events common to both
        common = reference.set_index(KEY).index.intersection(varied.set_index(KEY).index)
        reference = reference[reference.set_index(KEY).index.isin(common)]
        varied = varied[varied.set_index(KEY).index.isin(common)]
    nom = np.histogram(reference.loc[reference.part == part, variable], edges)[0].astype(float)
    var = np.histogram(varied.loc[varied.part == part, variable], edges)[0].astype(float)
    # the samples have different POT (the DENT pair is restricted to common events, so equal)
    scale = 1.0 if syst == "DENT" else pot(ref_file, ref_tree) / pot(DETVAR_FILE, var_tree)
    ratio = np.where(nom > 0, scale * var / np.maximum(nom, 1), 1.0)
    sigma = np.where((nom > 0) & (var > 0),
                     ratio * np.sqrt(1 / np.maximum(var, 1) + 1 / np.maximum(nom, 1)), 1.0)
    return ratio, sigma


def fit_range(study: str) -> tuple[float, float]:
    """The observable range of the fit, from the first <bins> element of the study XML."""
    text = closure.study_file(study, ".xml").read_text()
    line = re.search(r'<bins unit="[^"]*"[^>]*min="([^"]*)" max="([^"]*)"', text)
    if line is None:                        # an explicit-edges axis
        edges = re.search(r'<bins unit="[^"]*"\s+edges="([^"]*)"', text).group(1).split()
        return float(edges[0]), float(edges[-1])
    return float(line.group(1)), float(line.group(2))


def smooth_in_range(ratio: np.ndarray, sigma: np.ndarray, edges: np.ndarray,
                    window: tuple[float, float], method: str) -> np.ndarray:
    """Smooth only the bins inside the fit window; the others are left as they are."""
    import smoother_bias_variance as bv

    centers = 0.5 * (edges[1:] + edges[:-1])
    inside = (centers > window[0]) & (centers < window[1])
    out = ratio.copy()
    out[inside] = bv.smooth_all(ratio[inside], np.maximum(sigma[inside], 1e-4))[method]
    return out


def write_ratio_file(study: str, systs: list[str], tag: str, part: int,
                     method: str = "none") -> Path:
    """The study's ratio file with each listed systematic replaced by its part-made ratio."""
    variable, source = STUDIES[study]
    original = uproot.open(REPO / source)
    replaced = {}
    for syst in systs:
        target = f"{syst}_sbnd_Run1"
        edges = original[target].to_numpy()[1]
        new, sigma = ratio_from_part(syst, variable, edges, part)
        if method != "none":
            new = smooth_in_range(new, sigma, edges, fit_range(study), method)
        replaced[target] = new
    text = closure.study_file(study, ".xml").read_text()
    allowlist = re.findall(r'<allowlist[^>]*type="hist1d"[^>]*\sname="([^"]+)"', text)
    available = set(original.keys(cycle=False))
    wanted = set(re.findall(r'histname="([^"]+)"', text)) | set(allowlist)
    names = sorted(wanted & available)      # allowlist names that are only section names: skip
    path = RUN / "ratios" / f"{tag}.root"
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as out:
        for name in names:
            values, hist_edges = original[name].to_numpy()
            out[name] = (replaced.get(name, values), hist_edges)
    return path


def mc_xml(study: str, ratio_path: Path, mc_part: int) -> str:
    """The Tier B XML as a closure MC: unweighted nominal, only the main block, one part."""
    text = closure.study_file(study, ".xml").read_text()
    blocks = re.findall(r"<MCFile .*?</MCFile>\s*", text, flags=re.S)
    main = [b for b in blocks if 'treename="events/sbnd_large/selected"' in b.split(">", 1)[0]]
    others = [b for b in blocks if b not in main]
    sample_pot = float(uproot.open(NOMINAL_FILE)["events/sbnd_large/POT"].values()[0])
    header = (f'<MCFile treename="events/sbnd_large/selected" filename="work/sbnd_run1_cvmc.root" '
              f'scale = "3.0" pot="{sample_pot:.6e}">')
    block = re.sub(r"<MCFile[^>]*>", header, main[0], count=1).replace(' * cvwgt"', '"')

    def split(match):
        weight = match.group(2)
        return f'{match.group(1)}"{weight[:-1]} && ({SPLIT} == {mc_part}))"'

    block = re.sub(r'(additional_weight\s*=\s*)"([^"]*)"', split, block)
    for other in others:
        text = text.replace(other, "")
    text = text.replace(main[0], block)
    text = re.sub(r'(<detector name="SBND"\s+pot=")[^"]*"', rf'\g<1>{DETECTOR_POT:.6e}"', text)
    return text.replace(STUDIES[study][1], f"ratios/{ratio_path.name}")


def data_xml(study: str, tree: str, path: Path, data_part: int) -> str:
    """The data XML turned into a pseudo-data block reading one sample's events of one part."""
    text = closure.study_file(study, "_data.xml").read_text()
    sample_pot = float(uproot.open(path)[f"events/{tree}/POT"].values()[0])
    header = (f'<MCFile treename="events/{tree}/selected" filename="work/{path.name}" '
              f'scale = "3.0" pot="{sample_pot:.6e}">')

    def new_weight(match):
        weight = match.group(2)
        return f'{match.group(1)}"({weight}) && {FLASH} && ({SPLIT} == {data_part})"'

    text = re.sub(r"<MCFile[^>]*>", header, text, count=1)
    text = re.sub(r'(additional_weight\s*=\s*)"([^"]*)"', new_weight, text)
    return re.sub(r'(<detector name="SBND"\s+pot=")[^"]*"', rf'\g<1>{DETECTOR_POT:.6e}"', text)


def link_for(path: Path) -> str:
    """Name of the work/ symlink for a data file (created on demand)."""
    link = REPO / "work" / f"split_{path.stem}.root"
    if not link.exists():
        link.symlink_to(path)
    return link.name


def configure(study: str, syst: str, perm: int, variant: str) -> dict:
    """Write the ratio file, MC and pseudo-data XMLs and the ini for one configuration."""
    tag = f"sm_{study}_{syst}_{variant}_p{perm}"
    ratio_part, mc_part, data_part = perm % 3, (perm + 1) % 3, (perm + 2) % 3
    ratio_path = write_ratio_file(study, [syst], tag, ratio_part,
                                  "none" if variant == "unsmoothed" else variant)
    ref_file, ref_tree, var_tree = PAIRS[syst]
    directory = RUN / "xml"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{tag}.xml").write_text(mc_xml(study, ratio_path, mc_part))
    paths = {}
    samples = (("varied", (DETVAR_FILE, var_tree)), ("null", (ref_file, ref_tree)))
    for label, (file, tree) in samples:
        link = link_for(file)
        text = data_xml(study, tree, file, data_part).replace(f"work/{file.name}", f"work/{link}")
        (directory / f"{tag}_{label}_data.xml").write_text(text)
        paths[label] = directory / f"{tag}_{label}_data.xml"
    ini = closure.closure_ini(study).read_text()
    ini = re.sub(r"^xml\s*=.*$", f"xml = {directory / (tag + '.xml')}", ini, flags=re.M)
    ini = re.sub(r"^tag\s*=.*$", f"tag = {tag}", ini, flags=re.M)
    ini = re.sub(r"^log\s*=.*$", f"log = {tag}.log", ini, flags=re.M)
    for label in ("varied", "null"):
        text = re.sub(r"^data\s*=.*$", f"data = {paths[label]}", ini, flags=re.M)
        (RUN / "inis" / f"{tag}_{label}.ini").write_text(text)
    return {"tag": tag, "study": study, "syst": syst}


def process_and_fit(config: dict, label: str, keep_mc: bool) -> tuple[float, dict]:
    """Process (only the data bin when the MC bins exist) and run the global fit."""
    tag = config["tag"]
    ini = RUN / "inis" / f"{tag}_{label}.ini"
    if not keep_mc:
        for path in RUN.glob(f"{tag}_*.bin"):
            path.unlink()
    else:
        (RUN / f"{tag}_data.bin").unlink(missing_ok=True)
    command = [str(closure.PROFIT), "--config", str(ini), "-o", label, "-n", "4",
               "-l", f"logs/{tag}_{label}.log"]
    subprocess.run([*command, "process"], cwd=RUN, capture_output=True, text=True)
    done = subprocess.run([*command, "-s", "1", "global"], cwd=RUN, capture_output=True, text=True)
    (RUN / "logs" / f"{tag}_{label}.out").write_text(done.stdout)
    return closure.parse_fit(done.stdout)


def run(study: str, syst: str, perm: int, variant: str = "unsmoothed") -> pd.DataFrame:
    """Varied pseudo-data and the null for one configuration."""
    config = configure(study, syst, perm, variant)
    shown = closure.plotnames(study)[f"{syst}_sbnd_Run1"]
    rows = []
    for label, keep in (("varied", False), ("null", True)):
        chi2, values = process_and_fit(config, label, keep)
        others = {k: v for k, v in values.items() if k != shown and "#" not in k}
        worst = max(others, key=lambda k: abs(others[k])) if others else ""
        rows.append({"study": study, "syst": syst, "variant": variant, "perm": perm,
                     "data": label, "chi2_total": chi2, "fitted": values.get(shown),
                     "max_other_pull": abs(others.get(worst, 0.0)), "max_other": worst})
    return pd.DataFrame(rows)


def configure_multi(study: str, perm: int, variant: str) -> dict:
    """One ratio file and MC, every systematic's ratio made from the ratio part; many data XMLs."""
    tag = f"smm_{study}_{variant}_p{perm}"
    ratio_part, mc_part, data_part = perm % 3, (perm + 1) % 3, (perm + 2) % 3
    ratio_path = write_ratio_file(study, list(PAIRS), tag, ratio_part,
                                  "none" if variant == "unsmoothed" else variant)
    directory = RUN / "xml"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{tag}.xml").write_text(mc_xml(study, ratio_path, mc_part))
    samples = {f"varied_{syst}": (DETVAR_FILE, PAIRS[syst][2]) for syst in PAIRS}
    samples["null_nominal"] = (NOMINAL_FILE, "sbnd_large")
    samples["null_dent"] = (DETVAR_FILE, PAIRS["DENT"][1])
    ini = closure.closure_ini(study).read_text()
    ini = re.sub(r"^xml\s*=.*$", f"xml = {directory / (tag + '.xml')}", ini, flags=re.M)
    ini = re.sub(r"^tag\s*=.*$", f"tag = {tag}", ini, flags=re.M)
    ini = re.sub(r"^log\s*=.*$", f"log = {tag}.log", ini, flags=re.M)
    for label, (file, tree) in samples.items():
        link = link_for(file)
        text = data_xml(study, tree, file, data_part).replace(f"work/{file.name}", f"work/{link}")
        path = directory / f"{tag}_{label}_data.xml"
        path.write_text(text)
        (RUN / "inis" / f"{tag}_{label}.ini").write_text(
            re.sub(r"^data\s*=.*$", f"data = {path}", ini, flags=re.M))
    return {"tag": tag, "study": study, "labels": list(samples)}


def run_multi(study: str, perm: int, variant: str) -> pd.DataFrame:
    """Process once, then fit every pseudo-data sample; records every fitted pull."""
    config = configure_multi(study, perm, variant)
    names = closure.plotnames(study)
    rows = []
    for index, label in enumerate(config["labels"]):
        chi2, values = process_and_fit(config, label, keep_mc=index > 0)
        row = {"study": study, "variant": variant, "perm": perm, "data": label,
               "chi2_total": chi2}
        for syst in PAIRS:
            row[f"fit_{syst}"] = values.get(names[f"{syst}_sbnd_Run1"])
        rows.append(row)
    return pd.DataFrame(rows)


def run_all_multi(study: str, jobs: int, variants: list[str], out: str) -> pd.DataFrame:
    """Every variant and role permutation, ``jobs`` processed configurations at a time."""
    import concurrent.futures as cf

    tasks = [(perm, variant) for variant in variants for perm in range(3)]
    with cf.ThreadPoolExecutor(jobs) as pool:
        futures = {pool.submit(run_multi, study, perm, variant): (perm, variant)
                   for perm, variant in tasks}
        frames = []
        for future in cf.as_completed(futures):
            frames.append(future.result())
            table = pd.concat(frames, ignore_index=True)
            table.to_csv(REPO / "docs" / "detvar_smoothing" / out, index=False)
            print(f"done {futures[future]}", flush=True)
    return table


def run_all(study: str, jobs: int, variants: list[str], out: str) -> pd.DataFrame:
    """Every systematic, permutation and variant, ``jobs`` configurations at a time."""
    import concurrent.futures as cf

    tasks = [(syst, perm, variant) for variant in variants for syst in PAIRS for perm in range(3)]
    with cf.ThreadPoolExecutor(jobs) as pool:
        futures = {pool.submit(run, study, syst, perm, variant): (syst, perm, variant)
                   for syst, perm, variant in tasks}
        frames = []
        for future in cf.as_completed(futures):
            frames.append(future.result())
            table = pd.concat(frames, ignore_index=True)
            table.to_csv(REPO / "docs" / "detvar_smoothing" / out, index=False)
            print(f"done {futures[future]}", flush=True)
    return table


def main() -> None:
    """Run one configuration, or all of them."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("study")
    parser.add_argument("syst", nargs="?", choices=sorted(PAIRS))
    parser.add_argument("--perm", type=int, default=0)
    parser.add_argument("--all", action="store_true", help="all systematics and permutations")
    parser.add_argument("--multi", action="store_true",
                        help="one processed MC per split serves every systematic")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--variants", nargs="+", default=["unsmoothed"],
                        help="unsmoothed or a smoother_bias_variance method name")
    parser.add_argument("--out", help="output CSV name under docs/detvar_smoothing")
    args = parser.parse_args()
    out = args.out or f"closure2_{args.study}.csv"
    if args.multi:
        table = run_all_multi(args.study, args.jobs, args.variants, out)
    elif args.all:
        table = run_all(args.study, args.jobs, args.variants, out)
    else:
        table = run(args.study, args.syst, args.perm, args.variants[0])
    print(table.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
