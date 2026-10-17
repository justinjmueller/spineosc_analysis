#!/usr/bin/env python3
"""Closure test 1 for the detector-systematic ratios: inject one systematic at +1 sigma, fit.

For a configuration (an ini produced by scripts/dent_improvement/build_xmls.py), runs PROfit
``global`` with fake data in which one systematic is injected at +1 sigma (``--inject-systs``), with
every systematic floating, in absolute normalization (no shape-only). Reports the global best-fit
chi-square (which includes the nuisance penalty, so a 1 sigma pull gives about 1), the fitted value
of the injected parameter, and the largest other pull. A no-injection fit is the baseline.

    python scripts/detvar_smoothing/closure.py test1 muon_score_inclusive_s2 \
        [--systs DENT ChargeScale15 ...]

The processed binaries must already exist in runs/smoothing (run PROfit ``process`` on the ini
first; see docs/detvar_smoothing/PLAN.md).
"""
import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "runs" / "smoothing"
PROFIT = REPO.parent / "profit-asymmetric" / "build" / "bin" / "PROfit"
XML = REPO / "xml" / "dent_improvement" / "sbnd_run1"
DEFAULT_SYSTS = ["DENT", "ChargeScale15", "WireModxThetaXW", "MattModYZ"]


def study_file(study: str, suffix: str) -> Path:
    """A study's XML or ini: a generated one under runs/smoothing/studies, else the tracked one."""
    generated = RUN / "studies" / f"{study}{suffix}"
    return generated if generated.exists() else XML / f"{study}{suffix}"


def closure_ini(study: str) -> Path:
    """An ini for ``study`` without shape-only or area-norm, written under runs/smoothing/inis."""
    text = study_file(study, "_plot.ini").read_text()
    text = re.sub(r"^(shapeonly|area-norm)\s*=.*\n", "", text, flags=re.M)
    path = RUN / "inis" / f"closure_{study}.ini"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def plotnames(study: str) -> dict[str, str]:
    """Map each active hist1d systematic name to the name PROfit prints for it."""
    text = re.sub(r"<!--.*?-->", "", study_file(study, ".xml").read_text(), flags=re.S)
    entries = re.findall(r'<allowlist([^>]*type="hist1d"[^>]*)>', text)
    out = {}
    for attrs in entries:
        name = re.search(r'\sname="([^"]*)"', attrs)
        plot = re.search(r'plotname="([^"]*)"', attrs)
        if name and plot:
            out[name.group(1)] = plot.group(1)
    return out


def parse_fit(output: str) -> tuple[float, dict[str, float]]:
    """Global best-fit chi-square and the fitted nuisance values from a PROfit log."""
    chi2 = float(re.search(r"Global Best Fit chi\^2:\s*([-\d.eE+]+)", output).group(1))
    block = output.split("print_global_fit_results || Global Best Fit chi^2:")[-1]
    values = {}
    for line in block.splitlines():
        match = re.match(
            r"\s*\d+ print_global_fit_results \|\| (.+?)\s+:\s+([-\d.eE+]+)\s*$", line)
        if match:
            values[match.group(1)] = float(match.group(2))
    return chi2, values


def run_global(study: str, label: str, inject: dict[str, float]) -> tuple[float, dict[str, float]]:
    """One global fit with fake data and optional injected systematic shifts."""
    command = [str(PROFIT), "--config", str(closure_ini(study)), "-o", label,
               "-l", f"logs/{study}_{label}.log", "-n", "4", "--use-fake-data"]
    for name, shift in inject.items():
        command += ["--inject-systs", name, str(shift)]
    done = subprocess.run([*command, "global"], cwd=RUN, capture_output=True, text=True)
    (RUN / "logs").mkdir(exist_ok=True)
    (RUN / "logs" / f"{study}_{label}.out").write_text(done.stdout)
    return parse_fit(done.stdout)


def test1(study: str, systs: list[str]) -> pd.DataFrame:
    """Inject each systematic at +1 sigma and report chi-square and pulls."""
    names = plotnames(study)
    rows = []
    for syst in [None, *systs]:
        key = f"{syst}_sbnd_Run1" if syst else None
        inject = {key: 1.0} if key else {}
        chi2, values = run_global(study, f"cl1_{syst or 'none'}", inject)
        shown = names.get(key) if key else None
        others = {k: v for k, v in values.items()
                  if k != shown and "#" not in k and "(log)" not in k}
        worst = max(others, key=lambda k: abs(others[k])) if others else ""
        rows.append({"study": study, "injected": syst or "none", "chi2_total": chi2,
                     "fitted_injected": values.get(shown) if shown else None,
                     "max_other_pull": abs(others[worst]) if worst else 0.0,
                     "max_other_name": worst})
    return pd.DataFrame(rows)


def main() -> None:
    """Run the requested test."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("test", choices=["test1"])
    parser.add_argument("study")
    parser.add_argument("--systs", nargs="+", default=DEFAULT_SYSTS)
    args = parser.parse_args()
    table = test1(args.study, args.systs)
    out = REPO / "docs" / "detvar_smoothing" / f"test1_{args.study}.csv"
    table.to_csv(out, index=False)
    print(table.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
