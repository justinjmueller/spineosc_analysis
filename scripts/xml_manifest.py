#!/usr/bin/env python3
"""
Classify and check the analysis XMLs against ``xml/MANIFEST.toml``.

Usage::

    python scripts/xml_manifest.py            summary per status and per family
    python scripts/xml_manifest.py --list S   list families whose status is S
    python scripts/xml_manifest.py --check    exit 1 if a file is unclassified, a family has no
                                              files, a producer is missing, or two patterns tie

Matching uses ``fnmatch`` on the repository-relative path; the most specific pattern (the longest
once ``*`` is removed) wins. Covers ``*.xml`` and ``*.ini``.
"""

from __future__ import annotations

import argparse
import fnmatch
import sys
import tomllib
from collections import Counter, defaultdict
from pathlib import Path

#: Repository root.
REPO = Path(__file__).resolve().parents[1]

#: File suffixes that are classified.
SUFFIXES = (".xml", ".ini")

#: Producer suffixes whose existence is checked.
PRODUCER_SUFFIXES = (".py", ".sh", ".C")


def load_families(repo: Path = REPO) -> list[dict]:
    """
    Read the ``[[family]]`` entries of ``xml/MANIFEST.toml``.

    Parameters
    ----------
    repo : Path, optional
        Repository root.

    Returns
    -------
    list[dict]
        One dictionary per family.
    """
    return tomllib.loads((repo / "xml" / "MANIFEST.toml").read_text())["family"]


def specificity(pattern: str) -> int:
    """Length of a pattern with the wildcard characters removed."""
    return len(pattern.replace("*", ""))


def match(path: str, families: list[dict]) -> tuple[dict | None, bool]:
    """
    Find the most specific family matching a path.

    Parameters
    ----------
    path : str
        Repository-relative file path.
    families : list[dict]
        Families from :func:`load_families`.

    Returns
    -------
    family : dict or None
        The winning family, or None if nothing matches.
    tie : bool
        True if two different patterns match with equal specificity.
    """
    hits = [f for f in families if fnmatch.fnmatch(path, f["pattern"])]
    if not hits:
        return None, False
    best = max(specificity(f["pattern"]) for f in hits)
    top = [f for f in hits if specificity(f["pattern"]) == best]
    return top[0], len(top) > 1


def classify(repo: Path = REPO) -> dict:
    """
    Classify every XML/ini file under ``xml/``.

    Returns
    -------
    dict
        ``families``, ``files``, ``by`` (pattern -> files), ``unclassified`` and ``ties``.
    """
    families = load_families(repo)
    files = sorted(
        str(p.relative_to(repo)) for p in (repo / "xml").rglob("*") if p.suffix in SUFFIXES
    )
    by: dict[str, list[str]] = defaultdict(list)
    unclassified, ties = [], []
    for path in files:
        family, tie = match(path, families)
        if family is None:
            unclassified.append(path)
            continue
        by[family["pattern"]].append(path)
        if tie:
            ties.append(path)
    return {"families": families, "files": files, "by": by,
            "unclassified": unclassified, "ties": ties}


def problems(result: dict, repo: Path = REPO) -> list[str]:
    """List the problems that make ``--check`` fail."""
    bad = [f"unclassified: {p}" for p in result["unclassified"]]
    bad += [f"pattern tie for {p}" for p in result["ties"]]
    for fam in result["families"]:
        if not result["by"][fam["pattern"]]:
            bad.append(f"family without files: {fam['pattern']}")
        for producer in str(fam.get("producer", "")).split(","):
            name = producer.strip().split(" ")[0]
            if name.endswith(PRODUCER_SUFFIXES) and not (repo / name).exists():
                bad.append(f"missing producer {name} for {fam['pattern']}")
    return bad


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit code."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--list", metavar="STATUS", help="list families with this status")
    parser.add_argument("--check", action="store_true", help="exit 1 on any problem")
    args = parser.parse_args(argv)

    result = classify()
    families, by = result["families"], result["by"]
    if args.list:
        for fam in families:
            if fam["status"] == args.list:
                print(f"{fam['pattern']}  ({len(by[fam['pattern']])} files)")
        return 0

    status = Counter()
    for fam in families:
        status[fam["status"]] += len(by[fam["pattern"]])
    n_files, n_un = len(result["files"]), len(result["unclassified"])
    print(f"{n_files} xml/ini files, {n_files - n_un} classified, {n_un} unclassified")
    print("files per status:", dict(status))
    for fam in families:
        print(f"  {fam['status']:12s} {len(by[fam['pattern']]):4d}  {fam['pattern']}")
    bad = problems(result)
    for line in bad:
        print("WARN", line)
    return 1 if args.check and bad else 0


if __name__ == "__main__":
    sys.exit(main())
