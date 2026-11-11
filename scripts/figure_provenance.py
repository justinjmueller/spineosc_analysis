#!/usr/bin/env python3
"""List the figures that carry a provenance sidecar and flag the stale ones.

    python scripts/figure_provenance.py list  figures/ docs/
    python scripts/figure_provenance.py check figures/ --technote ../spineosc_technote
    python scripts/figure_provenance.py coverage docs/ figures/ --technote ../spineosc_technote
    python scripts/figure_provenance.py manifest docs/ figures/ --technote ../spineosc_technote
    python scripts/figure_provenance.py show  figures/x.pdf.prov.json       # the steps that remake it
    python scripts/figure_provenance.py restore figures/x.pdf.prov.json DIR  # write the recorded ini and XML

``check`` exits non-zero if any figure is MODIFIED, STALE or MISSING (see plotting/studies/provenance.py).
With ``--technote`` it also finds the figure of the same name in the technote tree and says whether that
copy is the recorded one, an older or edited one, or absent.
"""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "plotting" / "studies"))
import provenance  # noqa: E402


def find_sidecars(paths):
    for path in paths:
        path = Path(path)
        yield from ([path] if path.is_file() else sorted(path.rglob("*" + provenance.SUFFIX)))


def technote_status(technote: Path, rec: dict) -> str:
    copies = [p for p in (technote / "figures").rglob(rec["figure"]["name"])] if technote else []
    if not copies:
        return "not in technote"
    same = [c for c in copies if provenance.sha256(c) == rec["figure"]["sha256"]]
    return "technote: current" if same else "technote: DIFFERS"


INCLUDE = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}")


def technote_figures(technote: Path) -> dict:
    """Figures the technote includes: LaTeX argument -> files that satisfy it (searched under figures/)."""
    found = {}
    for tex in sorted(technote.rglob("*.tex")):
        for arg in INCLUDE.findall(tex.read_text(errors="replace")):
            name = Path(arg).name
            pattern = name if Path(name).suffix else name + ".*"
            found[arg] = sorted((technote / "figures").rglob(pattern)) or sorted(technote.rglob(pattern))
    return found


def coverage(technote: Path, search: list[str]) -> int:
    """Report which of the technote's figures have a provenance sidecar, and whether the copies match."""
    by_name = {}
    for sidecar in find_sidecars(search):
        rec = json.loads(sidecar.read_text())
        by_name.setdefault(rec["figure"]["name"], (sidecar, rec))
    rows, used = [], set()
    for arg, copies in technote_figures(technote).items():
        name = Path(arg).name if Path(arg).suffix else (copies[0].name if copies else Path(arg).name)
        if name not in by_name:
            rows.append(("no provenance", arg))
            continue
        sidecar, rec = by_name[name]
        used.add(name)
        status, _ = provenance.check(sidecar)
        current = any(provenance.sha256(c) == rec["figure"]["sha256"] for c in copies)
        rows.append((f"{status}, " + ("technote copy is the recorded one" if current else "technote copy differs"), arg))
    summary = {}
    for state, _ in rows:
        summary[state] = summary.get(state, 0) + 1
    print(f"{len(rows)} technote figures; {len(used)} of the {len(by_name)} figures with provenance are used")
    for state, n in sorted(summary.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {state}")
    missing = [arg for state, arg in rows if state == "no provenance"]
    print("\nwithout provenance, by folder:")
    folders = {}
    for arg in missing:
        folders[str(Path(arg).parent)] = folders.get(str(Path(arg).parent), 0) + 1
    for folder, n in sorted(folders.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {folder}")
    return 0


def manifest_report(technote: Path | None, search: list[str]) -> int:
    """Per manifest family: its status and technote section, the figures with provenance that come from it."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import xml_manifest

    families = xml_manifest.load_families()
    counts = {}
    for sidecar in find_sidecars(search):
        rec = json.loads(sidecar.read_text())
        status, _ = provenance.check(sidecar)
        patterns = {e["xml"].get("manifest", {}).get("pattern") for e in rec.get("profit", [])} - {None}
        for pattern in patterns:  # once per figure, however many inis share the tag
            bucket = counts.setdefault(pattern, {"OK": 0, "other": 0})
            bucket["OK" if status == "OK" else "other"] += 1
    print(f"{'status':10s} {'technote section':30s} {'figs OK':>7s} {'not OK':>7s} {'section figs':>12s}  family")
    for fam in families:
        section = fam.get("technote", "")
        have = counts.get(fam["pattern"], {"OK": 0, "other": 0})
        total = ""
        folders = [re.sub(r"\s*\(.*\)", "", part).strip() for part in section.split(",")]
        if technote and any((technote / "figures" / f).is_dir() for f in folders if f):
            total = str(sum(len([p for p in (technote / "figures" / f).rglob("*") if p.suffix in (".pdf", ".png")])
                            for f in folders if f and (technote / "figures" / f).is_dir()))
        print(f"{fam['status']:10s} {section[:30]:30s} {have['OK']:7d} {have['other']:7d} {total:>12s}  {fam['pattern']}"
              + ("  (proposed)" if "(proposed)" in fam.get("reviewed", "") else ""))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("list", "check", "coverage", "manifest", "show", "restore"))
    parser.add_argument("paths", nargs="+", help="directories (searched recursively) or sidecar files; for restore: SIDECAR DIR")
    parser.add_argument("--technote", type=Path, help="technote repository, to compare its copies of the figures")
    args = parser.parse_args(argv)

    if args.command == "restore":
        for path in provenance.restore(Path(args.paths[0]), Path(args.paths[1])):
            print(f"wrote {path}")
        return 0
    if args.command == "manifest":
        return manifest_report(args.technote, args.paths)
    if args.command == "coverage":
        if not args.technote:
            parser.error("coverage needs --technote")
        return coverage(args.technote, args.paths)
    if args.command == "show":
        for sidecar in args.paths:
            print(f"{sidecar}:")
            for i, step in enumerate(provenance.remake_steps(Path(sidecar)), 1):
                print(f"  {i}. {step}")
        return 0

    bad = 0
    for sidecar in find_sidecars(args.paths):
        rec = json.loads(sidecar.read_text())
        status, reasons = provenance.check(sidecar)
        line = f"{status:9s} {rec['figure']['name']:62s} {rec['script']['path']}"
        if args.technote:
            line += f"  [{technote_status(args.technote, rec)}]"
        print(line)
        if args.command == "check":
            for reason in reasons:
                print(f"          - {reason}")
        elif args.command == "list":
            print(f"          inputs: {', '.join(i['path'].rsplit('/', 1)[-1] for i in rec['inputs'])};"
                  f" args: {json.dumps(rec['arguments'])}")
        bad += status != "OK"
    return 1 if args.command == "check" and bad else 0


if __name__ == "__main__":
    sys.exit(main())
