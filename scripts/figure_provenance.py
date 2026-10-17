#!/usr/bin/env python3
"""List the figures that carry a provenance sidecar and flag the stale ones.

    python scripts/figure_provenance.py list  figures/ docs/
    python scripts/figure_provenance.py check figures/ --technote ../spineosc_technote
    python scripts/figure_provenance.py show  figures/x.pdf.prov.json       # the steps that remake it
    python scripts/figure_provenance.py restore figures/x.pdf.prov.json DIR  # write the recorded ini and XML

``check`` exits non-zero if any figure is MODIFIED, STALE or MISSING (see plotting/studies/provenance.py).
With ``--technote`` it also finds the figure of the same name in the technote tree and says whether that
copy is the recorded one, an older or edited one, or absent.
"""

import argparse
import json
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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("list", "check", "show", "restore"))
    parser.add_argument("paths", nargs="+", help="directories (searched recursively) or sidecar files; for restore: SIDECAR DIR")
    parser.add_argument("--technote", type=Path, help="technote repository, to compare its copies of the figures")
    args = parser.parse_args(argv)

    if args.command == "restore":
        for path in provenance.restore(Path(args.paths[0]), Path(args.paths[1])):
            print(f"wrote {path}")
        return 0
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
