#!/usr/bin/env python3
"""Build PROfit XML configurations from blueprints.

A blueprint names a template and supplies the handful of values that differ
between one configuration and its siblings; the template is an existing XML
with those values replaced by ``{{placeholder}}`` markers. The coverage
studies are the motivating case: thirteen 1263-line XMLs that differ in four
values apiece.

The correctness criterion is exact reproduction. ``--check`` renders every
blueprint and compares it byte-for-byte with the committed XML, so the
generator is only ever as trustworthy as that comparison -- it cannot quietly
reformat, reorder or drop anything.

    python3 schematics/build.py --check          # verify against xml/
    python3 schematics/build.py --out /tmp/gen   # render to a directory

Nothing under xml/ is written to. Generating into the tree is a separate
decision from proving the generator reproduces it.

Standard library only (tomllib needs Python 3.11+), matching the rest of the
repo's tooling.
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCHEMATICS = REPO / "schematics"
PLACEHOLDER = re.compile(r"\{\{(\w+)\}\}")
INCLUDE = re.compile(r"^\s*\{\{>\s*(\w+)\s*\}\}\s*$")
IF_OPEN = re.compile(r"^\s*\{\{#if\s+(\w+)\s*\}\}\s*$")
IF_CLOSE = re.compile(r"^\s*\{\{/if\}\}\s*$")


def detector_flags(detectors: list[str]) -> set[str]:
    """Expand a detector list into the flags templates may test.

    ``ICARUS`` is derived rather than listed, for the prose and comments that
    belong to the far detector generally rather than to one of its runs.

    Parameters
    ----------
    detectors : list of str
        The detectors a configuration includes.

    Returns
    -------
    set of str
        The flags `{{#if}}` may name.
    """
    flags = set(detectors)
    if flags & {"ICARUSRun2", "ICARUSRun4"}:
        flags.add("ICARUS")
    return flags


def apply_conditionals(lines: list[str], flags: set[str], source: str) -> list[str]:
    """Drop `{{#if TOKEN}}` regions whose token is not active.

    The enclosed lines are emitted verbatim when kept, which is the point:
    the XMLs carry hand-tuned column alignment that no formatting rule
    reproduces, so the generator must never rewrite a line it keeps.

    Nesting is refused. The moment a configuration needs more than
    "is this detector present", it wants the raw escape hatch rather than a
    template language.

    Parameters
    ----------
    lines : list of str
        Template lines.
    flags : set of str
        Active flags.
    source : str
        Template name, for error messages.

    Returns
    -------
    list of str
        Lines with inactive regions removed and the markers stripped.
    """
    out: list[str] = []
    active: str | None = None
    keeping = True
    for number, line in enumerate(lines, 1):
        opened = IF_OPEN.match(line)
        if opened:
            if active is not None:
                sys.exit(f"{source}:{number}: nested {{{{#if}}}} is not supported")
            active = opened.group(1)
            keeping = active in flags
            continue
        if IF_CLOSE.match(line):
            if active is None:
                sys.exit(f"{source}:{number}: {{{{/if}}}} without a matching {{{{#if}}}}")
            active, keeping = None, True
            continue
        if keeping:
            out.append(line)
    if active is not None:
        sys.exit(f"{source}: unclosed {{{{#if {active}}}}}")
    return out


def expand_includes(lines: list[str], source: str, depth: int = 0) -> list[str]:
    """Splice `{{> partial}}` lines in from schematics/partials.

    Parameters
    ----------
    lines : list of str
        Template lines.
    source : str
        Template name, for error messages.
    depth : int, optional
        Recursion guard.

    Returns
    -------
    list of str
        Lines with includes spliced in.
    """
    if depth > 4:
        sys.exit(f"{source}: includes nested too deeply")
    out: list[str] = []
    for line in lines:
        found = INCLUDE.match(line)
        if not found:
            out.append(line)
            continue
        partial = SCHEMATICS / "partials" / f"{found.group(1)}.xml"
        if not partial.exists():
            sys.exit(f"{source}: no partial {partial}")
        out.extend(expand_includes(partial.read_text().splitlines(), partial.name, depth + 1))
    return out


def render(template: str, values: dict[str, str], source: str,
           detectors: list[str] | None = None) -> str:
    """Substitute every ``{{placeholder}}`` in a template.

    Parameters
    ----------
    template : str
        The template text.
    values : dict of str to str
        Values to substitute, keyed by placeholder name.
    source : str
        Blueprint name, for error messages.

    Returns
    -------
    str
        The rendered configuration.

    Raises
    ------
    SystemExit
        If the template names a placeholder the blueprint does not supply, or
        the blueprint supplies a value the template never uses. Both are
        mistakes worth failing on: the first would emit a literal ``{{...}}``
        into an XML, the second usually means a typo'd key that is silently
        doing nothing.
    """
    lines = expand_includes(template.splitlines(), source)
    lines = apply_conditionals(lines, detector_flags(detectors or []), source)
    template = "\n".join(lines) + "\n"

    wanted = set(PLACEHOLDER.findall(template))
    missing = wanted - set(values)
    if missing:
        sys.exit(f"{source}: template needs values for {sorted(missing)}")
    unused = set(values) - wanted
    if unused:
        sys.exit(f"{source}: blueprint sets {sorted(unused)}, which the template never uses")
    return PLACEHOLDER.sub(lambda m: values[m.group(1)], template)


def load_blueprint(path: Path) -> tuple[str, dict[str, str], Path]:
    """Read one blueprint.

    Parameters
    ----------
    path : pathlib.Path
        The blueprint TOML file.

    Returns
    -------
    template : str
        Template text, unrendered.
    values : dict of str to str
        The substitution values, defaults filled in.
    target : pathlib.Path
        The XML this blueprint reproduces.
    detectors : list of str
        The detectors this configuration includes.
    """
    with path.open("rb") as handle:
        blueprint = tomllib.load(handle)

    template_path = SCHEMATICS / "templates" / f"{blueprint['template']}.xml"
    if not template_path.exists():
        sys.exit(f"{path.name}: no template {template_path}")

    values = dict(blueprint["values"])
    detectors = values.pop("detectors", ["SBND", "ICARUSRun2", "ICARUSRun4"])
    # A study's own data companion is the overwhelmingly common case, so it is
    # a default rather than a line in every blueprint. Three committed XMLs
    # reference a different study here and override it; see their blueprints.
    values.setdefault("detector_block_reference", f"{values['study']}_data.xml")
    # Each coverage study has its own detsys histogram file, named after it.
    values.setdefault("detsys", f"detsys_{values['study']}")

    target = REPO / blueprint.get("target", f"xml/coverage/{values['study']}.xml")
    return template_path.read_text(), values, target, detectors


def main() -> int:
    """Render every blueprint, checking against or writing out the XMLs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="compare against the committed XMLs and report differences")
    parser.add_argument("--out", type=Path,
                        help="render into this directory instead of checking")
    parser.add_argument("--blueprints", type=Path, default=SCHEMATICS / "blueprints",
                        help="blueprint directory (default: schematics/blueprints)")
    args = parser.parse_args()
    if not args.check and args.out is None:
        parser.error("give --check or --out")

    blueprints = sorted(args.blueprints.rglob("*.toml"))
    if not blueprints:
        sys.exit(f"no blueprints under {args.blueprints}")

    failures = 0
    for path in blueprints:
        template, values, target, detectors = load_blueprint(path)
        rendered = render(template, values, path.name, detectors)

        if args.out is not None:
            destination = args.out / target.relative_to(REPO)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(rendered)
            print(f"wrote {destination}")
            continue

        if not target.exists():
            print(f"MISSING  {target.relative_to(REPO)} (nothing to compare against)")
            failures += 1
            continue

        original = target.read_text()
        if rendered == original:
            print(f"ok       {target.relative_to(REPO)}")
            continue

        failures += 1
        diff = list(difflib.unified_diff(
            rendered.splitlines(), original.splitlines(),
            "generated", str(target.relative_to(REPO)), n=0, lineterm=""))
        print(f"DIFFERS  {target.relative_to(REPO)}")
        for line in diff[2:12]:
            print(f"         {line}")

    if args.check:
        print(f"\n{len(blueprints) - failures}/{len(blueprints)} reproduced exactly")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
