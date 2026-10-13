#!/usr/bin/env python3
"""Classify and check the analysis XMLs against xml/MANIFEST.toml.
  python scripts/xml_manifest.py            summary per status and per family
  python scripts/xml_manifest.py --list S   list files whose status is S
  python scripts/xml_manifest.py --check    exit 1 if any XML is unclassified or a family has no files / a missing producer
Matching: fnmatch on the repo-relative path, most specific (longest) pattern wins. Covers *.xml and *.ini."""
import fnmatch, pathlib, sys, tomllib
from collections import Counter, defaultdict
R = pathlib.Path(__file__).resolve().parents[1]
fam = tomllib.loads((R / "xml" / "MANIFEST.toml").read_text())["family"]
files = sorted(str(p.relative_to(R)) for p in (R / "xml").rglob("*") if p.suffix in (".xml", ".ini"))
def match(f):
    hits = [x for x in fam if fnmatch.fnmatch(f, x["pattern"])]
    return max(hits, key=lambda x: len(x["pattern"].replace("*", ""))) if hits else None
by = defaultdict(list); un = []
for f in files:
    m = match(f)
    (by[m["pattern"]].append(f) if m else un.append(f))
if "--list" in sys.argv:
    s = sys.argv[sys.argv.index("--list") + 1]
    for x in fam:
        if x["status"] == s: print(f"{x['pattern']}  ({len(by[x['pattern']])} files)")
    sys.exit(0)
bad = []
st = Counter()
for x in fam:
    n = len(by[x["pattern"]]); st[x["status"]] += n
    if n == 0: bad.append(f"family without files: {x['pattern']}")
    for p in str(x.get("producer", "")).split(", "):
        p = p.split(" ")[0]
        if p.endswith((".py", ".sh", ".C")) and not (R / p).exists(): bad.append(f"missing producer {p} for {x['pattern']}")
    if x["status"] == "needs-review": pass
print(f"{len(files)} xml/ini files, {len(files) - len(un)} classified, {len(un)} unclassified")
print("files per status:", dict(st))
for x in fam: print(f"  {x['status']:12s} {len(by[x['pattern']]):4d}  {x['pattern']}")
if un: print("unclassified (first 10):", un[:10])
for b in bad: print("WARN", b)
if "--check" in sys.argv and (un or bad): sys.exit(1)
