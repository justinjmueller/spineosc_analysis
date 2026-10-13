#!/usr/bin/env python3
"""Git clean filter: strips outputs and execution counts from a notebook on `git add`, leaving the working-tree notebook (with its outputs) untouched.
Install once per clone with scripts/setup_git_filters.sh (needs only the standard library)."""
import json, sys
nb = json.load(sys.stdin)
for c in nb.get("cells", []):
    if c.get("cell_type") == "code": c["outputs"] = []; c["execution_count"] = None
    c.get("metadata", {}).pop("execution", None)
nb.get("metadata", {}).pop("widgets", None)
json.dump(nb, sys.stdout, indent=1, ensure_ascii=False); sys.stdout.write("\n")
