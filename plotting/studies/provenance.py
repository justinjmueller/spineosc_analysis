"""
Figure provenance: a sidecar JSON next to each figure that records how it was made.

``record(figure, inputs, script, arguments)`` writes ``<figure>.prov.json`` holding the figure's
hash, the producing script (path, hash, git commit), the command line, the input files (size,
modification time and, for small files, a hash) and the selection and code versions.
``scripts/figure_provenance.py`` reads these back to list the figures and to flag the stale ones.

Large inputs (the release ntuples run to GB) are identified by size and modification time; files
up to ``HASH_LIMIT`` bytes are also hashed, so touching a small file without changing it is not
reported as a change.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

#: Inputs up to this size (bytes) are hashed; larger ones are identified by size and mtime.
HASH_LIMIT = 64 * 1024 * 1024

#: Suffix of the sidecar file (appended to the full figure file name).
SUFFIX = ".prov.json"

SCHEMA = 1


def sha256(path: Path) -> str:
    """Hex SHA-256 of a file's contents."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def describe_input(path: Path) -> dict:
    """Identify an input file: path, size, mtime and, when small, the hash."""
    path = Path(path)
    info = {"path": _relative(path)}
    if not path.exists():
        return {**info, "missing": True}
    stat = path.stat()
    info.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if stat.st_size <= HASH_LIMIT:
        info["sha256"] = sha256(path)
    return info


def _relative(path: Path) -> str:
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _git(*args: str) -> str:
    try:
        out = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return ""
    return out.stdout.strip()


def sidecar_path(figure: Path) -> Path:
    """Sidecar file of a figure."""
    figure = Path(figure)
    return figure.with_name(figure.name + SUFFIX)


def record(figure, inputs, script, arguments=None, versions=None, note: str = "") -> Path:
    """
    Write the provenance sidecar of a figure that has just been saved.

    Parameters
    ----------
    figure : path-like
        The figure file (must exist).
    inputs : iterable of path-like
        Every file the figure's numbers were read from.
    script : path-like
        The producing script (normally ``__file__``).
    arguments : dict, optional
        The parameters that select this figure (detector, variable, labels, ...).
    versions : dict, optional
        Selection and code versions stamped on the figure.
    note : str, optional
        Free text, e.g. "external input, not versioned".

    Returns
    -------
    Path
        The sidecar file.
    """
    figure, script = Path(figure), Path(script).resolve()
    rel_script = _relative(script)
    record_ = {
        "schema": SCHEMA,
        "figure": {"name": figure.name, "sha256": sha256(figure), "size": figure.stat().st_size},
        "script": {
            "path": rel_script,
            "sha256": sha256(script),
            "git_commit": _git("rev-parse", "--short", "HEAD"),
            "git_clean": _git("status", "--porcelain", "--", rel_script) == "",
        },
        "command": " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
        "arguments": arguments or {},
        "versions": versions or {},
        "inputs": [describe_input(p) for p in inputs],
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "note": note,
    }
    out = sidecar_path(figure)
    out.write_text(json.dumps(record_, indent=1) + "\n")
    return out


def check(sidecar: Path) -> tuple[str, list[str]]:
    """
    Compare a sidecar with the files it describes.

    Returns
    -------
    (status, reasons)
        ``OK``, or ``MODIFIED`` (the figure file no longer matches), ``STALE`` (the script or an
        input changed after the figure was made), ``MISSING`` (figure or an input is gone).
    """
    sidecar = Path(sidecar)
    rec = json.loads(sidecar.read_text())
    figure = sidecar.with_name(rec["figure"]["name"])
    reasons, status = [], "OK"

    def worse(new: str) -> None:
        nonlocal status
        order = ["OK", "STALE", "MODIFIED", "MISSING"]
        status = max(status, new, key=order.index)

    if not figure.exists():
        return "MISSING", ["figure file is gone"]
    if sha256(figure) != rec["figure"]["sha256"]:
        worse("MODIFIED")
        reasons.append("figure differs from the recorded one")
    script = (REPO / rec["script"]["path"]) if not Path(rec["script"]["path"]).is_absolute() else Path(rec["script"]["path"])
    if not script.exists():
        worse("MISSING")
        reasons.append(f"script {rec['script']['path']} is gone")
    elif sha256(script) != rec["script"]["sha256"]:
        worse("STALE")
        reasons.append(f"script {rec['script']['path']} changed since the figure was made")
    for item in rec["inputs"]:
        path = Path(item["path"]) if Path(item["path"]).is_absolute() else REPO / item["path"]
        if item.get("missing") or not path.exists():
            worse("MISSING")
            reasons.append(f"input {item['path']} is missing")
            continue
        stat = path.stat()
        if "sha256" in item:
            changed = sha256(path) != item["sha256"]
        else:
            changed = (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"])
        if changed:
            worse("STALE")
            reasons.append(f"input {item['path']} changed since the figure was made")
    return status, reasons
