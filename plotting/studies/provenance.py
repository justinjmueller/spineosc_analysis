"""
Figure provenance: a sidecar JSON next to each figure that records how it was made.

``record(figure, inputs, script, arguments)`` writes ``<figure>.prov.json`` holding the figure's
hash, the producing script (path, hash, git commit), the command line, the input files (size,
modification time and, for small files, a hash) and the selection and code versions.
``scripts/figure_provenance.py`` reads these back to list the figures and to flag the stale ones.

Inputs that are PROfit outputs (``<tag>_v1_PROplot.root``, ``<tag>_v1_surf.root``, ...) get a ``profit``
entry as well, found from the file name: the ini with ``tag = <tag>``, the XML it names (hash, git commit,
and the content itself when it is not committed), the ``work/<key>.root`` files the XML reads resolved to
real paths, the cache files, the PROfit version from the logs and the run commands noted in the ini header.
``restore()`` writes the recorded ini and XML back out so the figure can be remade from the sidecar alone.

Large inputs (the release ntuples run to GB) are identified by size and modification time; files
up to ``HASH_LIMIT`` bytes are also hashed, so touching a small file without changing it is not
reported as a change.
"""

from __future__ import annotations

import base64
import datetime
import functools
import gzip
import hashlib
import json
import re
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


def describe_input(path: Path, repo: Path = REPO) -> dict:
    """Identify an input: a file (size, mtime and, when small, hash) or a directory (count, size, newest mtime)."""
    path = Path(path)
    info = {"path": _relative(path, repo)}
    if not path.exists():
        return {**info, "missing": True}
    if path.is_dir():
        return {**info, **_directory_state(path)}
    stat = path.stat()
    info.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if stat.st_size <= HASH_LIMIT:
        info["sha256"] = sha256(path)
    return info


def _directory_state(path: Path) -> dict:
    """A directory input (e.g. hundreds of scan files): file count, total size and newest mtime."""
    files = [p for p in path.rglob("*") if p.is_file()]
    stats = [p.stat() for p in files]
    return {"directory": True, "n_files": len(files), "size": sum(s.st_size for s in stats),
            "mtime_ns": max((s.st_mtime_ns for s in stats), default=0)}


def _relative(path: Path, repo: Path = REPO) -> str:
    path = Path(path).resolve()
    try:
        return str(path.relative_to(Path(repo).resolve()))
    except ValueError:
        return str(path)


def _git(*args: str, repo: Path = REPO) -> str:
    try:
        out = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return ""
    return out.stdout.strip()



# ----------------------------------------------------------------------
# PROfit recipe: ini, XML, work/ file map, caches
# ----------------------------------------------------------------------
_SUFFIX = re.compile(r"_(?:PROplot|surf|fractional|ratio)")
_WORK = re.compile(r"work/([\w.+-]+\.root)")


@functools.lru_cache(maxsize=4)
def _ini_index(repo: Path) -> dict:
    """Map tag -> [ini paths] over jobs/ and xml/ (archives excluded)."""
    index: dict[str, list[Path]] = {}
    for top in ("jobs", "xml"):
        for ini in sorted((repo / top).rglob("*.ini")):
            if "archive" in ini.parts:
                continue
            for line in ini.read_text(errors="replace").splitlines():
                if line.lstrip().startswith("["):
                    break
                match = re.match(r"\s*tag\s*=\s*(\S+)", line)
                if match:
                    index.setdefault(match.group(1), []).append(ini)
    return index


def _ini_key(text: str, key: str) -> str:
    for line in text.splitlines():
        if line.lstrip().startswith("["):
            break
        match = re.match(rf"\s*{key}\s*=\s*(\S+)", line)
        if match:
            return match.group(1)
    return ""



# ----------------------------------------------------------------------
# xml/MANIFEST.toml: the family (status, purpose, technote section) of an XML or ini
# ----------------------------------------------------------------------
#: Manifest fields copied into a sidecar.
MANIFEST_FIELDS = ("pattern", "status", "purpose", "consumer", "technote", "regenerate", "superseded_by", "reviewed")


@functools.lru_cache(maxsize=4)
def _manifest(repo: Path):
    """(``match`` function, families) of ``xml/MANIFEST.toml``, or ``None`` if the repository has none."""
    manifest = Path(repo) / "xml" / "MANIFEST.toml"
    if not manifest.is_file():
        return None
    import importlib.util

    spec = importlib.util.spec_from_file_location("xml_manifest", REPO / "scripts" / "xml_manifest.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.match, module.load_families(Path(repo))


def manifest_family(path: str, repo: Path = REPO) -> dict:
    """The manifest family of a repository-relative path (only the fields in ``MANIFEST_FIELDS``), or ``{}``."""
    loaded = _manifest(Path(repo))
    if loaded is None:
        return {}
    family, _tie = loaded[0](path, loaded[1])
    return {key: family[key] for key in MANIFEST_FIELDS if family and key in family}


def _file_state(path: Path, repo: Path, embed_if_unversioned: bool) -> dict:
    """A configuration file: hash, git commit, and the content itself when git cannot restore it."""
    rel = _relative(path, repo)
    info = describe_input(path, repo)
    tracked = _git("ls-files", "--error-unmatch", rel, repo=repo) != ""
    clean = tracked and _git("status", "--porcelain", "--", rel, repo=repo) == ""
    info.update(git_commit=_git("log", "-1", "--format=%h", "--", rel, repo=repo) if tracked else "", git_clean=clean,
                manifest=manifest_family(rel, repo))
    if embed_if_unversioned and not clean:
        info["content_gz_b64"] = base64.b64encode(gzip.compress(path.read_bytes())).decode()
    return info


def _select_inis(inis: list[Path], remainder: str) -> list[Path]:
    """
    Pick the ini(s) of a tag that fit the rest of the output file name.

    Several inis can share a tag (a plain and a ``phase2`` plot ini, for example). The file name then carries the
    output tag (``..._cvrw_phase2_PROplot...``); keep the inis whose distinguishing name parts all occur in it, and of
    those the most specific.
    """
    if len(inis) <= 1:
        return list(inis)
    parts = {ini: set(ini.stem.split("_")) for ini in inis}
    common = set.intersection(*parts.values())
    wanted = set(remainder.split("_"))
    fits = {ini: parts[ini] - common for ini in inis if parts[ini] - common <= wanted}
    if not fits:
        return list(inis)
    best = max(len(extra) for extra in fits.values())
    return [ini for ini, extra in fits.items() if len(extra) == best]


def profit_recipe(output_file: Path, repo: Path = REPO) -> list[dict]:
    """
    Describe how a PROfit output file was made, or ``[]`` if it does not look like one.

    PROfit names its outputs ``<tag>_<output tag>_PROplot*.root`` (or ``_surf``); the longest ini ``tag`` that
    prefixes the file name is used, and the inis in ``jobs/`` and ``xml/`` with that tag are narrowed by the
    output tag (a ``phase2`` ini is only used for a ``..._phase2_...`` file).
    """
    repo, output_file = Path(repo), Path(output_file)
    name = output_file.name
    tags = [t for t in _ini_index(repo) if name.startswith(t + "_")]
    if not tags:
        return []
    tag = max(tags, key=len)
    inis = _select_inis(_ini_index(repo)[tag], _SUFFIX.split(name[len(tag) + 1:])[0])
    entries = []
    for ini in inis:
        text = ini.read_text(errors="replace")
        xml = repo / _ini_key(text, "xml")
        if not xml.is_file():
            continue
        keys = sorted(set(_WORK.findall(xml.read_text(errors="replace"))))
        work = []
        for key in keys:
            link = repo / "work" / key
            real = link.resolve() if link.exists() else link
            item = describe_input(real, repo)
            item["work_key"] = key
            upstream = real.with_name(real.name + SUFFIX)
            if upstream.is_file():
                item["upstream"] = json.loads(upstream.read_text())
            work.append(item)
        caches = [describe_input(output_file.parent / f"{tag}_{kind}.bin", repo) for kind in ("data", "prop", "syst")]
        logs = sorted(output_file.parent.glob(f"{tag}*.log"))
        version = ""
        for log in logs:
            found = re.search(r"PROfit version (v[\d.]+)", log.read_text(errors="replace")[:20000])
            version = found.group(1) if found else version
        hints = [m.group(0) for m in re.finditer(r"(?m)^#\s+\S*(?:PROfit|profit_run)\S*.*$", text)]
        entries.append({
            "tag": tag,
            "output": _relative(output_file, repo),
            "ini": {**_file_state(ini, repo, embed_if_unversioned=False), "content": text},
            "xml": _file_state(xml, repo, embed_if_unversioned=True),
            "work_files": work,
            "caches": caches,
            "profit_version": version,
            "run_hints": [h.lstrip("# ").strip() for h in hints],
        })
    return entries


def restore(sidecar: Path, destination: Path, repo: Path = REPO) -> list[Path]:
    """
    Write the recorded ini and XML of a figure under ``destination`` (at their repository paths).

    The ini comes from the sidecar; the XML from the embedded copy, or from git at the recorded commit.
    """
    rec = json.loads(Path(sidecar).read_text())
    written = []
    for entry in rec.get("profit", []):
        for kind in ("ini", "xml"):
            state = entry[kind]
            target = Path(destination) / state["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            if kind == "ini":
                target.write_text(state["content"])
            elif "content_gz_b64" in state:
                target.write_bytes(gzip.decompress(base64.b64decode(state["content_gz_b64"])))
            else:
                blob = subprocess.run(["git", "-C", str(repo), "show", f"{state['git_commit']}:{state['path']}"],
                                      capture_output=True, check=True).stdout
                target.write_bytes(blob)
            written.append(target)
    return written


def sidecar_path(figure: Path) -> Path:
    """Sidecar file of a figure."""
    figure = Path(figure)
    return figure.with_name(figure.name + SUFFIX)


def record(figure, inputs, script, arguments=None, versions=None, note: str = "", repo: Path = REPO, recipe: bool = True) -> Path:
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
    repo : path-like, optional
        Repository root that relative paths refer to.
    recipe : bool, optional
        Trace PROfit outputs among the inputs back to their ini and XML; turn off for figures that read
        many PROfit files, whose recipes are in the per-panel sidecars.

    Returns
    -------
    Path
        The sidecar file.
    """
    figure, script = Path(figure), Path(script).resolve()
    rel_script = _relative(script, repo)
    record_ = {
        "schema": SCHEMA,
        "figure": {"name": figure.name, "sha256": sha256(figure), "size": figure.stat().st_size},
        "script": {
            "path": rel_script,
            "sha256": sha256(script),
            "git_commit": _git("rev-parse", "--short", "HEAD", repo=repo),
            "git_clean": _git("status", "--porcelain", "--", rel_script, repo=repo) == "",
        },
        "command": " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
        "arguments": arguments or {},
        "versions": versions or {},
        "inputs": [describe_input(p, repo) for p in inputs],
        "profit": [entry for p in inputs for entry in profit_recipe(Path(p), repo)] if recipe else [],
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "note": note,
    }
    out = sidecar_path(figure)
    out.write_text(json.dumps(record_, indent=1) + "\n")
    return out


def check(sidecar: Path, repo: Path = REPO) -> tuple[str, list[str]]:
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
    script = (repo / rec["script"]["path"]) if not Path(rec["script"]["path"]).is_absolute() else Path(rec["script"]["path"])
    if not script.exists():
        worse("MISSING")
        reasons.append(f"script {rec['script']['path']} is gone")
    elif sha256(script) != rec["script"]["sha256"]:
        worse("STALE")
        reasons.append(f"script {rec['script']['path']} changed since the figure was made")
    for item in rec["inputs"]:
        path = Path(item["path"]) if Path(item["path"]).is_absolute() else repo / item["path"]
        if item.get("missing") or not path.exists():
            worse("MISSING")
            reasons.append(f"input {item['path']} is missing")
            continue
        stat = path.stat()
        if item.get("directory"):
            now = _directory_state(path)
            changed = (now["n_files"], now["size"], now["mtime_ns"]) != (item["n_files"], item["size"], item["mtime_ns"])
        elif "sha256" in item:
            changed = sha256(path) != item["sha256"]
        else:
            changed = (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"])
        if changed:
            worse("STALE")
            reasons.append(f"input {item['path']} changed since the figure was made")
    for entry in rec.get("profit", []):
        for kind in ("ini", "xml"):
            state = entry[kind]
            path = repo / state["path"]
            if not path.exists():
                worse("MISSING")
                reasons.append(f"{kind} {state['path']} is gone (restore it from the sidecar)")
            elif sha256(path) != state["sha256"]:
                worse("STALE")
                reasons.append(f"{kind} {state['path']} changed since the figure was made")
        family = manifest_family(entry["xml"]["path"], repo)
        if family.get("status") in ("archived", "deprecated"):
            worse("STALE")
            why = f"; superseded by {family['superseded_by']}" if "superseded_by" in family else ""
            reasons.append(f"manifest: {entry['xml']['path']} is now {family['status']} ({family['pattern']}){why}")
        recorded = entry["xml"].get("manifest", {})
        if recorded and family and recorded.get("status") != family.get("status"):
            reasons.append(f"note: manifest status of {entry['xml']['path']} was {recorded.get('status')}, now {family.get('status')}")
        consumers = [c.strip() for c in family.get("consumer", "").split(",") if c.strip()]
        if consumers and rec["script"]["path"] not in consumers:
            reasons.append(f"note: manifest consumer of {family['pattern']} is {', '.join(consumers)}, "
                           f"but this figure was made by {rec['script']['path']}")
        for item in entry["work_files"]:
            path = Path(item["path"]) if Path(item["path"]).is_absolute() else repo / item["path"]
            if not path.exists():
                worse("MISSING")
                reasons.append(f"work/{item['work_key']} -> {item['path']} is missing")
            elif "sha256" in item and sha256(path) != item["sha256"]:
                worse("STALE")
                reasons.append(f"work/{item['work_key']} ({item['path']}) changed")
            elif "sha256" not in item and (path.stat().st_size, path.stat().st_mtime_ns) != (item["size"], item["mtime_ns"]):
                worse("STALE")
                reasons.append(f"work/{item['work_key']} ({item['path']}) changed")
        config_mtime = max((repo / entry[k]["path"]).stat().st_mtime_ns for k in ("ini", "xml") if (repo / entry[k]["path"]).exists())
        for cache in entry["caches"]:
            path = Path(cache["path"]) if Path(cache["path"]).is_absolute() else repo / cache["path"]
            if path.exists() and path.stat().st_mtime_ns < config_mtime:
                # PROfit's cache key covers the XML but not the ini flags, so this is only a hint.
                reasons.append(f"note: cache {cache['path']} is older than the ini or XML (check that the cache matches)")
    return status, list(dict.fromkeys(reasons))


def remake_steps(sidecar: Path) -> list[str]:
    """The steps, in order, that remake a figure, as the sidecar knows them."""
    rec = json.loads(Path(sidecar).read_text())
    steps = [f"use the repository at commit {rec['script']['git_commit']}   # script and helpers as they were (script clean when run: {rec['script']['git_clean']})"]
    for entry in rec.get("profit", []):
        steps.append(f"scripts/figure_provenance.py restore {Path(sidecar).name} DIR   # writes {entry['ini']['path']} and {entry['xml']['path']}")
        regenerate = entry["xml"].get("manifest", {}).get("regenerate")
        if regenerate:
            steps.append(f"or regenerate the XML (manifest {entry['xml']['manifest']['pattern']}): {regenerate}")
        keys = ", ".join(w["work_key"] for w in entry["work_files"])
        steps.append(f"scripts/deploy.sh   # work/ must map: {keys}")
        for hint in entry["run_hints"] or [f"PROfit --config {entry['ini']['path']} process; PROfit --config {entry['ini']['path']} plot"]:
            steps.append(f"{hint}   # tag {entry['tag']}, PROfit {entry['profit_version'] or '?'}")
    steps.append(f"python {rec['script']['path']} {rec['command'].split(' ', 1)[1] if ' ' in rec['command'] else ''}".rstrip())
    return steps
