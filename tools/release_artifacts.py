#!/usr/bin/env python3
"""Build once, promote the same bytes: checks on frozen release files (ADR-0040).

Three subcommands, used by ``.github/workflows/release.yml`` and by the manual
``twine`` path in ``.agents/dev-contract.md`` alike::

    python tools/release_artifacts.py check-dist --version X.Y.Z DIST_DIR
    python tools/release_artifacts.py sums DIST_DIR > SHA256SUMS
    python tools/release_artifacts.py verify-index --index pypi --version X.Y.Z \\
        --sums SHA256SUMS [--download-dir DIR] [--attempts 20]

``check-dist`` asserts that ``DIST_DIR`` holds exactly the wheel
``chemthermo-X.Y.Z-py3-none-any.whl`` and the sdist ``chemthermo-X.Y.Z.tar.gz``,
that both declare version ``X.Y.Z`` in their metadata, and that neither ships
``tests/``, ``.agents/``, ``benchmarks/``, ``notebooks/``, ``database/`` or
``tools/`` (MANIFEST.in; the wheel carries only the ``chemthermo`` package).

``sums`` prints ``sha256sum``-format lines for the two files, in name order.

``verify-index`` reads the index's JSON API for ``chemthermo X.Y.Z`` and
requires the release to hold exactly the files named in ``SHA256SUMS`` with
exactly those digests. With ``--download-dir`` it also downloads each file and
re-hashes the downloaded bytes. It retries while the version is not yet
visible (new uploads take a moment to appear) and fails on any mismatch at
once. It needs no credentials.
"""

from __future__ import annotations

import argparse
import email.parser
import hashlib
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

PROJECT = "chemthermo"
INDEXES = {"pypi": "https://pypi.org", "testpypi": "https://test.pypi.org"}
#: Top-level source directories that must never reach a distribution.
EXCLUDED = ("tests", ".agents", "benchmarks", "notebooks", "database", "tools")


def expected_names(version: str) -> tuple[str, str]:
    return f"{PROJECT}-{version}-py3-none-any.whl", f"{PROJECT}-{version}.tar.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _metadata_version(text: str) -> str | None:
    return email.parser.HeaderParser().parsestr(text).get("Version")


def check_dist(version: str, dist: Path) -> list[str]:
    wheel_name, sdist_name = expected_names(version)
    present = sorted(p.name for p in dist.iterdir())
    if present != sorted([wheel_name, sdist_name]):
        return [f"{dist} holds {present}, expected exactly {[sdist_name, wheel_name]}"]
    problems: list[str] = []

    with zipfile.ZipFile(dist / wheel_name) as wheel:
        names = wheel.namelist()
        tops = {name.split("/", 1)[0] for name in names}
        allowed = {PROJECT, f"{PROJECT}-{version}.dist-info"}
        if tops != allowed:
            problems.append(f"wheel top-level entries {sorted(tops)}, expected {sorted(allowed)}")
        metadata = wheel.read(f"{PROJECT}-{version}.dist-info/METADATA").decode("utf-8")
        if _metadata_version(metadata) != version:
            problems.append(f"wheel METADATA Version {_metadata_version(metadata)!r}")

    root = f"{PROJECT}-{version}"
    with tarfile.open(dist / sdist_name, "r:gz") as sdist:
        members = sdist.getnames()
        for excluded in EXCLUDED:
            leaked = [m for m in members if m.startswith(f"{root}/{excluded}/")]
            if leaked:
                problems.append(f"sdist ships {excluded}/ ({len(leaked)} files, e.g. {leaked[0]})")
        pkg_info = sdist.extractfile(f"{root}/PKG-INFO")
        if pkg_info is None:
            problems.append("sdist has no PKG-INFO")
        elif _metadata_version(pkg_info.read().decode("utf-8")) != version:
            problems.append("sdist PKG-INFO Version does not match")
    return problems


def _fetch_release(index: str, version: str) -> dict[str, str] | None:
    """{filename: sha256} the index lists for the version; None while absent."""
    url = f"{INDEXES[index]}/pypi/{PROJECT}/{version}/json"
    request = urllib.request.Request(url, headers={"Cache-Control": "no-cache"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise
    return {item["filename"]: item["digests"]["sha256"] for item in payload["urls"]}


def _file_urls(index: str, version: str) -> dict[str, str]:
    url = f"{INDEXES[index]}/pypi/{PROJECT}/{version}/json"
    with urllib.request.urlopen(url, timeout=30) as response:
        payload = json.load(response)
    return {item["filename"]: item["url"] for item in payload["urls"]}


def read_sums(path: Path) -> dict[str, str]:
    sums: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            sums[Path(name.lstrip("*")).name] = digest
    return sums


def verify_index(
    index: str, version: str, sums: dict[str, str], download_dir: Path | None, attempts: int
) -> list[str]:
    listed: dict[str, str] | None = None
    for attempt in range(attempts):
        listed = _fetch_release(index, version)
        if listed is not None and set(listed) >= set(sums):
            break
        print(f"{index}: {PROJECT} {version} not complete yet (attempt {attempt + 1})")
        time.sleep(15)
    if listed is None:
        return [f"{index} has no {PROJECT} {version}"]
    problems: list[str] = []
    if set(listed) != set(sums):
        problems.append(f"{index} lists {sorted(listed)}, the frozen build is {sorted(sums)}")
    for name, digest in sorted(sums.items()):
        if name in listed and listed[name] != digest:
            problems.append(f"{index} {name} sha256 {listed[name]} != frozen {digest}")
        elif name in listed:
            print(f"{index}: {name} sha256 {digest} matches")
    if download_dir is not None and not problems:
        download_dir.mkdir(parents=True, exist_ok=True)
        for name, url in sorted(_file_urls(index, version).items()):
            target = download_dir / name
            urllib.request.urlretrieve(url, target)
            if sha256(target) != sums[name]:
                problems.append(f"downloaded {name} hashes to {sha256(target)}, not {sums[name]}")
            else:
                print(f"{index}: downloaded {name}, bytes match")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description="Checks on frozen chemthermo release files.")
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check-dist")
    check.add_argument("--version", required=True)
    check.add_argument("dist", type=Path)
    sums = commands.add_parser("sums")
    sums.add_argument("dist", type=Path)
    verify = commands.add_parser("verify-index")
    verify.add_argument("--index", choices=sorted(INDEXES), required=True)
    verify.add_argument("--version", required=True)
    verify.add_argument("--sums", type=Path, required=True)
    verify.add_argument("--download-dir", type=Path, default=None)
    verify.add_argument("--attempts", type=int, default=20)
    args = parser.parse_args()

    if args.command == "sums":
        for path in sorted(args.dist.iterdir()):
            print(f"{sha256(path)}  {path.name}")
        return 0
    if args.command == "check-dist":
        problems = check_dist(args.version, args.dist)
    else:
        problems = verify_index(
            args.index, args.version, read_sums(args.sums), args.download_dir, args.attempts
        )
    for problem in problems:
        print(f"FAIL: {problem}", file=sys.stderr)
    if not problems:
        print(f"{args.command}: ok")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
