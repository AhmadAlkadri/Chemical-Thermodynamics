#!/usr/bin/env python3
"""Refuse a release whose tag, commit and version strings disagree (ADR-0040).

Run from the root of a checkout of the release commit, before anything is
built or uploaded::

    python tools/release_preflight.py --tag v0.4.1 [--sha <full commit id>] [--require-on main]

It checks, and exits non-zero on the first class of failure it finds:

- the tag is exactly ``vX.Y.Z`` (a final release; prereleases ``X.Y.ZbN``
  stay GitHub-only by ADR-0031 and are not published by this path);
- the tag exists locally as an **annotated** tag and peels to ``HEAD``, and to
  ``--sha`` when given (the commit the release event or the operator names);
- ``--require-on BRANCH``: that commit is reachable from ``origin/BRANCH``;
- ``pyproject.toml`` ``[project].version`` and ``chemthermo.__version__``
  (read from ``src/chemthermo/__init__.py`` without importing it) are both
  ``X.Y.Z``;
- ``CHANGELOG.md`` has a ``## X.Y.Z (`` section;
- the working tree has no tracked changes.

It prints ``version=X.Y.Z`` on success, in the ``key=value`` form GitHub
Actions reads from ``$GITHUB_OUTPUT``. It does not contact any package index.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
import tomllib
from pathlib import Path

TAG_PATTERN = re.compile(r"^v(?P<version>(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*))$")


def _git(*args: str) -> str:
    completed = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {completed.stderr.strip()}")
    return completed.stdout.strip()


def _package_version(init_file: Path) -> str | None:
    """``__version__ = "..."`` at module level, read without importing."""
    for node in ast.parse(init_file.read_text(encoding="utf-8")).body:
        if (
            isinstance(node, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "__version__" for t in node.targets)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            return node.value.value
    return None


def check(tag: str, sha: str | None, require_on: str | None, root: Path) -> list[str]:
    """Every reason not to release ``tag`` from ``root``; empty means go."""
    match = TAG_PATTERN.match(tag)
    if match is None:
        return [f"tag {tag!r} is not vX.Y.Z (final releases only; see ADR-0031)"]
    version = match["version"]
    problems: list[str] = []

    try:
        kind = _git("cat-file", "-t", f"refs/tags/{tag}")
        tagged = _git("rev-parse", f"refs/tags/{tag}^{{commit}}")
        head = _git("rev-parse", "HEAD")
    except RuntimeError as error:
        return [f"tag {tag} cannot be resolved: {error}"]
    if kind != "tag":
        problems.append(f"tag {tag} is a {kind}, not an annotated tag")
    if tagged != head:
        problems.append(f"tag {tag} peels to {tagged}, but HEAD is {head}")
    if sha is not None and tagged != sha:
        problems.append(f"tag {tag} peels to {tagged}, not the expected commit {sha}")
    if require_on is not None:
        reachable = subprocess.run(
            ["git", "merge-base", "--is-ancestor", tagged, f"origin/{require_on}"],
            capture_output=True,
            check=False,
        )
        if reachable.returncode != 0:
            problems.append(f"{tagged} is not reachable from origin/{require_on}")

    with (root / "pyproject.toml").open("rb") as handle:
        project_version = tomllib.load(handle)["project"]["version"]
    if project_version != version:
        problems.append(f"pyproject.toml version is {project_version!r}, the tag says {version!r}")
    package_version = _package_version(root / "src" / "chemthermo" / "__init__.py")
    if package_version != version:
        problems.append(f"chemthermo.__version__ is {package_version!r}, the tag says {version!r}")

    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    if not re.search(rf"^## {re.escape(version)} \(", changelog, flags=re.MULTILINE):
        problems.append(f"CHANGELOG.md has no '## {version} (' section")

    if _git("status", "--porcelain", "--untracked-files=no"):
        problems.append("the working tree has tracked changes")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Refuse a release whose tag and versions disagree."
    )
    parser.add_argument("--tag", required=True, help="the release tag, vX.Y.Z")
    parser.add_argument("--sha", default=None, help="the full commit id the tag must peel to")
    parser.add_argument("--require-on", default=None, help="a branch the commit must be on")
    args = parser.parse_args()
    problems = check(args.tag, args.sha, args.require_on, Path.cwd())
    for problem in problems:
        print(f"FAIL: {problem}", file=sys.stderr)
    if problems:
        return 1
    print(f"version={args.tag[1:]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
