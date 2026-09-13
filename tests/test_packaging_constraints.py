"""Guard packaging constraints that keep a fresh install resolvable and importable."""

from __future__ import annotations

import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"


def _bibtexparser_specifier() -> str:
    payload = tomllib.loads(PYPROJECT_PATH.read_text(encoding="utf-8"))
    dependencies: list[str] = payload["project"]["dependencies"]
    for dependency in dependencies:
        name = dependency.split(";")[0].strip()
        for separator in ("==", ">=", "<=", "!=", "~=", ">", "<"):
            index = name.find(separator)
            if index != -1:
                name = name[:index].strip()
                break
        if name == "bibtexparser":
            return dependency
    raise AssertionError("bibtexparser dependency not found in pyproject.toml")


def test_bibtexparser_dependency_excludes_major_version_2() -> None:
    """chemthermo.citations imports the bibtexparser 1.x API (bparser, customization),
    which was removed in bibtexparser 2.x. The dependency specifier must exclude
    major version 2 so a fresh install does not resolve an incompatible release."""
    specifier = _bibtexparser_specifier()

    assert "<2" in specifier.replace(" ", ""), (
        f"bibtexparser specifier {specifier!r} must exclude 2.x releases via an upper bound (e.g. '<2')"
    )
    assert ">=1.4.0" in specifier.replace(" ", ""), (
        f"bibtexparser specifier {specifier!r} must still allow the 1.x API chemthermo uses"
    )
