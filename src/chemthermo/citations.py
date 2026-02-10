"""Bibliographic database management."""

from __future__ import annotations

from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Mapping

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode


@lru_cache(maxsize=1)
def load_references(path: Path | None = None) -> Mapping[str, dict[str, str]]:
    """Load the BibTeX database from package resources or an override path."""
    handle = None
    if path is None:
        resource = resources.files("chemthermo.data") / "references.bib"
        if resource.is_file():
            handle = resource.open("r", encoding="utf-8")
    else:
        resolved = Path(path)
        if resolved.exists():
            handle = resolved.open("r", encoding="utf-8")

    if handle is None:
        return {}

    with handle as stream:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bib_database = bibtexparser.load(stream, parser=parser)

    return {entry["ID"]: entry for entry in bib_database.entries}


class BibliographicDatabase:
    """Interface to the bibliographic data."""

    def __init__(self, path: Path | None = None):
        self._entries = load_references(path)

    def get_entry(self, key: str) -> dict[str, str] | None:
        """Return the raw BibTeX entry dict for a key."""
        return self._entries.get(key)

    def get_citation_text(self, key: str) -> str:
        """Format a human-readable citation string."""
        entry = self.get_entry(key)
        if not entry:
            return f"[Unknown citation: {key}]"

        authors = entry.get("author", "Unknown Author").replace("\n", " ")
        title = entry.get("title", "Untitled").replace("{", "").replace("}", "")
        year = entry.get("year", "n.d.")

        return f"{authors} ({year}). {title}."
