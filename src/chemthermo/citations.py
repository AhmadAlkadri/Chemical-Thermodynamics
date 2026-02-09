"""Bibliographic database management."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Mapping

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode, author


@lru_cache(maxsize=1)
def load_references(path: Path | None = None) -> Mapping[str, dict[str, str]]:
    """Load the BibTeX database from the default or specified path.
    
    Returns:
        A dictionary mapping citation keys to BibTeX entry dictionaries.
    """
    if path is None:
        # Default to database/references.bib relative to the package root or known location
        # For now, we assume it's in the standard location relative to this file
        # chemthermo/src/chemthermo/citations.py -> ../../../database/references.bib
        # This might need adjustment based on installation structure.
        # Ideally, resources should be used or a configured path.
        # Fallback to local dev path for now.
        path = Path(__file__).resolve().parents[2] / "database" / "references.bib"
    
    if not path.exists():
        raise FileNotFoundError(f"BibTeX database not found at {path}")

    with open(path, encoding="utf-8") as f:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bib_database = bibtexparser.load(f, parser=parser)
    
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
        
        # Simple formatting logic suitable for plain text output
        authors = entry.get("author", "Unknown Author").replace("\n", " ")
        title = entry.get("title", "Untitled").replace("{", "").replace("}", "")
        year = entry.get("year", "n.d.")
        
        return f"{authors} ({year}). {title}."
