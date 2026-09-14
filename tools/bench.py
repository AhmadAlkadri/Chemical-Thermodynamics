#!/usr/bin/env python3
"""Thin wrapper so the benchmark harness runs from a source checkout.

``python tools/bench.py --out benchmarks/after_<sha>.json`` is exactly
``python -m chemthermo.bench --out ...``; it exists so the command works
without the package being importable from the current working directory.
See ``benchmarks/README.md`` and ADR-0023.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from chemthermo.bench import main  # noqa: E402  (after the sys.path insert)

if __name__ == "__main__":
    sys.exit(main())
