"""``python -m chemthermo.bench`` entry point (ADR-0023)."""

from __future__ import annotations

import sys

from . import main

if __name__ == "__main__":
    sys.exit(main())
