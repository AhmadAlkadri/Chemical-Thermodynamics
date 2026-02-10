"""Backward-compatible shim for phase-boundary VLE utilities.

`chemthermo.vle` is kept for compatibility. Prefer importing these functions from
`chemthermo.phase_boundary` (or from top-level `chemthermo`).
"""

from __future__ import annotations

from .phase_boundary import bubble_pressure, bubble_temperature, dew_pressure, dew_temperature

__all__ = ["bubble_temperature", "dew_temperature", "bubble_pressure", "dew_pressure"]
