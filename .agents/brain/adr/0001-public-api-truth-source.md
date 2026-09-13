# ADR-0001: Public API truth source

Status: accepted
Date: 2026-02-03
Amended by: ADR-0013 (2026-09-13) - `chemthermo.vlle` is now DEPRECATED, not
merely "documented public"; see decision 2 below and ADR-0013.

## Context
Agents and contributors need a stable definition of "public API" to avoid
accidental breaking changes. The repo exposes a top-level package plus a
couple of explicitly documented subpackages, and it has no CLI entry points.

## Decision
Public API is defined by the following sources (in order):
1) `src/chemthermo/__init__.py` and its `__all__` list.
2) Subpackages explicitly documented in `README.md` and their `__all__` lists
   (currently `chemthermo.eos` and `chemthermo.vlle`). **`chemthermo.vlle` is
   DEPRECATED as of ADR-0013 (2026-09-13)**: it remains a public, importable
   subpackage for one deprecation cycle, but importing it or calling
   `get_vlle_engine()` emits `DeprecationWarning`, and it is no longer the
   documented way to reach multiphase equilibrium - see ADR-0011 and
   ADR-0013.
3) CLI entry points only if they are declared in `[project.scripts]` in
   `pyproject.toml` (none today).

Anything outside these sources is internal and may change without ADRs.

## Alternatives considered
- Everything under `src/chemthermo` is public (rejected: too broad).
- Documentation alone defines public API (rejected: docs can drift).
- Tests define public API (rejected: tests are not a user contract).

## Consequences
- Any change to public exports or documented public subpackages requires a new
  ADR (or an update that supersedes this one).
- If new subpackages are documented as public, they must be added to this ADR.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
