# ADR-0037: Remove the `chemthermo.vlle` plugin package before the first PyPI release

Status: accepted
Date: 2026-09-25

## Context
ADR-0013 deprecated `chemthermo.vlle` (the boundary to an external
`chemthermo_vlle` engine) because three-phase equilibrium is discovered
in-tree (ADR-0011, ADR-0020) and the only candidate plugin was a scaffold
with no physics. It kept the package importable, with a `DeprecationWarning`,
"for one deprecation cycle". 0.4.0 is the first release published to PyPI;
the owner decided (2026-09-25) not to carry the deprecated package into it.

## Decision
1. Delete `src/chemthermo/vlle/` and its tests. `import chemthermo.vlle` raises
   `ModuleNotFoundError` (pinned by a test).
2. `flash_tp(..., flash_mode="vlle")` keeps raising `ModelError`, now naming
   the in-tree route and this ADR.
3. The deprecated `flash_mode="gamma-phi"` stays (the CLI v1 contract of
   ADR-0004 exposes it); it keeps its `DeprecationWarning`.

## Alternatives considered
- Keep it deprecated through 0.4.0 (rejected by the owner: no reason to put a
  known-dead API on PyPI).
- Remove gamma-phi too (rejected: breaks the CLI v1 contract).

## Consequences
- Breaking for anyone importing `chemthermo.vlle` from a git install; allowed
  pre-1.0 with a minor bump (ADR-0031 item 2); called out in CHANGELOG.
- ADR-0013's "one deprecation cycle" is closed.

## Supersedes (optional)
Completes ADR-0013; amends ADR-0001's list of documented subpackages.

## Superseded by (optional)
None.
