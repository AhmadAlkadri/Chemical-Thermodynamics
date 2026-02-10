# ADR-0003: Public CLI entrypoint for TP flash

Status: accepted
Date: 2026-02-10

## Context
The package exposes public Python APIs but no public command-line entrypoint.
Users currently need Python scripts for simple TP flash runs.

## Decision
Add a public CLI entrypoint in `[project.scripts]`:

- `chemthermo = "chemthermo.cli:main"`

Add module execution support:

- `python -m chemthermo`

Initial CLI surface (v1) is intentionally minimal:

- `chemthermo tp-flash` using Peng-Robinson EOS in phi-phi mode.
- JSON output includes `cli_schema_version` and solver/diagnostic metadata.
- Exit code contract:
  - `0` success
  - `1` validation/runtime errors
  - `2` usage/argparse errors
  - `3` solver nonconvergence

## Alternatives considered
- Keep script-only examples and no public CLI (rejected: higher adoption friction).
- Add broad multi-command CLI in one slice (rejected: violates thin-slice scope).

## Consequences
- CLI now becomes part of the public API contract (ADR-0001 rule for scripts).
- Future CLI expansions should preserve backward compatibility for v1 fields or bump `cli_schema_version`.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
