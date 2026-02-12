# ADR-0004: Extend `tp-flash` CLI with gamma-phi mode

Status: accepted
Date: 2026-02-12

## Context
`flash_tp` already supports both `phi-phi` and `gamma-phi`, but the public CLI
surface only exposed phi-phi. Users had to write Python to run gamma-phi
despite the existing solver capability and NRTL activity model support.

## Decision
Extend `chemthermo tp-flash` with:

- `--flash-mode {phi-phi,gamma-phi}` (default: `phi-phi`).
- Gamma-phi wiring that uses `NRTL()` as the liquid activity model and
  `PengRobinsonEOS()` as the EOS.

Compatibility decisions:

- Keep `cli_schema_version = 1`.
- Keep existing top-level JSON shape and keys.
- Keep exit-code contract unchanged:
  - `0` success
  - `1` validation/runtime/model error
  - `2` usage/argparse error
  - `3` solver nonconvergence

## Alternatives considered
- Keep CLI phi-phi only and require Python for gamma-phi (rejected: avoidable
  usability gap).
- Add activity-model selection flags in this slice (rejected: broader than
  thin-slice scope).
- Bump CLI schema to v2 (rejected: no structural break required).

## Consequences
- Positive: gamma-phi becomes available through the same public CLI entrypoint.
- Positive: deterministic CLI output contract is preserved without schema churn.
- Tradeoff: gamma-phi pair coverage is limited by packaged NRTL parameter data.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
