# Steering Brief

## What changed since last brief (files + bullets)
- `tools/build_database.py`
  - Rebuilt canonical DB generation to match runtime schema, added deterministic `--check` against `src/chemthermo/data/components.json`, and added optional non-runtime mirror output at `database/components.mirror.json`.
- `src/chemthermo/cli.py`, `src/chemthermo/__main__.py`, `pyproject.toml`
  - Added public `chemthermo` CLI entrypoint and module execution with `tp-flash` command.
- `tests/test_database_build_tool.py`, `tests/test_cli_tp_flash.py`, `tests/validation/test_reference_case.py`, `tests/fixtures/cli/tp_flash_v1.json`
  - Added focused tests for DB tooling, CLI contract/exit codes, and optional external validation reference case.
- `examples/database/01_component_database_demo.py`, `examples/validation/00_reference_case.py`
  - Added golden-path scripts for DB visibility and deterministic external validation.
- `README.md`, `database/README.md`, `examples/README.md`, `.gitignore`
  - Clarified canonical DB path policy and documented CLI/validation workflows.
- `.agents/brain/adr/0003-cli-entrypoint.md`, `.agents/brain/brain.md`
  - Recorded CLI API decision and updated architecture/public API status.

## Current architecture (8-12 lines)
- Canonical runtime DB path is `src/chemthermo/data/components.json`.
- Runtime DB loading is package-resource based via `chemthermo.data` helpers.
- Raw DB source tables stay in `database/organics.txt` and `database/inorganics.txt`.
- `tools/build_database.py` is the canonical regeneration/check tool for packaged DB sync.
- Optional mirror output path is `database/components.mirror.json` (non-runtime, generated-only).
- Core thermodynamic flow remains `Component/Composition/Mixture -> models -> flash_tp -> FlashResult`.
- New CLI layer (`chemthermo tp-flash`) maps user inputs to `Mixture + PengRobinsonEOS + flash_tp`.
- CLI supports deterministic text/json outputs with schema version and diagnostics.
- External validation remains optional and script-driven (`examples/validation/00_reference_case.py`).

## Public API status (stable vs experimental)
- Stable:
  - `chemthermo` Python exports in `src/chemthermo/__init__.py`.
  - `chemthermo.eos` and `chemthermo.vlle` documented public subpackages.
  - CLI script `chemthermo` with `tp-flash` subcommand and defined exit-code contract.
- Experimental/placeholder:
  - PC-SAFT residual Helmholtz implementation remains placeholder in open-source build.
  - Validation policy promotion to required CI remains deferred.

## Risks / unknowns
- CLI JSON contract must remain backward-compatible or version-bumped.
- DB builder currently assumes raw table column conventions stay unchanged.
- Optional validation still depends on external `thermo` package availability.

## Next 3 recommended actions
- Add gamma-phi mode support to CLI `tp-flash`.
- Add richer provenance metadata in DB tooling with compatibility guardrails.
- Revisit external validation CI promotion when package maturity allows.

## One simplification / deletion candidate
- Remove or archive deprecated `database/components.json` once all docs/tooling users are migrated.

## Assumptions I'm making
- Runtime DB source of truth remains `src/chemthermo/data/components.json`.
- `database/components.mirror.json` is generated-only and should not be committed.
- CLI v1 scope remains TP flash with Peng-Robinson phi-phi only.

## How to validate quickly
- `python tools/build_database.py --check`
- `python examples/database/01_component_database_demo.py`
- `python -m chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- `python examples/validation/00_reference_case.py` (with `pip install -e ".[validation]"`)
