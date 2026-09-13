# Steering Brief

## What changed since last brief (files + bullets)
- `src/chemthermo/stability/` (`__init__.py`, `tp.py`, `results.py`, `settings.py`), `src/chemthermo/__init__.py`
  - Added public `stability_tp` (Michelsen tangent-plane stability at fixed T, P, z) with `StabilityResult`, `StabilitySettings`, `StabilityTrial`; min-Gibbs compressibility root selected generically by minimizing `sum_i w_i ln phi_i(w)` over the existing `EquationOfState` interface.
- `tests/test_stability_tp.py`, `tests/validation/test_stability_vs_thermo.py`
  - Invariant tests (tangent-plane identity, stationarity plus finite-difference gradient of tm, tm/tpd/sum(W) relations, marginal stability of converged equilibrium phases, permutation invariance, degenerate feeds) and a 7-state cross-check against `thermo` 0.6.0.
- `examples/basic/stability_tp_peng_robinson_demo.py`, `examples/validation/06_stability_vs_thermo.py`, `examples/README.md`, `README.md`
  - Golden path for the stability verdict and a deterministic external validation script.
- `.agents/brain/adr/0005-stability-tp-public-api.md`, `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Recorded the public stability API decision and opened the validation-case ledger ("thermodynamics exam").
- `src/chemthermo/cli.py`, `tests/test_cli_tp_flash.py`, `tests/fixtures/cli/tp_flash_gamma_phi_v1.json`
  - Added `--flash-mode {phi-phi,gamma-phi}` to `tp-flash`, wired gamma-phi to `NRTL()` activity model, and expanded deterministic CLI contract tests/fixture coverage.
- `.agents/brain/adr/0004-cli-tp-flash-gamma-phi.md`, `.agents/brain/brain.md`
  - Recorded CLI gamma-phi public-contract decision and advanced roadmap slices.
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
- Phase stability (`chemthermo.stability`) is a solver-independent sibling of `chemthermo.flash`; `flash_tp` does NOT consume it yet.
- Canonical runtime DB path is `src/chemthermo/data/components.json`.
- Runtime DB loading is package-resource based via `chemthermo.data` helpers.
- Raw DB source tables stay in `database/organics.txt` and `database/inorganics.txt`.
- `tools/build_database.py` is the canonical regeneration/check tool for packaged DB sync.
- Optional mirror output path is `database/components.mirror.json` (non-runtime, generated-only).
- Core thermodynamic flow remains `Component/Composition/Mixture -> models -> flash_tp -> FlashResult`.
- CLI `tp-flash` maps inputs to `Mixture + PengRobinsonEOS + flash_tp`, with optional `NRTL()` activity model when `--flash-mode gamma-phi`.
- CLI supports deterministic text/json outputs with schema version and diagnostics.
- External validation remains optional and script-driven (`examples/validation/00_reference_case.py`).

## Public API status (stable vs experimental)
- Stable:
  - `chemthermo` Python exports in `src/chemthermo/__init__.py`, including `stability_tp` / `StabilityResult` / `StabilitySettings` / `StabilityTrial` (ADR-0005).
  - `chemthermo.eos` and `chemthermo.vlle` documented public subpackages.
  - CLI script `chemthermo` with `tp-flash` subcommand and defined exit-code contract.
- Experimental/placeholder:
  - PC-SAFT residual Helmholtz implementation remains placeholder in open-source build.
  - Validation policy promotion to required CI remains deferred.

## Risks / unknowns
- `stability_tp` reporting "stable" is bounded by its deterministic trial set; it is not a global proof, and the docs must keep saying so.
- `PengRobinsonEOS.kij` is a scalar applied to the diagonal of `aij` as well, which is wrong; only `kij = 0.0` is validated today. Fix in the `pr-kij-matrix` slice.
- chemthermo uses the rounded PR constants 0.45724 / 0.07780 while `thermo` uses the exact roots; this bounds external agreement at ~2e-4 in ln(phi).
- CLI JSON contract must remain backward-compatible or version-bumped.
- DB builder currently assumes raw table column conventions stay unchanged.
- Optional validation still depends on external `thermo` package availability.

## Next 3 recommended actions
- `pr-kij-matrix`: fix the diagonal-kij bug and accept a per-pair kij matrix.
- `stability-tpd-nrtl`: activity-model (liquid-liquid) tangent-plane stability.
- `flash-auto-phase-detection`: let `flash_tp` consume `stability_tp` instead of the K-bound heuristic.

## One simplification / deletion candidate
- Remove or archive deprecated `database/components.json` once all docs/tooling users are migrated.

## Assumptions I'm making
- Runtime DB source of truth remains `src/chemthermo/data/components.json`.
- `database/components.mirror.json` is generated-only and should not be committed.
- CLI v1 scope remains TP flash with Peng-Robinson EOS and `--flash-mode {phi-phi,gamma-phi}`.

## How to validate quickly
- `python examples/basic/stability_tp_peng_robinson_demo.py`
- `python examples/validation/06_stability_vs_thermo.py` (with `pip install -e ".[validation]"`)
- `python tools/build_database.py --check`
- `python examples/database/01_component_database_demo.py`
- `python -m chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- `python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- `python examples/validation/00_reference_case.py` (with `pip install -e ".[validation]"`)
