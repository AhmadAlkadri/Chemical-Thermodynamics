# Steering Brief

## What changed since last brief (files + bullets)
- `src/chemthermo/models/peng_robinson.py`
  - Fixed the diagonal-kij bug: `aij`'s diagonal is now always unaffected by `kij` (pure-component `a_ii` never corrupted). Added per-pair `kij` support (`float` or `Mapping[tuple[str, str], float]` keyed by normalized component names, canonicalized to a sorted tuple in `__post_init__`), plus two private helpers (`_kij_matrix`, `_mixture_parameters`) so `fugacity_coefficients`/`compressibility_factor` no longer duplicate the mixing-rule code. Public method signatures unchanged; `kij=0.0` (default) is bit-identical to before.
- `tests/test_pr_eos.py`, `tests/validation/test_pr_kij_vs_thermo.py`
  - Unit coverage (pure-component invariance, diagonal-not-corrupted, scalar-vs-mapping equivalence, name normalization/symmetry, conflicting/identical-pair `ModelError`, unknown-pair no-op, permutation invariance, flash regression unchanged) and a `thermo` 0.6.0 PRMIX cross-check (binary + synthetic 3-component kij matrix, pure-component limit, end-to-end `flash_tp`/`stability_tp` vs `FlashVL`/`stability_test_Michelsen`).
- `examples/basic/tp_flash_pr_kij_demo.py`, `examples/README.md`, `README.md`
  - Golden path contrasting `kij=0.0` against a per-pair mapping on the same feed; README "Binary interaction parameters (kij)" section.
- `.agents/brain/adr/0006-pr-kij-matrix.md`, `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Recorded the per-pair kij API decision (ADR-0006) and added validation Case K-1 (pre-fix vs post-fix discrepancy against `thermo`: ~0.365 max |d ln phi| before, ~3.4e-4 after, at one representative state).
- `pyproject.toml`, `tests/test_packaging_constraints.py`, `README.md`
  - Constrained the `bibtexparser` dependency to `>=1.4.0,<2`: a fresh (non-editable) install previously resolved `bibtexparser` 2.x, whose removed `bparser`/`customization` modules `src/chemthermo/citations.py` imports at module load time, breaking `import chemthermo` entirely; added a test that reads the specifier from `pyproject.toml` and a README "Common issues" bullet.
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
- `PengRobinsonEOS.kij` is a scalar (off-diagonal only) or a name-keyed per-pair `Mapping`; `flash_tp` and `stability_tp` results for nonzero kij are now trustworthy (ADR-0006).
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
  - `PengRobinsonEOS(kij=...)` scalar-or-mapping constructor contract (ADR-0006).
  - `chemthermo.eos` and `chemthermo.vlle` documented public subpackages.
  - CLI script `chemthermo` with `tp-flash` subcommand and defined exit-code contract.
- Experimental/placeholder:
  - PC-SAFT residual Helmholtz implementation remains placeholder in open-source build.
  - Validation policy promotion to required CI remains deferred.

## Risks / unknowns
- `stability_tp` reporting "stable" is bounded by its deterministic trial set; it is not a global proof, and the docs must keep saying so.
- chemthermo uses the rounded PR constants 0.45724 / 0.07780 while `thermo` uses the exact roots; this bounds external agreement at a few times 1e-4 in ln(phi) (was ~2e-4 at the stability states, ~4-6e-4 at some nonzero-kij states -- see Case K-1).
- The illustrative kij value used in docs/examples/tests (0.0411, Methane/n-Decane) is explicitly NOT a validated literature parameter -- do not let it drift into being read as one.
- `thermo`'s own `CEOSLiquid`/`CEOSGas` root solver was observed to be numerically order-sensitive (not permutation-invariant) at some near-critical-locus states during `pr-kij-matrix` development; avoid cross-checking permutation invariance against `thermo` at such states (use chemthermo-internal invariance checks instead, as `tests/validation/test_pr_kij_vs_thermo.py` now does).
- CLI JSON contract must remain backward-compatible or version-bumped.
- DB builder currently assumes raw table column conventions stay unchanged.
- Optional validation still depends on external `thermo` package availability.
- (Removed) Unpinned `bibtexparser>=1.4.0` allowed a fresh/non-editable install to resolve `bibtexparser` 2.x, whose removed `bparser`/`customization` modules broke `import chemthermo`; this was invisible locally because the dev venv already had 1.4.4 installed. Now pinned to `>=1.4.0,<2` and covered by `tests/test_packaging_constraints.py`.

## Next 3 recommended actions
- `stability-tpd-nrtl`: activity-model (liquid-liquid) tangent-plane stability.
- `flash-auto-phase-detection`: let `flash_tp` consume `stability_tp` instead of the K-bound heuristic.
- (not yet scoped): define after `flash-auto-phase-detection` lands.

## One simplification / deletion candidate
- Remove or archive deprecated `database/components.json` once all docs/tooling users are migrated.

## Assumptions I'm making
- Runtime DB source of truth remains `src/chemthermo/data/components.json`.
- `database/components.mirror.json` is generated-only and should not be committed.
- CLI v1 scope remains TP flash with Peng-Robinson EOS and `--flash-mode {phi-phi,gamma-phi}`.

## How to validate quickly
- `python examples/basic/tp_flash_pr_kij_demo.py`
- `python -m pytest tests/validation/test_pr_kij_vs_thermo.py -q` (with `pip install -e ".[validation]"`)
- `python examples/basic/stability_tp_peng_robinson_demo.py`
- `python examples/validation/06_stability_vs_thermo.py` (with `pip install -e ".[validation]"`)
- `python tools/build_database.py --check`
- `python examples/database/01_component_database_demo.py`
- `python -m chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- `python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- `python examples/validation/00_reference_case.py` (with `pip install -e ".[validation]"`)
