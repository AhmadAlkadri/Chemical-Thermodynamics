# Project Brain: chemthermo

How to use this document
- Read this first for repo invariants, public API, and CI contract.
- Follow evidence tags before changing load-bearing behavior.
- Keep updates short; link to existing docs instead of duplicating.

**Agent Contract**
- Read first: `.agents/brain/brain.md`, `.agents/brain/steering-brief.md`, relevant ADRs in `.agents/brain/adr/`.
- Do-not-touch list (initial): public API boundaries (ADR-0001), CI contract (`.github/workflows/ci.yml`), compatibility rules (public API + schema versions), numerics invariants (SI units, composition tolerance, flash determinism).
- Definition of done: run canonical bootstrap, CI gates, and golden path from `.agents/dev-contract.md`; update docs; keep no TODOs in critical path.
- Stop if >3 plausible root causes: write an experiment/instrumentation plan first.
- If uncertain after 2 iterations, produce a minimal repro and document findings.
- Complexity receipts rule: any new abstraction must state why, bug prevented, cost, and what breaks if omitted.
- ADR rule: any public API or architectural change requires a new ADR (or update + supersede).
- **Thin Vertical Slice Rule** (ADR-0002): Every exposed capability must be end-to-end usable. No scaffolded/partial features.
- **Golden Path Rule** (ADR-0002): Every Thin Vertical Slice must ship with at least one golden path (example/test) that executes successfully from a clean environment.
- Keep claims factual; add evidence tags or links for load-bearing statements.

**Skills**
- Local skill catalog: `.agents/skills/README.md`
- Preferred repo workflow skill: `.agents/skills/chemthermo-change-loop/SKILL.md`

## 0) Repo at a glance
- Purpose: chemical engineering thermodynamics utilities packaged as `chemthermo`, SI units throughout. (source: README.md)
- Primary language/toolchain: Python 3.11+, setuptools build via `pyproject.toml`. (source: pyproject.toml)
- Primary entry points: `chemthermo` top-level API, `chemthermo.eos` registry, `chemthermo.vlle` plugin boundary, and CLI script `chemthermo`. (source: src/chemthermo/__init__.py, src/chemthermo/eos/__init__.py, src/chemthermo/vlle/__init__.py, src/chemthermo/cli.py, pyproject.toml)
- Tests: `pytest` (plus lint/type checks in CI). (source: .github/workflows/ci.yml)
- Golden Path command is maintained in `.agents/dev-contract.md`. (source: .agents/dev-contract.md)

## 1) Purpose and scope policy
- Purpose: provide core thermodynamics utilities (components, mixtures, EOS/activities, TP flash) in SI units. (source: README.md, src/chemthermo/flash/tp.py, src/chemthermo/models/base.py)
- Scope policy: VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## 1.1) Thin Vertical Slice Enforcement
- Handoff blocker: Missing explicit slice declaration, "After this change, user can X by running Y."
- Handoff blocker: No runnable golden-path command/test proving the claimed user flow.
- Handoff blocker: No focused test evidence for touched behavior.
- Handoff blocker: Missing CI-equivalent and install-smoke evidence from `.agents/dev-contract.md`.
- Handoff blocker: User-visible behavior changed without corresponding docs updates.

## 2) Public API surface (current)
Definition of public API follows ADR-0001 (source of truth rules in `.agents/brain/adr/0001-public-api-truth-source.md`).

Stable (public) entry points
- `chemthermo` top-level exports in `__all__` (core types, flash API, phase-stability API, models, parameters, exceptions, units, validation helpers). (source: src/chemthermo/__init__.py)
- Phase stability: `stability_tp`, `StabilityResult`, `StabilitySettings`, `StabilityTrial` (ADR-0005, ADR-0007). `stability_tp` reports `status` in {"stable","unstable","inconclusive"}; "stable" means "no negative tangent-plane distance was found from the deterministic trial set", not a global proof. (source: src/chemthermo/stability/, .agents/brain/adr/0005-stability-tp-public-api.md)
- `stability_tp(mixture, *, temperature_K, pressure_Pa, eos=None, activity_model=None, settings=None)` (ADR-0007). **Exactly one** of `eos` / `activity_model` is required; passing both raises `ModelError` because combined gamma-phi stability (activity liquid vs EOS vapor) is out of scope. With `activity_model` the test is liquid-liquid: `ln gamma_i` replaces `ln phi_i`, trials are pure-component-dominant only (no Wilson estimates), and `feed_branch` / `phase_branch` are `None`. `pressure_Pa` stays required and validated but is inert for an activity model (`diagnostics["pressure_dependent"]`, `diagnostics["model_family"]`). `StabilityTrial` gained `ssi_iterations`, `second_order_iterations`, `converged_stage`; `StabilitySettings` gained `second_order`, `ssi_iterations`, `second_order_max_iter`, `second_order_max_step`. Peng-Robinson results are bit-identical to the pre-slice values. (source: src/chemthermo/stability/, .agents/brain/adr/0007-stability-tangent-plane-evaluator.md)
- `NRTL.activity_coefficients` implements the standard Renon-Prausnitz
  equation with column sums. Public signature unchanged, but **returned values
  changed** for asymmetric parameters: before the `nrtl-gibbs-duhem-fix` slice
  the implementation used row sums and per-term denominators and violated
  Gibbs-Duhem (residuals ~1e-1, up to 0.89 off in `ln gamma` versus
  `thermo.NRTL`). Gamma-phi flash results with the packaged synthetic pairs
  shifted slightly (CLI Methane/Ethane vapor fraction 0.767092 -> 0.764835).
  See validation Cases N-1..N-3. (source: src/chemthermo/models/nrtl.py)
- `PengRobinsonEOS.kij` accepts a scalar (off-diagonal only; diagonal always unaffected) or a `Mapping[tuple[str, str], float]` keyed by normalized component-name pairs, default per-pair value `0.0` (ADR-0006). `flash_tp` and `stability_tp` results for nonzero `kij` are now trustworthy (previously the diagonal was silently corrupted; see ADR-0006). (source: src/chemthermo/models/peng_robinson.py, .agents/brain/adr/0006-pr-kij-matrix.md)
- EOS registry module (`chemthermo.eos`: `EOSProtocol`, `PCSAFTEOS`, `get_eos`, `list_eos`, `register_eos`). (source: src/chemthermo/eos/__init__.py)
- VLLE plugin boundary (`chemthermo.vlle`: `get_vlle_engine`, `VLLEEngine`, `VLLEResult`, and related types/errors). (source: src/chemthermo/vlle/__init__.py, README.md)

Implementation status notes
- `chemthermo.vlle`: Public plugin boundary for optional VLLE engines.
- `chemthermo.eos.pcsaft`: Public EOS registry entry and interface with implementation details evolving over time.

CLI entry points
- `chemthermo = "chemthermo.cli:main"` in `[project.scripts]` and module execution via `python -m chemthermo`. (source: pyproject.toml, src/chemthermo/__main__.py, src/chemthermo/cli.py)

## 3) Architecture
Text-only diagram
```
src/chemthermo/data/components.json -> data loaders -> Component/Composition/Mixture -> models (PR/NRTL) -> flash_tp -> FlashResult
src/chemthermo/data/components.json -> Component/Composition/Mixture -> models (PR | NRTL) -> stability_tp -> internal _TangentPlaneEvaluator (ln phi on the min-Gibbs root | ln gamma) -> Michelsen TPD (SSI + Newton) -> StabilityResult
chemthermo CLI -> parser -> Mixture + PengRobinsonEOS -> flash_tp -> text/json output

```

Key modules and flow
- Component databank lives in `src/chemthermo/data/components.json`, loaded via `chemthermo.data` helpers. Optional non-runtime mirror path is `database/components.mirror.json` and is never loaded by runtime code. (source: src/chemthermo/data/__init__.py, src/chemthermo/data/components.json, tools/build_database.py)
- Core domain objects: `Component`, `Composition`, `Mixture`. (source: src/chemthermo/core/component.py, src/chemthermo/core/composition.py, src/chemthermo/core/mixture.py)
- Flash solver (`flash_tp`) orchestrates models and returns `FlashResult`. Its single-phase decision is still a Wilson K-bound / Rachford-Rice heuristic, NOT a stability analysis. (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/results.py)
- Phase stability (`stability_tp`) implements Michelsen's tangent-plane test independently of the flash solver and returns `StabilityResult`. For an EOS it selects the lowest-Gibbs compressibility root generically by minimizing `sum_i w_i ln phi_i(w)` over the `phase="vapor"`/`phase="liquid"` calls of the existing `EquationOfState` interface. The two model families are isolated behind the **internal** (not exported) `_TangentPlaneEvaluator` protocol in `src/chemthermo/stability/_evaluator.py` with adapters `_EOSTangentPlane` and `_ActivityTangentPlane`; the solver, trivial detection, summary and result types never branch on the family (ADR-0007). Each trial is successive substitution followed, if needed, by a damped Newton stage on the stationarity condition in `ln W` (required near plait points, dormant for every validated PR state). (source: src/chemthermo/stability/tp.py, src/chemthermo/stability/_evaluator.py, src/chemthermo/stability/results.py)
- EOS registry provides named EOS factories. (source: src/chemthermo/eos/registry.py)
- Deeper usage docs: `README.md`, `examples/README.md`. (source: README.md, examples/README.md)

Key entry points (top paths)
- `README.md`, `pyproject.toml`, `.github/workflows/ci.yml`, `src/chemthermo/__init__.py`, `src/chemthermo/cli.py`, `src/chemthermo/flash/tp.py`, `src/chemthermo/models/peng_robinson.py`, `src/chemthermo/models/nrtl.py`, `src/chemthermo/eos/registry.py`, `src/chemthermo/parameters/nrtl.py`, `src/chemthermo/data/components.json`, `examples/basic/flash_tp_peng_robinson_demo.py`, `examples/validation/00_reference_case.py`, `tests/test_flash_tp.py`, `tests/test_cli_tp_flash.py`. (source: README.md, pyproject.toml, .github/workflows/ci.yml, src/chemthermo/__init__.py, src/chemthermo/cli.py, src/chemthermo/flash/tp.py, src/chemthermo/models/peng_robinson.py, src/chemthermo/models/nrtl.py, src/chemthermo/eos/registry.py, src/chemthermo/parameters/nrtl.py, src/chemthermo/data/components.json, examples/basic/flash_tp_peng_robinson_demo.py, examples/validation/00_reference_case.py, tests/test_flash_tp.py, tests/test_cli_tp_flash.py)

## 4) Key invariants and assumptions
- SI units everywhere: temperature in K, pressure in Pa. Enforced via validation helpers and documented usage. (source: README.md, src/chemthermo/validation.py, src/chemthermo/flash/tp.py)
- Composition fractions are non-negative and sum to 1 within `COMPOSITION_SUM_TOL` (1e-8) unless `normalize=True`. Enforced by `validate_fractions` and `Composition`. (source: src/chemthermo/validation.py, src/chemthermo/core/composition.py)
- `Mixture` length matches composition length and has at least one component. Enforced in `Mixture.__post_init__`. (source: src/chemthermo/core/mixture.py)
- `flash_tp` requires mole-fraction compositions and an EOS; gamma-phi requires an activity model. Enforced by runtime checks. (source: src/chemthermo/flash/tp.py)
- Flash solver determinism for fixed inputs/settings. Stated in docs. (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/settings.py)
- `FlashResult` phases non-empty; phase fractions in [0,1] sum to 1 within tolerance. Enforced in `FlashResult.__post_init__`. (source: src/chemthermo/flash/results.py)
- **Activity models must be thermodynamically consistent**: `ln gamma` must be
  the composition derivative of a single reduced excess Gibbs energy, so
  `sum_i x_i d ln gamma_i = 0` at fixed T, P (Gibbs-Duhem). For NRTL this is
  enforced by test, not by construction; any future activity model must carry
  the same check. Consistency must be tested with **asymmetric** parameters:
  symmetric binaries are blind to row/column mix-ups. (source:
  src/chemthermo/models/nrtl.py, tests/test_activity_nrtl.py, validation Case N-1)
- NRTL index convention: `tau[i, j] = tau_ij`, `alpha[i, j] = alpha_ij`,
  `G_ij = exp(-alpha_ij tau_ij)`, and every internal sum `S_j = sum_k G_kj x_k`,
  `C_j = sum_k tau_kj G_kj x_k` runs down a **column**. (source:
  src/chemthermo/models/nrtl.py module docstring)
- Packaged NRTL pair parameters are synthetic demo placeholders, not fitted or
  published data; published parameter sets live in `tests/fixtures/` with a
  citation and are never loaded by default. The Tessier (2000) Problem 2 set
  (`tests/fixtures/nrtl/tessier2000_problem2.json`) is third-party fitted data
  regressed from the DECHEMA Chemistry Data Series; its redistribution here is
  limited to that cited test fixture and it must never become packaged runtime
  data. (source: src/chemthermo/parameters/data/activity/nrtl.json,
  tests/fixtures/nrtl/tessier2000_problem1.json,
  tests/fixtures/nrtl/tessier2000_problem2.json)
- Tangent-plane stability is one criterion for both model families: with an
  activity model, `ln gamma_i` replaces `ln phi_i` in `tpd` and nothing else
  changes (the pure-liquid reference cancels between two liquid phases). Any
  new model family must enter through
  `chemthermo.stability._evaluator._TangentPlaneEvaluator`, never by branching
  inside the solver. (source: src/chemthermo/stability/_evaluator.py,
  .agents/brain/adr/0007-stability-tangent-plane-evaluator.md)

## 5) Error handling & validation policy
- Validation helpers raise `InputRangeError` for invalid temperatures/pressures; `CompositionError` for invalid fractions. (source: src/chemthermo/validation.py, src/chemthermo/exceptions.py)
- Model misuse or invalid model outputs raise `ModelError`; flash non-convergence raises `ConvergenceError`. (source: src/chemthermo/exceptions.py, src/chemthermo/flash/tp.py)
- Missing component properties raise `PropertyNotFoundError`. (source: src/chemthermo/core/component.py, src/chemthermo/exceptions.py)


## 6) Configuration & defaults
- Flash defaults: `FlashSettings(max_iter=100, tol=1e-8, damping=None)`. (source: src/chemthermo/flash/settings.py)
- Composition sum tolerance: `COMPOSITION_SUM_TOL = 1e-8`. (source: src/chemthermo/validation.py)
- Unit constants: `R_J_PER_MOL_K`, `STANDARD_T_K`, `STANDARD_P_PA`, pressure conversions. (source: src/chemthermo/units.py)
- NRTL parameters load from packaged JSON (`src/chemthermo/parameters/data/activity/nrtl.json`). (source: src/chemthermo/parameters/nrtl.py, src/chemthermo/parameters/data/activity/nrtl.json)

- No environment-variable configuration is documented in README or `pyproject.toml`. (source: README.md, pyproject.toml)

Cheap checks
- Canonical cheap checks are maintained in `.agents/dev-contract.md` to avoid command drift.
- Use the cheap-check subset first for fast feedback, then run full CI gates from `.agents/dev-contract.md`.

## 7) Testing & CI contract
- CI runs: ruff format check, ruff lint, pyright, pytest on Python 3.11. (source: .github/workflows/ci.yml)
- Optional validation tests compare against the `thermo` library and are skipped if not installed. (source: tests/validation/test_flash_vs_thermo.py, tests/validation/test_stability_vs_thermo.py, pyproject.toml)
- Validation case ledger ("thermodynamics exam"): `.agents/brain/validation-cases.md`. One entry per case with source, location, assumptions, parameters and provenance, expected outcome, tolerance achieved, independent route, and test path. Never record an expected value that was not read from a source or produced by an independent route.

## 8) Decisions log (index)
- ADR folder: `.agents/brain/adr/`
- Accepted ADRs:
  - `.agents/brain/adr/0001-public-api-truth-source.md`
  - `.agents/brain/adr/0002-thin-vertical-slices.md` (Adopted 2026-02-10)
  - `.agents/brain/adr/0003-cli-entrypoint.md` (Adopted 2026-02-10)
  - `.agents/brain/adr/0004-cli-tp-flash-gamma-phi.md` (Adopted 2026-02-12)
  - `.agents/brain/adr/0005-stability-tp-public-api.md` (Adopted 2026-09-13)
  - `.agents/brain/adr/0006-pr-kij-matrix.md` (Adopted 2026-09-13)
  - `.agents/brain/adr/0007-stability-tangent-plane-evaluator.md` (Adopted 2026-09-13)
- ADR rules: one decision per ADR; keep under 1 page; include status and supersedes fields.

## 9) Roadmap: next 3 increments (vertical slices)
- **Recently completed**
  - `stability-tpd-pr`: public `stability_tp` (Michelsen tangent-plane stability) with Peng-Robinson, min-Gibbs root selection, deterministic trial set, golden path and thermo cross-check.
  - CLI gamma-phi extension for `chemthermo tp-flash` via `--flash-mode`.
  - `pr-kij-matrix`: fixed the diagonal-kij bug and added per-pair `kij` support (`float` or name-keyed `Mapping`) to `PengRobinsonEOS`; `flash_tp` and `stability_tp` results for nonzero kij are now trustworthy. See ADR-0006 and validation Case K-1.
  - `stability-tpd-nrtl`: `stability_tp` now accepts `activity_model=` for liquid-liquid tangent-plane stability, behind an internal `_TangentPlaneEvaluator` contract (ADR-0007) that also serves the Peng-Robinson path unchanged; added a damped-Newton second stage (required near plait points), the cited Tessier (2000) Problem 2 fixture, and golden paths `examples/basic/stability_tp_nrtl_lle_demo.py` and `examples/validation/08_stability_nrtl_tessier2000.py`. Reproduces the published tangent-plane global minima of Problems 1 and 2. See validation Cases S-6, S-7, S-8.
  - `nrtl-gibbs-duhem-fix`: corrected the NRTL activity-coefficient equation (column sums, single first-term denominator); added Gibbs-Duhem / binary-reduction / permutation / regression tests, a tight `thermo` cross-check with asymmetric parameters, the cited Tessier (2000) Problem 1 fixture, and the Table 2 reproduction golden path `examples/validation/07_nrtl_tessier_stationary_points.py`. Packaged NRTL pairs are now labelled synthetic. No ADR (public signature unchanged). See validation Cases N-1, N-2, N-3.
- **Slice 1: `flash-auto-phase-detection`**
  - Capability: `flash_tp` decides 1-vs-2 phases from `stability_tp` instead of K-bound heuristics, and seeds K-values from the converged stationary point.
  - Requirements: keep the `FlashResult` shape and CLI JSON contract, add diagnostics for the stability verdict, and prove behavior change only where the heuristic was wrong.
- **Slice 2: (not yet scoped)**
  - To be defined after `flash-auto-phase-detection` lands.

## 10) Open questions / risks
