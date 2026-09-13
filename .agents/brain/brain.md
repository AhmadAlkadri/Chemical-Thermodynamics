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
- `FlashSettings` gained `phase_detection` (`"tangent-plane"` default, or `"wilson-heuristic"`) and `stability_settings: StabilitySettings | None` (ADR-0008). In phi-phi mode `flash_tp` decides one phase versus two from `stability_tp`, returns a single phase only when the feed is found stable, seeds the split from the tangent-plane minimizer, and raises `ConvergenceError` when stability is inconclusive. New diagnostics keys: `phase_detection`, `stability_status`, `tpd_min`, `feed_branch`, `stability_trials`, `k_seed`, `incipient_phase`, `mass_balance_residual`, `fugacity_residual`, `delta_g_split_rt`; a single-phase tangent-plane result no longer carries `k_min` / `k_max` / `max_delta_k` / `rr_*`. Converged two-phase results are unchanged to 8.6e-7 relative in vapor fraction. Gamma-phi is unchanged and reports `phase_detection = "wilson-heuristic"`. (source: src/chemthermo/flash/, .agents/brain/adr/0008-flash-tangent-plane-phase-detection.md)
- `flash_tp(mixture, *, temperature_K, pressure_Pa, eos=None, activity_model=None, flash_mode=None, settings=None)` (ADR-0009). `eos` is now optional and `flash_mode` defaults to `None` = infer: an activity model with no EOS is the new `"gamma-gamma"` (liquid-liquid) mode, anything else is `"phi-phi"`. Naming a mode explicitly still wins, so `flash_mode="phi-phi"` without an `eos` is a `ModelError`, as is `"gamma-gamma"` with one. A gamma-gamma result names its phases `"liquid1"` (feed-like) / `"liquid2"` (incipient-like) - **roles, not identities; they can swap between feeds on one tie-line** - carries `vapor_fraction=None`, and reports `equilibrium_residual` where phi-phi reports `fugacity_residual`, plus `ssi_iterations`, `second_order_iterations`, `converged_stage`. A stable feed returns one phase named `"liquid"`. `FlashSettings` gained `post_split_stability` (True), `second_order` (True), `ssi_iterations` (50), `second_order_max_iter` (100), `second_order_tol` (1e-12); the second-order stage applies to the liquid-liquid split only, so every phi-phi number is unchanged. Every two-phase result from the tangent-plane phi-phi path and the gamma-gamma path is re-tested phase by phase and reports `post_split_checked`, `post_split_stable`, `post_split_status`, `post_split_tpd_min`, `phase_stability_<name>`, `phase_stability_tpd_min_<name>`; an unstable phase raises `ConvergenceError` ("a third phase is required") unless `post_split_stability=False`. Gamma-phi and the legacy `wilson-heuristic` path cannot run the check and report `post_split_checked=False` with a reason. (source: src/chemthermo/flash/, .agents/brain/adr/0009-flash-liquid-liquid-activity.md)
- `flash_tp(..., flash_mode="modified-raoult")` and `stability_tp(..., activity_model=..., vapor="ideal")` (ADR-0010). A low-pressure gamma-phi model: an activity-coefficient liquid with Antoine pure-liquid reference fugacities (`f_i^0 = Psat_i(T)`, `phi^sat = 1`, Poynting = 1) against an **ideal-gas** vapor, both expressed against the common reference `ln(f_i/(x_i P))` (liquid: `ln gamma_i + ln(Psat_i/P)`; vapor: `0`), so one tangent-plane test detects vapor-liquid *or* liquid-liquid *or* neither and the candidate label says which. `stability_tp` gained `vapor: Literal["none","ideal"] = "none"` (only valid with `activity_model`; `"ideal"` with an `eos` is a `ModelError`), `feed_branch`/`phase_branch` then carry `"liquid"`/`"vapor"`, `model_family` is `"modified-raoult"` and diagnostics gain `antoine_valid_Tmin_K`/`antoine_valid_Tmax_K`. `flash_tp`'s new mode requires `activity_model`, forbids `eos`, and is **never inferred** (an activity model alone still means `gamma-gamma`). VLE returns `"liquid"`/`"vapor"` with a real `vapor_fraction`; LLE returns `"liquid1"`/`"liquid2"` with `vapor_fraction=None`; a stable feed returns one phase named by the feed candidate. Both phases are post-split-tested against **both** candidates. Antoine ranges are enforced (`InputRangeError`), never extrapolated. `flash_mode="gamma-phi"` is **DEPRECATED** in docs only - unchanged, not removed, removal would need its own ADR. (source: src/chemthermo/stability/_evaluator.py, src/chemthermo/flash/_detect.py, src/chemthermo/models/_antoine.py, .agents/brain/adr/0010-stability-phase-candidates-modified-raoult.md)
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
src/chemthermo/data/components.json -> data loaders -> Component/Composition/Mixture -> models (PR/NRTL) -> flash_tp -> stability_tp -> seeded Rachford-Rice/SSI split (+ second-order stage for LLE) -> post-split stability of each phase -> FlashResult
src/chemthermo/data/components.json -> Component/Composition/Mixture -> models (PR | NRTL | NRTL+Antoine+ideal gas) -> stability_tp -> internal _TangentPlaneEvaluator (a set of _PhaseCandidates; the min-Gibbs one wins) -> Michelsen TPD (SSI + Newton) -> StabilityResult
chemthermo CLI -> parser -> Mixture + PengRobinsonEOS -> flash_tp -> text/json output

```

Key modules and flow
- Component databank lives in `src/chemthermo/data/components.json`, loaded via `chemthermo.data` helpers. Optional non-runtime mirror path is `database/components.mirror.json` and is never loaded by runtime code. (source: src/chemthermo/data/__init__.py, src/chemthermo/data/components.json, tools/build_database.py)
- Core domain objects: `Component`, `Composition`, `Mixture`. (source: src/chemthermo/core/component.py, src/chemthermo/core/composition.py, src/chemthermo/core/mixture.py)
- Flash solver (`flash_tp`) orchestrates models and returns `FlashResult`. Reference flow: `flash_tp -> stability_tp -> split -> post-split stability` (ADR-0008, ADR-0009, ADR-0010), for phi-phi (EOS, vapor-liquid), gamma-gamma (activity model, liquid-liquid) and modified-raoult (activity liquid + ideal vapor, VLE *or* LLE) alike. The modified-Raoult split gives each phase the candidate the stability test assigned to it, so `K_i = exp(t_i^x(x) - t_i^y(y))` is `gamma_i Psat_i / P` for a liquid/vapor pair and `gamma_i^I / gamma_i^II` for two liquids - one update rule, two regimes. In phi-phi mode its single-phase decision is Michelsen's tangent-plane test (`flash_tp -> stability_tp -> split`, ADR-0008); the Wilson K-bound / Rachford-Rice heuristic remains reachable through `FlashSettings(phase_detection="wilson-heuristic")` and is what gamma-phi still uses. The split loop (Rachford-Rice plus fixed-point K updates) is shared by both paths and unchanged; only its seed differs. At most two phases are returned; the converged phases **are** now re-tested for stability and a phase set that needs a third phase raises instead of being returned (ADR-0009). The liquid-liquid split adds a second-order stage - a damped Newton minimization of the two-phase Gibbs energy whose gradient is the equal-activity residual - because successive substitution needs 536-3922 iterations near a plait point. `src/chemthermo/flash/tp.py` is a thin public orchestrator (validates inputs, resolves the mode, dispatches); the implementation lives in internal modules `_detect.py` (phase detection + split seeding), `_split.py` (the shared K-loop), `_second_order.py` (the liquid-liquid Newton stage), `_verify.py` (residuals + post-split stability), `_assemble.py` (`FlashResult` construction) and `_legacy.py` (the `wilson-heuristic` path) - internal only per ADR-0001, no public names moved (slice `flash-module-split`). (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/_detect.py, src/chemthermo/flash/_split.py, src/chemthermo/flash/_second_order.py, src/chemthermo/flash/_verify.py, src/chemthermo/flash/_assemble.py, src/chemthermo/flash/_legacy.py, src/chemthermo/flash/settings.py, src/chemthermo/flash/results.py)
- Phase stability (`stability_tp`) implements Michelsen's tangent-plane test independently of the flash solver and returns `StabilityResult`. For an EOS it selects the lowest-Gibbs compressibility root generically by minimizing `sum_i w_i ln phi_i(w)` over the `phase="vapor"`/`phase="liquid"` calls of the existing `EquationOfState` interface. The model families are isolated behind the **internal** (not exported) `_TangentPlaneEvaluator` protocol in `src/chemthermo/stability/_evaluator.py`; the solver, trivial detection, summary and result types never branch on the family (ADR-0007). An evaluator holds one or more `_PhaseCandidate`s and `_select_min_gibbs` keeps the one minimizing `sum_i w_i term_i(w)` - the only candidate-dependent part of `G/RT` - which is one rule for three families: the cubic's two compressibility roots (`_EOSTangentPlane`), a single activity liquid (`_ActivityTangentPlane`, label `None` because nothing was selected) and the modified-Raoult pair (`_ModifiedRaoultTangentPlane`: `ln gamma_i + ln(Psat_i/P)` versus `0`) (ADR-0010). A candidate's `optional` flag separates "this root may not exist here" (recorded and skipped) from "this model failed" (re-raised). Each trial is successive substitution followed, if needed, by a damped Newton stage on the stationarity condition in `ln W` (required near plait points, dormant for every validated PR state). (source: src/chemthermo/stability/tp.py, src/chemthermo/stability/_evaluator.py, src/chemthermo/stability/results.py)
- EOS registry provides named EOS factories. (source: src/chemthermo/eos/registry.py)
- Deeper usage docs: `README.md`, `examples/README.md`. (source: README.md, examples/README.md)

Key entry points (top paths)
- `README.md`, `pyproject.toml`, `.github/workflows/ci.yml`, `src/chemthermo/__init__.py`, `src/chemthermo/cli.py`, `src/chemthermo/flash/tp.py`, `src/chemthermo/models/peng_robinson.py`, `src/chemthermo/models/nrtl.py`, `src/chemthermo/eos/registry.py`, `src/chemthermo/parameters/nrtl.py`, `src/chemthermo/data/components.json`, `examples/basic/flash_tp_peng_robinson_demo.py`, `examples/validation/00_reference_case.py`, `tests/test_flash_tp.py`, `tests/test_cli_tp_flash.py`. (source: README.md, pyproject.toml, .github/workflows/ci.yml, src/chemthermo/__init__.py, src/chemthermo/cli.py, src/chemthermo/flash/tp.py, src/chemthermo/models/peng_robinson.py, src/chemthermo/models/nrtl.py, src/chemthermo/eos/registry.py, src/chemthermo/parameters/nrtl.py, src/chemthermo/data/components.json, examples/basic/flash_tp_peng_robinson_demo.py, examples/validation/00_reference_case.py, tests/test_flash_tp.py, tests/test_cli_tp_flash.py)

## 4) Key invariants and assumptions
- SI units everywhere: temperature in K, pressure in Pa. Enforced via validation helpers and documented usage. (source: README.md, src/chemthermo/validation.py, src/chemthermo/flash/tp.py)
- Composition fractions are non-negative and sum to 1 within `COMPOSITION_SUM_TOL` (1e-8) unless `normalize=True`. Enforced by `validate_fractions` and `Composition`. (source: src/chemthermo/validation.py, src/chemthermo/core/composition.py)
- `Mixture` length matches composition length and has at least one component. Enforced in `Mixture.__post_init__`. (source: src/chemthermo/core/mixture.py)
- `flash_tp` requires mole-fraction compositions and at least one model: an EOS for phi-phi, both for gamma-phi, and an activity model *alone* for gamma-gamma. Enforced by runtime checks. (source: src/chemthermo/flash/tp.py)
- Flash solver determinism for fixed inputs/settings, and permutation invariance under reordering the components (measured: |d beta| = 0.0, worst |d composition| = 2.2e-16). Stated in docs. (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/settings.py, tests/test_flash_phase_detection.py)
- **A phase split must be verified, not just converged.** Every two-phase `flash_tp` result reports `mass_balance_residual`, `fugacity_residual` and `delta_g_split_rt`, and `delta_g_split_rt < 0` is what makes the split an answer rather than a fixed point. Any future split solver must carry the same three. (source: src/chemthermo/flash/tp.py `_verify_split`, validation Case F-3)
- **Vapor/liquid naming is a convention where the model cannot tell.** A single-root (dense or supercritical) single-phase feed gets the min-Gibbs tie-break label; a two-phase result is named by the Wilson volatility *ranking*. `EquationOfState` exposes no molar volume. The ranking must never influence a verdict, a composition or a vapor fraction. (source: src/chemthermo/flash/tp.py, .agents/brain/adr/0008-flash-tangent-plane-phase-detection.md)
- **A two-phase answer is only an answer if both phases are stable.** Every converged two-phase result on a tangent-plane path is re-tested with `stability_tp`; a trial that converges onto the *partner* phase is marginal (one tangent plane, two phases - Case S-3), anything else negative means a third phase and raises. Measured worst post-split `tpd_min` on the in-repo PR grid: -7.0e-09, inside `tpd_tol = 1e-8` by a factor of 1.4. (source: src/chemthermo/flash/tp.py `_post_split_stability`, validation Case L-4)
- **Liquid phase names carry no identity.** `liquid1`/`liquid2` are roles assigned by the split seed, and two feeds on one tie-line may return the same pair of compositions under swapped labels. Any test or caller must compare the phase *set*. There is no volatility ordering to lean on as there is for vapor/liquid. (source: src/chemthermo/flash/tp.py, .agents/brain/adr/0009-flash-liquid-liquid-activity.md)
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
- Tangent-plane stability is one criterion for every model family: only the
  per-component "fugacity term" and the deterministic trial set differ. With an
  activity model alone, `ln gamma_i` replaces `ln phi_i` (the pure-liquid
  reference cancels between two liquid phases); with `vapor="ideal"` the liquid
  carries `ln gamma_i + ln(Psat_i/P)` and the vapor `0`. Any new model family
  must enter through
  `chemthermo.stability._evaluator._TangentPlaneEvaluator` (as a new evaluator
  or as a new `_PhaseCandidate`), never by branching inside the solver.
  (source: src/chemthermo/stability/_evaluator.py,
  .agents/brain/adr/0007-stability-tangent-plane-evaluator.md,
  .agents/brain/adr/0010-stability-phase-candidates-modified-raoult.md)
- **Where several phase descriptions compete, the lowest-Gibbs one is the
  physical one at that composition.** `sum_i w_i term_i(w)` is the only
  candidate-dependent part of `G/RT`, so minimizing it selects the right
  candidate - the min-Gibbs cubic root, or liquid-versus-vapor for the
  modified-Raoult pair. The same rule must serve any future candidate set
  (PC-SAFT density roots, a Gibbs-energy phase model). (source:
  src/chemthermo/stability/_evaluator.py `_select_min_gibbs`, validation
  Case R-2)
- **A vapor-pressure correlation is never extrapolated silently.** Antoine
  evaluation outside a component's stated `[Tmin_K, Tmax_K]` raises
  `InputRangeError`; the mixture's validity window is reported in diagnostics.
  (source: src/chemthermo/models/_antoine.py, validation Case R-2)

## 5) Error handling & validation policy
- Validation helpers raise `InputRangeError` for invalid temperatures/pressures; `CompositionError` for invalid fractions. (source: src/chemthermo/validation.py, src/chemthermo/exceptions.py)
- Model misuse or invalid model outputs raise `ModelError`; flash non-convergence raises `ConvergenceError`. (source: src/chemthermo/exceptions.py, src/chemthermo/flash/tp.py)
- Missing component properties raise `PropertyNotFoundError`, including a component with no Antoine record in `modified-raoult` mode. (source: src/chemthermo/core/component.py, src/chemthermo/models/_antoine.py, src/chemthermo/exceptions.py)


## 6) Configuration & defaults
- Flash defaults: `FlashSettings(max_iter=100, tol=1e-8, damping=None, phase_detection="tangent-plane", stability_settings=None, post_split_stability=True, second_order=True, ssi_iterations=50, second_order_max_iter=100, second_order_tol=1e-12)`. (source: src/chemthermo/flash/settings.py)
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
  - `.agents/brain/adr/0008-flash-tangent-plane-phase-detection.md` (Adopted 2026-09-13)
  - `.agents/brain/adr/0009-flash-liquid-liquid-activity.md` (Adopted 2026-09-13)
  - `.agents/brain/adr/0010-stability-phase-candidates-modified-raoult.md` (Adopted 2026-09-13)
- ADR rules: one decision per ADR; keep under 1 page; include status and supersedes fields.

## 9) Roadmap: next 3 increments (vertical slices)
- **Recently completed**
  - `flash-modified-raoult`: `flash_tp(..., flash_mode="modified-raoult")` and `stability_tp(..., activity_model=..., vapor="ideal")` give thermodynamically consistent **low-pressure vapor-liquid AND liquid-liquid** behavior from one tangent plane, by generalizing the internal evaluator to hold several **phase candidates** and keep the lowest-Gibbs one (the same rule that already selected the min-Gibbs cubic root). The EOS and activity-only paths are unchanged (flash bit-identity fixture and the pinned Peng-Robinson stability numbers untouched). `gamma-phi` is documented as DEPRECATED. See ADR-0010 and validation Cases R-1..R-4. Golden paths `examples/basic/flash_tp_modified_raoult_demo.py` and `examples/validation/10_modified_raoult_water_butanol.py`.
  - `flash-lle-activity`: `flash_tp` accepts `activity_model=` with no `eos` and returns a verified liquid-liquid tie-line discovered (not assumed) by the tangent-plane test, with a second-order stage that minimizes the two-phase Gibbs energy; and every two-phase result on the tangent-plane paths is now re-tested phase by phase, so a state that needs a third phase raises instead of being returned. See ADR-0009 and validation Cases L-1..L-4. Golden paths `examples/basic/flash_tp_nrtl_lle_demo.py` and `examples/validation/09_lle_tessier2000_tie_lines.py`.
  - `flash-auto-phase-detection`: phi-phi `flash_tp` now decides 1-vs-2 phases from `stability_tp` and seeds the split from the tangent-plane minimizer; the legacy heuristic stays reachable as `FlashSettings(phase_detection="wilson-heuristic")`. Verdict agreement with `thermo`'s `FlashVL` over a 175-state PR grid rose from 166/175 to 175/175; 75 states in a wider 1144-state scan moved from `ConvergenceError` to a clean single-phase answer and 2 from a wrong single-phase verdict to a two-phase split, with 0 regressions. See ADR-0008 and validation Cases F-1, F-2, F-3. Golden path `examples/basic/flash_tp_auto_phase_demo.py`.
  - `stability-tpd-pr`: public `stability_tp` (Michelsen tangent-plane stability) with Peng-Robinson, min-Gibbs root selection, deterministic trial set, golden path and thermo cross-check.
  - CLI gamma-phi extension for `chemthermo tp-flash` via `--flash-mode`.
  - `pr-kij-matrix`: fixed the diagonal-kij bug and added per-pair `kij` support (`float` or name-keyed `Mapping`) to `PengRobinsonEOS`; `flash_tp` and `stability_tp` results for nonzero kij are now trustworthy. See ADR-0006 and validation Case K-1.
  - `stability-tpd-nrtl`: `stability_tp` now accepts `activity_model=` for liquid-liquid tangent-plane stability, behind an internal `_TangentPlaneEvaluator` contract (ADR-0007) that also serves the Peng-Robinson path unchanged; added a damped-Newton second stage (required near plait points), the cited Tessier (2000) Problem 2 fixture, and golden paths `examples/basic/stability_tp_nrtl_lle_demo.py` and `examples/validation/08_stability_nrtl_tessier2000.py`. Reproduces the published tangent-plane global minima of Problems 1 and 2. See validation Cases S-6, S-7, S-8.
  - `nrtl-gibbs-duhem-fix`: corrected the NRTL activity-coefficient equation (column sums, single first-term denominator); added Gibbs-Duhem / binary-reduction / permutation / regression tests, a tight `thermo` cross-check with asymmetric parameters, the cited Tessier (2000) Problem 1 fixture, and the Table 2 reproduction golden path `examples/validation/07_nrtl_tessier_stationary_points.py`. Packaged NRTL pairs are now labelled synthetic. No ADR (public signature unchanged). See validation Cases N-1, N-2, N-3.
- **Slice 1: `flash-vlle-phase-addition`**
  - Capability: multiphase flash with phase **addition and removal** (multiphase Rachford-Rice / Michelsen multiphase split) so `FlashResult` can carry more than two phases and the post-split failures of ADR-0009 and ADR-0010 become answers rather than a `ConvergenceError`.
  - Concrete target: water / 1-butanol at 101325 Pa and T3 = 366.2138 K must return three phases, and the ~0.135 K window just below T3 (where the modified-Raoult path currently refuses, because it seeds from the deepest tangent-plane minimum, which is the vapor) must return the two-liquid pair - which needs *removal* of the vapor, not only addition of a liquid. See validation Case R-3.
  - Requirements: extend the `FlashResult` phase-naming contract beyond `liquid`/`vapor`/`liquid1`/`liquid2`, keep the two-phase results of `flash-auto-phase-detection`, `flash-lle-activity` and `flash-modified-raoult` unchanged where no third phase exists, and give the vapor/liquid label a basis better than volatility ordering (ADR-0008 decision 3) if more than two phases can appear. The modified-Raoult path already has that basis (candidate labels) and can be the first consumer.
- **Slice 2: VLLE**
  - After phase addition lands, wire the multi-phase split through the `chemthermo.vlle` boundary.
- **Also open**
  - Full gamma-phi (an EOS vapor against an activity liquid) phase detection, still blocked on a reference fugacity carrying `phi^sat` and a Poynting correction (ADR-0007). ADR-0010 discharges the low-pressure case only. The legacy `flash_mode="gamma-phi"` is deprecated and its removal needs its own ADR.
  - An accelerated / second-order **phi-phi** split: 1 state in the 1144-state scan is weakly unstable and near-critical (`tpd_min = -1.2e-3`) and still exhausts the iteration limit. The second-order machinery now exists (ADR-0009) but is wired to the liquid-liquid path only, deliberately, so that no phi-phi number moved in that slice.

## 10) Open questions / risks
