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
- **Thin Vertical Slice Rule**: Every exposed capability must be end-to-end usable. No scaffolded/partial features.
- **Golden Path Rule**: Every Thin Vertical Slice must ship with at least one golden path (example/test) that executes successfully from a clean environment.
- Keep claims factual; add evidence tags or links for load-bearing statements.

**Skills**
- Local skill catalog: `.agents/skills/README.md`
- Preferred repo workflow skill: `.agents/skills/chemthermo-change-loop/SKILL.md`

## 0) Repo at a glance
- Purpose: chemical engineering thermodynamics utilities packaged as `chemthermo`, SI units throughout. (source: README.md)
- Primary language/toolchain: Python 3.11+, setuptools build via `pyproject.toml`. (source: pyproject.toml)
- Primary entry points: `chemthermo` top-level API, `chemthermo.eos` registry, `chemthermo.vlle` plugin boundary. (source: src/chemthermo/__init__.py, src/chemthermo/eos/__init__.py, src/chemthermo/vlle/__init__.py)
- Tests: `pytest` (plus lint/type checks in CI). (source: .github/workflows/ci.yml)
- Golden Path command is maintained in `.agents/dev-contract.md`. (source: .agents/dev-contract.md)

## 1) Purpose and scope policy
- Purpose: provide core thermodynamics utilities (components, mixtures, EOS/activities, TP flash) in SI units. (source: README.md, src/chemthermo/flash/tp.py, src/chemthermo/models/base.py)
- Scope policy: VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## 2) Public API surface (current)
Definition of public API follows ADR-0001 (source of truth rules in `.agents/brain/adr/0001-public-api-truth-source.md`).

Stable (public) entry points
- `chemthermo` top-level exports in `__all__` (core types, flash API, models, parameters, exceptions, units, validation helpers). (source: src/chemthermo/__init__.py)
- EOS registry module (`chemthermo.eos`: `EOSProtocol`, `PCSAFTEOS`, `get_eos`, `list_eos`, `register_eos`). (source: src/chemthermo/eos/__init__.py)
- VLLE plugin boundary (`chemthermo.vlle`: `get_vlle_engine`, `VLLEEngine`, `VLLEResult`, and related types/errors). (source: src/chemthermo/vlle/__init__.py, README.md)

Implementation status notes
- `chemthermo.vlle`: Public plugin boundary for optional VLLE engines.
- `chemthermo.eos.pcsaft`: Public EOS registry entry and interface with implementation details evolving over time.

CLI entry points
- No CLI scripts are defined in `pyproject.toml` (no `[project.scripts]` section). (source: pyproject.toml)

## 3) Architecture
Text-only diagram
```
components.json -> data loaders -> Component/Composition/Mixture -> models (PR/NRTL) -> flash_tp -> FlashResult

```

Key modules and flow
- Component databank lives in `src/chemthermo/data/components.json`, loaded via `chemthermo.data` helpers. (source: src/chemthermo/data/__init__.py, src/chemthermo/data/components.json)
- Core domain objects: `Component`, `Composition`, `Mixture`. (source: src/chemthermo/core/component.py, src/chemthermo/core/composition.py, src/chemthermo/core/mixture.py)
- Flash solver (`flash_tp`) orchestrates models and returns `FlashResult`. (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/results.py)
- EOS registry provides named EOS factories. (source: src/chemthermo/eos/registry.py)
- Deeper usage docs: `README.md`, `examples/README.md`. (source: README.md, examples/README.md)

Key entry points (top paths)
- `README.md`, `pyproject.toml`, `.github/workflows/ci.yml`, `src/chemthermo/__init__.py`, `src/chemthermo/flash/tp.py`, `src/chemthermo/models/peng_robinson.py`, `src/chemthermo/models/nrtl.py`, `src/chemthermo/eos/registry.py`, `src/chemthermo/parameters/nrtl.py`, `src/chemthermo/data/components.json`, `examples/basic/flash_tp_peng_robinson_demo.py`, `tests/test_flash_tp.py`. (source: README.md, pyproject.toml, .github/workflows/ci.yml, src/chemthermo/__init__.py, src/chemthermo/flash/tp.py, src/chemthermo/models/peng_robinson.py, src/chemthermo/models/nrtl.py, src/chemthermo/eos/registry.py, src/chemthermo/parameters/nrtl.py, src/chemthermo/data/components.json, examples/basic/flash_tp_peng_robinson_demo.py, tests/test_flash_tp.py)

## 4) Key invariants and assumptions
- SI units everywhere: temperature in K, pressure in Pa. Enforced via validation helpers and documented usage. (source: README.md, src/chemthermo/validation.py, src/chemthermo/flash/tp.py)
- Composition fractions are non-negative and sum to 1 within `COMPOSITION_SUM_TOL` (1e-8) unless `normalize=True`. Enforced by `validate_fractions` and `Composition`. (source: src/chemthermo/validation.py, src/chemthermo/core/composition.py)
- `Mixture` length matches composition length and has at least one component. Enforced in `Mixture.__post_init__`. (source: src/chemthermo/core/mixture.py)
- `flash_tp` requires mole-fraction compositions and an EOS; gamma-phi requires an activity model. Enforced by runtime checks. (source: src/chemthermo/flash/tp.py)
- Flash solver determinism for fixed inputs/settings. Stated in docs. (source: src/chemthermo/flash/tp.py, src/chemthermo/flash/settings.py)
- `FlashResult` phases non-empty; phase fractions in [0,1] sum to 1 within tolerance. Enforced in `FlashResult.__post_init__`. (source: src/chemthermo/flash/results.py)

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
- Optional validation tests compare against the `thermo` library and are skipped if not installed. (source: tests/validation/test_flash_vs_thermo.py, pyproject.toml)

## 8) Decisions log (index)
- ADR folder: `.agents/brain/adr/`
- Accepted ADRs:
  - `.agents/brain/adr/0001-public-api-truth-source.md`
- ADR rules: one decision per ADR; keep under 1 page; include status and supersedes fields.

## 9) Roadmap: next 3 increments (vertical slices)
- **Slice 1: Component Database Expansion (End-to-End)**
  - Capability: Users can access a broad, realistic set of components at runtime.
  - Requirements: Regenerate `components.json` from `organics.txt`, add representative tests, add Golden Path example.
- **Slice 2: CLI TP Flash Tool**
  - Capability: Users can run a TP flash from the command line without writing Python.
  - Requirements: Minimal CLI entrypoint, integration test, Golden Path (README/script).
- **Slice 3: Minimal External Validation Slice**
  - Capability: Users can trust numerical correctness for at least one EOS + mixture combination.
  - Requirements: Deterministic validation against external reference, CI or documented optional run, Golden Path validation example.

## 10) Open questions / risks
