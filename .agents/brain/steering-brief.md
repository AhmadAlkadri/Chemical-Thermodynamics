# Steering Brief

## What changed since last brief (files + bullets)
- `src/chemthermo/stability/_evaluator.py` (new), `src/chemthermo/stability/tp.py`, `settings.py`, `results.py`
  - `stability_tp` now takes `eos=None, activity_model=None` (exactly one required; both raises `ModelError` because combined gamma-phi stability is out of scope, ADR-0007). Introduced the INTERNAL `_TangentPlaneEvaluator` protocol (`ln_fugacity_terms(w) -> (ndarray, str | None)`, `initial_estimates(z, active)`) with `_EOSTangentPlane` (ln phi on the min-Gibbs root, Wilson + pure trials) and `_ActivityTangentPlane` (ln gamma, pure-component trials only); the solver, trivial detection, summary and result types no longer know the model family. Added a damped-Newton second stage in `ln W` (FD Jacobian, step cap, backtracking line search) after `settings.ssi_iterations`; new settings `second_order`, `ssi_iterations` (50), `second_order_max_iter`, `second_order_max_step`; new trial fields `ssi_iterations`, `second_order_iterations`, `converged_stage`; `feed_branch`/`phase_branch` are `None` for activity models. Peng-Robinson results are bit-identical (pinned to 1e-12 with per-trial iteration counts).
- `tests/fixtures/nrtl/tessier2000_problem2.json` (new), `tests/conftest.py`
  - Cited Tessier (2000) Table 4 parameters for n-propanol / n-butanol / benzene / water. Table 4 prints G and tau, not alpha; recovered alpha is symmetric to 2.597e-05 and the rounded symmetric average reproduces the printed G to 4.413e-06. Third-party DECHEMA-regressed data: test fixture only, never packaged runtime data.
- `tests/test_stability_activity.py`, `tests/validation/test_stability_nrtl_tessier2000.py`, `tests/test_stability_tp.py`
  - Activity-path invariants (tangent-plane identity, stationarity + FD gradient of tm, tm/tpd/sum(W), permutation invariance, determinism, degenerate feeds), negative controls (pure component, ideal solution, n-butanol/water LLE with an independent binodal solve), a proof that SSI alone cannot solve the near-plait feed, the Tessier Problem 1 / Problem 2 reproduction, a `thermo` 0.6.0 ln-gamma cross-check (max |dD| = 4.7e-16), and a bit-level PR regression pin.
- `examples/basic/stability_tp_nrtl_lle_demo.py`, `examples/validation/08_stability_nrtl_tessier2000.py`, `README.md`, `examples/README.md`
  - Golden paths for the activity-model stability verdict and for the published-minima reproduction (neither needs `thermo`), plus README coverage of the activity usage and the gamma-phi "not supported" note.
- `.agents/brain/adr/0007-stability-tangent-plane-evaluator.md`, `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Recorded the evaluator contract (and why it stays internal), and added validation Cases S-6 (Problem 1 minima), S-7 (Problem 2 minima + the stable control + one disputed printed D) and S-8 (n-butanol/water LLE control).
- `src/chemthermo/models/nrtl.py`
  - Corrected the NRTL activity-coefficient equation. The previous code computed `S = G @ x` (row sums) and divided the first term term-by-term; the standard Renon-Prausnitz form needs the column sums `S_j = sum_k G_kj x_k`, `C_j = sum_k tau_kj G_kj x_k` and a single denominator `S_i` in the first term. The implementation is now vectorized (`G.T @ x`, `(tau * G).T @ x`) and the module docstring derives the equation and states the index convention (`tau[i, j] = tau_ij`, `G_ij = exp(-alpha_ij tau_ij)`). Public signature and the single-component `[1.0]` shortcut unchanged; returned values change for asymmetric parameters (pre-fix: Gibbs-Duhem residuals ~1e-1, up to 0.89 off in `ln gamma` versus `thermo.NRTL`).
- `src/chemthermo/parameters/data/activity/nrtl.json`, `src/chemthermo/parameters/nrtl.py`
  - Labelled the two packaged pairs (Methane/Ethane, Benzene/Water) honestly: new top-level `provenance` block and per-pair `"source": "synthetic-demo"`. They are illustrative placeholders added ad hoc in commits e5ccd8f/ecd476e, not fitted or published data. The loader already ignored unknown keys, so no `schema_version` change; that tolerance is now documented in the `NRTLParameters` docstring together with the index convention.
- `tests/fixtures/nrtl/tessier2000_problem1.json`, `tests/conftest.py`
  - Cited published parameter set (Tessier, Brennecke & Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1; parameters originally McDonald & Floudas, AIChE J. 41 (1995) 1798) plus Table 2 stationary points, deliberately kept out of the packaged defaults and loaded via `NRTLParameters.from_pairs`. New `tests/conftest.py` exposes it as session fixtures.
- `tests/test_activity_nrtl.py`, `tests/validation/test_nrtl_tessier2000.py`
  - Gibbs-Duhem (4 compositions x 4 simplex directions, worst residual 2.07e-10), binary reduction against independently written two-component formulas (1.1e-16), permutation invariance, symmetric-ternary and zero-tau limits, a hard-coded regression guard against the row-sum bug, a `thermo` 0.6.0 cross-check with the asymmetric Tessier parameters (max |d ln gamma| 8.88e-16), and the Table 2 reproduction with a from-scratch tangent-plane distance and damped-Newton stationary-point solve.
- `examples/validation/07_nrtl_tessier_stationary_points.py`, `examples/README.md`, `README.md`
  - Golden path reproducing Table 2 (no optional dependency), plus a README "NRTL activity coefficients" section carrying the equation, the correctness note and the synthetic-parameter disclosure.
- `tests/fixtures/cli/tp_flash_gamma_phi_v1.json`
  - Regenerated: the corrected gammas move the CLI gamma-phi Methane/Ethane result (vapor fraction 0.7670920475829917 -> 0.7648352545438684). Schema, phase names, iteration count and termination reason unchanged.
- `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Added validation Cases N-1 (Gibbs-Duhem + binary reduction), N-2 (`thermo` cross-check) and N-3 (Tessier Table 2, including two printed D values documented as typographical errors rather than accommodated), and a new brain invariant that activity models must satisfy Gibbs-Duhem and must be tested with asymmetric parameters.
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
- `NRTL` implements the standard Renon-Prausnitz equation (column sums); it satisfies Gibbs-Duhem to ~2e-10 and matches `thermo` to ~9e-16 with asymmetric parameters. Packaged pair parameters are synthetic placeholders; published sets live in `tests/fixtures/`.
- `PengRobinsonEOS.kij` is a scalar (off-diagonal only) or a name-keyed per-pair `Mapping`; `flash_tp` and `stability_tp` results for nonzero kij are now trustworthy (ADR-0006).
- Phase stability (`chemthermo.stability`) is a solver-independent sibling of `chemthermo.flash`; `flash_tp` does NOT consume it yet. It serves both an EOS and an activity model through the internal `_TangentPlaneEvaluator` contract, and each trial runs successive substitution then an optional damped-Newton stage.
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
  - `chemthermo` Python exports in `src/chemthermo/__init__.py`, including `stability_tp` / `StabilityResult` / `StabilitySettings` / `StabilityTrial` (ADR-0005, ADR-0007). `stability_tp` takes exactly one of `eos=` / `activity_model=`.
  - `PengRobinsonEOS(kij=...)` scalar-or-mapping constructor contract (ADR-0006).
  - `chemthermo.eos` and `chemthermo.vlle` documented public subpackages.
  - CLI script `chemthermo` with `tp-flash` subcommand and defined exit-code contract.
- Experimental/placeholder:
  - PC-SAFT residual Helmholtz implementation remains placeholder in open-source build.
  - Validation policy promotion to required CI remains deferred.

## Risks / unknowns
- `stability_tp` reporting "stable" is bounded by its deterministic trial set; it is not a global proof, and the docs must keep saying so. Measured limit: on Tessier (2000) Problem 2 the trial set reaches 7 of the 10 non-trivial published stationary points; the three it misses all have D > 0 (Cases S-6, S-7).
- Combined gamma-phi stability (activity liquid vs EOS vapor) raises `ModelError`. It stays unsupported until there is a consistent pure-liquid reference fugacity (ADR-0007); the current gamma-phi flash does not carry one correctly.
- The second-order stage uses a central-difference Jacobian of `ln phi`/`ln gamma`. For a cubic EOS the min-Gibbs branch can switch between finite-difference probes, which would make that Jacobian noisy; it has not been observed because every validated PR state converges inside the 50-iteration SSI budget and never enters the stage. Watch for it if PR states that need the stage ever appear.
- chemthermo uses the rounded PR constants 0.45724 / 0.07780 while `thermo` uses the exact roots; this bounds external agreement at a few times 1e-4 in ln(phi) (was ~2e-4 at the stability states, ~4-6e-4 at some nonzero-kij states -- see Case K-1).
- The illustrative kij value used in docs/examples/tests (0.0411, Methane/n-Decane) is explicitly NOT a validated literature parameter -- do not let it drift into being read as one.
- `thermo`'s own `CEOSLiquid`/`CEOSGas` root solver was observed to be numerically order-sensitive (not permutation-invariant) at some near-critical-locus states during `pr-kij-matrix` development; avoid cross-checking permutation invariance against `thermo` at such states (use chemthermo-internal invariance checks instead, as `tests/validation/test_pr_kij_vs_thermo.py` now does).
- CLI JSON contract must remain backward-compatible or version-bumped.
- DB builder currently assumes raw table column conventions stay unchanged.
- Optional validation still depends on external `thermo` package availability.
- (Removed) Unpinned `bibtexparser>=1.4.0` allowed a fresh/non-editable install to resolve `bibtexparser` 2.x, whose removed `bparser`/`customization` modules broke `import chemthermo`; this was invisible locally because the dev venv already had 1.4.4 installed. Now pinned to `>=1.4.0,<2` and covered by `tests/test_packaging_constraints.py`.

## Next 3 recommended actions
- `flash-auto-phase-detection`: let `flash_tp` consume `stability_tp` instead of the K-bound heuristic, and seed K-values from the converged stationary point.
- (not yet scoped): define after `flash-auto-phase-detection` lands.

## One simplification / deletion candidate
- Remove or archive deprecated `database/components.json` once all docs/tooling users are migrated.

## Assumptions I'm making
- Runtime DB source of truth remains `src/chemthermo/data/components.json`.
- `database/components.mirror.json` is generated-only and should not be committed.
- CLI v1 scope remains TP flash with Peng-Robinson EOS and `--flash-mode {phi-phi,gamma-phi}`.

## How to validate quickly
- `python examples/basic/stability_tp_nrtl_lle_demo.py`
- `python examples/validation/08_stability_nrtl_tessier2000.py`
- `python examples/basic/tp_flash_pr_kij_demo.py`
- `python -m pytest tests/validation/test_pr_kij_vs_thermo.py -q` (with `pip install -e ".[validation]"`)
- `python examples/basic/stability_tp_peng_robinson_demo.py`
- `python examples/validation/06_stability_vs_thermo.py` (with `pip install -e ".[validation]"`)
- `python tools/build_database.py --check`
- `python examples/database/01_component_database_demo.py`
- `python -m chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- `python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- `python examples/validation/00_reference_case.py` (with `pip install -e ".[validation]"`)
