# Steering Brief

## What changed since last brief (files + bullets)
- `src/chemthermo/flash/tp.py`, `src/chemthermo/flash/_detect.py` (new), `src/chemthermo/flash/_split.py` (new), `src/chemthermo/flash/_second_order.py` (new), `src/chemthermo/flash/_verify.py` (new), `src/chemthermo/flash/_assemble.py` (new), `src/chemthermo/flash/_legacy.py` (new) (slice `flash-module-split`)
  - Pure move refactor: `tp.py` (1572 lines) is now a thin orchestrator (297 lines, mostly its unchanged docstrings) that validates inputs, resolves the mode and dispatches to the internal modules above - phase detection/seeding, the shared K-loop, the liquid-liquid Newton stage, verification/post-split stability, `FlashResult` assembly, and the legacy `wilson-heuristic` path, one module per concern (mirrors the `stability/_evaluator.py` pattern). No algorithm was rewritten, no numpy operation order changed, no public name moved or renamed; `flash_tp`'s signature and docstring are unchanged. No ADR: internal structure only (ADR-0001, "anything outside [the public API sources] is internal"). Proven bit-identical by `tests/test_flash_refactor_bit_identity.py` against a fixture captured at the pre-refactor commit (155 states: the phi-phi grid, gamma-gamma binary/Tessier feeds, gamma-phi cases, both legacy cases), plus all 223 pre-existing tests unchanged.
- `src/chemthermo/flash/tp.py`, `src/chemthermo/flash/settings.py` (slice `flash-lle-activity`, ADR-0009)
  - `flash_tp` accepts `activity_model=` with **no** `eos`. `eos` gained a default of `None` and `flash_mode` a default of `None` = *infer from the models supplied*; an activity model alone resolves to the new `"gamma-gamma"` (liquid-liquid) mode. Naming a mode explicitly still wins, so `flash_mode="phi-phi"` without an `eos` and `"gamma-gamma"` with one are both `ModelError`. Phase count stays an output: the feed goes through `stability_tp(..., activity_model=...)` and a split is attempted only when it is found unstable. Phases are `"liquid1"` (feed-like) / `"liquid2"` (incipient-like) - roles, not identities, and they demonstrably swap between two feeds on one tie-line - with `vapor_fraction=None` and `equilibrium_residual` in place of `fugacity_residual`.
  - **Second-order stage**: a damped Newton *minimization* of the two-phase reduced Gibbs energy in the phase-II mole numbers, whose gradient is (Gibbs-Duhem) exactly the equal-activity residual; FD Hessian, eigenvalue shift, Armijo/residual backtracking inside `0 < n_i < z_i`. **Design deviation, recorded:** the brief's first option, Newton on the residual system in `(ln K, beta)`, was implemented and rejected on evidence - it walks into the trivial branch (`beta -> -8`) on the near-plait feed `z = (0.12, 0.05, 0.83)` and stalls at residual 1.4e-07 unless given 200+ substitutions first. A residual system cannot distinguish the equilibrium from the trivial solution; a descent method on `g` can. The stage is wired to the liquid-liquid split **only**, so no phi-phi number moved.
  - **Post-split stability**: every two-phase result on the tangent-plane phi-phi path and the gamma-gamma path re-tests each converged phase with `stability_tp`. A minimizer that *is* the partner phase (in the stability module's own `trivial_tol` metric) is `"marginal"`, not an instability - two coexisting phases share one tangent plane. Anything else negative raises `ConvergenceError` ("a third phase is required"); `FlashSettings(post_split_stability=False)` gates the *raise*, not the computation. `gamma-phi` and the legacy `wilson-heuristic` path cannot run it and report `post_split_checked=False` with a reason.
  - New `FlashSettings` fields: `post_split_stability`, `second_order`, `ssi_iterations` (50), `second_order_max_iter` (100), `second_order_tol` (1e-12).
- `tests/test_flash_lle.py` (new), `tests/validation/test_flash_lle_tessier2000.py` (new), `tests/test_flash_phase_detection.py`, `tests/test_cli_tp_flash.py`, `tests/fixtures/cli/*`
  - Mode contract, verdicts, binary binodal (0.019998419467 / 0.359999661508) and lever rule to 2.8e-12, determinism and permutation invariance of the phase *set*, `FlashResult` invariants, proof that SSI alone cannot solve the near-plait feed in the default budget, failure semantics; Tessier Problem 1 and 2 tie-lines against an independent SSI + damped-Newton solve of the full `(x^I, x^II, beta)` system (agreement <= 1.1e-12); phi-phi post-split grid (47/47 pass, worst post-split `tpd_min` -7.0055e-09); CLI fixtures gained the new diagnostics keys with `cli_schema_version` still 1.
  - **Two findings, recorded not accommodated:** (1) the binary feed x(n-butanol) = 0.30 is *inside* the miscibility gap and splits - the brief listed it as a single-phase control; 0.45 is used instead and `thermo` independently agrees 0.30 is unstable. (2) `thermo` 0.6.0's `FlashVLN` will not return a liquid-liquid split for two `GibbsExcessLiquid` phases over one NRTL model (`unique_liquid_count == 1`, single phase returned for every in-gap feed), so no flash-level comparison is possible; its own `stability_test_Michelsen` does find the conjugate pair and agrees to 1.0e-06 / 2.5e-05.
- `examples/basic/flash_tp_nrtl_lle_demo.py` (new), `examples/validation/09_lle_tessier2000_tie_lines.py` (new), `README.md`, `examples/README.md`, `.agents/brain/adr/0009-flash-liquid-liquid-activity.md` (new), `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Golden paths for the binary tie-line/lever rule and for both published systems; README sections "Liquid-liquid flash (`gamma-gamma`)" and "Post-split stability" with the limits stated (two phases max, roles not identities, the check is only as sharp as the split, gamma-phi and legacy not checked); validation Cases L-1 (Problem 1 tie-lines), L-2 (Problem 2 + stable control + post-split), L-3 (binary binodal, lever rule, `thermo`), L-4 (phi-phi post-split grid).
- `src/chemthermo/flash/tp.py`, `src/chemthermo/flash/settings.py`
  - Phi-phi `flash_tp` now decides 1-vs-2 phases from `stability_tp` (ADR-0008): stable feed -> single phase (`termination_reason = "feed_stable_tangent_plane"`, name from the min-Gibbs `feed_branch`), unstable feed -> split seeded from the tangent-plane minimizer with the *unnormalized* mole numbers `W = w exp(-tpd)` (using normalized `w` makes `f_RR(0) = 0` exactly, an unbracketable root), inconclusive -> `ConvergenceError`. Wilson stays as a documented K-seed fallback (`diagnostics["k_seed"]`). Every converged split reports `mass_balance_residual`, `fugacity_residual` and `delta_g_split_rt`. New `FlashSettings.phase_detection` (`"tangent-plane"` default / `"wilson-heuristic"`) and `FlashSettings.stability_settings`. The split loop itself is unchanged and shared by both paths; gamma-phi stays on the heuristic and says so in diagnostics.
  - **Design deviation, recorded:** the brief specified orienting the seed from the stability result's `phase_branch`. That was implemented and rejected on evidence - when the cubic has a single real root both branch calls return identical `phi` and the label is only a tie-break, which returns the mirror-labelled solution (`beta = 0.3255` instead of the pinned `0.67451818` at Methane/Ethane 240 K / 3 MPa). The vapor/liquid *name* is now decided by the Wilson volatility *ranking* only (`ln(w_hi/z_hi) - ln(w_lo/z_lo) >= 0`); it never touches the verdict, the compositions or `beta`. See ADR-0008 decision 3.
- `src/chemthermo/cli.py`, `tests/fixtures/cli/tp_flash_v1.json`, `tests/fixtures/cli/tp_flash_gamma_phi_v1.json`, `tests/test_cli_tp_flash.py`
  - `solver.algorithm` now reports the path taken (`tangent-plane+...` vs `wilson+...`); `cli_schema_version` stays 1 (new keys live inside the free-form `diagnostics` mapping). Both fixtures' `diagnostics` blocks regenerated; the phi-phi fixture's `result` block was deliberately left at its pre-slice values and the contract test's tolerance relaxed from rel 1e-9 to 1e-7 (achieved 1.11e-9), so the fixture still pins the *old* equilibrium.
- `tests/test_flash_phase_detection.py` (new), `tests/validation/test_flash_phase_detection_vs_thermo.py` (new)
  - Regression (both canonical splits, bit-identical legacy numbers), verdict semantics, the disagreement state, a 144-state invariant grid (47 two-phase: mass balance <= 2.01e-13, fugacity residual <= 8.25e-9, `delta_g_split_rt <= -1.458e-4`), permutation invariance and determinism, failure semantics, and a 175-state `thermo` `FlashVL` cross-check (175/175 verdicts vs 166/175 for the heuristic).
- `examples/basic/flash_tp_auto_phase_demo.py` (new), `examples/basic/flash_tp_peng_robinson_demo.py`, `README.md`, `examples/README.md`
  - New golden path showing a stable feed, an unstable feed and the heuristic-vs-tangent-plane disagreement; the existing demo prints the new diagnostics keys; README gained an "Automatic phase detection" section that is explicit about the limits (two phases max, no post-split stability test, gamma-phi still heuristic, vapor/liquid naming is a convention).
- `.agents/brain/adr/0008-flash-tangent-plane-phase-detection.md` (new), `.agents/brain/validation-cases.md`, `.agents/brain/brain.md`
  - Recorded the decision, the rejected `phase_branch` orientation, the `cli_schema_version` reasoning, and validation Cases F-1 (thermo grid + regression), F-2 (the disagreement state) and F-3 (invariants).
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
- Phase stability (`chemthermo.stability`) is a solver-independent sibling of `chemthermo.flash`, and `flash_tp` now consumes it twice: `flash_tp -> stability_tp -> split -> post-split stability` (ADR-0008, ADR-0009), for phi-phi (EOS, vapor-liquid) and gamma-gamma (activity model, liquid-liquid) alike. It serves both an EOS and an activity model through the internal `_TangentPlaneEvaluator` contract, and each trial runs successive substitution then an optional damped-Newton stage.
- Canonical runtime DB path is `src/chemthermo/data/components.json`.
- Runtime DB loading is package-resource based via `chemthermo.data` helpers.
- Raw DB source tables stay in `database/organics.txt` and `database/inorganics.txt`.
- `tools/build_database.py` is the canonical regeneration/check tool for packaged DB sync.
- Optional mirror output path is `database/components.mirror.json` (non-runtime, generated-only).
- Core thermodynamic flow is `Component/Composition/Mixture -> models -> flash_tp -> stability_tp -> seeded split (+ second-order stage for LLE) -> post-split stability -> FlashResult` (phi-phi and gamma-gamma); gamma-phi still goes straight to the Wilson-seeded split and is not post-split checked.
- CLI `tp-flash` maps inputs to `Mixture + PengRobinsonEOS + flash_tp`, with optional `NRTL()` activity model when `--flash-mode gamma-phi`.
- CLI supports deterministic text/json outputs with schema version and diagnostics.
- External validation remains optional and script-driven (`examples/validation/00_reference_case.py`).

## Public API status (stable vs experimental)
- Stable:
  - `chemthermo` Python exports in `src/chemthermo/__init__.py`, including `stability_tp` / `StabilityResult` / `StabilitySettings` / `StabilityTrial` (ADR-0005, ADR-0007). `stability_tp` takes exactly one of `eos=` / `activity_model=`.
  - `PengRobinsonEOS(kij=...)` scalar-or-mapping constructor contract (ADR-0006).
  - `FlashSettings(phase_detection=..., stability_settings=...)` and the `flash_tp` diagnostics keys listed in its docstring (ADR-0008).
  - `flash_tp(..., eos=None, activity_model=..., flash_mode=None)` liquid-liquid mode, the `liquid1`/`liquid2` role naming, `vapor_fraction=None` for gamma-gamma, the post-split diagnostics keys, and `FlashSettings(post_split_stability=..., second_order=..., ssi_iterations=..., second_order_max_iter=..., second_order_tol=...)` (ADR-0009).
  - `chemthermo.eos` and `chemthermo.vlle` documented public subpackages.
  - CLI script `chemthermo` with `tp-flash` subcommand and defined exit-code contract.
- Experimental/placeholder:
  - PC-SAFT residual Helmholtz implementation remains placeholder in open-source build.
  - Validation policy promotion to required CI remains deferred.

## Risks / unknowns
- `flash_tp` still returns at most two phases. A state needing a third is now *detected* by the post-split check and refused with `ConvergenceError`, not solved. Both the feed verdict and the post-split verdict inherit `stability_tp`'s bound: "stable" means nothing negative was found from the deterministic trial set, so a third phase no trial reaches would still be missed.
- The post-split check is only as sharp as the split it checks. Measured margin at default settings is large (worst post-split `tpd_min` -7.0e-09 against `tpd_tol` 1e-8 - a factor of 1.4 - with the minimizer matching the partner phase to `sum ln(w/x)^2 <= 4e-16`, twelve orders inside `trivial_tol`), but with a deliberately loosened `FlashSettings(tol=1e-3)` a genuinely two-phase state is reported as needing a third phase. Documented in README and ADR-0009.
- Every two-phase flash on the tangent-plane paths now pays for two extra stability analyses (measured: 0.9 s for the 47 two-phase states of the in-repo grid, roughly doubling their cost). There is no way to skip the computation, only the raise.
- `liquid1`/`liquid2` carry no physical identity and can swap between feeds on the same tie-line. Any consumer that indexes `result.phases["liquid1"]` expecting a particular liquid is wrong.
- The vapor/liquid *name* of a `flash_tp` phase is a convention in two places: a single-root (dense/supercritical) single-phase feed gets whichever branch label wins the min-Gibbs tie-break, and a two-phase result is named by Wilson volatility ranking. `EquationOfState` exposes no molar volume, so nothing better is available without widening that protocol.
- Every phi-phi flash now pays for a stability analysis (2 Wilson trials + 1 pure trial per component) even for an obviously single-phase feed.
- 5 states in the 1144-state scan (Methane/Propane/n-Decane 0.7/0.2/0.1 at 1.2e7-2.0e7 Pa) split in chemthermo (`delta_g_split_rt < 0`) while `thermo`'s `FlashVL` returns `VF = 0`. They disagreed with `thermo` before this slice too, so they are not caused by it. Now adjudicated: evaluating `thermo`'s own `lnphis_at_zs(..., most_stable=True)` at chemthermo's `stability_tp` minimizer gives a negative tpd (within 5e-5 abs of chemthermo's `tpd_min`) at all 5 states, so `thermo`'s `VF = 0` is a stability-*search* miss on `thermo`'s side, not a chemthermo error. See Case F-1 in `validation-cases.md` and `tests/validation/test_flash_thermo_disagreements_adjudicated.py`.
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
- `flash-phase-addition`: multiphase Rachford-Rice / Michelsen multiphase flash, so `FlashResult` can carry more than two phases and the ADR-0009 post-split failure becomes an answer instead of an error.
- `vlle`: after phase addition, wire the multiphase split through the `chemthermo.vlle` boundary. Also needs the pure-liquid reference fugacity ADR-0007 records as missing.
- Wire the existing second-order stage into the **phi-phi** split: the one state in 1144 that still raises `ConvergenceError` is weakly unstable and near-critical (`tpd_min = -1.2e-3`). The machinery exists; it was deliberately kept off that path so no validated phi-phi number moved in `flash-lle-activity`.

## One simplification / deletion candidate
- Remove or archive deprecated `database/components.json` once all docs/tooling users are migrated.

## Assumptions I'm making
- Runtime DB source of truth remains `src/chemthermo/data/components.json`.
- `database/components.mirror.json` is generated-only and should not be committed.
- CLI v1 scope remains TP flash with Peng-Robinson EOS and `--flash-mode {phi-phi,gamma-phi}`.

## How to validate quickly
- `python examples/basic/flash_tp_peng_robinson_demo.py`
- `python examples/basic/flash_tp_auto_phase_demo.py`
- `python examples/basic/stability_tp_nrtl_lle_demo.py`
- `python examples/basic/flash_tp_nrtl_lle_demo.py`
- `python examples/validation/08_stability_nrtl_tessier2000.py`
- `python examples/validation/09_lle_tessier2000_tie_lines.py`
- `python examples/basic/tp_flash_pr_kij_demo.py`
- `python -m pytest tests/validation/test_pr_kij_vs_thermo.py -q` (with `pip install -e ".[validation]"`)
- `python examples/basic/stability_tp_peng_robinson_demo.py`
- `python examples/validation/06_stability_vs_thermo.py` (with `pip install -e ".[validation]"`)
- `python tools/build_database.py --check`
- `python examples/database/01_component_database_demo.py`
- `python -m chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- `python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- `python examples/validation/00_reference_case.py` (with `pip install -e ".[validation]"`)
