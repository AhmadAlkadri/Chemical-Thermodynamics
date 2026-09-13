# ADR-0008: `flash_tp` decides 1-vs-2 phases from a tangent-plane stability test

Status: accepted
Date: 2026-09-13

## Context
`flash_tp` declared "single phase" from two heuristics on an *initial estimate*:

1. all Wilson K-values `<= 1` (call it liquid) or all `>= 1` (call it vapor);
2. no Rachford-Rice root for the Wilson K-values (call it liquid if `f(0) < 0`,
   vapor otherwise).

Neither is a thermodynamic criterion. Both can mis-classify a feed whose Wilson
estimate is poor, and a converged two-phase solution was never checked for
anything at all. ADR-0005 recorded this as the motivation for building
`stability_tp`, and named this slice as the one that would consume it.

Measured at HEAD `9a714c9` over 1144 Peng-Robinson states (8 databank mixtures,
T 150-450 K, P 1e5-3.2e7 Pa, `kij = 0`):

- **76 states** where the heuristic starts a split that never converges and
  `flash_tp` raises `ConvergenceError`. On **75** of them the feed is in fact
  stable, and the tangent-plane path returns a clean single-phase result;
  `thermo` `FlashVL` with the same constants agrees single-phase on all 75. The
  76th feed is genuinely unstable and near-critical, and still fails.
- **2 states** where the heuristic returns a *single phase* (`rr_no_root`) for a
  feed that is demonstrably two-phase (validation Case F-2).

## Decision

### 1. The phi-phi reference path runs stability first
`flash_tp -> stability_tp(feed) -> single phase | seeded split`.

- `status == "stable"` -> single-phase `FlashResult`, phase name taken from the
  stability result's `feed_branch` (the minimum-Gibbs compressibility root of
  the feed), `vapor_fraction` 1.0 / 0.0 to keep the `FlashResult` contract, and
  `termination_reason = "feed_stable_tangent_plane"`.
- `status == "unstable"` -> the converged stationary point seeds the K-values
  and the existing successive-substitution / Rachford-Rice loop (Michelsen's
  recommended first-stage phase split) runs unchanged.
- `status == "inconclusive"` -> `ConvergenceError`. A stability search that
  could not converge must not silently produce a single-phase answer.

### 2. The seed uses the *unnormalized* trial mole numbers
With the normalized incipient composition `w`, `K_i = w_i / z_i` gives
`f_RR(0) = sum_i w_i - 1 = 0` exactly - a degenerate root that bisection cannot
bracket. Michelsen's unnormalized mole numbers follow from the stationary-point
identity `tpd = -ln(sum_i W_i)`, so `W = w * exp(-tpd)` and
`f_RR(0) = sum_i W_i - 1 > 0` for an unstable feed. Both seeds have the same
fixed point; only the bracket differs. If the seeded K still yields no root the
Wilson estimate is tried as a documented fallback and `diagnostics["k_seed"]`
records which seed was used; if neither brackets a root for an unstable feed,
`ConvergenceError` is raised rather than a single-phase answer returned.

### 3. Which converged phase is *named* "vapor" is decided by volatility ordering
The successive-substitution loop is symmetric under swapping the two phases
(`x <-> y`, `beta <-> 1 - beta`, `K <-> 1/K`), so the seed's orientation decides
the labels. The design brief for this slice specified using the stability
result's `phase_branch` for that. **That rule was implemented, measured and
rejected**: whenever the cubic has a single real root at the trial composition -
which is common - both `phase="vapor"` and `phase="liquid"` return identical
fugacity coefficients and the branch label is only the tie-break order inside
`_ln_phi_min_gibbs`. Concretely, at the canonical Methane/Ethane 240 K / 3 MPa
feed the stability minimizer is `w = (0.150, 0.850)`, an ethane-rich *liquid*
with a single root at `Z = 0.0914`, yet `phase_branch` reads `"vapor"`. Seeding
from that label returns the mirror-labelled solution, `beta = 0.3255` instead of
the pinned `0.67451818`.

The label is therefore decided by the only volatility ordering the package has.
Writing `hi` for the feed component with the largest Wilson K and `lo` for the
smallest, the incipient phase is named vapor-like when
`ln(w_hi / z_hi) - ln(w_lo / z_lo) >= 0`. Only the *ranking* is used, never the
magnitudes, and it decides the name only - never the verdict, the compositions
or the vapor fraction. `EquationOfState` exposes no molar volume, so no
density-based identification is available. `diagnostics["incipient_phase"]`
records the outcome.

### 4. The converged split is verified and the residuals are reported
`mass_balance_residual` (`max_i |z_i - (beta y_i + (1-beta) x_i)|`),
`fugacity_residual` (`max_i |ln(x_i phi_i^L) - ln(y_i phi_i^V)|`) and
`delta_g_split_rt` are computed and recorded. The Gibbs term uses the
minimum-Gibbs branch for the feed and the branch the solver converged on for
each phase; since the min-Gibbs branch is never higher in Gibbs energy, that
combination is a conservative (upper-bound) estimate of the reduction, so
`delta_g_split_rt < 0` remains a valid proof of improvement. The phases
themselves are **not** re-tested for stability in this slice.

### 5. `FlashSettings` gains two fields
`phase_detection: str = "tangent-plane"` (or `"wilson-heuristic"`) and
`stability_settings: StabilitySettings | None = None`. Both are keyword fields
appended after `damping`, so positional construction is unaffected.

### 6. Gamma-phi stays on the heuristic
Its diagnostics report `phase_detection = "wilson-heuristic"` and its numbers
are bit-identical to the pre-slice values. A gamma-phi stability test needs a
consistent pure-liquid reference fugacity that this package does not carry
(ADR-0007); refusing is better than a silently wrong tangent plane.

### 7. `cli_schema_version` stays 1
ADR-0003/ADR-0004 require a bump for a *structural* break. `diagnostics` is a
free-form mapping whose keys are documented as implementation details, so adding
keys inside it removes nothing and changes no type. The one v1 value that did
change is `solver.algorithm`, which now reports the path actually taken
(`"tangent-plane+rachford-rice+fixed-point"` for phi-phi,
`"wilson+rachford-rice+fixed-point"` for gamma-phi and for the legacy mode); the
key, its type and every other field are unchanged. Asserted in
`tests/test_cli_tp_flash.py::test_cli_tp_flash_json_diagnostics_carry_the_phase_detection_keys`.

## What changes for users
- **Diagnostics.** New keys `phase_detection`, `stability_status`, `tpd_min`,
  `feed_branch`, `stability_trials`, `k_seed`, `incipient_phase`,
  `mass_balance_residual`, `fugacity_residual`, `delta_g_split_rt`. A
  single-phase tangent-plane result no longer carries `k_min`, `k_max`,
  `max_delta_k`, `rr_f0`, `rr_f1` or `rr_status`: no K-values are computed on
  that path, and reporting zeros would be fiction.
- **Verdicts change where the heuristic was wrong.** States that used to raise
  `ConvergenceError` now return a single phase; states the heuristic called
  single-phase may now split. Converged two-phase results are unchanged to
  8.6e-7 relative in vapor fraction (they are the same equilibrium reached from
  a different seed, so they differ only within the K-update tolerance).
- **A new failure mode.** `phase_detection="tangent-plane"` raises
  `ConvergenceError` when stability is inconclusive. Pass
  `phase_detection="wilson-heuristic"` to get the old behavior.
- **Single-phase naming is a convention for single-root fluids.** For a dense or
  supercritical feed with one real compressibility root the two branches
  coincide and the reported `"vapor"` / `"liquid"` name is a tie-break, not a
  phase identification. Genuinely subcritical liquids have three roots and are
  named correctly by the min-Gibbs rule.

## Alternatives considered
- **Keep the heuristic.** Rejected: it is measurably wrong (counts above), and
  `stability_tp` already exists precisely to replace it.
- **Always run a full Gibbs minimization.** Rejected: a global minimization is a
  much larger dependency and a much larger claim than this package can honestly
  support; Michelsen's local test with a deterministic trial set is what is
  implemented and what is documented as bounded.
- **Test both converged phases for stability now (phase addition).** Rejected
  for this slice: it changes the number of phases a `FlashResult` may contain,
  which is a separate contract change. Deferred to `flash-phase-addition-lle`.
  `tpd_min` of the *feed* is recorded so the groundwork is visible.
- **Use `phase_branch` for the vapor/liquid orientation** (the original design).
  Rejected on evidence; see decision 3.
- **Use a compressibility factor or molar volume to identify the vapor.**
  Rejected: `EquationOfState` exposes only `fugacity_coefficients`, and widening
  that protocol is a public-API change this slice does not need. It is the right
  fix if labelling ever has to be more than a convention.
- **Bump `cli_schema_version` to 2.** Rejected: no structural break.

## Consequences
- Positive: the single-phase verdict is now a thermodynamic statement with a
  recorded `tpd_min`, and every two-phase answer ships its own verification
  residuals.
- Positive: verdict agreement with `thermo`'s `FlashVL` over a 175-state
  Peng-Robinson grid rises from 166/175 (heuristic) to 175/175.
- Tradeoff: every phi-phi flash now pays for a stability analysis (one trial per
  component plus two Wilson trials) even when the feed is obviously one phase.
- Tradeoff: two more `FlashSettings` fields to keep compatible, and a second
  code path (the legacy one) that must stay alive.
- Tradeoff: the vapor/liquid *name* still leans on the Wilson volatility
  ordering. That is honest for VLE and meaningless for a liquid-liquid split -
  which this slice cannot produce anyway, since it returns at most two phases
  named `liquid`/`vapor`.
- Known limitation: near-critical states where successive substitution is too
  slow still raise `ConvergenceError` (1 state in 1144 on the scanned grid);
  this slice adds no acceleration or second-order stage to the *flash*.

## Next slices
- `flash-phase-addition-lle`: test each converged phase for stability, add or
  remove a phase, and let `FlashResult` carry more than two phases (LLE, then
  VLLE).
- Gamma-phi phase detection once a consistent pure-liquid reference fugacity
  exists (ADR-0007).

## Supersedes (optional)
None. Discharges ADR-0005's "what future slices will change" note about
`flash-auto-phase-detection`.

## Superseded by (optional)
None.
