# ADR-0021: An equation-of-state stability trial runs on one fixed density root

Status: accepted
Date: 2026-09-13

## Context
ADR-0012 pinned each *modified-Raoult* tangent-plane trial to one phase
candidate, because re-selecting the lowest-Gibbs candidate inside the iteration
makes the successive-substitution map

    ln W_i <- d_i - min-Gibbs term_i(w)

discontinuous where the candidates cross, and a vapor-like trial started on the
far side of that crossing gets dragged onto the liquid candidate and collapses
to the trivial solution. It explicitly *declined* to do the same for the two
compressibility roots of an equation of state, and gave a reason:

> A missing compressibility root is the *same* model failing to be evaluable at
> that composition: there is nothing else to iterate on there, and minimum-Gibbs
> root selection at every evaluation is what Michelsen and Mollerup require.

The first half of that is right and the second half does not follow. Whether a
root is "the same model" or "a different model" says nothing about whether the
fixed-point map is continuous, and the map is *not* continuous across the
composition where the two roots exchange Gibbs energy. Validation Case P-9 (iv)
is the failure, measured on the same repository, in the same shape ADR-0012
described:

PC-SAFT (Water 2B / n-Hexane, `k_ij = 0`), 101325 Pa, `T = 335 K` - above the
binary three-phase temperature `T3 = 334.807826336 K` - feed `z_water = 0.7`:

- `flash_tp` returned **two liquids**, `G/RT = -1.1618107137`, against the
  vapor-liquid pair's `-1.1643059308`: the answer was metastable by
  **2.495e-03 RT**;
- the hexane-rich liquid `(0.02273, 0.97727)` was reported **stable**
  (`tpd_min = -3.0e-09`): all four trials reached the trivial solution or the
  partner liquid;
- yet the vapor stationary point at `y = (0.21356, 0.78644)` has
  **`tpd = -6.5237e-03`**.

From the Wilson vapor-like start the *liquid* root is the lower-Gibbs one at
the intermediate compositions, so the substitution used liquid `ln phi` and the
iterate slid onto the liquid root, exactly as the 363 K modified-Raoult feed of
ADR-0012 did. The phase-addition search of ADR-0020 was never at fault: a phase
count is never better than the stability test that produced it.

The reason ADR-0012 could not simply pin both families at once is real, and it
is what this ADR has to answer: the modified-Raoult candidates are two
different models that both exist at every composition, whereas a density root
can genuinely **not exist** at an iterate, and "iterate on the vapor root" then
has no meaning.

## Decision

### 1. A multicomponent EOS trial names the root its start estimates
`_EOSTangentPlane.initial_estimates` returns, for `n` active components, the
same `n + 2` trials it always did, each now pinned:

| label            | surface | initial estimate `W0`       |
| ---------------- | ------- | --------------------------- |
| `wilson-vapor`   | vapor   | `z K^Wilson`                |
| `wilson-liquid`  | liquid  | `z / K^Wilson`              |
| `pure-<name>`    | liquid  | component `<name>` dominant |

The pure-component-dominant estimates are liquid-like by construction and are
Michelsen's recommendation for finding a second *liquid*. Running them on the
vapor root as well was **measured** over the 144-state Peng-Robinson grid, the
188-state PC-SAFT grid of validation Case F-4 and the two water / n-hexane
states above: 334 states, **zero verdict changes** and not one state whose
`tpd_min` fell by more than `2.6e-13` - the same stationary point reached by a
different trial, never a new one. They are therefore not run, and the
measurement is recorded rather than the trials added (ADR-0012 decision 3 made
the same call on the ideal-gas surface, for a different reason: there the
argument was structural, here it is empirical).

### 2. A feed with one active component keeps minimum-Gibbs selection
There is no composition degree of freedom, both Wilson estimates *are* the
feed, and pinning them would test a root the feed is not on. This is ADR-0012
decision 3's exception, unchanged.

### 3. A pinned trial evaluates one root per iterate
`ln_terms_on_surface` asks the model for the named branch only. That is half
the model calls the minimum-Gibbs rule made at the same iterate, and it is why
this slice made the equation-of-state stability test *faster*: the 188-state
PC-SAFT grid of Case F-4 goes from 76.6 s to 44.3 s of `stability_tp` alone,
and the default test suite from 244.2 s to 231.7 s before the new golden path
is added (247.3 s for 750 tests with it, against 244.2 s for 738).

### 4. The fallback: a root that does not exist is not an error
Where a model has a single admissible root, `PengRobinsonEOS` and `PCSAFTEOS`
both answer `fugacity_coefficients(..., phase="vapor")` and `phase="liquid"`
with the **same** numbers - documented behaviour of both, not an accident. So a
pinned trial walks the only surface the model has there, without an exception
and without a wrong phase. Where the named branch *raises*, the lowest-Gibbs
candidate is used instead, which is ADR-0012 decision 5 unchanged.

Both are reported, in `StabilityTrial.surface_fallback` (bool) and the new
`surface_fallback_count` (int), and summed into
`diagnostics["surface_fallback_evaluation_count"]`. They are **counted where
the solver compares the branches**, which for a pinned EOS trial is its
stopping point: `_reported_terms` has to evaluate the lower envelope there
anyway (decision 5), so the comparison is free. It is not counted at every
intermediate iterate, and that is a deliberate cost decision: detecting it
there would mean evaluating both branches at every iterate - doubling the model
calls, and undoing decision 3 - for a number that cannot change any result,
because where the model has one root the pinned surface and the lowest-Gibbs
surface *are* the same surface.

What the un-instrumented frequency is, measured once with a both-branches build
and recorded here rather than paid for at run time: over the 144-state
Peng-Robinson grid the build recorded **3889 single-root evaluations against
4378 trial iterations**, and over the 188-state PC-SAFT grid **14975 against
16898**. (The two are not the same denominator - the second-order stage makes
several model evaluations per iteration - so read it as a ratio, about nine in
ten, not as an exact percentage.) The single-root case is the *common* one, and
the fixed surface therefore changes the iteration only in the remaining tenth -
which is where Case P-9 (iv) lived.

### 5. The reported distance and branch are still the lowest-Gibbs ones
Unchanged from ADR-0012 decision 4. The stationarity residual is measured on
the trial's own root, because that is the equation it solves; the
tangent-plane distance is evaluated with the **minimum-Gibbs** root at the
converged composition, because the distance is measured to the lower envelope;
and `StabilityTrial.phase_branch` / `StabilityResult.phase_branch` still name
that minimum-Gibbs root, relabelled by `EquationOfState.phase_identity`
(ADR-0017) exactly as before. `chemthermo.flash._detect` reads `phase_branch`
and `feed_branch` to seed and pin phases, so **its semantics are untouched by
this ADR**: the incipient phase is still pinned to the branch that is lowest in
Gibbs at the stationary point.

### 6. The evaluator contract grows two small things
`ln_terms_on_surface` now returns a named `_SurfaceTerms(terms, label,
fell_back, min_gibbs)`; `min_gibbs` carries the lower-envelope terms when the
evaluation already produced them, so the reporting step does not repeat a pass
over the model (for `surface is None` that keeps the unpinned families
bit-identical *and* spares the call they never used to make). A second method,
`ln_report_terms(w, surface)`, is what a pinned trial calls at its stopping
point; it is where decision 4's counting happens.

## Alternatives considered
- **Evaluate both roots at every iterate and count every fallback exactly.**
  Implemented and measured first. It is bit-identical to what shipped - the
  iterate walks the same root either way - and it costs what the minimum-Gibbs
  rule cost, which put the default suite at 311 s against a 244 s baseline.
  Rejected: it buys a diagnostic, not a result, and it forgoes decision 3's
  halving.
- **Leave the EOS family unpinned and fix Case P-9 (iv) in the search.**
  Rejected: the search is not wrong. `flash_tp` returned a genuine stationary
  state of the model; it was never told a lower one existed.
- **Add more initial estimates instead.** Rejected for ADR-0012's reason: the
  repaired `wilson-vapor` trial uses the same starting point it always did.
  Only the surface it walks changed.
- **Report the pinned root's own distance as `tpd`.** Rejected: it is not the
  tangent-plane distance, and it would over-report instability wherever the
  pinned root lies above the other.

## Consequences
- Validation Case P-9 (iv)'s pinned miss is **retired**. At `z_water = 0.7` the
  41-point scan across `[T3 - 1 K, T3 + 1 K]` now switches `LLE -> VLE` exactly
  once, with zero `ConvergenceError`, and the boundary bisects to
  `2.68e-07 K` above the independently computed `T3` (`z_water = 0.3`:
  `2.09e-07 K` below it; the shipped example bisects a smaller bracket fewer
  times and reports `4.88e-07 K`, its own resolution). FeOs's chemical
  potentials at the new vapour-liquid phases agree to `9.0e-12`, `3.6e-11` and
  `2.8e-11` with matched universal constants, against an asserted 1e-08.
- **No verdict, tie line or pinned number moved.** Over the 144-state
  Peng-Robinson grid and the 188-state PC-SAFT grid: zero verdict changes, zero
  branch-label changes, worst `|delta tpd_min|` `1.8e-15` and `1.2e-12`, worst
  minimizing-composition move `7.3e-12` and `2.2e-15`. What did change is that
  18 and 38 states respectively now *find* a stationary point where no
  non-trivial one was reachable before - every one of them **above** the
  tangent plane, so every one of them still `"stable"`.
- The bit-identity fixture was regenerated to
  `tests/fixtures/flash/refactor_bit_identity_v3.json` after a state-by-state
  audit: **122 of 155 states bit-identical**, 18 changed only
  `diagnostics["tpd_min"]` from `0.0` to a positive number (the 18 above), 15
  changed only in the last bits (worst composition move `1.1e-16`, worst
  diagnostics move `5.3e-15`). No phase name, phase set, composition or phase
  fraction moved beyond those last bits, and no diagnostics key appeared or
  disappeared. `v1` and `v2` are kept, unreferenced, for history.
- `StabilityResult.diagnostics` gains `trial_surfaces`,
  `surface_fallback_trial_count`, `surface_fallback_evaluation_count` and
  `minimizing_trial_surface` for the EOS family, which previously had none of
  them. `chemthermo.flash` forwards a fixed whitelist of stability diagnostics,
  so none of these reaches a `FlashResult`.
- The honesty note on `StabilityResult` still stands, and ADR-0012's wording of
  it applies here word for word: fixing the surfaces enlarges the set of
  reachable stationary points; it does not turn a local search into a global
  proof.

## Supersedes (optional)
None. Extends ADR-0012 from the modified-Raoult candidate pair to the density
roots of an equation of state, and reverses its "Pin every family to a surface,
cubics included - rejected" alternative on measured evidence. Refines ADR-0005
(minimum-Gibbs root selection, still the rule for *reporting*) and ADR-0007
(evaluator contract).
