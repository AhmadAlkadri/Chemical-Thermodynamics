# ADR-0028: The map's 36 refusals are three wrong starting points

Status: accepted
Date: 2026-09-14

## Context

ADR-0027 gave this repository a coverage instrument and then a number:
**2074 of 2110 states converge, 36 refuse, 0 converge and violate an
invariant**, and every one of the 36 is polyethylene / n-pentane at 453 K
(ledger Case R-MAP-1). The roadmap's first item became "fix what the map
ranked first", which was a `_phi_phi_second_order` band of 27 states; the
second and fourth were a dilute-feed class of 7 and a stability stall of 2.

They are three different defects and they have the same shape. Not one of them
is a stage that steps wrongly: in each, a stage is handed a **starting point**
that is not near anything, and then behaves correctly all the way into a wall.
That is why this ADR changes no equation, no tolerance, no step rule, no
acceptance test and no public signature; what it changes is what each stage is
started from when the previous attempt has already failed.

Three candidate framings were ruled out before any of this was written, and
they are recorded because each was the obvious one:

- *"the polymer needs looser tolerances"* - no. Every state below converges to
  an equal-fugacity residual between `2.5e-14` and `3.2e-12` against a `1e-08`
  acceptance. Nothing was relaxed.
- *"successive substitution needs a bigger budget"* - no, in the split. At
  8.7 MPa the K-loop's `max_delta_k` is `2.69e+128` after 100 iterations: it is
  diverging, not converging slowly. (It **is** the answer in the *stability*
  test; see decision 3, where that was measured rather than assumed.)
- *"the second-order stage needs a better Newton"* - no. ADR-0026 already put
  a Gill-Murray modified-Newton direction into the log-space stage, and it was
  measured then to leave the 5.2 MPa band alone. It repairs those states here
  the moment the stage is started somewhere else (decision 1), which is the
  measurement that says the direction was never what was missing.

### Class 1: 27 states where the stage continues from a divergence

`split-non-convergence` / `phi-phi`, 27 states: `Mw = 53 000` at 1 wt%
(3.6-5.4 and 7.5-8.7 MPa) and 15 wt% (3.9-7.5 and 8.1 MPa), and `Mw = 16 400`
at 5 wt%, 7.5 MPa.

Traced at 3.6, 5.1 and 8.7 MPa. Successive substitution runs its full 100
iterations and **diverges**: `max_delta_k` = `1.54e+89`, `5.46e+105`,
`2.69e+128`, growing monotonically with pressure. What it hands over at
3.6 MPa is a pair of phases, one at `6.4e-90` polymer and one at `0.9925`
polymer, with `beta = 0.99999` - a state whose phase-II mole numbers
`beta * y` exceed the feed by four orders of magnitude, i.e. not a split at
all. The linear second-order stage pulls that back into its box by shrinking
`beta` to `6.9e-06` (correctly: the box is `0 < n_i < z_i`) and then minimizes
from there; the log-space retry that follows it is seeded from the *same*
iterate through `seed_from_iterate`. Both fail, and the state raises with a
residual of `3.897e-05`.

The tangent-plane stationary point at the very same state is
`w = (5.61e-09, 1.0)` with `tpd_min = -3.889e-05` - a phase the stability test
converged on, already in `_flash_tp_tangent_plane`'s hand, and never passed to
either stage. Started from it, the ADR-0024 stage with the ADR-0026 safeguard
reaches `4.55e-13` in **11** iterations.

The 15 wt% half of the band needs less than that: from the stationary point the
*unsafeguarded* stage converges in 5 to 8 iterations. So the band is two
sub-classes of one defect, and the seed alone settles fourteen of the twenty-six
`Mw = 53 000` states.

### Class 2: 7 dilute states that converge on the trivial solution

`beta-outside-window` (4 states, `Mw = 16 400`, 1 wt%, 0.3-1.2 MPa) and
`log-space` (3 states, `Mw = 53 000`, 1 wt%, 0.3-0.9 MPa). Both reach
`_phi_phi_log_space` through the ADR-0024 stationary-point route.

At 0.3 MPa on the `Mw = 16 400` chain the stage **converges**, to a residual of
`8.283e-13`, on `x^I = x^II = z` with `beta` an exact `0.0`. That is the
*trivial* solution, and the stage's acceptance rule is why it is an attractor:
a step is kept when it lowers the two-phase Gibbs energy **or** when it lowers
the residual, and the residual at the trivial solution is identically zero.
`_flash_tp_tangent_plane` then refuses the result, correctly, because a single
phase contradicts a `tpd_min` of `-455.03`.

What puts it on that side is the seed's phase fraction. `log_space_seed` builds
`n_i = beta K_i z_i / ((1 - beta) + beta K_i)` at a fixed `beta = 0.5`, chosen
in ADR-0024 as the neutral value because the stage is a descent method from
wherever it starts. The equilibrium here is nowhere near neutral: the melt
holds `4.074e-04` of a 1 wt% feed. Measured, at 0.3 MPa: `beta = 0.5` gives the
trivial solution in 17 iterations with or without the curvature safeguard;
`beta = 0.99` and above converge in 6 to 8 iterations to
`beta = 0.9995926`, i.e. a melt fraction of **`4.074163e-04`** - the number
Case R-MAP-1 had already derived from the lever rule on the 5 wt% tie line
(`4.074164e-04`), without running anything.

Class 3 in the map's ranking (`log-space`, the `Mw = 53 000` chain) is the same
seed placement expressed through the other exit: the stage does not converge on
the trivial solution there, it *crawls* next to it - `beta` reaching `6.1e-02`,
`1.6e-02` and `0.0` at the three pressures, residual `8.0e-02` to `1.7e-06` -
and the state raises instead of being refused. Same feed, same shape, same fix.

### Class 4: 2 states where the Newton stage has no admissible step

`stability-inconclusive`, `Mw = 53 000`, 15 wt%, 10.5 and 10.8 MPa. All four
trials end `second_order_no_progress` at stationarity residuals of `0.685` and
`1.746`, so `stability_tp` returns `inconclusive` and `flash_tp` can neither
split the feed nor call it one phase.

This is a genuinely marginal state and the neighbours say so: 9.9 and 10.2 MPa
converge to non-trivial stationary points at `tpd = +1.013e-04` and
`+2.290e-04`, and 11.1 and 11.4 MPa converge to the **trivial** one. Between
them the non-trivial stationary point is merging with the trivial one and the
stationarity system is ill-conditioned along the merging direction.

Four things were measured before choosing, since ADR-0025's clamp had already
been ruled out in Case R-MAP-1 (`|ln W| <= 11.09`):

| what was changed | 10.5 MPa | 10.8 MPa |
| --- | --- | --- |
| nothing (shipped) | inconclusive, residual 0.685 | inconclusive, residual 1.746 |
| `second_order_max_iter` 100 -> 1000 | unchanged, still 65 iterations | unchanged |
| `second_order_max_step` 4.0 -> 0.5 | unchanged | unchanged |
| substitution only, 300 / 3000 iterations | residual **53.4**, `w -> (7e-27, 1)` | residual **36.2** |
| `ssi_iterations` 50 -> 300, then the same Newton stage | **stable**, 3 of 4 trials converge in 307 | **stable**, same |

So the Newton stage is not broken and successive substitution is not simply
slow - on its own it walks away. The **handover** at iteration 50 is early: the
substitution iterate is still travelling, and the Newton stage started from
where it is at 50 has no admissible step, while the same stage started from
where it is at 300 converges in seven.

## Decision

### 1. A seed ladder from the stationary point, walked only after a refusal

`chemthermo.flash._log_space.stability_seed_ladder` returns three ordered
`(seed, curvature_safeguard)` attempts:

1. the ADR-0024 seed at the neutral phase fraction, unsafeguarded;
2. the same seed with the ADR-0026 curvature safeguard;
3. the same seed at `lever_rule_phase_fraction`, safeguarded.

`chemthermo.flash._detect._walk_stability_seed_ladder` runs them in order and
keeps the **first** that returns a residual at or under `FlashSettings.tol`
*and* a phase fraction strictly inside `(0, 1)`. Both halves are load-bearing:
the trivial solution satisfies the first exactly.

Two callers, two gates, both of which end the flash today:

- `_phi_phi_second_order` walks the whole ladder when its two stages have left
  a residual above `tol` - the line that raises class 1's message.
- `_phi_phi_log_space` runs entries 1 and 2 itself already (they are literally
  the two `log_space_split` calls ADR-0024 and ADR-0026 put there), so it walks
  the ladder from entry 3, when those two left a residual above `tol` **or** a
  phase fraction outside `(0, 1)`. The second disjunct is class 2; without it
  the trivial solution passes for converged and the caller refuses it one
  frame up.

`_phi_phi_second_order`'s gate is deliberately **not** widened with the same
disjunct. No state measured anywhere reaches that line with a converged
residual and an unusable `beta`, and shipping a branch nothing exercises is
what ADR-0002 forbids; the case is carried by the other entry point, which
does exercise it.

If no ladder entry converges, nothing is adopted: the caller is left with
exactly the state it had and raises the message it raised before, so a state
that still refuses refuses identically.

### 2. `lever_rule_phase_fraction`, the seed's phase fraction from a bound

If the stationary phase `w` holds a fraction `lambda` of one mole of feed, the
complement is `(z - lambda w) / (1 - lambda)` and has to be non-negative, so

    lambda <= min_i z_i / w_i

and the bound is attained exactly when the complement is empty of the component
attaining it. At a polymer melt stationary point - `w` an essentially pure
polymer, the solvent phase holding `exp(-450)` of it - the bound is tight to a
factor of a few, and being a bound it can never put the seed on the wrong side
of the trivial solution. `w` is phase II when the stationary point is the
vapour-like one and phase I otherwise, the same convention `log_space_seed`
already uses, so the bound is returned as `lambda` or as `1 - lambda`
accordingly, clamped a hair inside `(0, 1)`.

`log_space_seed` gains one keyword, `phase_fraction`, whose default is the
`_SEED_PHASE_FRACTION` it has always used. Every seed formed before this ADR is
formed by the same expression.

### 3. One stability retry with the full substitution budget

`stability_tp` re-runs the whole analysis, from the same trial compositions,
with `ssi_iterations` raised to `max_iter` - the budget successive substitution
already has when the second-order stage is switched off - **and only when the
first pass returned `inconclusive`**. If the retry is also inconclusive the
first pass's result is returned unchanged, so the refusal message does not
move either. It is reported as
`diagnostics["substitution_budget_retry"]` with
`diagnostics["substitution_budget_retry_from"]`.

Nothing inside `_run_trial` or `_second_order_stage` changes: not the
finite-difference Jacobian, not the step cap, not the line search, not a
tolerance. What changes is which iterate the Newton stage is handed, on a
second attempt, for a state that had no verdict at all.

### 4. Dormancy is by construction, and the record says so

Every gate above is on the *failure* path of a state that ended in a raised
`ConvergenceError` or an `inconclusive` verdict. No converged answer can reach
any of them. The claim is checked rather than only argued: the full 2110-state
sweep is re-run and every previously converged state is identical on every
field the record carries, and `refactor_bit_identity_v3.json` and all nine
ADR-0023 benchmark result hashes are unchanged.

## Consequences

- **The map refuses nothing.** All 36 states converge or return a verdict; the
  `polymer` family joins the other five. See the regenerated record and ledger
  Case P-17.
- `tests/test_robustness_map.py`'s `PINNED_REFUSALS` is **empty**, which is now
  the assertion a regression has to break.
- Three new diagnostics values, all conditional on the route that produced
  them: `log_space_seed` gains `"stability-w-lever-rule"`, and
  `stability_tp`'s diagnostics gain `substitution_budget_retry` and
  `substitution_budget_retry_from`.
- **Cost.** A refused state used to be cheap. A state that now converges
  through the ladder pays for the attempts before the one that worked (up to
  three log-space stages), and a stability retry pays for a second full set of
  trials. Both are bounded by the existing budgets and neither is on a path a
  converging state takes, so the sweep's converged states cost what they cost.
- **Not claimed.** That these are the right answers because the solver says so:
  every one of the 34 new splits is checked against a solver sharing no code
  with `flash_tp`, and the two stability verdicts against a trial-free scan of
  equation (2). And not that the ladder is complete - it is three entries
  chosen because three defects were measured, and a fourth geometry may well
  need a fourth.

## Alternatives considered

- **Route to the log-space stage whenever any mole fraction in the iterate is
  below ~1e-12 instead of `TRACE_MOLE_FRACTION = 1e-30`.** This was the first
  candidate and it is the wrong lever: the class 1 states are already *in* the
  log-space stage (through `_phi_phi_second_order`'s ADR-0024 retry) and fail
  there, because of the seed and not the parametrization. Lowering the trace
  threshold would also re-route states that converge today, which is precisely
  what ADR-0024 declined to do.
- **Change `_SEED_PHASE_FRACTION` from 0.5 to the lever rule everywhere.** It
  repairs class 2 and it moves every state ADR-0024, ADR-0025 and ADR-0026
  converged, since they all start from that expression. Measured on the band it
  is also *worse* on 12 of the 27 class-1 states. Rejected on bit-identity, and
  it is the same rejection ADR-0026 recorded for the same constant.
- **A modified-Newton stage in the linear (non-log) second-order stage.** The
  linear stage cannot represent a phase at `exp(-300)` at all - ADR-0024's whole
  premise - so a better direction in it would still be minimizing over a box
  that excludes the answer.
- **Raise `StabilitySettings.ssi_iterations` from 50 to 300 by default.** It
  repairs class 4 and it changes the handover point for every trial in the
  repository, which moves converged stationary points and therefore converged
  splits. Rejected on bit-identity; the retry gets the same iterate on the
  states that need it and nobody else.
- **A trust region, or letting the Newton stage continue as substitution after
  a stalled step.** Both are changes to `_second_order_stage` itself, i.e. to
  the stability solver's mathematics, for a defect measured to be in the
  handover. The retry is the smaller change and it is the one the measurement
  points at.
- **Accept `beta` outside `(0, 1)` when the residual converges.** Considered
  for class 2 and rejected outright: `beta = 0` there is the trivial solution,
  not a small phase. Returning it would be returning "one phase" for a feed the
  tangent-plane test proves is two.

## References

- ADR-0024 (the log-space split stage and its seed), ADR-0025 (`ln W` summed in
  logs), ADR-0026 (the curvature safeguard and the convex denominators),
  ADR-0027 (the robustness map), ADR-0002 (no unexercised paths).
- Validation Case P-17 in `.agents/brain/validation-cases.md`, with an
  amendment to Case R-MAP-1.
- Michelsen, *Fluid Phase Equilibria* **9** (1982) 1 (the stationary-point
  search and the two-stage recommendation); Gill, Murray & Wright, *Practical
  Optimization*, sec. 4.4.2 (the modified-Newton direction ADR-0026 uses).
- `benchmarks/robustness_9adf390.json` / `.md`, `benchmarks/README.md`,
  `tests/test_robustness_map.py`, `tests/test_pcsaft_polymer.py`,
  `tests/test_flash_log_space_stage.py`,
  `tests/validation/test_pcsaft_polymer_vs_feos.py`.
