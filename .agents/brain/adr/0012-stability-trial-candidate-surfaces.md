# ADR-0012: A tangent-plane trial runs on one fixed phase-candidate surface

Status: accepted
Date: 2026-09-13

## Context
ADR-0007 introduced the internal tangent-plane evaluator and ADR-0010 gave it
*phase candidates*: at each composition `w` the evaluator returns the terms of
the candidate with the lowest Gibbs energy, `ln phi_i` for a cubic's two
compressibility branches, `ln gamma_i + ln(Psat_i/P)` versus `0` for the
modified-Raoult pair. One rule served both, and the solver never branched on
the model family.

That rule was applied *inside* every trial iteration, and for the
modified-Raoult pair it is wrong. Measured, at T = 363.0 K, P = 101325 Pa, for
the 1-propanol / n-butanol / water feed of validation Case V-2,
`z = (0.15493061, 0.04905728, 0.79601210)`:

- the feed's lowest-Gibbs candidate is the liquid;
- the reduced tangent-plane distance is `+8.07e-04` at the tie-triangle's
  liquid I, `+1.33e-03` at liquid II, and **`-9.92e-03` at the equilibrium
  vapor**, so the feed is provably unstable;
- yet `stability_tp(..., vapor="ideal")` reported **stable**, with every trial
  collapsed onto the trivial solution.

The cause is the switching itself, not the initial estimates. From the
Raoult-vapor start the *liquid* candidate has the lower Gibbs energy at the
intermediate compositions, so the successive-substitution update used the
liquid terms and dragged the iterate onto the liquid surface, where the only
reachable stationary point is the feed. The fixed-point map

    ln W_i <- d_i - min-Gibbs term_i(w)

is discontinuous across the surface where the two candidates cross, and a
discontinuous map is not the map Michelsen's first stage assumes.

A cubic's two roots do **not** have this problem, and the difference is not one
of degree. A missing compressibility root is the *same* model failing to be
evaluable at that composition: there is nothing else to iterate on there, and
minimum-Gibbs root selection at every evaluation is what Michelsen and Mollerup
require and what Cases S-4 and F-1 validate against `thermo`. The
modified-Raoult candidates are two *different* models - an activity-coefficient
liquid and an ideal gas - whose surfaces both exist, and are smooth, at every
composition. Choosing between them at every iterate is therefore not a
robustness measure; it is a different, non-monotone iteration over a
non-smooth envelope, and it can hide the very stationary point the test exists
to find.

## Decision

### 1. An initial estimate may name the candidate surface its trial runs on
`initial_estimates(z, active)` returns `_InitialEstimate(label, w0, surface)`.
When `surface` is a candidate label, the solver evaluates *that* candidate at
every successive-substitution and Newton iterate. When it is `None` the solver
calls `ln_fugacity_terms` exactly as before, so an evaluator that names no
surface is bit-for-bit unchanged; this is the guard
`tests/test_flash_refactor_bit_identity.py` holds.

### 2. Cubic roots and the activity-only liquid name no surface
The equation-of-state evaluator keeps minimum-Gibbs root selection at every
iterate (ADR-0005, validated in Cases S-4 / F-1), and the activity-only
evaluator has a single candidate, so there is nothing to pin. Both implement
`ln_terms_on_surface` for contract completeness only; the solver never calls it
for them.

### 3. The modified-Raoult trial set, one surface per trial
For `n` active components, `n + 2` trials - unchanged in count:

| label            | surface | initial estimate `W0`          |
| ---------------- | ------- | ------------------------------ |
| `raoult-vapor`   | vapor   | `z K^Raoult`                   |
| `raoult-liquid`  | liquid  | `z / K^Raoult`                 |
| `pure-<name>`    | liquid  | component `<name>` dominant    |

The vapor surface needs exactly one trial and its starting point is
irrelevant: the ideal-gas term is identically zero, so the
successive-substitution map on that surface is the *constant* map
`ln W_i <- d_i`, which lands on the surface's unique stationary point in one
step with an exactly zero residual at the next evaluation. Running the
pure-component estimates on the vapor surface as well would add `n` trials that
all return the `raoult-vapor` point; they are deliberately not run, and the
reason is recorded rather than the trials added. The liquid surface has no such
structure and keeps its `n + 1` starts.

A feed with a single active component is the exception: there is no composition
degree of freedom, the feed itself is the only admissible trial, and naming a
surface would pin it to the wrong one for half the feeds (a pure vapor tested
on the liquid surface can never meet the stationarity condition). That
degenerate trial keeps minimum-Gibbs selection.

### 4. The reported distance is still the lowest-Gibbs one
The stationarity residual (equation 5) is measured on the trial's own surface,
because that is the equation the trial solves. The tangent-plane distance
(equation 2) is evaluated with the **minimum-Gibbs** candidate at the converged
composition, because the distance is measured to the lower envelope of the
candidates, not to one sheet of it. A trial therefore records both `surface`
(what it iterated on) and `phase_branch` (the lowest-Gibbs candidate where it
stopped). They normally agree.

### 5. An unavailable optional surface falls back, and says so
If a named candidate is `optional` and not evaluable at an iterate, there is no
surface to walk and the solver uses the lowest-Gibbs candidate there, setting
`StabilityTrial.surface_fallback`. Not reachable through the modified-Raoult
pair, whose candidates are both evaluable everywhere; the path is kept generic
for the future multi-root models ADR-0007 anticipated. A *mandatory* candidate
that fails still raises.

### 6. Diagnostics
`StabilityTrial` gains `surface: str | None` and `surface_fallback: bool`.
`StabilityResult.diagnostics` gains `trial_surfaces` (per-surface trial counts
as a deterministic `"<label>:<count>"` string, in trial order),
`surface_fallback_trial_count`, and `minimizing_trial_surface`. All three are
absent when no trial named a surface, so the diagnostics of the EOS and
activity-only paths are unchanged.

## Alternatives considered
- **More initial estimates, same switching rule.** Rejected: the miss is not a
  reachability problem. The repaired `raoult-vapor` trial uses the *same*
  starting point it always did; only the surface it walks changed. Adding
  starts would have left the discontinuous map in place and made the trial set
  larger for no proof.
- **Pin every family to a surface, cubics included.** Rejected: it would
  contradict ADR-0005 and Michelsen and Mollerup, and it has no defined meaning
  where a root does not exist. The asymmetry is the point of this ADR.
- **Report the trial's own surface distance as `tpd`.** Rejected: it is not the
  tangent-plane distance. It would over-report instability wherever the pinned
  surface lies above the other candidate.
- **Detect the crossing and damp across it.** Rejected as unearned complexity:
  it needs a new tolerance, it does not make the map smooth, and the fixed
  surface already removes the crossing from the iteration entirely.

## Consequences
- The 363 K feed of validation Case V-2 is now found unstable
  (`tpd_min = -1.1680295426e-02`, minimizing trial `raoult-vapor` on the vapor
  surface, incipient candidate `vapor`) and `flash_tp` returns the verified
  three-phase answer. Case V-2's pinned miss is retired.
- Validation Case V-5 (new) maps the verdict over 75-76 feeds per temperature
  at 363 / 364 / 365 K against an independent lowest-Gibbs classifier: **zero
  disagreements** at all three temperatures.
- phi-phi, gamma-gamma and gamma-phi results are unchanged, bit-for-bit.
- The honesty note on `StabilityResult` still stands. Fixing the surfaces
  enlarges the set of reachable stationary points; it does not turn a local
  search into a global proof. One trial of the repaired set - `pure-Water` at
  the Case V-2 feed - still does not converge, because the stationarity
  Jacobian there has an eigenvalue of about `3.3e-09` (near-plait, condition
  number `~8.2e+08`); it is recorded, not accommodated, and the verdict does
  not depend on it.

## Supersedes (optional)
None. Refines ADR-0007 (evaluator contract) and ADR-0010 (phase candidates).
