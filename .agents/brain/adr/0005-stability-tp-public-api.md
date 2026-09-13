# ADR-0005: Public `stability_tp` tangent-plane phase stability API

Status: accepted
Date: 2026-09-13

## Context
`flash_tp` decides "single phase" from Wilson K-bounds and a Rachford-Rice root
check. Those are heuristics on an *initial estimate*, not a thermodynamic
criterion: they can call a genuinely two-phase feed single-phase (and the
reverse) near critical or dense states, and they say nothing about
liquid-liquid splits.

The phase-equilibrium campaign (activity-model stability, automatic 1-vs-2 phase
detection inside `flash_tp`, LLE, VLLE) needs a real stability criterion as a
first-class, independently testable capability rather than as a private helper
inside the flash solver.

## Decision
Add a public `chemthermo.stability_tp` implementing Michelsen's tangent-plane
stability analysis (Michelsen, Fluid Phase Equilibria 9 (1982) 1-19; Michelsen
and Mollerup, "Thermodynamic Models: Fundamentals and Computational Aspects").

Public exports added to `src/chemthermo/__init__.py` `__all__`:
`stability_tp`, `StabilityResult`, `StabilitySettings`, `StabilityTrial`.
Implementation lives in the internal `chemthermo.stability` package
(`tp.py`, `results.py`, `settings.py`); only the four names above are public.

Specific decisions:

1. **Result contract.** `StabilityResult` carries `stable`, a three-valued
   `status` in `{"stable", "unstable", "inconclusive"}`, `tpd_min` (reduced
   tangent-plane distance, dimensionless in units of RT), the minimizing trial
   composition `w`, the implied incipient K-values `w_i / z_i`, the selected
   min-Gibbs root branch, a per-trial record tuple, and a diagnostics mapping.
   `"stable"` is documented as "no negative tangent-plane distance was found
   from this deterministic trial set", never as a global proof.

2. **Status semantics.** `"unstable"` requires a converged, non-trivial
   stationary point with `tpd < settings.tpd_tol` in magnitude and negative in
   sign. `"inconclusive"` is reserved for the case where *no* trial converged at
   all. Convergence to the trivial solution `w -> z` is a successful trial
   outcome (it is a genuine stationary point with `tpd = 0`), so an all-trivial
   outcome is reported as `"stable"`, not `"inconclusive"`. Reporting a clearly
   single-phase state as `"inconclusive"` would be the misleading answer.

3. **Minimum-Gibbs root selection is done generically over the existing
   fugacity interface.** The module calls
   `EquationOfState.fugacity_coefficients(..., phase="vapor")` and
   `phase="liquid"` and keeps the branch minimizing `sum_i w_i ln phi_i(w)`,
   which is the reduced residual Gibbs energy and therefore selects the
   lowest-Gibbs compressibility root at fixed `(T, P, w)`. No new model hook and
   no model-specific code is introduced. Cubic EOS implementations that return
   the same values on both labels (single real root) are handled for free.

4. **Deterministic, small trial set.** Two Wilson-based estimates (vapor-like
   `W = K_wilson * z`, liquid-like `W = z / K_wilson`) plus one
   pure-component-dominant estimate per component. Every trial's outcome is
   recorded in the result. No randomized or adaptive restarts.

5. **No new abstraction.** No generic "phase model" or "phase thermodynamics"
   protocol is introduced by this slice; the existing `EquationOfState`
   interface is sufficient for one model family.

## Alternatives considered
- **Global optimization of tpd** (tunneling, interval, or stochastic methods).
  Rejected: adds a solver dependency and large complexity for a guarantee this
  slice does not need, and would be dishonest to claim without a proof.
- **Require the caller to specify which phase to test.** Rejected: it pushes the
  root-selection question onto users, who would then get physically wrong
  tangent planes whenever they guessed the branch wrong.
- **Fold stability directly into `flash_tp`.** Rejected for this slice: it would
  change `flash_tp`'s behavior, diagnostics and CLI output in the same change
  that introduces the new numerics, making regressions hard to attribute.
  `flash_tp` will consume `stability_tp` in a later, separate slice.
- **Flat `src/chemthermo/stability.py` module.** Rejected for symmetry with
  `chemthermo.flash` (settings / results / solver split), which keeps the
  frozen-dataclass contracts readable.

## Consequences
- Positive: phase stability becomes testable on its own, with invariants
  (tangent-plane identity at the feed, stationarity, marginal stability of
  converged equilibrium phases) that do not depend on the flash solver.
- Positive: min-Gibbs root selection is now implemented once, correctly, and can
  be reused by the flash solver later.
- Tradeoff: four more public names to keep compatible.
- Tradeoff: `"stable"` is a bounded claim, so callers must read the docstring.
- Known limitation inherited from this slice: `PengRobinsonEOS.kij` is a scalar
  that is (incorrectly) applied to the diagonal as well, so only `kij = 0.0` is
  used and validated here.

## What future slices will change
- `pr-kij-matrix`: fix the diagonal-kij bug and accept a per-pair kij matrix.
  Stability results for non-zero kij only become trustworthy after that.
- `stability-tpd-nrtl`: liquid-liquid stability driven by an activity model. A
  narrow phase-thermodynamics contract (something that returns `ln phi` or
  `ln gamma`-derived tangent-plane intercepts) becomes *earned* at that point,
  because a second model family will need it. It is deliberately not introduced
  now.
- `flash-auto-phase-detection`: `flash_tp` consumes `stability_tp` in place of
  the K-bound heuristic, and the converged stationary point seeds the K-value
  iteration.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
