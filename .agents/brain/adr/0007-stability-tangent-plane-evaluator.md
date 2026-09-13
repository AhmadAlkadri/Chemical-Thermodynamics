# ADR-0007: Internal tangent-plane evaluator contract, and `stability_tp` with an activity model

Status: accepted
Date: 2026-09-13

## Context
ADR-0005 shipped `stability_tp` for equations of state and deliberately
introduced **no** abstraction: "no generic phase model or phase thermodynamics
protocol is introduced by this slice; the existing `EquationOfState` interface
is sufficient for one model family". It also recorded what would change that:
"`stability-tpd-nrtl`: a narrow phase-thermodynamics contract becomes *earned*
at that point, because a second model family will need it."

That point is here. Liquid-liquid stability with an activity-coefficient model
uses exactly the same mathematics: with two liquid phases at the same
pure-liquid reference state, `mu_i^0` cancels from the tangent-plane distance
and `ln gamma_i(w)` takes the place of `ln phi_i(w)` in

    d_i    = ln z_i + ln f_i(z)
    tpd(w) = sum_i w_i [ ln w_i + ln f_i(w) - d_i ]

with `f = phi` for an EOS and `f = gamma` for an activity model. Michelsen's
`tm(W)`, the stationarity condition `ln W_i + ln f_i(w) - d_i = 0`, and the
relations `tm* = 1 - sum W` and `tpd = -ln(sum W)` all carry over unchanged.

Only two things differ between the families, and both are small:

1. the fugacity term itself (an EOS must additionally pick the minimum-Gibbs
   compressibility root; an activity model has a single branch), and
2. the set of deterministic initial estimates (Wilson K-values are built from
   `Tc`, `Pc` and `omega` and say nothing about a liquid-liquid split).

## Decision

### 1. Public signature
`stability_tp(mixture, *, temperature_K, pressure_Pa, eos=None,
activity_model=None, settings=None)`.

- Exactly one of `eos` and `activity_model` must be given; neither or both
  raises `ModelError`.
- The **combined gamma-phi case** (activity-coefficient liquid tested against an
  equation-of-state vapor) is explicitly out of scope and is the error message
  attached to passing both. Reason, recorded here rather than in a code comment:
  a gamma-phi tangent plane needs a *consistent pure-liquid reference fugacity*
  `f_i^0(T, P)` to put both phases on one Gibbs surface, including a Poynting
  correction and a real pure-liquid saturation fugacity. `chemthermo`'s current
  gamma-phi flash does not carry such a reference correctly, so a gamma-phi
  stability test built on it would return a tangent plane that is not the
  physical one - and it would do so silently, which is worse than refusing.
- `pressure_Pa` stays **required** for API uniformity and is always validated.
  For an activity-only model it does not affect the result;
  `diagnostics["pressure_dependent"]` records `False` for the activity family
  and `True` for the EOS family, and `diagnostics["model_family"]` records
  `"eos"` / `"activity"`.
- `StabilityResult.feed_branch` and `StabilityResult.phase_branch` become
  `str | None` and are `None` for the activity family (no root selection is
  performed). `feed_branch`'s dataclass default changes from `"vapor"` to
  `None`. EOS behavior is unchanged.
- `StabilityTrial` gains `ssi_iterations`, `second_order_iterations` and
  `converged_stage`; `iterations` is now their sum (identical to before
  whenever the second-order stage does not run).

### 2. An INTERNAL evaluator contract
`chemthermo/stability/_evaluator.py` defines a private `Protocol`:

```python
class _TangentPlaneEvaluator(Protocol):
    model_family: str          # "eos" | "activity"
    pressure_dependent: bool

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]: ...
    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[tuple[str, np.ndarray]]: ...
```

with two adapters:

- `_EOSTangentPlane` - returns `ln phi` on the minimum-Gibbs root (the ADR-0005
  selection rule, unchanged) plus the branch label; estimates are the two Wilson
  estimates plus one pure-component-dominant estimate per component.
- `_ActivityTangentPlane` - returns `ln gamma` and `None`; estimates are
  pure-component-dominant only (Michelsen's recommendation for liquid-liquid
  stability), with the feed itself as the single admissible trial when only one
  component is present.

The solver (`tp.py`), trivial-solution detection, the summary and the result
types take an evaluator and never branch on the model family.

**It is internal, not public**, because the honest case for it is "two
implementations exist", not "users need to plug in a third". Nothing in
`_evaluator.py` is exported from `chemthermo` or `chemthermo.stability`; the
only new public surface in this slice is the two keyword arguments on
`stability_tp` and the three new `StabilityTrial` fields.

What would make it public: a third family that a *user* supplies (a
Gibbs-energy phase model, a PC-SAFT-style model that needs its own root
selection, or an electrolyte model), or a second consumer inside the package
that needs the same abstraction (`flash_tp` consuming stability in the
`flash-auto-phase-detection` slice will consume `stability_tp`, not the
evaluator). At that point the contract would need: a stable name, a documented
composition/units convention, a statement of which exceptions an implementation
may raise, and a conformance test other implementations can run. None of that
is worth writing for two in-tree adapters.

Model-specific fast paths stay possible: an adapter may cache, may short-circuit
the two-branch root search when it knows there is one root, or may supply extra
deterministic initial estimates - all behind the same two methods, with no
solver change.

### 3. A second-order stage
`StabilitySettings` gains `second_order=True`, `ssi_iterations=50`,
`second_order_max_iter=100` and `second_order_max_step=4.0`. After
`ssi_iterations` successive substitutions that have not met `tol`, the trial
switches to a damped Newton solve of the stationarity condition in the variables
`u_i = ln W_i` (which keeps `W_i > 0` by construction), with a
central-difference Jacobian, a cap on `|delta u|` and a backtracking line search
on `max_i |g_i|`.

It is needed, not speculative: at the near-plait-point feeds of Tessier et al.
(2000) Table 2 the fixed-point map has a contraction ratio close to one. With
`second_order=False` and a 1000-iteration budget, *no* trial at feed
(0.148, 0.052, 0.80) reaches `tol = 1e-10`; with the second-order stage every
trial converges in 5-8 Newton steps. This is asserted in
`tests/test_stability_activity.py::test_successive_substitution_alone_cannot_solve_the_near_plait_feed`.

The default `ssi_iterations = 50` is chosen so that the stage stays **dormant**
for every validated Peng-Robinson state (the largest iteration count there is
17), which is what keeps the EOS results bit-identical.

## Alternatives considered
- **Two separate public functions** (`stability_tp` and `stability_tp_activity`).
  Rejected: the physics is one criterion, and duplicating the solver would let
  the two copies drift - exactly the failure the tangent-plane identities are
  meant to catch.
- **A public `PhaseThermodynamics` protocol.** Rejected as premature: see above.
  It would freeze a contract designed from two examples.
- **Dispatch on `isinstance(model, ActivityModel)` inside `tp.py`.** Rejected:
  it spreads family knowledge through the solver (root selection, trial sets,
  branch labels, diagnostics) instead of isolating it in one place, and it makes
  the third family a diff across the whole module.
- **Second-order stage always on from iteration 1.** Rejected for this slice: it
  would change every existing Peng-Robinson number in ways unrelated to the
  capability being added, making a regression impossible to attribute. The
  Michelsen/Mollerup recommendation is in any case "a few substitutions first",
  which also gives the Newton stage a better starting point.
- **Adding `scipy.optimize`.** Rejected: a new runtime dependency for ~60 lines
  of damped Newton.

## Consequences
- Positive: liquid-liquid stability with an activity model is a first-class
  capability; the published global tangent-plane minima of Tessier et al. (2000)
  Problems 1 and 2 are reproduced (validation Cases S-6 and S-7).
- Positive: the EOS path is provably untouched - the canonical 240 K / 3 MPa
  ternary is pinned to 1e-12 including per-trial iteration counts
  (`tests/test_stability_tp.py::test_canonical_ternary_numbers_are_pinned_to_the_pre_activity_slice`).
- Positive: near-plait-point feeds are now solvable at all.
- Tradeoff: `eos` becomes optional in the signature, so a caller who passes
  neither model now gets a `ModelError` at runtime rather than a `TypeError`.
  Callers who passed `eos=` by keyword (the only supported form) are unaffected.
- Tradeoff: `feed_branch` / `phase_branch` may now be `None`; consumers that
  assumed `str` must handle it.
- Tradeoff: two more solver settings to keep compatible, and a second code path
  whose numerics have to be validated (they are: Cases S-6, S-7, S-8).
- Known limitation: gamma-phi stability is unsupported, so `stability_tp` cannot
  yet be used to test a *vapor-liquid* split described by an activity model.

## Supersedes (optional)
None. Discharges ADR-0005's "no new abstraction" decision and its "what future
slices will change" note about `stability-tpd-nrtl`.

## Superseded by (optional)
None.
