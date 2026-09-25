# ADR-0010: Phase candidates in the tangent-plane evaluator, and modified-Raoult `flash_tp`

Status: accepted
Date: 2026-09-13

## Context
ADR-0007 isolated the two model families behind the internal
`_TangentPlaneEvaluator` contract and recorded the one thing it could not do:

> a gamma-phi tangent plane needs a *consistent pure-liquid reference fugacity*
> `f_i^0(T, P)` to put both phases on one Gibbs surface [...] `chemthermo`'s
> current gamma-phi flash does not carry such a reference correctly, so a
> gamma-phi stability test built on it would return a tangent plane that is not
> the physical one - and it would do so silently, which is worse than refusing.

That reference is writable at low pressure and nowhere else cheaply: with
`phi_i^sat = 1`, Poynting = 1 and an ideal-gas vapor, `f_i^0 = Psat_i(T)` and
the equilibrium condition is modified Raoult's law,
`y_i P = x_i gamma_i(x) Psat_i(T)`. The databank already carries Antoine
coefficients for every component and nothing used them.

ADR-0009 left a second gap: `flash_tp` with an activity model could return two
liquids but not a liquid and a vapor, so the whole VLLE campaign had no
vapor-liquid path that a stability test could reach.

The two gaps share one mechanism. Michelsen's test needs, at each composition,
the *lowest* Gibbs energy the mixture can have there. For a cubic EOS that is
the minimum-Gibbs compressibility root (ADR-0005 decision 3). For an activity
liquid against an ideal vapor it is whichever of the two is lower. Those are
the same question asked of a different set of alternatives.

## Decision

### 1. An evaluator holds phase *candidates*, and picks the lowest-Gibbs one
`chemthermo/stability/_evaluator.py` gains an internal `_PhaseCandidate`
protocol - a `label`, an `optional` flag, and
`ln_fugacity_terms(w) -> ndarray` returning that candidate's contribution to the
shared reference `ln( f_i / (x_i P) )` - plus one selection function
`_select_min_gibbs`, which keeps the candidate minimizing `sum_i w_i term_i(w)`.
That sum is the only candidate-dependent part of

    G/RT = sum_i w_i [ g_i^0/RT + ln(w_i P / P^0) ] + sum_i w_i term_i(w)

so it selects the lowest-Gibbs candidate. `_TangentPlaneEvaluator` is unchanged
in shape: `ln_fugacity_terms(w)` still returns `(terms, label)`, where `label`
now names the selected *candidate* rather than specifically a compressibility
root, and is `None` when there is only one candidate and therefore no choice
was made.

The three families are now three candidate sets:

| family | candidates | term |
|---|---|---|
| `"eos"` | the cubic's vapor and liquid roots | `ln phi_i(w)` |
| `"activity"` | one activity-model liquid | `ln gamma_i(w)` |
| `"modified-raoult"` | an activity liquid and an ideal gas | `ln gamma_i(w) + ln(Psat_i/P)` and `0` |

`optional` distinguishes "this candidate may not exist here" (a cubic root,
recorded and skipped) from "this candidate failed" (the activity model, which
re-raises). Answering with the surviving candidate in the second case would
report the wrong *phase*, not a missing branch.

**The EOS and activity-only adapters are behaviourally unchanged**: they are the
two-candidate and one-candidate special cases of the same rule, with the same
floating-point operations in the same order and the same error strings. Proven
by `tests/test_flash_refactor_bit_identity.py` and the pinned canonical
Peng-Robinson stability numbers, both unchanged.

A future PC-SAFT (several density roots) or a user-supplied Gibbs-energy phase
model plugs in as further candidates, with no solver change. That is the second
reason this ADR exists: the abstraction is now *earned twice* rather than once,
which is what ADR-0007 said would be needed.

### 2. `stability_tp(..., activity_model=..., vapor="ideal")`
A new keyword `vapor: Literal["none", "ideal"] = "none"`. `"none"` is the
ADR-0007 liquid-liquid test, unchanged. `"ideal"` adds the ideal-gas candidate
and is valid only with `activity_model` (`ModelError` otherwise; an equation of
state already supplies its own vapor branch). With it:

- `feed_branch` and `phase_branch` report `"liquid"` / `"vapor"`, so the caller
  learns *what* the feed is and *what* the incipient phase is;
- `diagnostics["model_family"] == "modified-raoult"` and
  `pressure_dependent` is True (`Psat/P` depends on P);
- `diagnostics["antoine_valid_Tmin_K"]` / `["antoine_valid_Tmax_K"]` record the
  intersection of the components' Antoine ranges.

Combined gamma-phi stability against an **equation-of-state** vapor stays
unsupported. This decision narrows ADR-0007's refusal to the case where the
reference fugacity is still missing; it does not overturn it.

### 3. `flash_mode="modified-raoult"`
Requires `activity_model`, forbids `eos`, and is **never inferred** - an
activity model on its own still means `"gamma-gamma"`. Which vapor model applies
at a given pressure is the caller's physical judgement, and inferring it would
silently change every existing liquid-liquid call.

The candidate labels decide the regime:

| feed candidate | incipient candidate | result |
|---|---|---|
| liquid | vapor | VLE (bubble side), phases `liquid`/`vapor` |
| vapor | liquid | VLE (dew side), phases `vapor`/`liquid` |
| liquid | liquid | LLE, phases `liquid1`/`liquid2`, `vapor_fraction = None` |

so no volatility-ordering heuristic is needed to name the phases (contrast
ADR-0008 decision 3): the tangent-plane test says which candidate each phase is.
`vapor_fraction` is set for a vapor-liquid result and `None` for a
liquid-liquid one, as in ADR-0009. A stable feed returns one phase named by the
feed candidate.

The split gives each phase the candidate the stability test assigned to it, and
the K update becomes one rule for both regimes:

    K_i = y_i / x_i = exp( t_i^x(x) - t_i^y(y) )

which is `gamma_i Psat_i / P` for a liquid/vapor pair and
`gamma_i^I / gamma_i^II` for two liquids. The ADR-0009 second-order stage is
generalized to two different term functions; its derivation is unaffected
(equation (2) there holds phase by phase, and an additive
composition-independent reference is a constant of the split).

Post-split stability re-tests **both** phases against **both** candidates. An
unstable phase raises `ConvergenceError` ("a third phase is required"), exactly
as in ADR-0009 decision 4.

### 4. Antoine ranges are enforced, not extrapolated
`Psat_i(T)` comes from the databank record in the form it is stored in,
`ln( P^sat / bar ) = A - B / (T/K + C)` (Koretsky 2012; verified by water
returning 1.0131 bar at 373.15 K). A temperature outside a component's
`[Tmin_K, Tmax_K]` raises `InputRangeError`. A vapor-pressure fit extrapolated
past its range is wrong by an unbounded and silent amount, and every number
downstream of it would inherit that silently.

### 5. `flash_mode="gamma-phi"` is DEPRECATED
Not removed, not changed, not warned about at runtime; marked deprecated in the
`flash_tp` docstring, in `FLASH_MODES`, and in the README, in favour of
`"modified-raoult"`. The physical reason: gamma-phi sets
`K_i = gamma_i phi_i^L / phi_i^V` with `gamma` from the activity model *and*
`phi^L` from the equation of state evaluated on the liquid mixture, so the
liquid's nonideality is counted twice; and it carries no pure-liquid reference
fugacity at all (no `Psat_i`, no `phi_i^sat`, no Poynting), so its two phases are
not on one Gibbs surface. That is also why it has no stability test and no
post-split check. Removal, if it happens, gets its own ADR; the CLI is
untouched by this slice.

## Alternatives considered
- **A separate `flash_vle_gamma(...)` entry point.** Rejected for ADR-0009's
  reason: the phase *count* and now the phase *kind* are outputs of this
  package, not inputs. The whole point is that one call returns VLE, LLE or a
  single phase.
- **Inferring `"modified-raoult"` from the models supplied.** Rejected: it would
  silently reinterpret every existing `flash_tp(..., activity_model=...)` call
  as vapor-liquid, and it would assume a pressure regime the package cannot
  check.
- **Making `_PhaseCandidate` public.** Rejected, for ADR-0007's reason and with
  the same trigger recorded: a user-supplied third family, or a second consumer
  needing the abstraction. `chemthermo.flash` does now import
  `modified_raoult_candidates`, but as an internal collaborator inside one
  package boundary, not as a contract offered to users.
- **Extrapolating Antoine with a warning.** Rejected: a warning is not a
  number, and the number would still be wrong.
- **A Poynting factor and `phi^sat` from the EOS.** Deliberately out of scope:
  that is the full gamma-phi model, it needs liquid molar volumes the
  `EquationOfState` protocol does not expose, and shipping a half-corrected
  reference would be exactly the silent wrongness ADR-0007 refused.
- **Retrying the split from the second-deepest stationary point when the
  post-split check fails.** Rejected for this slice: see Consequences.

## Consequences
- Positive: one call now covers low-pressure vapor-liquid *and* liquid-liquid
  behavior, discovered rather than assumed, with the same verification every
  other tangent-plane path carries. Validation Cases R-1 to R-4.
- Positive: the liquid-liquid answer is bit-comparable with the `gamma-gamma`
  path (worst deviation 1.7e-15 over four feeds), because the reference offset
  cancels between two liquids - a property the tests assert rather than assume.
- Positive: the pure-liquid reference fugacity ADR-0007 recorded as missing now
  exists for the regime where it is honest, and the `vapor="ideal"` keyword
  leaves room for `"eos"` later without another signature change.
- Tradeoff: two new public keyword values (`flash_mode="modified-raoult"`,
  `vapor="ideal"`) and two new diagnostics keys.
- Tradeoff: `feed_branch` / `phase_branch` now carry candidate labels for a
  third family; consumers that assumed "EOS root or None" must handle it.
- **Known limitation, measured.** In a narrow window just below the three-phase
  temperature (about 0.135 K wide for water / 1-butanol at z = 0.2, 1 atm) the
  deepest tangent-plane minimum is the *vapor*, the split converges a
  vapor-liquid pair whose liquid lies inside the miscibility gap, and the
  post-split test refuses it with "a third phase is required". The refusal is
  correct - that pair is not the equilibrium - but the *diagnosis* is only half
  right: below T3 the resolution is not a third phase but a different pair of
  two, reached by adding the second liquid and then removing the vapor. Phase
  addition **and removal** is `flash-vlle-phase-addition`. Recorded in
  validation Case R-3 and printed by
  `examples/validation/10_modified_raoult_water_butanol.py`.
- Known limitation: ideal vapor, no Poynting, no `phi^sat`, and a temperature
  range bounded by the Antoine fits. The mode is for low pressure, and says so.
- Known limitation: `"stable"` remains bounded by the deterministic trial set
  (now the two Raoult estimates plus one pure-component estimate per
  component).

## Next slices
- `flash-vlle-phase-addition`: multiphase Rachford-Rice with phase addition and
  removal, so the three-phase state at T3 and the window below it become
  answers instead of errors.
- `vlle`: wire the multiphase split through the `chemthermo.vlle` boundary.

## Supersedes (optional)
None. Discharges ADR-0007's "gamma-phi stability needs a consistent pure-liquid
reference fugacity" note for the low-pressure case only, and ADR-0009's note
that a VLLE state "needs the pure-liquid reference fugacity that ADR-0007
records as missing".

## Superseded by (optional)
None.
