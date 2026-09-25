# ADR-0019: Each phi-phi phase on its own density root, and `liquid1` / `liquid2` for an EOS

Status: accepted
Date: 2026-09-13

## Context

`chemthermo.flash._split._solve_k_loop` evaluated the phi-phi split's two
phases on **fixed** density/compressibility branches: `x` always with
`phase="liquid"`, `y` always with `phase="vapor"`, at every iterate, from
ADR-0008 through ADR-0018. The same fixed assignment was written into the
second-order stage (ADR-0016's `_ln_phi_function(..., "liquid")` /
`(..., "vapor")`) and into the two-phase naming rule (ADR-0017's
`_orient_two_phase_labels`, which asked `phase_identity` on those same two
labels).

That makes a **liquid-liquid split from an equation of state inexpressible at
any state where a vapour root still exists**, which is the defect validation
Case P-7(iii) recorded, measured, and deliberately left standing:

- PC-SAFT (2B water, `k_ij = 0`), water / n-hexane, `z = (0.5, 0.5)`,
  298.15 K, 101,325 Pa. The tangent-plane test is **right**: the feed is
  `unstable` with `tpd_min = -9.281926e-01`, `feed_branch = "liquid"` and
  `phase_branch = "liquid"` - both branches are liquids, and the test says so.
- The isotherm nevertheless has two density roots at the feed composition
  (43.01 and 13,419.37 mol/m^3), so the fixed pairing had a vapour root to land
  on. It converged on a water-rich liquid `(0.99992175, 0.00007826)` against a
  hexane-rich **vapour** at 43.3 mol/m^3, with
  `delta_g_split_rt = +0.2594` - a "split" whose Gibbs energy is *above* the
  feed's - and the post-split stability test refused it (`ConvergenceError`,
  "a third phase is required"), with both phases individually unstable
  (`tpd = -1.53` and `-0.78`).
- The true answer is two liquids - water's and n-hexane's vapour pressures sum
  to about 23 kPa, far below 1 atm - and **FeOs's own `State.tp_flash` at the
  same state returns exactly that**.
- At 1 MPa the same tie line *was* reachable, because above about 0.6 MPa the
  isotherm has a single root and `"vapor"` and `"liquid"` name the same root.
  But both converged phases then measured as liquids, ADR-0017's two-phase
  naming rule fell through to its Wilson-ranking fallback, and the pair came
  back named `"liquid"` / `"vapor"` with `vapor_fraction = 0.503164` - which
  was really the hexane-rich **liquid**'s fraction.

Two separate things were therefore wrong: *which root each phase may sit on*,
and *what the two phases are called when both are liquids*. ADR-0018 named the
fix `flash-phase-addition-eos` and said it needed "the phi-phi path to put two
phases on the same density branch" plus "the `liquid1` / `liquid2` vocabulary
ADR-0017's measurement already justifies". This ADR does those two things and
**not** the third-phase search, which stays the next slice.

## Decision

### 1. Each phase sits on the branch the stability test found *that phase* on

`chemthermo.flash._split._PhaseRoot` holds one phase's branch. The
tangent-plane path builds two of them (`_detect._phi_phi_roots`) and hands them
to `_solve_k_loop`, to the second-order stage and to the naming rule, so every
evaluation of that phase - successive substitution, the Newton stage's energy,
gradient and finite-difference Hessian, and the final verification - is on the
same root.

Which branch goes to which phase follows the seed's own orientation.
`_stability_k_seed` sets `K = y / x` with `y` the phase the Wilson ranking
calls vapour-like, so:

| `incipient_phase` | `x` is | `y` is | branch of `x` | branch of `y` |
|---|---|---|---|---|
| `"vapor"` | feed-like | incipient | `feed_branch` | `phase_branch` |
| `"liquid"` | incipient | feed-like | `phase_branch` | `feed_branch` |

`feed_branch` and `phase_branch` are `StabilityResult`'s own labels, already
corrected by ADR-0017's compressibility measurement. For water / n-hexane at
1 atm both are `"liquid"`, which is precisely the information the split was
throwing away. When the documented Rachford-Rice fallback replaced the
stability seed with Wilson K-values (`diagnostics["k_seed"] == "wilson"`) those
roles do not exist and the historical `("liquid", "vapor")` pair is used, so
that path is unchanged.

### 2. The branch is *held* for the split; lowest-Gibbs is the fallback and the check

Michelsen & Mollerup's rule is that each phase sits on the root of lowest Gibbs
energy - the branch minimizing `sum_i w_i ln phi_i(w)`, which is what
`chemthermo.stability._evaluator._ln_phi_min_gibbs` already applies inside the
stability test. Re-selecting that root at **every iterate** of the split was
implemented first and is **rejected** (see "Alternatives considered"): it lets
successive substitution walk a phase onto its partner's branch and collapse to
the trivial solution.

What is shipped instead is ADR-0012's rule for the stability module's own
trials - "a trial pinned to one candidate iterates on that candidate's Gibbs
surface" - applied to the split, for the same reason ADR-0012 gave. The
lowest-Gibbs rule then appears in two places rather than inside the loop:

- as `_PhaseRoot`'s **fallback**, when the pinned branch is not evaluable at a
  composition (`_PhaseRoot.fallbacks` counts it);
- as the **check**, because `_post_split_stability` re-tests every converged
  phase with `stability_tp`, which uses `_ln_phi_min_gibbs`. A split that
  converged onto a higher-Gibbs root is refused there rather than returned.

`delta_g_split_rt` is exact rather than an upper bound as a result: the feed
term was already on the min-Gibbs branch, and both phase terms are now on the
roots the phases actually sit on.

### 3. Naming: `liquid1` / `liquid2` when both phases measure as liquids

`_detect._name_two_phase_result` asks `EquationOfState.phase_identity`
(ADR-0017) for each phase **on the root it converged on**, which since decision
1 need not be `("liquid", "vapor")`:

| measured identities | names (positional with `(x, y)`) | `vapor_fraction` | `phase_regime` | `phase_label_method` |
|---|---|---|---|---|
| one liquid, one vapour | `"liquid"` / `"vapor"`, the vapour-identified phase carrying the fraction | `beta` or `1 - beta` | `"VLE"` | `"compressibility"` |
| both liquid | `"liquid1"` / `"liquid2"` | `None` | `"LLE"` | `"compressibility"` |
| both vapour, or either identity unavailable | `"liquid"` / `"vapor"` in the historical `(x, y)` orientation | `beta` | `"VLE"` | `"wilson-ranking"` |

`vapor_fraction = None` for a liquid-liquid result is the rule ADR-0011 already
set for every phase set with no vapour in it: reporting a number there would be
fiction.

**The `liquid1` / `liquid2` order is composition-based**, unlike the
gamma-gamma path's, whose two names are *roles assigned by the seed* and may
swap between two feeds on one tie line (ADR-0009). `_liquid_order` names
`liquid1` the phase with the larger mole fraction of the **first** component,
ties broken by the second and so on, and finally by position. Two feeds on one
tie line therefore come back with the same labels on the same phases, which is
what makes a lever-rule comparison across feeds meaningful (validation Case
P-8(ii) checks three feeds). The order is relative to the mixture's component
order, so permuting the components permutes which phase is `liquid1`; the phase
*set* is unchanged, and a test asserts exactly that.

The last row keeps ADR-0008 decision 3's Wilson-ranking orientation as the
documented last resort: a model that does not implement `phase_identity` (the
default returns `None`), or a near-critical split whose two phases both measure
`"vapor"`. Ranking two same-side phases by the **magnitude** of `kappa` was
considered and is deliberately not done here - `EquationOfState` exposes the
verdict, not the number, and putting `kappa` itself on the public interface is
a separate decision. No state in any validated grid reaches that row for a
two-phase result; a test reaches it with a `PengRobinsonEOS` subclass whose
`phase_identity` returns `None`, so the fallback is exercised rather than only
asserted to exist.

### 4. `phase_i_branch` / `phase_ii_branch`, conditional keys

A two-phase phi-phi result reports the identity of each converged root as
`diagnostics["phase_i_branch"]` and `["phase_ii_branch"]`, **only when the pair
is not the historical `("liquid", "vapor")`**. Absent means `("liquid",
"vapor")`. This is the same principle ADR-0016 decision 6 used for its four
stage keys: a result the pre-slice code could also have produced carries the
diagnostics mapping it carried before, bit for bit, and a new key appears
exactly where there is something new to say.

What is reported is the ADR-0017 **identity** measured on the converged root,
not the raw candidate label. On a single-real-root state `"liquid"` and
`"vapor"` name the same root, so the raw label there is a tie-break between two
names for one thing; reporting it would make these keys fire on states where no
number moved (measured: Methane / n-Pentane at 360 K / 8 MPa on the ADR-0017
grid). The identity is also the same quantity `feed_branch` / `phase_branch`
report for the stability test, so the three keys mean one thing.

### 5. What is deliberately untouched

The legacy `phase_detection="wilson-heuristic"` path and the deprecated
`gamma-phi` mode both call `_solve_k_loop` **without** the two `_PhaseRoot`
holders and therefore take the unchanged fixed-branch code path. `gamma-gamma`,
`modified-raoult`, the multiphase Rachford-Rice, the phase addition/removal
search, every stability-solver equation, NRTL, the Peng-Robinson and PC-SAFT
physics and `chemthermo.vlle` are not touched. No public name is added:
`phase_i_branch` / `phase_ii_branch` are new keys inside the existing free-form
`diagnostics` mapping, which ADR-0008 decision 7 already settled needs no
`cli_schema_version` bump, and `solver.algorithm` in the CLI output is
unchanged.

## Bit-identity evidence

Two independent checks, because "the outputs did not move" and "the solver did
the same arithmetic" are different claims.

**Outputs.** `tests/test_flash_refactor_bit_identity.py` pins the entire
observable surface of `flash_tp` - phase names, every phase's composition,
phase fractions, `vapor_fraction` and the full `diagnostics` mapping - for 155
states, floats compared with `==`. It passes **unchanged against
`refactor_bit_identity_v2.json`**: no fixture regeneration was needed and none
was done. The audit policy for this slice was that a v3 fixture would be
written only with a per-state ledger entry; no state required one.

**Arithmetic.**
`tests/test_flash_eos_lle.py::test_the_pinned_root_is_the_historical_branch_on_the_whole_phi_phi_grid`
instruments `_PhaseRoot` and replays all 144 Peng-Robinson phi-phi states of
that fixture. On all 47 two-phase states, at **every** iterate of the split,
the fugacity coefficients each phase was evaluated with are exactly (`==`) the
ones `phase="liquid"` for phase I and `phase="vapor"` for phase II would have
produced. The 97 single-phase states build no `_PhaseRoot` at all. So on that
whole grid the branch each phase is pinned to *is* the historical one, and the
new code path is arithmetically the old one.

Two smaller facts found while proving this, both fixed rather than tolerated:

- `_PhaseRoot.fugacity_coefficients` takes the composition the loop already
  normalized and passes it to the model **unchanged**. Normalizing again - a
  vector whose sum is one only to the last bit - perturbs the model's argument
  and moved compositions in the fifteenth digit on one grid state
  (Ethane / n-Heptane, 240 K, 2.0e5 Pa). `ln_fugacity_terms`, the entry point
  the second-order stage uses, does normalize, because the pre-slice callable
  it replaces did.
- The two model calls keep their historical order (`y` first, then `x`) and the
  `K = phi^x / phi^y` update keeps its floating-point spelling. Returning `ln
  phi` from the holder and recombining as `exp(t^x - t^y)` - the shape the
  modified-Raoult path uses - would have been mathematically identical and
  numerically different.

Beyond the fixture: the full suite (712 passing tests before this slice's own
were added), the PC-SAFT Cases P-1 to P-7, the ADR-0016 reference state and its
two siblings, the 16-state F-4 subset, the 188-state slow grid and the teqp
cross-checks all pass with their pinned numbers unchanged.

## Alternatives considered

- **Re-select the lowest-Gibbs root at every iterate** (each phase, every
  evaluation - the literal Michelsen & Mollerup rule). **Implemented first,
  measured, rejected.** On the ADR-0016 reference state - PC-SAFT carbon
  dioxide / n-decane, `z = (0.9, 0.1)`, 240 K, 1.0 MPa - successive
  substitution walks the liquid phase to `x = (0.99066, 0.00934)` at the
  fifteenth iterate, where the *vapour* root has the lower Gibbs energy. Both
  phases then sit on the same root at nearly the same composition, the
  iteration collapses onto the trivial solution, and the split returns
  `beta = -3.247203e+09`, which `_detect`'s guard correctly refuses. Three
  ADR-0016 states, three states of the F-4 grid, the 16-state subset and both
  teqp cross-checks failed with it. The rule is right *at the solution* and
  wrong *as an iteration map*, which is exactly why ADR-0012 pins a stability
  trial to one candidate surface; decision 2 applies the same reasoning, and
  the post-split stability test still enforces the rule at the solution.
- **Keep the branch fixed but seed only the first evaluation**, re-selecting
  afterwards. This was the intermediate the measurement above ruled out too:
  the collapse happens at iterate 15, not iterate 1.
- **Emit `phase_i_branch` / `phase_ii_branch` on every two-phase phi-phi
  result.** Rejected: it changes the pinned `diagnostics` mapping of all 47
  two-phase fixture states without changing a single number, so it would have
  forced a v3 fixture whose entire ledger read "two keys added". The
  conditional form says the same thing and leaves the evidence intact.
- **Report the raw candidate label rather than the measured identity** in those
  two keys. Rejected; see decision 4.
- **A `kappa`-magnitude tie-break for two phases that both measure `"vapor"`.**
  Rejected for now: `EquationOfState.phase_identity` returns a verdict, not the
  number, and adding a `kappa` accessor is a public-interface decision of its
  own (ADR-0001) for a row no validated state reaches.
- **Name the two liquids by role, as `gamma-gamma` does** (`liquid1` = the
  feed-like phase). Rejected: the phi-phi path has a composition-based order
  available and role names would make the three-feed lever-rule comparison of
  Case P-8(ii) meaningless. The two conventions now differ, which is documented
  in `flash_tp`'s docstring rather than smoothed over.
- **Solve the 1 atm state by adding a third phase instead.** That is the
  correct fix for a genuinely three-phase state and it is the next slice; it is
  *not* the fix here, because 298.15 K / 1 atm water / n-hexane is a genuinely
  **two**-phase state that the split could not express.

## Consequences

- Positive: `flash_tp(mixture, temperature_K=298.15, pressure_Pa=101325.0,
  eos=PCSAFTEOS())` on water / n-hexane `z = 0.5/0.5` returns `liquid1` /
  `liquid2` whose compositions, densities and phase amounts match FeOs's own
  `tp_flash` to 1e-8, 1e-6 relative and 1e-8 (validation Case P-8). Equation-of-
  state liquid-liquid equilibrium is representable at all, at any pressure.
- Positive: the 1 MPa tie line keeps its numbers and loses its mislabel - the
  `wilson-ranking` fallback no longer fires there, and `vapor_fraction` is
  `None` instead of a liquid's fraction reported as a vapour's.
- Positive: Peng-Robinson gains the same capability with no model change. At
  298.15 K / 1 atm with `k_ij = 0` the same binary now returns a verified
  liquid-liquid split (`delta_g_split_rt = -1.1217`, post-split stable) where
  it previously raised. What that model *predicts* there - essentially zero
  hexane in the water-rich phase - is reported, not asserted.
- Tradeoff: a phi-phi split now makes one extra `phase_identity` call per phase
  only in the cases that already made them; the loop itself makes exactly the
  same number of `fugacity_coefficients` calls as before, because the pinned
  branch is tried first and succeeds. The measured suite time is unchanged
  within noise.
- Tradeoff: `liquid1` / `liquid2` mean different things on the phi-phi path
  (composition order) and on the gamma-gamma / modified-Raoult paths (seed
  roles). Documented in `flash_tp`; unifying them would change pinned
  gamma-gamma numbers and is not this slice's business.
- Known limitation, unchanged: a phi-phi state that needs **three** phases still
  raises. That is now the *only* remaining phi-phi phase-count gap, and it is
  reachable on the same binary - see below.

## Next slice

`flash-phase-addition-eos`: give the phi-phi path the ADR-0011 phase
addition/removal search that `modified-raoult` already has. **Target state:
water / n-hexane at 1 atm, `z = (0.5, 0.5)`, between 328 K and 335 K.** Measured
with this slice's code: 322, 324, 326 and 328 K return a stable two-liquid
result; 335, 336, 340 and 350 K return a stable vapour-liquid result; 330 K and
334 K raise `ConvergenceError` because the post-split stability test finds a
converged phase unstable - i.e. the three-phase window is bracketed, not
guessed. ADR-0011 supplies the multiphase Rachford-Rice and the search;
ADR-0019 supplies the per-phase root bookkeeping a three-phase EOS set needs
(two liquids and a vapour cannot share one branch assignment either).

## Supersedes (optional)

Amends **ADR-0008 decision 3** further: the Wilson-ranking orientation, already
demoted to a near-critical fallback by ADR-0017, is now also bypassed whenever
both converged phases measure as liquids. It remains the documented last resort
for a model without `phase_identity` and for a both-vapour split.

Amends **ADR-0016**'s statement that the phi-phi second-order stage evaluates
"`ln phi` on the liquid root branch for phase I and on the vapor root branch
for phase II": it now evaluates each phase on that phase's own branch. Every
ADR-0016 number is unchanged (see "Bit-identity evidence").

Amends **ADR-0018**'s "Next slice" note, which expected the two-phases-on-one-
branch fix and the `liquid1` / `liquid2` vocabulary to arrive inside
`flash-phase-addition-eos`. They arrived first, as their own slice; the
phase-addition work is unchanged and still next.

Resolves the limitation recorded in validation Case **P-7(iii)**, which is
amended in `.agents/brain/validation-cases.md` and superseded there by Case
**P-8**.

## Superseded by (optional)

None.
