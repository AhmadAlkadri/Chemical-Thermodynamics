# TP flash (`flash_tp`)

Isothermal-isobaric flash: automatic phase-count detection, the flash modes, three-phase equilibrium, binary interaction parameters and the NRTL model. Units are SI throughout (K, Pa, mole fractions).

# Current API usage

SI units are used throughout (temperature in K, pressure in Pa).
```python
from chemthermo import Component, Composition, Mixture, PengRobinsonEOS, flash_tp

component_names = ("Methane", "Ethane", "Propane")
components = tuple(Component.from_database(name) for name in component_names)
z = (0.50, 0.30, 0.20)  # Mole fractions.

composition = Composition(fractions=z, basis="mole", normalize=False)
mixture = Mixture(components=components, composition=composition)
eos = PengRobinsonEOS()

temperature_K = 240.0
pressure_Pa = 3.0e6

result = flash_tp(
    mixture,
    temperature_K=temperature_K,
    pressure_Pa=pressure_Pa,
    eos=eos,
)
print(result.vapor_fraction)
print(result.phases["liquid"].composition.fractions)
print(result.phases["vapor"].composition.fractions)
```

See `examples/basic/flash_tp_peng_robinson_demo.py` for a runnable script that prints a
table-style summary.

## Automatic phase detection

In phi-phi mode `flash_tp` decides one phase versus two with Michelsen's
tangent-plane stability test, not with Wilson K-value heuristics (ADR-0008).
`eos=` accepts any `EquationOfState`, which since ADR-0015 means
`PengRobinsonEOS()` **or** `PCSAFTEOS()` (see
[PC-SAFT](#phase-equilibrium-with-pc-saft)):

```
flash_tp -> stability_tp(feed) -> single phase | split seeded from the minimizer
```

- The feed is tested first. A **single-phase** result means the feed was *found
  stable*; `termination_reason` is `"feed_stable_tangent_plane"` and
  `diagnostics["tpd_min"]` records the smallest tangent-plane distance found.
- An **unstable** feed seeds the K-values from the converged stationary point
  (`diagnostics["k_seed"] == "stability"`; the Wilson estimate stays as a
  fallback and is reported as `"wilson"` when it is used), and the
  Rachford-Rice / successive-substitution split runs from there.
- An **inconclusive** stability analysis raises `ConvergenceError` rather than
  quietly returning one phase.
- Every converged two-phase result reports its own verification residuals:
  `mass_balance_residual`, `fugacity_residual`
  (`max_i |ln(x_i phi_i^L) - ln(y_i phi_i^V)|`) and `delta_g_split_rt`, which
  must be negative for the split to be an improvement on the feed.

```python
from chemthermo import Mixture, PengRobinsonEOS, flash_tp

mixture = Mixture.from_database(("Methane", "n-Pentane"), (0.60, 0.40))
result = flash_tp(mixture, temperature_K=175.0, pressure_Pa=1.778e6, eos=PengRobinsonEOS())

print(result.phase_names())                        # ['liquid', 'vapor']
print(result.vapor_fraction)                       # 0.0792360...
print(result.diagnostics["stability_status"])      # 'unstable'
print(result.diagnostics["tpd_min"])               # -0.04265697...
print(result.diagnostics["delta_g_split_rt"])      # -0.00176661...
```

The legacy behavior stays reachable:

```python
from chemthermo import FlashSettings

legacy = FlashSettings(phase_detection="wilson-heuristic")
```

Runnable demo:

```bash
python examples/basic/flash_tp_auto_phase_demo.py
```

**What this does and does not give you.**

- Converged two-phase results are unchanged: the two paths reach the same
  equilibrium from different seeds and agree to 8.6e-7 relative in vapor
  fraction over the validated grid. Verdicts change only where the heuristic was
  wrong. Over a 175-state Peng-Robinson grid, verdict agreement with `thermo`'s
  `FlashVL` (same `Tc`/`Pc`/`omega`, `kij = 0`) rises from 166/175 to 175/175;
  see validation Cases F-1 and F-2.
- **Two phases at most.** A state that needs a third phase is now *reported* -
  the post-split check below raises rather than returning a two-phase answer the
  package has itself proved wrong - but it is not solved. Multiphase flash is
  the next slice. The two phases may be **two liquids**; see "Liquid-liquid
  equilibrium from an equation of state" below.
- **The converged phases are re-tested for stability** (ADR-0009); see
  "Post-split stability" below.
- **A split that successive substitution cannot finish is finished by a
  second-order stage** (ADR-0016); see "When successive substitution
  oscillates" below.
- `"stable"` means no negative tangent-plane distance was found from the
  deterministic trial set, not a global proof (same bound as `stability_tp`).
- **Liquid-liquid splits have their own mode**, `gamma-gamma`; see
  "Liquid-liquid flash" below. **Low-pressure vapor-liquid equilibrium from an
  activity model** is `modified-raoult`; see below.
- **Gamma-phi is still heuristic.** Its diagnostics say
  `phase_detection == "wilson-heuristic"` and its numbers are unchanged. A
  gamma-phi stability test needs a consistent pure-liquid reference fugacity
  that this package does not yet carry; see ADR-0007 for why building one on the
  current gamma-phi flash would produce a silently wrong tangent plane.
- **Vapor/liquid naming is measured by compressibility, not guessed** (ADR-0017).
  `EquationOfState.phase_identity` computes a dimensionless isothermal
  compressibility ratio (`kappa = P / (rho * dP/drho)`, 1 for an ideal gas and
  well below 1 for a liquid) at the root a phase actually converged on, and
  the model that has one implements it (Peng-Robinson and PC-SAFT both do).
  For a single-phase result this replaces the pre-ADR-0017 min-Gibbs
  tie-break - the case where the cubic has a single real root and both
  branch labels would otherwise coincide is exactly what `kappa` was added to
  settle, e.g. a compressed CO2-rich liquid now reports `"liquid"` with
  `vapor_fraction = 0.0` instead of a tie-broken `"vapor"`. For a two-phase
  result the phase measured as a vapor is named `"vapor"` and the other
  `"liquid"`; when **both** measure as liquids the pair is named
  `"liquid1"` / `"liquid2"` instead (ADR-0019, below); and when both measure as
  vapors (near-critical states, where any label is a convention) or the model
  does not implement `phase_identity`, the historical Wilson volatility ranking
  is kept. `diagnostics["phase_label_method"]` records which rule decided
  (`"compressibility"`, `"wilson-ranking"`, or `"tie-break"` for a model
  without `phase_identity`). None of this ever decides the verdict or the
  compositions - only which already-converged phase (or `1 - vapor_fraction`)
  each name attaches to.

## Liquid-liquid equilibrium from an equation of state

`flash_tp(..., eos=...)` returns **two liquids** when that is what the state
is, at any pressure (ADR-0019):

```python
import chemthermo as ct

mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5], normalize=True)
result = ct.flash_tp(
    mixture, temperature_K=298.15, pressure_Pa=101325.0, eos=ct.PCSAFTEOS()
)

print(result.phase_names())                 # ['liquid1', 'liquid2']
print(result.vapor_fraction)                # None - neither phase is a vapor
print(result.diagnostics["phase_regime"])   # 'LLE'
print(result.phases["liquid1"].composition.fractions)
# (0.999983257..., 1.674279...e-05)     the water-rich liquid
print(result.phases["liquid2"].composition.fractions)
# (0.006312223..., 0.993687776...)      the hexane-rich liquid
```

```bash
python examples/basic/flash_tp_pcsaft_lle_demo.py
```

Until ADR-0019 that call raised. The split evaluated one phase on the model's
`"liquid"` density root and the other on its `"vapor"` root, always - so
wherever a vapor root still exists (water / n-hexane at 1 atm has one) the only
pair it could offer was a vapor-liquid one, its Gibbs energy came out *above*
the feed's, and the post-split stability test refused it. Each phase now sits
on the branch the tangent-plane stability test found **that phase** on, so both
may be liquids.

What to expect from the result:

- **Naming.** Both phases measured as liquids -> `"liquid1"` / `"liquid2"`,
  `vapor_fraction = None` (reporting a vapor fraction for a set with no vapor
  in it would be fiction), `diagnostics["phase_regime"] == "LLE"`. One of each
  -> `"liquid"` / `"vapor"` with a real `vapor_fraction`, exactly as before.
- **`liquid1` is the phase with the larger mole fraction of the first
  component**, ties broken by the next component. That order is deterministic
  and composition-based, so two feeds on one tie line come back with the same
  labels on the same phases - unlike the `gamma-gamma` path's `liquid1` /
  `liquid2`, which are roles assigned by the seed and may swap. Permuting the
  mixture's components permutes which phase is `liquid1`; the phase *set* does
  not change.
- **`diagnostics["phase_i_branch"]` / `["phase_ii_branch"]`** report the
  measured identity of the density root each phase converged on, and are
  present **only** when that pair is not the historical `("liquid", "vapor")` -
  read them with `.get()`.
- Every verification residual is reported as usual: the split must be a Gibbs
  decrease against the feed, must satisfy equal fugacities and the mass
  balance, and every converged phase is re-tested for stability.

Nothing about vapor-liquid results changed: the 155-state bit-identity fixture
(`tests/test_flash_refactor_bit_identity.py`) passes unchanged, and a separate
test replays all 144 Peng-Robinson phi-phi states in it to show that on every
iterate each phase was evaluated on exactly the root the old fixed assignment
would have used. A phi-phi state that needs **three** phases still raises; that
is the next slice (validation Case P-8 brackets one, water / n-hexane at 1 atm
between 328 K and 335 K).

## When successive substitution oscillates

Successive substitution on the K-values is only linearly convergent, and on
some feeds it oscillates instead: the K-values cross 1 back and forth and the
implied vapor fraction leaves `[0, 1]`, or there is momentarily no vapor
fraction at all. Before ADR-0016 that ended the flash with
`ConvergenceError("Rachford-Rice failed to bracket a vapor fraction.")`, for
feeds the tangent-plane test had already *proved* to be two-phase.

Two things now keep it going, and neither changes any answer that converged
before:

- the vapor fraction may leave `[0, 1]` **during** iteration (the "negative
  flash" of Whitson & Michelsen 1989), solved on the window
  `1 / (1 - K_max) < beta < 1 / (1 - K_min)` - precisely the range over which
  every phase mole fraction is non-negative. The in-window solver is called
  first and its answer is returned unchanged when it has one, so an iterate
  that worked before is bit-for-bit unchanged;
- the second-order Gibbs-energy minimization already used by the liquid-liquid
  split finishes the job when successive substitution has spent its whole
  `max_iter` budget.

A result that needed the second stage says so, and **only then** carries the
extra keys `converged_stage`, `ssi_iterations`, `second_order_iterations` and
`negative_flash_steps` - read them with `.get()`. A converged vapor fraction
outside `(0, 1)` is still an error after the stage has run: it contradicts the
stability verdict that started the split. (Since ADR-0024 a first stage that
converges *on* such a vapor fraction - the trivial solution - is handed to the
second-order stage before that error is raised, rather than refused outright.)

Since ADR-0024 there is a **third** stage for splits whose compositions leave
machine range - a polymer/solvent vapour-liquid state, where the equilibrium
vapour holds `exp(-450)` of polymer and the tangent-plane minimizer's K-values
bracket no vapor fraction at all. It is the same Gibbs minimization written in
`u = ln n`, it runs only where the first two cannot, and it reports
`converged_stage = "second-order-log"` plus `log_space_*` keys. See
[Polymers](#polymers-adr-0022) and ADR-0024.

`FlashSettings(second_order=False)` turns the stage off;
`phase_detection="wilson-heuristic"` is the full pre-ADR-0016 behavior and
still fails on these feeds, deliberately.

Measured (validation Case F-4): over a 188-state PC-SAFT grid (CO2 / n-decane
and methane / n-hexane) the `ConvergenceError` count goes **4 to 0**, all 184
states that converged before are bit-identical, and over a 1248-state
Peng-Robinson grid 1247 states are unchanged while the one near-critical state
that used to exhaust the iteration limit now converges.

```bash
python examples/validation/15_flash_split_robustness.py
```

## Liquid-liquid flash (`gamma-gamma`)

Pass `activity_model=` with **no** `eos` to get a liquid-liquid split. Both
phases are liquids described by the same model, so the pure-liquid reference
cancels and the equilibrium condition is equality of activities,
`x_i^I gamma_i^I = x_i^II gamma_i^II`. The mode is inferred from the models you
supply and reported as `flash_mode == "gamma-gamma"`; you can also name it
explicitly.

```python
from chemthermo import Mixture, NRTL, NRTLParameters, flash_tp

# n-butanol / water, NRTL parameters from Tessier, Brennecke & Stadtherr,
# Chem. Eng. Sci. 55 (2000) 1785, Table 1 (pair 2-3).
parameters = NRTLParameters.from_pairs(
    [("n-Butanol", "Water", 0.90047, 3.51307, 0.48, 0.48)]
)
mixture = Mixture.from_database(("n-Butanol", "Water"), (0.10, 0.90))

result = flash_tp(
    mixture,
    temperature_K=298.15,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
)

print(result.phase_names())        # ['liquid1', 'liquid2']
print(result.vapor_fraction)       # None - neither phase is a vapor
print(result.phase_fractions)      # {'liquid1': 0.7647..., 'liquid2': 0.2352...}
print(result.phases["liquid1"].composition.fractions)  # (0.019998..., 0.980001...)
print(result.phases["liquid2"].composition.fractions)  # (0.359999..., 0.640000...)
print(result.diagnostics["equilibrium_residual"])      # 8.4e-14
print(result.diagnostics["delta_g_split_rt"])          # -0.01053...
```

Runnable demos:

```bash
python examples/basic/flash_tp_nrtl_lle_demo.py
python examples/validation/09_lle_tessier2000_tie_lines.py
```

Notes and limits:

- **The number of liquid phases is an output.** The feed is tested with
  `stability_tp(..., activity_model=...)` first; a stable feed returns a single
  phase named `"liquid"`. There is no "assume two liquids" mode.
- **`liquid1` / `liquid2` are roles, not identities, on this path.** `liquid1`
  is the phase the split started from as feed-like, `liquid2` the one started
  from the tangent-plane minimizer, so two feeds on the same tie-line can come
  back with the same two compositions under swapped labels. Compare the phase
  *set*, not `result.phases["liquid1"]`. (The **phi-phi** path's `liquid1` /
  `liquid2` are ordered by composition instead and do not swap; see
  "Liquid-liquid equilibrium from an equation of state" above.)
- `vapor_fraction` is always `None` for a gamma-gamma result; use
  `result.phase_fractions`.
- **Verification.** Every split reports `mass_balance_residual`,
  `equilibrium_residual` (`max_i |ln(x_i^I gamma_i^I) - ln(x_i^II gamma_i^II)|`)
  and `delta_g_split_rt`, which must be negative.
- **Two stages.** Successive substitution converges linearly with a ratio close
  to one near a plait point - 536 to 3922 iterations on the Tessier et al.
  (2000) Problem 1 feeds - so after `FlashSettings.ssi_iterations` (default 50)
  the solver switches to a damped Newton *minimization* of the two-phase Gibbs
  energy, whose gradient is exactly the equal-activity residual.
  `ssi_iterations`, `second_order_iterations` and `converged_stage` in
  `diagnostics` record what happened. Measured over the nine validated feeds of
  Tessier Problems 1 and 2, the final residual is at round-off (worst 1.8e-14).
- **Not vapor-liquid.** Passing both an `eos` and an `activity_model` to
  `gamma-gamma` is a `ModelError`. For vapor-liquid equilibrium from an activity
  model use `flash_mode="modified-raoult"` (below); the older
  `flash_mode="gamma-phi"` is unchanged, still heuristic, and deprecated.
- The scope is the same as `stability_tp`'s: `"stable"` means no negative
  tangent-plane distance was found from the deterministic trial set, not a
  global proof.

## Low-pressure vapor-liquid and liquid-liquid (`modified-raoult`)

`flash_mode="modified-raoult"` gives an activity-coefficient liquid a vapor to
be in equilibrium with, at the one place where the pure-liquid reference
fugacity can be written down honestly: low pressure, where `phi_i^V = 1`,
`phi_i^sat = 1`, the Poynting factor is 1, and `f_i^0 = Psat_i(T)`. The
equilibrium condition is modified Raoult's law,

```
y_i P = x_i gamma_i(x) Psat_i(T)
```

Both phases are put on **one** Gibbs surface by measuring them against the same
reference `ln( f_i / (x_i P) )`:

```
liquid candidate:  ln gamma_i(w) + ln( Psat_i(T) / P )
vapor candidate:   0
```

Michelsen's test is then run on whichever candidate has the **lower Gibbs
energy** at each composition - the same rule that picks the minimum-Gibbs root
of a cubic EOS. So a single call decides *one phase or two*, and if two,
*vapor-liquid or liquid-liquid*, and says which (ADR-0010).

```python
from chemthermo import Mixture, NRTL, NRTLParameters, flash_tp

# 1-propanol / water, NRTL parameters from Tessier, Brennecke & Stadtherr,
# Chem. Eng. Sci. 55 (2000) 1785, Table 1 (pair 1-3). LLE-fitted and
# temperature independent: an illustration of the method, not a correlation.
parameters = NRTLParameters.from_pairs(
    [("1-Propanol", "Water", -0.07149, 2.7425, 0.3, 0.3)]
)
mixture = Mixture.from_database(("1-Propanol", "Water"), (0.50, 0.50))

result = flash_tp(
    mixture,
    temperature_K=361.0,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
    flash_mode="modified-raoult",
)

print(result.phase_names())                     # ['liquid', 'vapor']
print(result.diagnostics["phase_regime"])       # 'VLE'
print(result.diagnostics["incipient_phase"])    # 'vapor'
print(result.vapor_fraction)                    # 0.0869266874...
print(result.phases["vapor"].composition.fractions)  # (0.4423209..., 0.5576790...)
print(result.diagnostics["equilibrium_residual"])    # 1.6e-15
```

The stability test on its own answers the same question without solving a
split:

```python
from chemthermo import stability_tp

result = stability_tp(
    mixture,
    temperature_K=380.0,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
    vapor="ideal",          # add an ideal-gas candidate to the liquid
)
print(result.status)        # 'stable'
print(result.feed_branch)   # 'vapor'  - the feed is a superheated vapor
```

Runnable demos:

```bash
python examples/basic/flash_tp_modified_raoult_demo.py
python examples/validation/10_modified_raoult_water_butanol.py
```

Notes and limits:

- **The mode is never inferred.** An `activity_model` with no `eos` still means
  `gamma-gamma`; you must name `flash_mode="modified-raoult"`. Which vapor model
  applies at a given pressure is your physical judgement, not something the
  package should guess. Passing an `eos` to this mode is a `ModelError`.
- **Phase naming comes from the candidates, not from a heuristic.** VLE returns
  `liquid` / `vapor` with a real `vapor_fraction`; LLE returns
  `liquid1` / `liquid2` with `vapor_fraction = None` (the same role-not-identity
  caveat as `gamma-gamma`); a single phase is named by the feed's candidate.
  There is no volatility-ordering convention here, unlike phi-phi.
- **Every phase is re-tested against both candidates.** A vapor-liquid answer
  whose liquid is inside a miscibility gap, or a liquid-liquid answer that
  should be boiling, is not returned: it is the starting point of the phase
  search below.
- **This is the mode that can return three phases.** See
  [Three phases, and how many there are](#three-phases-and-how-many-there-are).
- **Honest limits of the model itself:** ideal vapor, so low pressure only; no
  Poynting correction; no `phi^sat`; and the temperature must lie inside every
  component's Antoine validity range, which is **enforced** - outside it the
  call raises `InputRangeError` instead of extrapolating a vapor-pressure fit.
  `diagnostics["antoine_valid_Tmin_K"]` and `["antoine_valid_Tmax_K"]` report
  the window. Antoine coefficients come from the packaged databank in the form
  `ln( P^sat / bar ) = A - B / (T/K + C)` (Koretsky 2012).
- `"stable"` still means no negative tangent-plane distance was found from the
  deterministic trial set (two Raoult estimates plus one pure-component estimate
  per component), not a global proof.

## `gamma-phi` is deprecated

`flash_mode="gamma-phi"` still works and its numbers are unchanged, but it is
**deprecated in favour of `"modified-raoult"`** and should not be used for new
work. The reason is physical, not stylistic:

- it sets `K_i = gamma_i phi_i^L / phi_i^V` with `gamma_i` from the activity
  model *and* `phi_i^L` from the equation of state evaluated on the liquid
  mixture, so the liquid's nonideality is counted **twice**; and
- it carries **no pure-liquid reference fugacity at all** - no `Psat_i`, no
  `phi_i^sat`, no Poynting - so its two phases are not on one Gibbs surface.

That is also why it has no stability test, no automatic phase detection and no
post-split check (ADR-0007, ADR-0009). Use `"modified-raoult"` at low pressure
and `"phi-phi"` at high pressure. Removal, if it happens, will get its own ADR;
nothing about `gamma-phi` or the CLI changed in this release.

## Post-split stability

Every two-phase result from the tangent-plane phi-phi path, the liquid-liquid
path and the modified-Raoult path is re-tested: each converged phase is fed back
into `stability_tp` with the same model (and, for modified Raoult, against
**both** phase candidates).

```python
print(result.diagnostics["post_split_status"])            # 'stable'
print(result.diagnostics["post_split_tpd_min"])           # 1.3e-14
print(result.diagnostics["phase_stability_liquid1"])      # 'stable'
print(result.diagnostics["phase_stability_tpd_min_liquid1"])
```

- If a phase is genuinely unstable, the phase set is not the answer. On the
  `modified-raoult` path the incipient phase found there is **added** and the
  set re-solved (see below); on the phi-phi and `gamma-gamma` paths `flash_tp`
  raises `ConvergenceError` saying that the solution is **not a stable phase
  set and a third phase is required**. Either way it never returns an answer it
  has proved wrong.
- `FlashSettings(post_split_stability=False)` returns the result anyway; the
  flag gates the *raise*, not the computation, so the diagnostics show the
  failure either way.
- **Converging onto the partner phase is not an instability.** Two coexisting
  phases share one tangent plane, so a stability test on either finds the other
  with `tpd = 0` up to the split's own convergence tolerance. Such a minimizer
  is reported as `"marginal"`. Measured over the in-repo Peng-Robinson grid, the
  most negative post-split `tpd_min` is -7.0e-09, inside the default
  `tpd_tol = 1e-8`.
- **The check is only as sharp as the split.** With a deliberately loosened
  `FlashSettings(tol=...)` a genuinely two-phase state can be reported as
  needing a third phase, because its phases are then too inaccurate for the
  minimizer to be recognised as the partner.
- **Two paths cannot run it and say so.** `gamma-phi` has no stability test
  (ADR-0007) and the legacy `phase_detection="wilson-heuristic"` path exists to
  reproduce pre-ADR-0008 behavior unchanged; both report
  `diagnostics["post_split_checked"] == False` with a
  `post_split_skipped_reason`.

## Three phases, and how many there are

`flash_tp` discovers the number of equilibrium phases, up to
`FlashSettings.max_phases` (default 3), on the `modified-raoult` path
(ADR-0011) and - since ADR-0020 - on the **equation-of-state** (phi-phi) path
as well. Nothing is told how many phases there are:

```
stability of the feed  ->  a two-phase split  ->  stability of every phase
                       ->  add the phase that was found, re-solve
                       ->  a phase fraction goes to zero or below? remove it
                       ->  every phase stable and every fraction positive: done
```

The inner solve is the multiphase Rachford-Rice written as the constrained
convex minimization of `F(beta) = -sum_i z_i ln(t_i)` (Okuno, Johns &
Sepehrnoori, SPE J 15 (2010) 313) over a feasible region built from
`x_i^j >= 0`, which contains no pole. That region does **not** constrain the
sign of the phase fractions, so a phase that should not be there converges to a
non-positive fraction - the "negative flash" - instead of making the solve
fail, and that sign is the removal signal. Successive substitution hands over
to a Newton minimization of the total Gibbs energy in the non-reference phases'
mole numbers (ADR-0011).

```python
from chemthermo import Mixture, NRTL, NRTLParameters, flash_tp

# 1-propanol / n-butanol / water, NRTL parameters from Tessier, Brennecke &
# Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1. LLE-fitted and
# temperature independent: an illustration of the method, not a correlation.
parameters = NRTLParameters.from_pairs(
    [
        ("1-Propanol", "n-Butanol", -0.61259, 0.7164, 0.3, 0.3),
        ("1-Propanol", "Water", -0.07149, 2.7425, 0.3, 0.3),
        ("n-Butanol", "Water", 0.90047, 3.51307, 0.48, 0.48),
    ]
)
mixture = Mixture.from_database(
    ("1-Propanol", "n-Butanol", "Water"), (0.13418838, 0.08427618, 0.78153544)
)

result = flash_tp(
    mixture,
    temperature_K=364.0,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
    flash_mode="modified-raoult",
)

print(result.phase_names())                          # ['liquid1', 'liquid2', 'vapor']
print(result.diagnostics["phase_regime"])            # 'VLLE'
print(result.diagnostics["phase_set_history"])       # 'L -> LV -> LLV'
print(result.phase_fractions)                        # ~1/3 each
print(result.vapor_fraction)                         # 0.3333333...
print(result.diagnostics["delta_g_vs_two_phase_rt"]) # -3.69e-04  (< 0)
```

```bash
python examples/basic/flash_tp_vlle_demo.py
python examples/validation/11_vlle_water_propanol_butanol.py
```

Notes and limits:

- **`max_phases` caps the search, it does not choose the answer.** The
  one-versus-two decision is Michelsen's stability test, not a setting, so
  `max_phases=1` behaves like `max_phases=2`. `max_phases=2` reproduces the
  pre-ADR-0011 behavior exactly: a state needing a third phase raises.
- **`FlashSettings(post_split_stability=False)` skips the search**, returning
  the converged two-phase answer with the failure in `diagnostics`. The flag
  means "do not police the phase set".
- **`gamma-gamma` does not search.** The activity-only path still stops at two
  phases whatever `max_phases` says, and raises as before: no state in this
  repository needs a third *liquid* there, and an unexercised path is not a
  shipped capability (ADR-0011 "What remains", narrowed by ADR-0020 decision 5).
  The phi-phi path does search; see the next section.
- **Diagnostics keys appear only on a result that entered the search**
  (`phase_set_history`, `phases_added`, `phases_removed`,
  `delta_g_vs_two_phase_rt`, `rachford_rice_iterations`). Every single- and
  two-phase result of the earlier releases is unchanged down to the last bit,
  diagnostics included.
- **Phase names.** `vapor` is the ideal-gas candidate; `liquid`, or
  `liquid1` / `liquid2` / ..., are the liquids, numbered in the order the search
  created them. On this path the liquid numbers are **roles, not identities** -
  compare the phase *set*. (On the phi-phi path below they are ordered by
  composition instead, and are comparable between feeds.) `vapor_fraction` is
  the `vapor` phase's fraction when there is one and `None` when there is not.
- **Each stability trial runs on one fixed phase candidate** (ADR-0012). The
  modified-Raoult tangent plane is the lower envelope of two *different* models
  - an activity-coefficient liquid and an ideal gas - and re-selecting the
  lower one at every iterate makes the successive-substitution map
  discontinuous where the two surfaces cross. A vapor-like trial could then be
  dragged onto the liquid surface and collapse onto the feed, hiding a real
  instability: that is what made one feed inside the 363 K tie-triangle come
  back a single liquid (the pinned miss of validation Case V-2, now fixed). A
  trial is therefore pinned to one candidate, reports the stationarity residual
  on that surface, and still reports its tangent-plane distance with the
  lowest-Gibbs candidate at the converged composition - the distance is to the
  envelope, not to one sheet. `StabilityTrial.surface` says what a trial
  iterated on, `StabilityTrial.phase_branch` where it stopped. **The density
  roots of an equation of state are pinned the same way** (ADR-0021): the
  vapour-like Wilson start iterates on the vapour root, the liquid-like and
  pure-component starts on the liquid root, and where the model has a single
  admissible root - both `phase` labels name it - the trial walks that one and
  records the fact in `StabilityTrial.surface_fallback`. Over a grid of 75-76 feeds
  per temperature at 363 / 364 / 365 K, the phase-count verdict now agrees with
  an independent lowest-Gibbs classifier at **every** feed (validation Case
  V-5, `python examples/validation/12_vlle_verdict_map.py`).
- **The phase count is never better than the stability test that produced it.**
  `"stable"` means no negative tangent-plane distance was found from the
  deterministic trial set. Pinning the surfaces enlarges the set of stationary
  points a trial can reach; it does not turn a local search into a global
  proof.
- **No performance work was done.** A three-phase solve costs a two-phase solve
  plus several stability tests plus roughly 50 successive substitutions, a
  handful of Newton steps and ~100-150 Rachford-Rice Newton iterations.
- **Multiphase equilibrium is in-tree.** `chemthermo.vlle` - a public plugin
  boundary (`get_vlle_engine`, `VLLEEngine`, `VLLEResult`) added before any of
  this existed - is now **deprecated** (ADR-0013): importing it, or calling
  `get_vlle_engine()`, emits a `DeprecationWarning` pointing here, and
  `flash_mode="vlle"` raises `ModelError` with the same pointer. The exported
  names are unchanged and stay importable for one deprecation cycle; nothing
  is removed yet.

## Three phases from an equation of state (ADR-0020)

The same search serves `flash_tp(..., eos=...)`. Each phase is pinned to its
own density root (ADR-0019), so a vapour and two liquids sit on three
independent roots, and a phase added by the search is pinned to the branch the
stability test found *it* on.

```python
import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

# water / ethanol / n-hexane, PC-SAFT with 2B water and ethanol, k_ij = 0.
mixture = ct.Mixture.from_database(
    ("Water", "Ethanol", "n-Hexane"), (0.4, 0.3, 0.3), normalize=True
)
result = ct.flash_tp(
    mixture, temperature_K=333.0, pressure_Pa=101325.0, eos=PCSAFTEOS()
)

print(result.phase_names())                      # ['liquid1', 'liquid2', 'vapor']
print(result.diagnostics["phase_regime"])        # 'VLLE'
print(result.diagnostics["phase_set_history"])   # 'L -> LL -> LLV'
print(result.vapor_fraction)                     # 0.0853317...
print(result.diagnostics["delta_g_vs_two_phase_rt"])  # -2.89e-04  (< 0)
```

```bash
python examples/basic/flash_tp_pcsaft_vlle_demo.py
python examples/validation/18_pcsaft_vlle_water_hexane.py
```

What is different from the `modified-raoult` path:

- **`liquid1` / `liquid2` / `liquid3` are ordered by the first component's mole
  fraction**, not by the order the search created them, so the names mean the
  same thing at every feed inside a tie triangle (the ADR-0019 rule, extended
  to more than two liquids). `vapor` is whichever phase `phase_identity`
  measures as one (ADR-0017), and `diagnostics["phase_label_method"]` records
  that the names were measured rather than conventional.
- **A phase set with no vapour is `phase_regime = "LLE"` with
  `vapor_fraction = None`**, however many liquids it has. Peng-Robinson with
  `k_ij = 0` on the ternary above at 280 K returns three liquids.
- **Removal is what resolves the interesting binary case.** Water / n-hexane at
  1 atm just below its three-phase temperature runs
  `V -> LV -> LLV -> LL`: the deepest tangent-plane minimum from the feed is a
  vapour, the vapour-liquid pair is unstable towards a second liquid, and the
  three-phase solve then drives the vapour amount negative. Those temperatures
  raised `ConvergenceError` before ADR-0020.
- **A binary at fixed pressure has no three-phase region.** Gibbs' phase rule
  leaves one degree of freedom, so three phases meet at a single temperature,
  and the three phase *amounts* are not determined by the mass balance there.
  `flash_tp` returns the two-phase answer on either side of it, and a
  three-phase answer needs a ternary.
- **A three-phase PC-SAFT flash is slow** - about 35 s on the development
  machine, nearly all of it in density-root solves. Peng-Robinson costs about
  0.1 s. No performance work has been done (`perf-profile-baseline` is the
  slice that will do it).
- **Still capped by the stability test.** At `z_water = 0.7` and above the
  three-phase temperature, water / n-hexane comes back as two liquids whose
  Gibbs energy is 2.5e-03 RT *above* the vapour-liquid pair's, because the
  deterministic trial set misses the vapour stationary point from the
  hexane-rich liquid. Measured and pinned as validation Case P-9 (iv).

## Binary interaction parameters (`kij`)

`PengRobinsonEOS.kij` accepts either a scalar (applied to every `i != j` pair,
never to the diagonal) or a mapping from an unordered pair of component names
to a per-pair value. Names are matched case/whitespace-insensitively; a pair
missing from the mapping defaults to `0.0`.

```python
from chemthermo import Mixture, PengRobinsonEOS, flash_tp

mixture = Mixture.from_database(("Methane", "n-Decane"), (0.50, 0.50))
eos = PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.0411})

result = flash_tp(mixture, temperature_K=350.0, pressure_Pa=3.0e6, eos=eos)
print(result.vapor_fraction)
```

Both orders of a pair (`("A", "B")` and `("B", "A")`) refer to the same value;
giving both with *different* values raises `ModelError`, as does a pair
naming the same component twice. The diagonal is always unaffected by `kij`,
so pure-component fugacities never change with it (this was a bug in earlier
versions -- see ADR-0006). Runnable demo:

```bash
python examples/basic/tp_flash_pr_kij_demo.py
```

## NRTL activity coefficients

`NRTL` implements the standard Renon-Prausnitz equation (AIChE J. 14 (1968)
135) with the index convention `tau[i, j] = tau_ij`, `alpha[i, j] = alpha_ij`
and `G_ij = exp(-alpha_ij tau_ij)`:

```text
S_j = sum_k G_kj x_k                      (column sums)
C_j = sum_k tau_kj G_kj x_k
ln gamma_i = C_i / S_i + sum_j x_j G_ij / S_j * (tau_ij - C_j / S_j)
```

**Correctness note.** Versions before the `nrtl-gibbs-duhem-fix` slice summed
`G` along rows instead of columns and divided the first term term-by-term. The
result was not the composition derivative of any excess Gibbs energy: with
asymmetric parameters it violated the Gibbs-Duhem relation
`sum_i x_i d ln gamma_i = 0` (residuals of order 1e-1) and differed from
`thermo.NRTL` by up to 0.89 in `ln gamma`. The current implementation agrees
with `thermo` to ~9e-16 and satisfies Gibbs-Duhem to ~2e-10 (central
differences, step 1e-6). Symmetric binaries were unaffected, so gamma-phi flash
results with the packaged pairs shifted only slightly (vapor fraction
0.767092 -> 0.764835 for the Methane/Ethane CLI example).

**Packaged parameters are synthetic.** The pairs shipped in
`src/chemthermo/parameters/data/activity/nrtl.json` (Methane/Ethane and
Benzene/Water) are illustrative placeholders so that the gamma-phi path,
`--flash-mode gamma-phi` and the examples run out of the box. They are **not
fitted to data and not taken from any publication**; do not use them for
engineering work. Supply your own:

```python
from chemthermo import NRTL, NRTLParameters

parameters = NRTLParameters.from_pairs(
    [("1-Propanol", "Water", -0.07149, 2.7425, 0.3, 0.3)]
)
model = NRTL(parameters=parameters)
```

A citation-backed published parameter set (Tessier, Brennecke & Stadtherr,
Chem. Eng. Sci. 55 (2000) 1785, Table 1) lives in
`tests/fixtures/nrtl/tessier2000_problem1.json`, deliberately outside the
packaged defaults. Reproduce the paper's Table 2 tangent-plane stationary
points with:

```bash
python examples/validation/07_nrtl_tessier_stationary_points.py
```
