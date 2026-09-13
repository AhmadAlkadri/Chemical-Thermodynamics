# Chemical-Thermodynamics
Chemical engineering thermodynamics utilities packaged as the `chemthermo`
library (src layout, SI units).

## Local install

macOS/Linux (development install):
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest -q
```

Windows PowerShell (development install):
```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest -q
```

Non-dev install (wheel/non-editable):
```bash
python -m pip install .
```

Tiny import and data check:
```bash
python - <<'PY'
from chemthermo import Component

methane = Component.from_database("Methane")
print(methane.name, methane.tc_k, methane.pc_pa)
PY
```

## Common issues

- Editable install fails with an old pip: run `python -m pip install --upgrade pip` inside the active `.venv`, then retry.
- Build backend errors mentioning setuptools/wheel: run `python -m pip install --upgrade setuptools wheel`, then retry `pip install -e ".[dev]"`.
- `ModuleNotFoundError: No module named 'bibtexparser.bparser'`: your environment resolved `bibtexparser` 2.x; chemthermo's citation loader uses the 1.x API. Reinstall with a `bibtexparser<2` constraint (the pinned `pyproject.toml` requirement already enforces this for new installs).

## Current API usage

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

### Automatic phase detection

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
  the next slice.
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
- **Vapor/liquid naming is a convention in two places.** For a single-phase
  result the name is the minimum-Gibbs compressibility root branch; when the
  cubic has a single real root (dense or supercritical fluids) both branches
  coincide and the name is a tie-break, not a phase identification. For a
  two-phase result the phase named `vapor` is the one enriched, relative to the
  feed, in the component with the largest Wilson K over the one with the
  smallest. `EquationOfState` exposes no molar volume, so no density-based
  identification is available; this decides the *name* only, never the verdict
  or the compositions.

### When successive substitution oscillates

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
outside `(0, 1)` is still an error: it contradicts the stability verdict that
started the split.

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

### Liquid-liquid flash (`gamma-gamma`)

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
- **`liquid1` / `liquid2` are roles, not identities.** `liquid1` is the phase
  the split started from as feed-like, `liquid2` the one started from the
  tangent-plane minimizer. Nothing distinguishes two liquids the way volatility
  distinguishes a vapor from a liquid, so no attempt is made to name them by
  composition: two feeds on the same tie-line can come back with the same two
  compositions under swapped labels. Compare the phase *set*, not
  `result.phases["liquid1"]`.
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

### Low-pressure vapor-liquid and liquid-liquid (`modified-raoult`)

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

### `gamma-phi` is deprecated

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

### Post-split stability

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

### Three phases, and how many there are

`flash_tp(..., flash_mode="modified-raoult")` discovers the number of
equilibrium phases, up to `FlashSettings.max_phases` (default 3). Nothing is
told how many phases there are:

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
- **Only `modified-raoult` searches.** The phi-phi and `gamma-gamma` paths still
  stop at two phases whatever `max_phases` says, and raise as before: no state
  in this repository needs a third phase on either, and an unexercised path is
  not a shipped capability. See ADR-0011 "What remains".
- **Diagnostics keys appear only on a result that entered the search**
  (`phase_set_history`, `phases_added`, `phases_removed`,
  `delta_g_vs_two_phase_rt`, `rachford_rice_iterations`). Every single- and
  two-phase result of the earlier releases is unchanged down to the last bit,
  diagnostics included.
- **Phase names.** `vapor` is the ideal-gas candidate; `liquid`, or
  `liquid1` / `liquid2` / ..., are the liquids, numbered in the order the search
  created them. The liquid numbers are **roles, not identities** - compare the
  phase *set*. `vapor_fraction` is the `vapor` phase's fraction when there is
  one and `None` when there is not.
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
  iterated on, `StabilityTrial.phase_branch` where it stopped. Cubic
  compressibility roots are deliberately *not* pinned: a missing root is the
  same model failing to exist at that composition, so minimum-Gibbs root
  selection at every iterate stays the rule there. Over a grid of 75-76 feeds
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

### Binary interaction parameters (`kij`)

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

### NRTL activity coefficients

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

## Phase stability (tangent-plane analysis)

`stability_tp` answers "is this feed one phase or more?" at fixed T, P and z
using Michelsen's tangent-plane-distance criterion, independently of `flash_tp`.

```python
from chemthermo import Mixture, PengRobinsonEOS, stability_tp

mixture = Mixture.from_database(("Methane", "Ethane", "Propane"), (0.50, 0.30, 0.20))

result = stability_tp(
    mixture,
    temperature_K=240.0,
    pressure_Pa=3.0e6,
    eos=PengRobinsonEOS(),
)

print(result.status)             # "unstable"
print(result.tpd_min)            # -0.3492770207  (dimensionless, units of RT)
print(result.trial_composition)  # incipient-phase mole fractions w
print(result.k_values)           # w_i / z_i for the incipient phase
```

Runnable demo:

```bash
python examples/basic/stability_tp_peng_robinson_demo.py
```

### Liquid-liquid stability with an activity model

Pass `activity_model=` instead of `eos=` to test a liquid feed for a
liquid-liquid split. Both phases are liquids with the same pure-liquid
reference state, so that reference cancels and `ln gamma_i` takes the place of
`ln phi_i` in the same tangent-plane distance - no other change.

```python
from chemthermo import Mixture, NRTL, NRTLParameters, stability_tp

# n-butanol / water, NRTL parameters from Tessier, Brennecke & Stadtherr,
# Chem. Eng. Sci. 55 (2000) 1785, Table 1 (pair 2-3).
parameters = NRTLParameters.from_pairs(
    [("n-Butanol", "Water", 0.90047, 3.51307, 0.48, 0.48)]
)
mixture = Mixture.from_database(("n-Butanol", "Water"), (0.10, 0.90))

result = stability_tp(
    mixture,
    temperature_K=298.15,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
)

print(result.status)             # "unstable"
print(result.tpd_min)            # -0.0299944888
print(result.trial_composition)  # (0.419473, 0.580527) - the incipient phase
```

Runnable demo:

```bash
python examples/basic/stability_tp_nrtl_lle_demo.py
```

### Adding an ideal-gas candidate (`vapor="ideal"`)

`stability_tp(..., activity_model=..., vapor="ideal")` adds a **second phase
candidate** - an ideal gas - to the activity-model liquid, using the databank
Antoine coefficients for the pure-liquid reference fugacity `f_i^0 = Psat_i(T)`
(ADR-0010). At each trial composition the candidate with the lower Gibbs energy
is used, so one test covers vapor-liquid *and* liquid-liquid behavior:

```python
result = stability_tp(
    Mixture.from_database(("1-Propanol", "Water"), (0.50, 0.50)),
    temperature_K=361.0,
    pressure_Pa=101325.0,
    activity_model=NRTL(parameters=parameters),
    vapor="ideal",
)
print(result.status)         # 'unstable'
print(result.feed_branch)    # 'liquid'  - what the feed is
print(result.phase_branch)   # 'vapor'   - what the incipient phase is
```

- `vapor="none"` (the default) is the liquid-liquid test above, unchanged.
- `vapor="ideal"` is only valid with `activity_model`; combining it with an
  `eos` is a `ModelError` (an equation of state supplies its own vapor branch).
- `feed_branch` / `phase_branch` become `"liquid"` / `"vapor"` instead of
  compressibility-root labels, and `pressure_Pa` now *does* affect the result.
- The temperature must lie inside every component's Antoine validity range;
  outside it the call raises `InputRangeError` rather than extrapolating.
  `diagnostics["antoine_valid_Tmin_K"]` / `["antoine_valid_Tmax_K"]` report it.
- Combined gamma-phi stability against an **equation-of-state** vapor is still
  unsupported: that needs a reference fugacity with `phi^sat` and a Poynting
  correction, which this package does not carry. ADR-0010 narrows ADR-0007's
  refusal to the low-pressure case; it does not overturn it.

Reproduce the published tangent-plane global minima of Tessier et al. (2000)
Problems 1 and 2:

```bash
python examples/validation/08_stability_nrtl_tessier2000.py
```

Notes:
- `status` is `"unstable"` when a trial converges to a stationary point with a
  negative tangent-plane distance, `"stable"` when trials converge and none
  does, and `"inconclusive"` when no trial converges.
- **What "stable" means here:** no negative tangent-plane distance was found
  from the deterministic trial set (for an EOS: two Wilson estimates plus one
  pure-component-dominant estimate per component; for an activity model: the
  pure-component-dominant estimates only). This is *not* a global proof of
  stability; Michelsen's test is a local stationary-point search and a
  stationary point that no initial estimate reaches can hide an instability.
- `tpd_min` is dimensionless (units of RT). At a stationary point it equals
  `-ln(sum_i W_i)`, so `sum(W) > 1` is the instability signal.
- Exactly one of `eos` and `activity_model` must be given. **The combined
  gamma-phi case with an equation-of-state vapor is still not supported** and
  raises `ModelError`: it needs a reference fugacity carrying `phi^sat` and a
  Poynting correction that this package does not have, and returning a
  plausible-looking wrong tangent plane would be worse than refusing. See
  ADR-0007. At low pressure use `activity_model=...` with `vapor="ideal"`
  (below), where the reference is `Psat_i(T)` and is exact for the model.
- `pressure_Pa` is required and validated for every family, but it does not
  affect an activity-only result (`vapor="none"`);
  `result.diagnostics["pressure_dependent"]` says which case applies, and
  `["model_family"]` is `"eos"`, `"activity"` or `"modified-raoult"`.
- `feed_branch` and `phase_branch` name the **lowest-Gibbs phase candidate** at
  the feed and at the minimizing trial: the compressibility root for an EOS,
  `"liquid"` / `"vapor"` for `vapor="ideal"`, and `None` for an activity model
  on its own (one candidate, so no selection was made). The selection rule is
  the same for all three: keep the candidate minimizing `sum_i w_i term_i(w)`,
  which is the only candidate-dependent part of `G/RT`.
- Each trial runs successive substitution for `settings.ssi_iterations` (default
  50) and then, if still above `settings.tol`, a damped Newton stage on the
  stationarity condition. Near a plait point successive substitution alone does
  not converge at all; `trial.ssi_iterations`, `trial.second_order_iterations`
  and `trial.converged_stage` record what actually happened.
- `flash_tp` **does** consume this now, twice: in phi-phi and gamma-gamma mode
  its single-phase decision is this test (see "Automatic phase detection" and
  "Liquid-liquid flash" above, ADR-0008 and ADR-0009), and every two-phase
  result is re-tested phase by phase ("Post-split stability"). To get the
  tie-line rather than only the verdict, call `flash_tp` with the same
  `activity_model`. Gamma-phi still uses the K-bound heuristic.

## CLI usage

Run TP flash without writing Python:

```bash
chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json
```

Run gamma-phi mode (NRTL liquid activity + Peng-Robinson EOS):

```bash
chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json
```

Module execution is also supported:

```bash
python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000
```

Notes:
- `--flash-mode` defaults to `phi-phi`; valid choices are `phi-phi` and `gamma-phi`.
- Gamma-phi mode currently uses `NRTL()` for the liquid activity model with the
  packaged **synthetic** pair parameters (see "NRTL activity coefficients").
- NRTL pair coverage is data-dependent; missing pair data returns a runtime validation/model error.

## PC-SAFT (`chemthermo.eos`)

`chemthermo.eos.PCSAFTEOS` implements the **non-associating** PC-SAFT equation
of state of Gross & Sadowski, *Ind. Eng. Chem. Res.* **40** (2001) 1244-1260
(hard-chain plus dispersion), with the packaged pure-component parameters of
that paper's Table 1 (ADR-0014).

```python
from chemthermo.eos import PCSAFTEOS

eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
x = [0.5, 0.5]

eos.residual_helmholtz(temperature_K=300.0, volume_m3=1 / 200.0, composition=x)
# -0.0908711764                     A^res / (R T), volume_m3 is MOLAR volume

eos.compressibility_factor(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# 0.9102033819

eos.pressure_Pa(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# 454071.12

eos.ln_fugacity_coefficients(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# [0.0365173954, -0.2096785687]     natural logs, one per component
```

The state of those four methods is always `(T, molar density, x)` (or
`(T, molar volume, x)` for `residual_helmholtz`). Binary interaction parameters
use the same contract as Peng-Robinson (ADR-0006) - a scalar applied to every
off-diagonal pair, or a mapping keyed by an unordered pair of component names:

```python
PCSAFTEOS(components=("Methane", "n-Decane"), kij={("Methane", "n-Decane"): 0.03})
```

### Phase equilibrium with PC-SAFT

Since ADR-0015 `PCSAFTEOS` also implements the pressure-based
`chemthermo.models.EquationOfState` interface, so it goes straight into the
stability test and the flash - **no solver changed for this**:

```python
import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.30, 0.70])
eos = PCSAFTEOS()                      # components come from the mixture

ct.stability_tp(mixture, temperature_K=300.0, pressure_Pa=3.0e6, eos=eos).status
# 'unstable'

result = ct.flash_tp(mixture, temperature_K=300.0, pressure_Pa=3.0e6, eos=eos)
result.vapor_fraction
# 0.16844005
result.phases["liquid"].composition.fractions
# (0.16086974, 0.83913026)
result.phases["vapor"].composition.fractions
# (0.98686249, 0.01313751)
```

Leave `components` out and the model takes its order from the `Mixture`; pass
it and a mixture whose names differ is rejected rather than silently reordered.

Under the hood, a call that names a pressure and a phase has to solve
`P_model(T, rho, x) = P` first. That root set is public:

```python
eos.density_roots(
    temperature_K=300.0, pressure_Pa=21858.084278856164, composition=[1.0],
    mixture=ct.Mixture.from_database(["n-Hexane"], [1.0]),
)
# (8.868596301913758, 7518.498733715524)     mol/m^3, ascending
```

Only mechanically stable roots (`dP/drho > 0`) are returned, ascending, so the
first is the vapour-like candidate and the last the liquid-like one -
`phase="vapor"` and `phase="liquid"` pick exactly those, which is the
Peng-Robinson convention verbatim. **One root is a normal answer**, not a
failure: at a dense or supercritical state both labels name the same state and
the phase name that comes back is the min-Gibbs tie-break convention. The
spinodal-branch root is found and discarded, never returned.
`PCSAFTEOS.molar_volume(..., phase=)` is the reciprocal of the selected root.

Two roots close enough together to fall inside one step of the scan grid
(`5e-4` in packing fraction) are not resolved, and the state is then reported
with one root fewer; that happens only where the isotherm is nearly tangent to
the target pressure, i.e. near a critical or spinodal state.

### Parameters

Eleven compounds ship with the package: Methane, Ethane, Propane, n-Butane,
n-Pentane, n-Hexane, n-Heptane, n-Octane, n-Decane, Nitrogen and Carbon
dioxide. Their provenance is recorded in
`src/chemthermo/parameters/data/eos/pcsaft.json`: the primary citation is the
paper above, but the paper is paywalled and **was not read directly** when this
file was written - the values were transcribed from two independent secondary
sources that cite it and agree digit for digit (FeOs' `gross2001.json` and
Clapeyron.jl's `PCSAFT_like.csv`), and the 42 universal constants were taken
from teqp's source and Wikipedia's PC-SAFT article, which likewise agree.

Supply your own with `PCSAFTParameters.from_records(...)`:

```python
from chemthermo import PCSAFTParameters
from chemthermo.eos import PCSAFTEOS

parameters = PCSAFTParameters.from_records(
    [{"name": "My fluid", "m": 2.5, "sigma_A": 3.6, "epsilon_k_K": 210.0}]
)
eos = PCSAFTEOS(components=("My fluid",), parameters=parameters)
```

A component with no record raises `PCSAFTParameterError`.

### Limits, stated plainly

- **Non-associating only.** The association term (Gross & Sadowski, *IECR* **41**
  (2002) 5510) and the polar terms are not implemented. Do not use this for
  water, alcohols, acids or amines: nothing in the code stops you, and the
  answer will be wrong.
- **Two phases at most.** The phi-phi flash decides one phase versus two and
  stops there, whatever `FlashSettings(max_phases=...)` says (ADR-0009,
  ADR-0011). Automatic phase *addition* exists only on the `modified-raoult`
  path. A three-phase PC-SAFT system is returned as one or two phases with no
  error.
- **The `(T, rho, x)` methods still choose nothing.** `compressibility_factor`,
  `pressure_Pa` and `ln_fugacity_coefficients` take the density you give them.
  Between the two spinodals `Z` is negative, `pressure_Pa` returns the (real)
  negative pressure, and `ln_fugacity_coefficients` raises `ModelError` instead
  of returning a `nan`. Use `density_roots` (or the `phase=` interface) when
  you want the model to choose.
- **`kij` is yours to justify.** The packaged parameter set is pure components
  only; no PC-SAFT binary-interaction table ships with this package, and any
  `kij` used in the docs or examples (0.03 for methane / n-decane) is
  **illustrative**, not a literature-validated value.
- **No temperature derivative**, so no residual enthalpy or entropy, and no
  phase densities in `FlashResult` (compute them with `density_roots` at the
  converged composition).
- Validated against [teqp](https://github.com/usnistgov/teqp) (NIST, MIT,
  automatic differentiation) to better than 3e-14 in `A^res/RT`, `Z` and
  `ln phi` over fourteen states; against its `pure_VLE_T` saturation solver for
  n-hexane at 300 K and 400 K (densities to 2.2e-16 relative); and, for phase
  equilibrium, against teqp's own traced methane / n-hexane 300 K isotherm -
  seven tie lines to `|dx1| <= 2.0e-9` and `|dy1| <= 3.7e-12`, phase densities
  to 1.3e-9 relative, and teqp's own fugacity coefficients evaluated at
  chemthermo's converged phases giving equal fugacities to 3.8e-9 relative.
  Bubble pressures located by bisecting `stability_tp`'s verdict match teqp's
  `mix_VLE_Tx` to 1.4e-8 relative. See validation Cases P-0 to P-5,
  `examples/validation/13_pcsaft_vs_teqp.py` and
  `examples/validation/14_pcsaft_flash_vs_teqp.py`.

## EOS extension points

The public repo defines a minimal residual-Helmholtz EOS protocol and registry
hooks in `chemthermo.eos`:

```python
from chemthermo.eos import EOSProtocol, list_eos

print(list_eos())      # ['pcsaft']
```

`EOSProtocol` requires:
- `num_components()`
- `residual_helmholtz(temperature_K, volume_m3, composition)` - reduced
  residual Helmholtz energy `A^res/(R T)`, with `volume_m3` the **molar**
  volume in m^3/mol

`get_eos("pcsaft", components=[...], kij=..., parameters=...)` builds the
model above through the registry.

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## Examples and notebooks

- Scripted demos: `examples/README.md`
- Jupyter notebooks: `notebooks/README.md`

## Optional validation dependencies

Install the reference libraries used by validation tests - `thermo` (cubic EOS
and activity models) and `teqp` (PC-SAFT, NIST):

```bash
pip install -e ".[validation]"
```

Deterministic single-case validation scripts:

```bash
python examples/validation/00_reference_case.py
python examples/validation/06_stability_vs_thermo.py
python examples/validation/13_pcsaft_vs_teqp.py
python examples/validation/14_pcsaft_flash_vs_teqp.py
```

Every validation test and script skips cleanly when its optional dependency is
missing. `examples/validation/15_flash_split_robustness.py` (phi-phi split
robustness, validation Case F-4) needs no optional dependency at all.

## Database source of truth

Canonical packaged runtime DB path:

- `src/chemthermo/data/components.json`

Regeneration/check command:

```bash
python tools/build_database.py --check
```

Optional non-runtime mirror path (generated, git-ignored):

- `database/components.mirror.json`

## Project brain

This repository maintains lightweight architectural context and decision history in `.agents/`:

- `.agents/brain/brain.md` - current architecture, invariants, and contributor/agent contract
- `.agents/brain/adr/` - Architecture Decision Records (why key design choices were made)
- `.agents/brain/steering-brief.md` - short summaries of recent changes and next steps

Contributors and AI agents should read `.agents/brain/brain.md` before making structural or API changes.
Internal agent workflow and vendored skills: see `.agents/`.
