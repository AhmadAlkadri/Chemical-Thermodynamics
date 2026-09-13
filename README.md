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
tangent-plane stability test, not with Wilson K-value heuristics (ADR-0008):

```
flash_tp -> stability_tp(feed) -> single phase | split seeded from the minimizer
```

- The feed is tested first. A **single-phase** result means the feed was *found
  stable*; `termination_reason` is `"feed_stable_tangent_plane"` and
  `diagnostics["tpd_min"]` records the smallest tangent-plane distance found.
- An **unstable** feed seeds the K-values from the converged stationary point
  (`diagnostics["k_seed"] == "stability"`; the Wilson estimate stays as a
  fallback and is reported as `"wilson"` when it is used), and the existing
  Rachford-Rice / successive-substitution split runs unchanged.
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
- **Both phases are re-tested against both candidates.** A vapor-liquid answer
  whose liquid is inside a miscibility gap, or a liquid-liquid answer that
  should be boiling, raises `ConvergenceError` ("a third phase is required")
  rather than being returned.
- **Near a three-phase state it refuses.** For water / 1-butanol at 1 atm and
  z(butanol) = 0.20 there is a window about **0.135 K wide just below** the
  three-phase temperature T3 = 366.2138 K where the deepest tangent-plane
  minimum is the vapor, the converged vapor-liquid pair is not the equilibrium,
  and `flash_tp` raises. The refusal is correct; the *resolution* below T3 is
  not a third phase but a different pair of two, which needs phase addition
  **and removal**. See validation Case R-3.
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

- If a phase is genuinely unstable, `flash_tp` raises `ConvergenceError` saying
  that the two-phase solution is **not a stable phase set and a third phase is
  required**. It refuses rather than returning an answer it has proved wrong.
  Multiphase flash is the next slice.
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

## EOS extension points

The public repo defines a minimal residual-Helmholtz EOS protocol and registry
hooks in `chemthermo.eos`:

```python
from chemthermo.eos import EOSProtocol, list_eos

print(list_eos())
```

`EOSProtocol` requires:
- `num_components()`
- `residual_helmholtz(temperature_K, volume_m3, composition)`

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## Examples and notebooks

- Scripted demos: `examples/README.md`
- Jupyter notebooks: `notebooks/README.md`

## Optional validation dependencies

Install the reference library used by validation tests:

```bash
pip install -e ".[validation]"
```

Deterministic single-case validation scripts:

```bash
python examples/validation/00_reference_case.py
python examples/validation/06_stability_vs_thermo.py
```

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
