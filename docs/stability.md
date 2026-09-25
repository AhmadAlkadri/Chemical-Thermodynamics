# Phase stability (`stability_tp`)

Michelsen tangent-plane stability analysis. **"Stable" always means: no negative tangent-plane distance was found from a deterministic trial set** - it is not a global proof.

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

Since ADR-0025 the iteration carries Michelsen's unnormalized mole numbers
without clamping them, so a stationary point outside the exponential's range is
found rather than missed. Where that happens - a polymer melt against a
solvent-vapour feed needs `ln W_polymer = 1452` - the normalized
`trial_composition` rounds to an exact `0.0` in some component and
`result.trial_ln_W` is what carries the magnitude; `diagnostics["ln_sum_W"]`
is `-tpd_min` (equation (7)) and `sum_W` itself is `inf`. Nothing else changed:
the diagnostics key `log_space_trial_count` appears **only** when the route
engaged, and it does not on any Peng-Robinson, PC-SAFT or activity-model state
in this repository. See [Polymers](#polymers-adr-0022) and ADR-0025.

Runnable demo:

```bash
python examples/basic/stability_tp_peng_robinson_demo.py
```

## Liquid-liquid stability with an activity model

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

## Adding an ideal-gas candidate (`vapor="ideal"`)

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
  "Liquid-liquid flash" in [flash.md](flash.md), ADR-0008 and ADR-0009), and
  every two-phase result is re-tested phase by phase ("Post-split stability"). To get the
  tie-line rather than only the verdict, call `flash_tp` with the same
  `activity_model`. Gamma-phi still uses the K-bound heuristic.
