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
  "Liquid-liquid flash" below.
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
- **Not gamma-phi.** Passing both an `eos` and an `activity_model` to
  `gamma-gamma` is a `ModelError`; combined vapor-liquid equilibrium with an
  activity liquid is `flash_mode="gamma-phi"`, which is unchanged and still
  heuristic.
- The scope is the same as `stability_tp`'s: `"stable"` means no negative
  tangent-plane distance was found from the deterministic trial set, not a
  global proof.

### Post-split stability

Every two-phase result from the tangent-plane phi-phi path and from the
liquid-liquid path is re-tested: each converged phase is fed back into
`stability_tp` with the same model.

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
  gamma-phi case (activity-coefficient liquid tested against an
  equation-of-state vapor) is not supported** and raises `ModelError`: it needs
  a consistent pure-liquid reference fugacity that this package does not yet
  carry, and returning a plausible-looking wrong tangent plane would be worse
  than refusing. See ADR-0007.
- `pressure_Pa` is required and validated for both families, but it does not
  affect an activity-model result; `result.diagnostics["pressure_dependent"]`
  says which case applies, and `["model_family"]` is `"eos"` or `"activity"`.
- `feed_branch` and `phase_branch` are the minimum-Gibbs compressibility root
  labels for an EOS and are `None` for an activity model (single branch).
  For an EOS, fugacity coefficients for both the feed and every trial are taken
  from the compressibility root with the lowest Gibbs energy at that state.
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
