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

Notes:
- `status` is `"unstable"` when a trial converges to a stationary point with a
  negative tangent-plane distance, `"stable"` when trials converge and none
  does, and `"inconclusive"` when no trial converges.
- **What "stable" means here:** no negative tangent-plane distance was found
  from the deterministic trial set (two Wilson estimates plus one
  pure-component-dominant estimate per component). This is *not* a global proof
  of stability; Michelsen's test is a local stationary-point search and a
  stationary point that no initial estimate reaches can hide an instability.
- `tpd_min` is dimensionless (units of RT). At a stationary point it equals
  `-ln(sum_i W_i)`, so `sum(W) > 1` is the instability signal.
- Fugacity coefficients for both the feed and every trial are taken from the
  compressibility root with the lowest Gibbs energy at that state.
- `flash_tp` does **not** consume this yet; its single-phase decision is still a
  K-bound heuristic.

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
