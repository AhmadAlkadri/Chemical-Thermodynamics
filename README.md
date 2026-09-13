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
- Gamma-phi mode currently uses `NRTL()` for the liquid activity model.
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
