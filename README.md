# Chemical-Thermodynamics
Chemical engineering thermodynamics utilities packaged as the `chemthermo`
library (src layout, SI units).

## Quickstart / Current API usage

Install from the repo (editable):
```bash
pip install -e .
```

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

See `examples/flash_tp_peng_robinson_demo.py` for a runnable script that prints a
table-style summary.

## EOS extension points (PC-SAFT scaffold)

The public repo defines a minimal residual-Helmholtz EOS protocol and registry
hooks in `chemthermo.eos`:

```python
from chemthermo.eos import EOSProtocol, get_eos, list_eos

print(list_eos())  # ["pcsaft"]
eos = get_eos("pcsaft", components=["Methane"])
```

`EOSProtocol` requires:
- `num_components()`
- `residual_helmholtz(temperature_K, volume_m3, composition)`

The PC-SAFT implementation is a placeholder. Parameter access is exposed via
`get_pcsaft_parameters()` and raises `PCSAFTParameterError` until a parameter
dataset is provided.

## Examples and notebooks

- Scripted demos: `examples/README.md`
- Jupyter notebooks: `notebooks/README.md`

## VLLE support

The public `chemthermo` repository defines a minimal VLLE API boundary only.
The working VLLE implementation lives in an optional private plugin package,
`chemthermo_vlle`.

Example usage:

```python
from chemthermo.vlle import get_vlle_engine

engine = get_vlle_engine()
result = engine.solve(
    temperature_K=300.0,
    pressure_Pa=101325.0,
    z=[0.5, 0.5],
)
print([phase.name for phase in result.phases])
```

If the plugin is not installed, `get_vlle_engine()` raises an error instructing
you to install `chemthermo_vlle`.

## Optional validation dependencies

Install the reference library used by validation tests:

```bash
pip install -e ".[validation]"
```

## Project brain

This repository maintains lightweight architectural context and decision history in `AGENT/`:

- `brain.md` — current architecture, invariants, and contributor/agent contract
- `adr/` — Architecture Decision Records (why key design choices were made)
- `steering-brief.md` — short summaries of recent changes and next steps

Contributors and AI agents should read `AGENT/brain.md` before making structural or API changes.
