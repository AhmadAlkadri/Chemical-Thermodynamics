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

> [!WARNING]
> **PC-SAFT is explicitly unsupported.**
> While some placeholder code for PC-SAFT exists in the codebase (`chemthermo.eos.pcsaft`), it is **non-functional and out of scope** for this repository. It will not be developed further.

## Examples and notebooks

- Scripted demos: `examples/README.md`
- Jupyter notebooks: `notebooks/README.md`

## VLLE support

> [!WARNING]
> **VLLE is explicitly out of scope.**
> This repository focuses on VLE (Vapor-Liquid Equilibrium) only. Any symbols or hooks related to VLLE (`chemthermo.vlle`) are legacy boundaries and are not supported.

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
