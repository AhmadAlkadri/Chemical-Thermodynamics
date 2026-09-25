# Extending: equation-of-state extension points

How to plug in a new equation of state.

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
