# Chemical-Thermodynamics
Thermodynamic flash calculation solvers written in Python. Implementing different equations of state (in order of difficulty).

## Peng–Robinson equation of state helper

The :mod:`flash.eos` module provides a lightweight Peng–Robinson
implementation that can evaluate compressibility factors, molar volumes, and
mass densities for pure components or mixtures. It relies on critical property
data stored in ``database/database.h5``. A minimal example:

```python
from flash import Mixture, PengRobinsonEOS

mixture = Mixture.from_names(["methane", "ethane"], [0.5, 0.5])
eos = PengRobinsonEOS()

result = eos.compute_volumes(mixture=mixture, temperature=300, pressure=10)
print(result.z_factors)
print(result.preferred_phase_volume("vapor"))
```

When attached to a :class:`flash.system.SimpleFlashCalculator`, the EOS adds
vapor-phase compressibility factor, molar volume, and density columns to the
reported phase dataframe.

## Antoine saturation pressure

Component metadata includes Antoine coefficients in the form
``ln(P_bar) = A - B / (T + C)`` with temperature in Kelvin. Each
:class:`flash.system.Component` exposes ``saturation_pressure(T, units="bar")`` to
evaluate vapor pressure in either bar or kPa:

```python
from flash import Component

methanol = Component.from_database("methanol")
print(methanol.saturation_pressure(338.0, units="bar"))
```

## Ideal-gas heat capacity

Heat capacity coefficients from ``database/Cp_gas.txt`` follow
``Cp/R = A + B*T + C*T^2 + D*T^-2 + E*T^3`` (T in Kelvin) with the tabulated
``B``, ``C``, ``D``, and ``E`` values scaled by 1e-3, 1e-6, 1e-5, and 1e-9,
respectively. Each :class:`flash.system.Component` can evaluate the correlation:

```python
from flash import Component

methane = Component.from_database("methane")
print(methane.heat_capacity(400.0))           # J/mol/K
print(methane.heat_capacity(400.0, units="R"))
```
