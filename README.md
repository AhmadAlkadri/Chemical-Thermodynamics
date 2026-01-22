# Chemical-Thermodynamics
Work in progress: this repository is being refactored into the `chemthermo`
package (src layout, SI units). Legacy modules and notebooks remain for now and
will be archived during the migration.

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

## Enthalpy, entropy, and internal energy

Each :class:`flash.system.Component` now provides ideal-gas enthalpy, entropy,
and internal energy relative to the reference state 298.15 K and 1 bar. The
phi–phi flash implementation combines those ideal contributions with
Peng–Robinson residuals to report ``h[J/mol]``, ``s[J/mol/K]``, and ``u[J/mol]``
for each phase.

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

## VLE helpers

With :class:`flash.eos.PengRobinsonEOS` attached to a
:class:`flash.system.SimpleFlashCalculator`, ``flash`` performs a phi–phi flash
and can trace pure-component P–T envelopes and binary T–x curves:

```python
from flash import Component, Mixture, PengRobinsonEOS, SimpleFlashCalculator

calc = SimpleFlashCalculator(eos_model=PengRobinsonEOS())
methanol = Component.from_database("methanol")
pt_curve = calc.pure_PT_envelope(methanol, temperatures=[320, 340, 360])

mixture = Mixture.from_names(["methanol", "water"], [0.5, 0.5])
tx = calc.binary_Tx_diagram(mixture.components, pressure=1.0, x1_grid=[0.1, 0.5, 0.9])
```
