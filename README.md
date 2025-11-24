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
