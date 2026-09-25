# chemthermo

Phase equilibrium for chemical engineering in Python: tangent-plane
stability, multiphase TP flash, Peng-Robinson, PC-SAFT (with association and
polymers) and NRTL, in SI units, with every numerical claim traceable to a
recorded cross-check.

**Status: beta (0.4.0).** The science is checked against independent
implementations, but the public API may still change before 1.0. Python
>= 3.11; tested on macOS arm64 and Linux x86_64. MIT licence.

## Install

```bash
pip install chemthermo
```

or pin a tagged commit from GitHub:

```bash
pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@v0.4.0"
```

Runtime dependencies: numpy, pydantic, bibtexparser (< 2).

## Quick start

A flash that decides for itself how many phases there are:

```python
import chemthermo as ct

mixture = ct.Mixture.from_database(["Methane", "Ethane", "Propane"], [0.5, 0.3, 0.2])
result = ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS())

print(result.phase_names())                 # ['liquid', 'vapor']
print(result.vapor_fraction)                # 0.468...
print(result.phases["vapor"].composition.fractions)
```

Stability of a feed, and the phase it would split off:

```python
check = ct.stability_tp(
    ct.Mixture.from_database(["Methane", "n-Hexane"], [0.5, 0.5]),
    temperature_K=300.0, pressure_Pa=2.0e6, eos=ct.PengRobinsonEOS(),
)
print(check.status, check.tpd_min)          # unstable -1.2558...
print(check.phase_branch, check.trial_composition)
```

Three liquids, PC-SAFT, and residual properties:

```python
from chemthermo.eos import PCSAFTEOS

three = ct.flash_tp(
    ct.Mixture.from_database(["Water", "Ethanol", "n-Hexane"], [0.2, 0.4, 0.4]),
    temperature_K=280.0, pressure_Pa=101325.0, eos=ct.PengRobinsonEOS(),
)
print(three.phase_names())                  # ['liquid1', 'liquid2', 'liquid3']

water = PCSAFTEOS(components=("Water",))
rho = max(water.density_roots(temperature_K=300.0, pressure_Pa=1.0e5, composition=[1.0]))
props = water.residual_properties(temperature_K=300.0, density_mol_m3=rho, composition=[1.0])
print(props["h_res"], props["s_res_tp"])    # H^res/RT and S^res/R (ideal gas at same T, P)
```

From the command line:

```bash
chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 \
    --temperature-k 240 --pressure-pa 3e6 --format json
chemthermo stability-tp --components Methane,n-Hexane --z 0.5,0.5 \
    --temperature-k 300 --pressure-pa 2e6 --eos pc-saft
```

## What it does

| Capability | Models | Checked against (ledger case) |
| --- | --- | --- |
| Tangent-plane stability, incipient-phase composition (`stability_tp`) | Peng-Robinson, PC-SAFT, NRTL (liquid-liquid), NRTL + ideal gas | `thermo` 0.6 Michelsen test, same verdicts (S-4); Tessier et al. (2000) published global minima to 2e-11 (S-6, S-7); teqp (P-4) |
| TP flash with automatic phase count, 1-3 phases (`flash_tp`) | Peng-Robinson and PC-SAFT (`phi-phi`), NRTL (`gamma-gamma`), NRTL + Raoult (`modified-raoult`) | `thermo` `FlashVL` (F-1); teqp's traced PC-SAFT isotherm, tie lines to 2e-9 (P-5); FeOs flash and chemical potentials (P-7, P-8); a published multiphase Rachford-Rice table (V-4); a ternary VLLE tie triangle (V-1) |
| Peng-Robinson with per-pair `kij` | cubic EOS | `thermo` PRMIX (K-1); agreement floored at ~1e-4 in `ln phi` by chemthermo's rounded constants (S-4) |
| PC-SAFT residual Helmholtz, `Z`, `ln phi`, density roots | Gross & Sadowski 2001 | teqp, better than 3e-14 (P-1, P-3) |
| PC-SAFT association (2B scheme) | Gross & Sadowski 2002 | FeOs, association term to 3e-15 (P-6) |
| PC-SAFT polymers (segments per mass, long chains in log space) | monodisperse chains | FeOs, 5e-12 with matched constants (P-12 - P-14) |
| PC-SAFT residual `H`, `S`, `U`, `G` and `d(A^res/RT)/dT` | incl. association | teqp `Ar10` to 5e-16 (P-19); FeOs residual entropy/enthalpy to 2e-15 (P-20) |
| NRTL activity coefficients | NRTL | `thermo` NRTL to 9e-16 (N-2) |
| Command line: `tp-flash`, `stability-tp` | PR, PC-SAFT | golden JSON fixtures; exit-code contract |

Packaged data: critical constants, acentric factors and Antoine coefficients
for **82 components** (from Koretsky, *Engineering and Chemical
Thermodynamics*); PC-SAFT parameters for **16** (11 non-associating, plus
water, methanol, ethanol, 1-propanol and n-butanol with 2B association).

**What "checked" means here.** The comparisons above are code against
independent implementations (teqp, FeOs, `thermo`) or published worked
problems. They show the equations are implemented correctly. They say nothing
about how well a model with these parameters matches experiment, and nothing
in this repository claims that.

## Limits, stated plainly

- **"Stable" is not a proof.** It means no negative tangent-plane distance was
  found from a deterministic set of trial compositions. A phase that no trial
  reaches can be missed, and a flash's phase count inherits that limit.
- **Bring your own interaction parameters.** No binary `kij` table ships.
  The two packaged NRTL pairs are **synthetic placeholders** for demos: real
  NRTL work needs your own parameters (`NRTLParameters.from_pairs`). PC-SAFT
  with `kij = 0` gets water / hydrocarbon mutual solubilities badly wrong.
- **Residual properties only.** No ideal-gas heat capacities are packaged, so
  no total enthalpy, entropy or `Cp`; no `Cp^res` yet.
- **PC-SAFT scope:** no polar terms; only the 2B association scheme is
  validated; no induced association; polymers are monodisperse; no polymer
  parameters are packaged.
- **Modified Raoult is a low-pressure model** (ideal vapour, no Poynting
  correction); Antoine ranges are enforced, never extrapolated.
- **Cost:** a Peng-Robinson or activity-model flash takes milliseconds; PC-SAFT
  flashes take up to a few seconds (polymer and associating states the
  longest), and a three-phase PC-SAFT flash about 35 s.
- **Deprecated:** `flash_mode="gamma-phi"` (kept for the CLI v1 contract;
  emits `DeprecationWarning`). `chemthermo.vlle` was removed in 0.4.0.

Longer, model-specific limits are in the docs below.

## Documentation

- [Installation and development setup](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/installation.md)
- [TP flash](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/flash.md) - phase detection, flash modes, three phases, `kij`, NRTL
- [Phase stability](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/stability.md)
- [PC-SAFT](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/pcsaft.md) - parameters, association, polymers, residual properties, limits and their validation
- [Command line](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/cli.md)
- [Extending with a new equation of state](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/extending.md)
- [Benchmarks and the robustness map](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/docs/benchmarks.md)
- Runnable scripts: [`examples/`](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/examples/README.md); notebooks: [`notebooks/`](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/notebooks/README.md)
- [Changelog](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/CHANGELOG.md)

## How it is validated and developed

Every number quoted above has an entry in the validation ledger,
[`.agents/brain/validation-cases.md`](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/.agents/brain/validation-cases.md), with
its source, tolerance, achieved value and test. Design decisions are
recorded as ADRs in [`.agents/brain/adr/`](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/tree/main/.agents/brain/adr/). The tests that
compare against teqp, FeOs and `thermo` run with `pip install -e
".[validation]"` and skip cleanly otherwise. A 2505-state robustness map
(`python -m chemthermo.bench robustness`) sweeps every model family; its only
refusals are 5 states of the deprecated gamma-phi path.

Contributions: see [CONTRIBUTING.md](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/CONTRIBUTING.md).

## Licence

MIT - see [LICENSE](https://github.com/AhmadAlkadri/Chemical-Thermodynamics/blob/main/LICENSE).
