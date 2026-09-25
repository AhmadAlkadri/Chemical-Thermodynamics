# PC-SAFT (`chemthermo.eos.PCSAFTEOS`)

The PC-SAFT equation of state (Gross & Sadowski 2001, association 2002): properties, phase equilibrium, parameters, polymers, association, residual properties, and its limits with the validation behind them.

`chemthermo.eos.PCSAFTEOS` implements the PC-SAFT equation of state of
Gross & Sadowski, *Ind. Eng. Chem. Res.* **40** (2001) 1244-1260 (hard chain
plus dispersion, ADR-0014), **plus the association term** of Gross & Sadowski,
*Ind. Eng. Chem. Res.* **41** (2002) 5510-5515 (ADR-0018), with the packaged
pure-component parameters of both papers' tables.

```python
from chemthermo.eos import PCSAFTEOS

eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
x = [0.5, 0.5]

eos.residual_helmholtz(temperature_K=300.0, volume_m3=1 / 200.0, composition=x)
# -0.0908711764                     A^res / (R T), volume_m3 is MOLAR volume

eos.compressibility_factor(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# 0.9102033819

eos.pressure_Pa(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# 454071.12

eos.ln_fugacity_coefficients(temperature_K=300.0, density_mol_m3=200.0, composition=x)
# [0.0365173954, -0.2096785687]     natural logs, one per component
```

The state of those four methods is always `(T, molar density, x)` (or
`(T, molar volume, x)` for `residual_helmholtz`). Binary interaction parameters
use the same contract as Peng-Robinson (ADR-0006) - a scalar applied to every
off-diagonal pair, or a mapping keyed by an unordered pair of component names:

```python
PCSAFTEOS(components=("Methane", "n-Decane"), kij={("Methane", "n-Decane"): 0.03})
```

## Phase equilibrium with PC-SAFT

Since ADR-0015 `PCSAFTEOS` also implements the pressure-based
`chemthermo.models.EquationOfState` interface, so it goes straight into the
stability test and the flash - **no solver changed for this**:

```python
import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.30, 0.70])
eos = PCSAFTEOS()                      # components come from the mixture

ct.stability_tp(mixture, temperature_K=300.0, pressure_Pa=3.0e6, eos=eos).status
# 'unstable'

result = ct.flash_tp(mixture, temperature_K=300.0, pressure_Pa=3.0e6, eos=eos)
result.vapor_fraction
# 0.16844005
result.phases["liquid"].composition.fractions
# (0.16086974, 0.83913026)
result.phases["vapor"].composition.fractions
# (0.98686249, 0.01313751)
```

Leave `components` out and the model takes its order from the `Mixture`; pass
it and a mixture whose names differ is rejected rather than silently reordered.

Under the hood, a call that names a pressure and a phase has to solve
`P_model(T, rho, x) = P` first. That root set is public:

```python
eos.density_roots(
    temperature_K=300.0, pressure_Pa=21858.084278856164, composition=[1.0],
    mixture=ct.Mixture.from_database(["n-Hexane"], [1.0]),
)
# (8.868596301913758, 7518.498733715524)     mol/m^3, ascending
```

Only mechanically stable roots (`dP/drho > 0`) are returned, ascending, so the
first is the vapour-like candidate and the last the liquid-like one -
`phase="vapor"` and `phase="liquid"` pick exactly those, which is the
Peng-Robinson convention verbatim. **One root is a normal answer**, not a
failure: at a dense or supercritical state both labels name the same state,
and the phase name that comes back is now a compressibility measurement, not
a tie-break (ADR-0017) - see `eos.phase_identity(...)` and the "Vapor/liquid
naming" note in [flash.md](flash.md). Which root a *split* phase sits on is decided per phase by
the stability test (ADR-0019), which is what lets two liquids coexist at a
pressure where a vapour root also exists. The spinodal-branch root is found and
discarded, never returned.
`PCSAFTEOS.molar_volume(..., phase=)` is the reciprocal of the selected root.

Two roots close enough together to fall inside one step of the scan grid
(`5e-4` in packing fraction) are not resolved, and the state is then reported
with one root fewer; that happens only where the isotherm is nearly tangent to
the target pressure, i.e. near a critical or spinodal state.

## Parameters

Eleven compounds ship with the package: Methane, Ethane, Propane, n-Butane,
n-Pentane, n-Hexane, n-Heptane, n-Octane, n-Decane, Nitrogen and Carbon
dioxide. Their provenance is recorded in
`src/chemthermo/parameters/data/eos/pcsaft.json`: the primary citation is the
paper above, but the paper is paywalled and **was not read directly** when this
file was written - the values were transcribed from two independent secondary
sources that cite it and agree digit for digit (FeOs' `gross2001.json` and
Clapeyron.jl's `PCSAFT_like.csv`), and the 42 universal constants were taken
from teqp's source and Wikipedia's PC-SAFT article, which likewise agree.
The five associating records (below) come from the 2002 paper, which is
paywalled too and was likewise **not read directly** (`pubs.acs.org` returns
HTTP 403); they were transcribed from FeOs' `gross2002.json` and Clapeyron.jl's
`PCSAFT_like.csv` / `PCSAFT_assoc.csv`, whose `source` column is that paper's
DOI, and the two agree digit for digit.

Supply your own with `PCSAFTParameters.from_records(...)`:

```python
from chemthermo import PCSAFTParameters
from chemthermo.eos import PCSAFTEOS

parameters = PCSAFTParameters.from_records(
    [{"name": "My fluid", "m": 2.5, "sigma_A": 3.6, "epsilon_k_K": 210.0}]
)
eos = PCSAFTEOS(components=("My fluid",), parameters=parameters)
```

A component with no record raises `PCSAFTParameterError`.

## Polymers (ADR-0022)

A polymer's PC-SAFT parameters are published per unit **mass** - `m/M` in
mol/g - because the chain length depends on the molar mass of the particular
sample, so one parameter set covers every molar mass of that polymer. Say that
with `segments_per_g` and `MW_g_mol`, and the segment number is derived for
you:

```python
import chemthermo as ct
from chemthermo import PCSAFTParameters, PCSAFTRecord
from chemthermo.eos import PCSAFTEOS

parameters = PCSAFTParameters.from_records(
    [
        PCSAFTRecord(                      # m = 0.0263 * 16400 = 431.32
            name="Polyethylene",
            segments_per_g=0.0263,         # mol/g, as published
            MW_g_mol=16400.0,              # this sample
            sigma_A=4.0217,
            epsilon_k_K=247.5,
            source="your citation here",
        ),
        PCSAFTRecord(name="n-Pentane", m=2.6896, sigma_A=3.7729, epsilon_k_K=231.20),
    ]
)

polymer = ct.Component.custom(          # no databank entry, no critical point
    "Polyethylene", mw_kg_per_mol=16.4, formula="(C2H4)n", volatile=False
)
mixture = ct.Mixture.from_components(
    [polymer, ct.Component.from_database("n-Pentane")],
    [2.314804e-04, 0.99976852],         # 5 wt% polymer, as mole fractions
)
result = ct.flash_tp(
    mixture,
    temperature_K=453.0,
    pressure_Pa=8.0e6,
    eos=PCSAFTEOS(parameters=parameters, kij=-0.006),
)
# -> liquid1 / liquid2, x_polymer = 9.75e-04 and 1.15e-05, vapor_fraction None
```

Exactly one of `m` and `segments_per_g` may be given; both, or neither, raises,
because they are two spellings of the same parameter.

`Component.custom(name, *, mw_kg_per_mol, formula=, tc_k=, pc_pa=, omega=,
volatile=, antoine=, source=)` builds a component the databank does not carry.
The molar mass is **required** and is in kg/mol like every other molar mass
here; the critical constants are **optional**, because a polymer has none, and
`component.tc_k` raises `PropertyNotFoundError` rather than returning a
placeholder if you ask for one you did not supply. `volatile=False` replaces
that component's Wilson K-value *estimate* by a fixed `1e-10` ("essentially
absent from the vapour-like trial"), which is what keeps the deterministic
stability trial set complete without critical constants. It is an initial
estimate only: no verdict, composition or fugacity is a function of it, and
nothing else in the package reads `volatile`. The packaged databank schema is
untouched - `ComponentData` still requires `Tc` / `Pc` / `omega` and
`schema_version` is still 1.

**What is and is not claimed.**

- **No polymer parameters are packaged, and none should be read from here.**
  The polyethylene numbers above are a *cited test fixture*,
  `tests/fixtures/pcsaft/martini2009_polymers.json`: as tabulated by Martini,
  Cismondi, Barbosa & Brignole, *Sep. Sci. Technol.* **44** (2009) (author
  manuscript, CONICET open repository), citing Gross & Sadowski, *IECR* **41**
  (2002) 1084. **That primary table was not read** - it is paywalled - and no
  second open source printing the same three values was found. Unlike the
  packaged 2001 records, which two independent sources confirm digit for digit,
  these rest on one. Supply and cite your own.
- **A polymer here is monodisperse**: one chain length, one component. Real
  samples are not (the ones these parameters describe have polydispersities of
  1.14 to 2.94), and representing a distribution as several pseudo-components
  is a capability this package does not have.
- **Nothing is compared against measurement.** The `k_ij = -0.006` above was
  fitted *by the cited source* to cloud-point data this package never touches.
  The qualitative behaviours that are checked - the split appearing on heating
  at fixed pressure, and the cloud-point pressure rising with `k_ij` - are
  compared with statements in that source's text, not with digitized figures.
- **A long chain can exceed the exponential's range.** A polyethylene of
  `Mw = 53000` (`m = 1393.9`) in n-pentane at 453 K has `ln phi = -1690.6`, and
  `exp` of that is an exact `0.0`. `EquationOfState.log_fugacity_coefficients`
  (optional, `None` by default; implemented by `PCSAFTEOS`) is the log-space
  route the flash falls back to **only** where `phi` does not exist as a
  double, so every previously converging number is unchanged. Measured: 0 of
  274 branch evaluations in log space for water / n-hexane, 0 of 600 for the
  `Mw = 16400` polymer, 1025 of 1025 for the `Mw = 53000` one.
- **Below the solvent's saturation pressure the equilibrium is vapour-liquid,
  and the split runs in log mole numbers** (ADR-0024). At 453 K and under about
  2.6 MPa n-pentane still has a vapour root, and the answer is a solvent vapour
  over a solvent-swollen melt: at 1 MPa the melt holds `x_pentane = 0.971468`
  (13 wt% solvent) and the vapour holds `ln y_polymer = -450.53`. Neither the
  seed nor the iteration fits in linear mole numbers - the tangent-plane
  minimizer is an essentially pure melt whose K-values span `1e+180` and
  bracket no Rachford-Rice root at all - so `flash_tp` seeds and finishes that
  split in `u = ln n` instead, and says so in
  `diagnostics["converged_stage"] == "second-order-log"`. Every state that
  converged before ADR-0024 still takes the ordinary path, bit for bit.

  ```python
  result = ct.flash_tp(mixture, temperature_K=453.0, pressure_Pa=1.0e6, eos=eos)
  result.phases["liquid"].composition.fractions   # (0.0285316, 0.9714684) - the melt
  result.vapor_fraction                           # 0.99189
  result.diagnostics["log_space_ln_x_min"]        # -450.53 = ln y_polymer
  ```

- **A mole fraction may come back as exactly `0.0`.** For a `Mw = 53000` chain
  the vapour's polymer mole fraction is `exp(-1315)`, and `0.0` is the nearest
  double there is. The number is not lost - it is
  `diagnostics["log_space_ln_x_min"]` for the component named in
  `["log_space_ln_x_min_component"]` - and the material balance is then
  *exact*, because the melt holds every mole of polymer the feed had. Read
  `diagnostics["log_space_residual"]` for that component's equal-fugacity
  residual: `fugacity_residual` is taken over the components present in both
  phases and cannot see it.

- **The longest chain is in reach since ADR-0025.** `Mw = 53000` (`m = 1393.9`)
  at 0.5 and 1 MPa used to raise, because the stability iteration clamped
  `ln W` to `[-700, 700]` and the melt's stationary point sits at
  `ln W_polymer = 1452`. The normalization is now done in logarithms where -
  and only where - that clamp would have engaged, so the melt is found in three
  successive substitutions (`tpd_min = -1452.21`) and `flash_tp` returns the
  verified vapour-liquid split: melt `x_solvent = 0.97881` (5.92 wt% solvent),
  `beta_vapor = 0.99662`, `ln y_polymer = -1507.98`. A stationary point whose
  mole numbers leave machine range is reported as
  `StabilityResult.trial_ln_W`, because the normalized `trial_composition`
  rounds to `(1.0, 0.0)` there; `diagnostics["ln_sum_W"]` is equation (7)'s
  `-tpd`, and `sum_W` itself is then `inf`.

- **The whole 0.3-3.6 MPa band converges since ADR-0026.** Six states of that
  68-state sweep still raised at the previous slice - 0.3, 2.8 and 2.9 MPa,
  where the log-space stage spent its budget crawling next to the trivial
  solution, and 3.0 to 3.2 MPa, where a `K` of `1e-18` made
  `1 + (K - 1)` cancel to an exact zero and Rachford-Rice reported a bracket it
  had. Both are repaired where, and only where, the flash would otherwise
  refuse: `diagnostics["log_space_curvature_safeguard"]` and
  `diagnostics["rachford_rice_convex_denominators"]` say which. No tolerance
  changed, and the 62 states that converged before are bit-identical.

- **Since ADR-0028 the robustness map refuses nothing.** The 36 states the map
  of ADR-0027 refused were all this system, and all three causes were a stage
  being started in the wrong place rather than stepping wrongly; the log-space
  stage is now retried from the tangent-plane stationary point, at a phase
  fraction the lever rule bounds, and the stability test is retried once with
  substitution's full budget. Each retry runs only where the flash was about to
  refuse, and all 2074 states that converged before are identical in the
  regenerated record.

See ADR-0022, ADR-0024, ADR-0025, ADR-0026 and ADR-0028, validation Cases P-12,
P-13, P-14, P-15, P-16 and P-17, `examples/basic/pcsaft_polymer_demo.py`,
`examples/validation/20_pcsaft_polymer_vs_feos.py`,
`examples/validation/21_pcsaft_polymer_vle.py` and
`examples/validation/22_stability_log_space.py`.

## Association (ADR-0018)

Five more compounds ship with association parameters from Gross & Sadowski,
*IECR* **41** (2002) 5510, all in the **2B** scheme (one proton-donor site and
one proton-acceptor site per molecule, bonding A-B only): **Water, Methanol,
Ethanol, 1-Propanol and n-Butanol**. Nothing about the call signatures changes
- a mixture containing one of them simply gets the extra term:

```python
from chemthermo.eos import PCSAFTEOS

eos = PCSAFTEOS(components=("Water",))
eos.associates()
# True

eos.residual_helmholtz_terms(
    temperature_K=300.0, volume_m3=1 / 55000.0, composition=[1.0]
)
# {'hard-chain': 5.0793036678, 'dispersion': -8.8028886238,
#  'association': -5.7039482251, 'total': -9.4275331811}

eos.site_fractions(temperature_K=300.0, density_mol_m3=55000.0, composition=[1.0])
# [0.0356448106, 0.0356448106]    fraction of sites NOT hydrogen bonded
```

`residual_helmholtz_terms` and `site_fractions` are the two new methods; both
are inspection helpers, and `site_fractions` returns `[]` for a non-associating
mixture (where `residual_helmholtz_terms` reports `'association': 0.0`).

Cross-association between two associating components is computed from the pure
parameters by the Wolbach-Sandler rules (arithmetic mean in
`epsilon^AB`, geometric mean in `kappa^AB` scaled by
`[sqrt(sigma_i sigma_j) / (0.5 (sigma_i + sigma_j))]^3`), so a water/ethanol
mixture needs no extra data. `k_ij` still defaults to 0 and still only enters
the dispersion term.

Supply your own sites with an `association` block:

```python
from chemthermo import PCSAFTParameters

parameters = PCSAFTParameters.from_records(
    [
        {
            "name": "My alcohol", "m": 2.5, "sigma_A": 3.6, "epsilon_k_K": 210.0,
            "association": {"scheme": "2B", "kappa_ab": 0.03, "epsilon_ab_k_K": 2500.0},
        }
    ]
)
```

`na` / `nb` default to `1` / `1`; `scheme` is a label and is checked against
them for the names `2B`, `3B` and `4C`. `PCSAFTAssociationRecord` is exported
from `chemthermo` if you prefer to build the block as an object.

## Limits, stated plainly

- **The 2B scheme is what is validated.** Other site counts (`na`, `nb`) are
  implemented by the same general equations and will run, but nothing
  cross-checks them, and **induced association** (a non-associating component
  solvating with an associating one) is not modelled at all. The **polar**
  (dipolar / quadrupolar) terms are still not implemented, so a strongly polar
  non-associating compound is still out of scope.
- **`kij = 0` is a poor model for water with a hydrocarbon.** The packaged set
  has no binary parameters at all, and PC-SAFT without a fitted one is known to
  get water / hydrocarbon mutual solubilities wrong by an order of magnitude.
  The cross-checks below are code checks against another implementation, not
  evidence about the model.
- **A phase count is only as good as the stability test behind it.** Up to
  `FlashSettings(max_phases=...)` phases (default 3) are found by phase
  addition/removal (ADR-0011, ADR-0020); each phase sits on the density root
  the stability test found it on, and two liquids are named `liquid1` /
  `liquid2` with `vapor_fraction = None` (ADR-0019). A three-phase PC-SAFT
  flash takes about 35 s. See "Three phases from an equation of state" in
  [flash.md](flash.md) and validation Cases P-9 and P-10.
- **The `(T, rho, x)` methods still choose nothing.** `compressibility_factor`,
  `pressure_Pa` and `ln_fugacity_coefficients` take the density you give them.
  Between the two spinodals `Z` is negative, `pressure_Pa` returns the (real)
  negative pressure, and `ln_fugacity_coefficients` raises `ModelError` instead
  of returning a `nan`. Use `density_roots` (or the `phase=` interface) when
  you want the model to choose.
- **`kij` is yours to justify.** The packaged parameter set is pure components
  only; no PC-SAFT binary-interaction table ships with this package, and any
  `kij` used in the docs or examples (0.03 for methane / n-decane) is
  **illustrative**, not a literature-validated value.
- **Temperature derivative and residual properties** (ADR-0034, including
  association since its C2 amendment): `residual_helmholtz_temperature_derivative` and
  `residual_properties` (`h_res`, `u_res`, `s_res_tv` / `s_res_tp`, `g_res_tv` /
  `g_res_tp`, reduced by `RT` or `R`; the suffix names the ideal-gas reference,
  same `T` and volume or same `T` and pressure). Checked against FeOs for water,
  ethanol and their mixtures (2.1e-15 with matched constants). **Residual only**: no total enthalpy, entropy or `Cp` (the
  databank has no ideal-gas heat capacities), no `Cp^res`, and no phase
  densities in `FlashResult` (compute them with `density_roots` at the
  converged composition).

  ```python
  from chemthermo.eos import PCSAFTEOS

  eos = PCSAFTEOS(components=("n-Hexane",))
  rho = max(eos.density_roots(temperature_K=300.0, pressure_Pa=1.0e6, composition=[1.0]))
  props = eos.residual_properties(temperature_K=300.0, density_mol_m3=rho, composition=[1.0])
  props["h_res"] * 8.314462618 * 300.0   # H^res in J/mol, about -31.5 kJ/mol
  ```
- Validated against [teqp](https://github.com/usnistgov/teqp) (NIST, MIT,
  automatic differentiation) to better than 3e-14 in `A^res/RT`, `Z` and
  `ln phi` over fourteen states; against its `pure_VLE_T` saturation solver for
  n-hexane at 300 K and 400 K (densities to 2.2e-16 relative); and, for phase
  equilibrium, against teqp's own traced methane / n-hexane 300 K isotherm -
  seven tie lines to `|dx1| <= 2.0e-9` and `|dy1| <= 3.7e-12`, phase densities
  to 1.3e-9 relative, and teqp's own fugacity coefficients evaluated at
  chemthermo's converged phases giving equal fugacities to 3.8e-9 relative.
  Bubble pressures located by bisecting `stability_tp`'s verdict match teqp's
  `mix_VLE_Tx` to 1.4e-8 relative. See validation Cases P-0 to P-5,
  `examples/validation/13_pcsaft_vs_teqp.py` and
  `examples/validation/14_pcsaft_flash_vs_teqp.py`.
- The **association** term is validated against
  [FeOs](https://github.com/feos-org/feos) (MIT OR Apache-2.0, automatic
  differentiation), which teqp cannot serve because its PC-SAFT has no
  association. Over eighteen states - pure water and pure ethanol at four each,
  plus water/ethanol, water/n-hexane and a ternary - the association
  contribution agrees to **3.1e-15** and, once FeOs's own universal constants
  are substituted in, `A^res/RT`, `Z` and `ln phi` agree to **1.3e-11**. (FeOs
  hard-codes the 2001 paper's 42 universal constants to fourteen figures where
  the paper prints ten; as shipped that floors the comparison at 3.9e-09 in `Z`
  and 1.7e-06 in `ln phi`, which is an input difference, not a model
  difference.) Pure-water saturation at 373.15 K matches FeOs's own solver to
  3.7e-10 relative, and FeOs's fugacities evaluated at chemthermo's converged
  phases give equal fugacities to 6.1e-10 (water/ethanol VLE) and 1.2e-11
  (water/n-hexane LLE). See validation Cases P-6 and P-7 and
  `examples/validation/16_pcsaft_association_vs_feos.py`.
- The **liquid-liquid** tie line is validated against FeOs's own two-phase
  flash as well as its chemical potentials: water / n-hexane at 298.15 K, at
  1 atm and at 1 MPa, compositions to 6.5e-09 and phase amounts to 3.3e-09
  absolute, densities to 5.6e-09 relative, and FeOs's chemical potentials at
  chemthermo's phases equal to 1.2e-11 with matched universal constants
  (2.6e-06 as shipped). See validation Case P-8 and
  `examples/validation/17_pcsaft_lle_vs_feos.py`.
- The **polymer/solvent** capability is validated against FeOs the same way:
  `A^res/RT`, `Z` and `ln phi` for polyethylene melts and polyethylene /
  n-pentane liquids over eleven states agree to 5.0e-12 with matched universal
  constants (1.4e-06 as shipped), and FeOs's chemical potentials at
  chemthermo's converged phases are equal to 4.5e-13. Two facts about the
  reference are recorded rather than hidden: FeOs's own `tp_flash` **raises**
  on this system at 5 and 8 MPa and returns a degenerate pair at 3 MPa
  (it converges at 10 MPa, where the tie line agrees to 1.2e-09 in the
  solvent-rich polymer mole fraction), and `k_ij` cannot be given to FeOs's
  PC-SAFT from Python in feos 0.10.1, so those comparisons run at `k_ij = 0` on
  both sides. The fitted-`k_ij` answer is checked against an independently
  written two-equation Newton instead, to 1.7e-15. See validation Cases P-12
  and P-13 and `examples/validation/20_pcsaft_polymer_vs_feos.py`.
- The polymer/solvent **vapour-liquid** split of ADR-0024 is validated the same
  way, one regime lower in pressure: at 0.5, 1 and 2 MPa FeOs's chemical
  potentials at chemthermo's converged phases are equal to **2.6e-12** with
  matched universal constants (2.5e-07 as shipped), and the melt's solvent
  mole fraction agrees with a one-dimensional equal-fugacity solve - the vapour
  taken as exactly pure solvent - to 1.2e-15 - 2.5e-14 absolute, at the fitted
  `k_ij = -0.006`. The polymer's own equal-fugacity condition, which is only
  writable in logarithms (`ln f_polymer` of -492 to -531), closes to 7.4e-13.
  See validation Case P-14 and
  `examples/validation/21_pcsaft_polymer_vle.py`.
