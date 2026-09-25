"""PC-SAFT residual properties of a methane / n-hexane mixture.

What this shows
---------------
``chemthermo.eos.PCSAFTEOS`` evaluates the non-associating PC-SAFT model of
Gross & Sadowski (2001) at a state given as ``(T, molar density, x)``:

  * ``A^res/(R T)``  - the reduced residual Helmholtz energy,
  * ``Z``            - the compressibility factor, ``1 + rho (d a_res/d rho)``,
  * ``P``            - ``Z rho R T``,
  * ``ln phi_i``     - one natural log per component.

The same mixture is shown at two densities: a gas-like one (200 mol/m^3, about
4.5 bar) and a liquid-like one (11000 mol/m^3, about 255 bar - a compressed
liquid, printed to show that the dense branch is reachable and well behaved).

What this deliberately does **not** do
--------------------------------------
It never solves for a density at a given pressure. This slice implements the
thermodynamics, not the root finder: the caller says which state they are on.
A density between the two spinodals has a negative ``Z`` and no fugacity
coefficients at all, which the script demonstrates at the end by catching the
refusal instead of printing a ``nan``. Density roots at ``(T, P)`` and the
wiring into ``flash_tp`` / ``stability_tp`` are the next slice
(``pcsaft-density-roots-flash``, ADR-0014).

Parameters and their provenance
-------------------------------
Packaged pure-component parameters from Table 1 of

    J. Gross and G. Sadowski, "Perturbed-Chain SAFT: An Equation of State Based
    on a Perturbation Theory for Chain Molecules", Ind. Eng. Chem. Res. 40
    (2001) 1244-1260,

verified against two independent secondary sources (see the ``provenance``
block in ``src/chemthermo/parameters/data/eos/pcsaft.json``). ``k_ij`` defaults
to zero; the last section shows a nonzero per-pair value on a different binary.

Needs no optional dependency. All inputs are SI; compositions are mole
fractions.
"""

from __future__ import annotations

import math

from chemthermo import ModelError
from chemthermo.eos import PCSAFTEOS

COMPONENTS = ("Methane", "n-Hexane")
COMPOSITION = [0.5, 0.5]
TEMPERATURE_K = 300.0


def _report(eos: PCSAFTEOS, temperature_K: float, density_mol_m3: float, x: list[float]) -> None:
    a_res = eos.residual_helmholtz(
        temperature_K=temperature_K,
        volume_m3=1.0 / density_mol_m3,
        composition=x,
    )
    z_factor = eos.compressibility_factor(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=x
    )
    pressure = eos.pressure_Pa(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=x
    )
    print(f"  molar density : {density_mol_m3:>14,.1f} mol/m^3")
    print(f"  molar volume  : {1.0 / density_mol_m3:>14.6e} m^3/mol")
    print(f"  A^res/(R T)   : {a_res:>14.10f}")
    print(f"  Z             : {z_factor:>14.10f}")
    print(f"  P             : {pressure:>14.6e} Pa  ({pressure / 1e5:,.4f} bar)")
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=x
    )
    for name, value in zip(eos.components, ln_phi):
        print(f"  ln phi[{name:<9}]: {value:>14.10f}   phi = {math.exp(value):.10f}")

    identity = sum(xi * li for xi, li in zip(x, ln_phi))
    expected = a_res + z_factor - 1.0 - math.log(z_factor)
    print(
        "  identity      : sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z  ->  "
        f"{identity:+.12f} vs {expected:+.12f}  (|d| = {abs(identity - expected):.2e})"
    )


def main() -> None:
    eos = PCSAFTEOS(components=COMPONENTS)
    print("PC-SAFT (Gross & Sadowski 2001, non-associating), kij = 0")
    print(f"components     : {', '.join(COMPONENTS)}")
    print(f"composition x  : {COMPOSITION}")
    print(f"temperature    : {TEMPERATURE_K:.2f} K")

    m, sigma_A, epsilon_k_K = eos.component_parameters()
    print("\nPure-component parameters (Gross & Sadowski 2001, Table 1):")
    print(f"  {'component':<12}{'m':>10}{'sigma / A':>12}{'eps/k / K':>12}")
    for name, m_i, sigma_i, eps_i in zip(COMPONENTS, m, sigma_A, epsilon_k_K):
        print(f"  {name:<12}{m_i:>10.4f}{sigma_i:>12.4f}{eps_i:>12.2f}")

    print("\n--- gas-like density ---")
    _report(eos, TEMPERATURE_K, 200.0, COMPOSITION)

    print("\n--- liquid-like density ---")
    _report(eos, TEMPERATURE_K, 11000.0, COMPOSITION)

    print("\n--- inside the spinodal: no fugacity coefficients exist there ---")
    density = 5000.0
    z_factor = eos.compressibility_factor(
        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=COMPOSITION
    )
    pressure = eos.pressure_Pa(
        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=COMPOSITION
    )
    print(f"  molar density : {density:>14,.1f} mol/m^3")
    print(f"  Z             : {z_factor:>14.10f}   (P = {pressure:.4e} Pa, negative)")
    try:
        eos.ln_fugacity_coefficients(
            temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=COMPOSITION
        )
    except ModelError as exc:
        print(f"  ln phi        : refused -- {exc}")
    else:  # pragma: no cover - defensive
        raise AssertionError("expected a ModelError at a negative compressibility factor")

    print("\n--- a nonzero per-pair kij (same contract as PengRobinsonEOS, ADR-0006) ---")
    pair = ("Methane", "n-Decane")
    x = [0.3, 0.7]
    for kij in (0.0, 0.03):
        model = PCSAFTEOS(components=pair, kij=kij)
        ln_phi = model.ln_fugacity_coefficients(
            temperature_K=350.0, density_mol_m3=100.0, composition=x
        )
        formatted = ", ".join(f"{value:+.10f}" for value in ln_phi)
        print(f"  kij = {kij:<5} ->  ln phi = ({formatted})")
    print(
        "  A pure fluid is never affected by kij: the diagonal of the kij matrix is\n"
        "  always zero, exactly as for Peng-Robinson."
    )


if __name__ == "__main__":
    main()
