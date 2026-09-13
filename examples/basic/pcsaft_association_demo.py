"""PC-SAFT with association: pure water, and a water / n-hexane liquid-liquid split.

What this shows
---------------
``chemthermo.eos.PCSAFTEOS`` now carries the association term of

    J. Gross and G. Sadowski, "Application of the Perturbed-Chain SAFT Equation
    of State to Associating Systems", Ind. Eng. Chem. Res. 41 (2002) 5510-5515,

so water and the C1-C4 1-alkanols are usable components rather than components
the model had to refuse (ADR-0018). Five packaged records gained an
``association`` block with the 2B scheme (one proton-donor site and one
proton-acceptor site per molecule): Water, Methanol, Ethanol, 1-Propanol and
n-Butanol.

The script prints four things:

1. **The site fractions.** The new inner unknown is ``X``, the fraction of
   association sites that are *not* hydrogen bonded. It is 1 in a dilute gas
   and a few per cent in liquid water, and it is what the association term is
   a function of.
2. **Pure water properties** at three states, with the association
   contribution shown next to the hard-chain and dispersion ones, so its size
   is visible rather than buried in a total.
3. **Water's saturation state at 373.15 K**, found by equal fugacity on the
   two density roots. PC-SAFT with these parameters puts it near 100.9 kPa;
   the experimental value at that temperature is 101.325 kPa by the definition
   of the normal boiling point. That line is a **remark about the model**, not
   a claim about this code.
4. **A water / n-hexane liquid-liquid split** at 298.15 K and 1 MPa, from the
   unchanged ``flash_tp``. Two honest caveats come with it, printed by the
   script:

   * ``k_ij = 0``. PC-SAFT with no binary interaction parameter is known to be
     poor for water / hydrocarbon mutual solubilities, so the numbers below
     are a demonstration that the machinery works, not a prediction anyone
     should use.
   * The pressure is 1 MPa, not 1 atm, and that is not cosmetic. At 1 atm the
     isotherm still has a vapour density root, and the ``phi-phi`` split can
     only pair a vapour-root phase with a liquid-root one - so it converges on
     a spurious vapour-liquid pair, the post-split stability test catches it,
     and ``flash_tp`` raises. The script shows that too. Above about 0.6 MPa
     the vapour root is gone, both phases sit on the single (liquid) root, and
     the same machinery produces the real liquid-liquid tie line.

Labels: both converged phases are liquids by the ADR-0017 compressibility
criterion, but the phi-phi path has no ``liquid1`` / ``liquid2`` naming, so it
falls back to the Wilson ranking and calls them ``"liquid"`` and ``"vapor"``.
``vapor_fraction`` is then the hexane-rich *liquid*'s fraction. The script
prints the labels it actually got; fixing them is the next slice's job.

Validated against FeOs in ``examples/validation/16_pcsaft_association_vs_feos.py``
(validation Cases P-6 and P-7).

Needs no optional dependency. All inputs are SI; compositions are mole
fractions.
"""

from __future__ import annotations

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

WATER = ("Water",)
WATER_HEXANE = ["Water", "n-Hexane"]
FEED = [0.5, 0.5]

#: (T/K, rho/(mol/m^3), what it is)
PURE_WATER_STATES = [
    (300.0, 55000.0, "compressed liquid"),
    (373.15, 48755.5163, "saturated liquid"),
    (373.15, 33.2726, "saturated vapour"),
]

LL_TEMPERATURE_K = 298.15
LL_PRESSURE_PA = 1.0e6
ATMOSPHERIC_PA = 101325.0


def _pure_water_properties() -> None:
    eos = PCSAFTEOS(components=WATER)
    print("1) Pure water, term by term (all dimensionless, per mole)")
    print(
        f"   {'state':>18}  {'a_hc':>10}  {'a_disp':>10}  {'a_assoc':>10}  "
        f"{'a_res':>10}  {'Z':>9}  {'X (site)':>9}"
    )
    for temperature, density, label in PURE_WATER_STATES:
        terms = eos.residual_helmholtz_terms(
            temperature_K=temperature, volume_m3=1.0 / density, composition=[1.0]
        )
        sites = eos.site_fractions(
            temperature_K=temperature, density_mol_m3=density, composition=[1.0]
        )
        z_factor = eos.compressibility_factor(
            temperature_K=temperature, density_mol_m3=density, composition=[1.0]
        )
        print(
            f"   {label:>18}  {terms['hard-chain']:>10.5f}  {terms['dispersion']:>10.5f}  "
            f"{terms['association']:>10.5f}  {terms['total']:>10.5f}  "
            f"{z_factor:>9.5f}  {sites[0]:>9.6f}"
        )
    print(
        "   X is the fraction of association sites that are NOT bonded: ~1 in the\n"
        "   vapour, a few per cent in the liquid. The association term is what that\n"
        "   hydrogen-bond network costs in free energy."
    )


def _saturation(eos: PCSAFTEOS, temperature: float, low: float, high: float) -> tuple[float, ...]:
    """Bisect on equal fugacity between the two density roots."""

    def gap(pressure: float) -> float:
        roots = eos.density_roots(
            temperature_K=temperature, pressure_Pa=pressure, composition=[1.0]
        )
        return (
            eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=roots[-1], composition=[1.0]
            )[0]
            - eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=roots[0], composition=[1.0]
            )[0]
        )

    low_value = gap(low)
    for _ in range(80):
        middle = 0.5 * (low + high)
        value = gap(middle)
        if low_value * value <= 0.0:
            high = middle
        else:
            low, low_value = middle, value
    pressure = 0.5 * (low + high)
    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=pressure, composition=[1.0])
    return pressure, roots[-1], roots[0]


def _water_saturation() -> None:
    eos = PCSAFTEOS(components=WATER)
    pressure, rho_liquid, rho_vapour = _saturation(eos, 373.15, 5.0e4, 2.0e5)
    print("\n2) Pure water saturation at 373.15 K (equal fugacity on the two roots)")
    print(f"   Psat      = {pressure:>14,.2f} Pa")
    print(f"   rho_L     = {rho_liquid:>14,.4f} mol/m^3")
    print(f"   rho_V     = {rho_vapour:>14,.4f} mol/m^3")
    print(
        f"   Model versus experiment (a remark, not a check): water's normal boiling\n"
        f"   point is 373.15 K at 101,325 Pa by definition, so this model is "
        f"{(101325.0 - pressure) / 101325.0 * 100.0:.2f} %\n"
        "   low there. Cross-checked against FeOs, not against measurement."
    )


def _atmospheric_attempt() -> None:
    mixture = ct.Mixture.from_database(WATER_HEXANE, FEED)
    eos = PCSAFTEOS()
    print("\n3) Water / n-hexane at 298.15 K and 1 atm: the feed is unstable ...")
    stability = ct.stability_tp(
        mixture, temperature_K=LL_TEMPERATURE_K, pressure_Pa=ATMOSPHERIC_PA, eos=eos
    )
    print(f"   stability_tp -> {stability.status!r}, tpd_min = {stability.tpd_min:.6e}")
    try:
        ct.flash_tp(mixture, temperature_K=LL_TEMPERATURE_K, pressure_Pa=ATMOSPHERIC_PA, eos=eos)
    except ct.ConvergenceError as error:
        first_sentence = str(error).split(". ")[0]
        print(f"   ... and flash_tp refuses rather than returning it:\n     {first_sentence}.")
    else:  # pragma: no cover - would mean the limitation below was fixed
        print("   ... and flash_tp returned a result (the documented limitation is gone).")
    print(
        "   Why: at 1 atm the isotherm still has a vapour density root, and the phi-phi\n"
        "   split pairs one vapour-root phase with one liquid-root phase - it has no way\n"
        "   to put both phases on the liquid root. The post-split stability test catches\n"
        "   the spurious pair. See ADR-0018; the fix is the flash-phase-addition-eos slice."
    )


def _liquid_liquid_split() -> None:
    mixture = ct.Mixture.from_database(WATER_HEXANE, FEED)
    eos = PCSAFTEOS()
    print("\n4) The same system at 1 MPa, where the vapour root no longer exists")
    roots = eos.density_roots(
        temperature_K=LL_TEMPERATURE_K,
        pressure_Pa=LL_PRESSURE_PA,
        composition=FEED,
        mixture=mixture,
    )
    print(f"   density roots at the feed composition: {len(roots)} -> {roots[0]:,.2f} mol/m^3")
    result = ct.flash_tp(
        mixture, temperature_K=LL_TEMPERATURE_K, pressure_Pa=LL_PRESSURE_PA, eos=eos
    )
    for name, phase in result.phases.items():
        fractions = list(phase.composition.fractions)
        identity = eos.phase_identity(
            mixture=mixture,
            temperature_K=LL_TEMPERATURE_K,
            pressure_Pa=LL_PRESSURE_PA,
            composition=fractions,
            phase="liquid",
        )
        density = eos.density_roots(
            temperature_K=LL_TEMPERATURE_K,
            pressure_Pa=LL_PRESSURE_PA,
            composition=fractions,
            mixture=mixture,
        )[0]
        amounts = ", ".join(
            f"{component} {value:.8f}"
            for component, value in zip(mixture.component_names, fractions)
        )
        print(
            f"   phase named {name!r:>9}: {amounts}\n"
            f"     measured identity = {identity!r}, rho = {density:,.2f} mol/m^3"
        )
    print(
        f"   phase_label_method = {result.diagnostics['phase_label_method']!r}, "
        f"vapor_fraction = {result.vapor_fraction:.6f}"
    )
    print(
        f"   dG/RT of the split = {result.diagnostics['delta_g_split_rt']:.6e} (negative, so it\n"
        f"   is an answer and not just a fixed point); fugacity residual "
        f"{result.diagnostics['fugacity_residual']:.2e};\n"
        f"   post-split stability {result.diagnostics['post_split_status']!r}."
    )
    print(
        "   Labels: both phases measure as liquids, but the phi-phi path cannot name them\n"
        "   'liquid1' / 'liquid2', so ADR-0017 falls back to the Wilson ranking and\n"
        "   'vapor_fraction' is really the hexane-rich liquid's fraction. Recorded, not hidden."
    )
    print(
        "   k_ij = 0 here. PC-SAFT without a fitted binary parameter is known to be poor\n"
        "   for water / hydrocarbon mutual solubilities; this is a check of the code, not\n"
        "   of the model. Commonly quoted experimental values are ~5e-4 mole fraction\n"
        "   water in hexane and ~2e-6 hexane in water at 298 K (tabulated figures, not\n"
        "   verified against a primary source here)."
    )


def main() -> None:
    print("PC-SAFT with association (Gross & Sadowski 2002), ADR-0018")
    print("=" * 78)
    _pure_water_properties()
    _water_saturation()
    _atmospheric_attempt()
    _liquid_liquid_split()
    print("\n" + "=" * 78)
    print(
        "Scope: the 2B scheme is what is packaged and validated. General (na, nb) site\n"
        "counts are implemented but not cross-checked, induced association is not\n"
        "modelled, and there is still no temperature derivative anywhere in PC-SAFT."
    )


if __name__ == "__main__":
    main()
