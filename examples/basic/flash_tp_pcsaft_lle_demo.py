"""Liquid-liquid equilibrium from an equation of state (ADR-0019).

What this shows
---------------
``flash_tp(..., eos=PCSAFTEOS())`` on water / n-hexane at 298.15 K and 1 atm
returns **two liquids**, named ``liquid1`` and ``liquid2``, with
``vapor_fraction = None``.

That state used to raise. The phi-phi split evaluated one phase on the model's
``"liquid"`` density root and the other on its ``"vapor"`` root, always - so at
any state where a vapour root still exists, the only pair it could offer was a
vapour-liquid one. At 298.15 K and 1 atm water / n-hexane *does* still have a
vapour root, but the answer is two liquids (the two pure vapour pressures sum
to about 23 kPa, far below 1 atm), so the split converged on a water-rich
liquid against a hexane-rich vapour whose Gibbs energy was **above** the
feed's, and the post-split stability test - correctly - refused it. Validation
Case P-7(iii) recorded that as a limitation of the flash path.

Since ADR-0019 each phase sits on the branch the tangent-plane stability test
found *that phase* on. Here the feed branch and the incipient branch are both
``"liquid"``, so both phases stay on the liquid root and the split is a
liquid-liquid one. The two converged phases are then named from
``EquationOfState.phase_identity`` (ADR-0017, a compressibility criterion):
both measure as liquids, so they are ``liquid1`` / ``liquid2`` rather than
``liquid`` / ``vapor``, and there is no vapour whose fraction could be
reported.

The script prints four sections:

1. **The state, and why it used to be hard.** The stability verdict, the two
   branches it reports (both ``"liquid"`` - that is the information the split
   now uses), and the two density roots the feed composition really has, so
   the premise of the old failure is shown rather than asserted.
2. **The liquid-liquid split at 1 atm**, with every verification residual and
   the measured identity, density and compressibility ratio of each phase.
3. **Three feeds, one tie line.** ``z = 0.5/0.5``, ``0.2/0.8`` and ``0.8/0.2``
   return the same two compositions with different amounts, checked against the
   lever rule.
4. **The same tie line at 1 MPa**, where the vapour root no longer exists. That
   split was already reachable before this slice, but it came back named
   ``liquid`` / ``vapor`` by the Wilson-ranking fallback with a
   ``vapor_fraction`` that was really the hexane-rich *liquid*'s fraction.

One honest caveat, printed by the script: ``k_ij = 0``. PC-SAFT without a
fitted binary interaction parameter is known to be poor for water /
hydrocarbon mutual solubilities, and it is about an order of magnitude out on
both here. This demo shows that the machinery works, not that the numbers are
a prediction anyone should use; the cross-check against an independent
implementation is ``examples/validation/17_pcsaft_lle_vs_feos.py``.

Needs no optional dependency. All inputs are SI; compositions are mole
fractions.

Run from the repo root::

    python examples/basic/flash_tp_pcsaft_lle_demo.py
"""

from __future__ import annotations

from typing import Sequence

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD

NAMES = ("Water", "n-Hexane")
TEMPERATURE_K = 298.15
ATMOSPHERE_PA = 101325.0
HIGH_PRESSURE_PA = 1.0e6
FEED = (0.5, 0.5)


def _mixture(z: Sequence[float] = FEED) -> ct.Mixture:
    return ct.Mixture.from_database(list(NAMES), list(z), normalize=True)


def _kappa(eos: PCSAFTEOS, pressure_Pa: float, composition: Sequence[float]) -> float:
    """``kappa = P / (rho dP/drho)`` at the liquid-like root, by finite difference.

    The same dimensionless isothermal-compressibility ratio ADR-0017 names
    phases by, recomputed here from the public ``pressure_Pa`` so the printed
    number is independent of ``phase_identity``'s own analytic derivative.
    """
    density = eos.density_roots(
        mixture=_mixture(),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        composition=list(composition),
    )[-1]
    bound = PCSAFTEOS(components=NAMES)
    step = 1e-4 * density
    slope = (
        bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density + step,
            composition=list(composition),
        )
        - bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density - step,
            composition=list(composition),
        )
    ) / (2.0 * step)
    return pressure_Pa / (density * slope)


def _print_phases(result: ct.FlashResult, eos: PCSAFTEOS, pressure_Pa: float) -> None:
    mixture = _mixture()
    for name, phase in result.phases.items():
        fractions = list(phase.composition.fractions)
        density = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions,
        )[-1]
        identity = eos.phase_identity(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions,
            phase="liquid",
        )
        amounts = ", ".join(
            f"{component} {value:.8f}" for component, value in zip(NAMES, fractions)
        )
        print(
            f"    {name:<8} {amounts}\n"
            f"             fraction {result.phase_fractions[name]:.8f}, "
            f"rho = {density:>12,.2f} mol/m^3\n"
            f"             identity {identity!r}, kappa = "
            f"{_kappa(eos, pressure_Pa, fractions):.3e} "
            f"(< {KAPPA_LIQUID_THRESHOLD} is a liquid)"
        )


def _residuals(result: ct.FlashResult) -> None:
    diagnostics = result.diagnostics
    print(
        f"    phase_regime = {diagnostics['phase_regime']!r}, "
        f"vapor_fraction = {result.vapor_fraction!r}, "
        f"phase_label_method = {diagnostics['phase_label_method']!r}"
    )
    print(
        f"    dG_split/RT = {float(diagnostics['delta_g_split_rt']):.6e} "
        "(negative: the split is an answer,\n"
        "    not just a fixed point); fugacity residual "
        f"{float(diagnostics['fugacity_residual']):.2e}; mass balance "
        f"{float(diagnostics['mass_balance_residual']):.2e};\n"
        f"    post-split stability {diagnostics['post_split_status']!r}."
    )


def the_state_and_the_old_failure() -> None:
    print("1) Water / n-hexane at 298.15 K and 1 atm: the state, and the old failure")
    eos = PCSAFTEOS()
    mixture = _mixture()
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
    )
    print(
        f"    stability_tp -> {stability.status!r}, tpd_min = {stability.tpd_min:.6e}, "
        f"feed branch {stability.feed_branch!r}, incipient branch {stability.phase_branch!r}"
    )
    roots = eos.density_roots(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        composition=list(FEED),
    )
    print(
        f"    the isotherm has {len(roots)} density roots at the feed: "
        + " and ".join(f"{value:,.2f}" for value in roots)
        + " mol/m^3,"
    )
    print(
        "    so a vapour root exists here and the pre-ADR-0019 split - which pinned one\n"
        "    phase to the liquid root and the other to the vapour root - had one to land\n"
        "    on. Both branches the stability test reports are 'liquid', which is the\n"
        "    information the split now uses."
    )


def the_liquid_liquid_split() -> None:
    print("\n2) The split at 1 atm")
    eos = PCSAFTEOS()
    result = ct.flash_tp(
        _mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
    )
    _print_phases(result, eos, ATMOSPHERE_PA)
    _residuals(result)
    print(
        f"    diagnostics['phase_i_branch'] = {result.diagnostics['phase_i_branch']!r}, "
        f"['phase_ii_branch'] = {result.diagnostics['phase_ii_branch']!r}\n"
        "    (those two keys are absent from a vapour-liquid result, where the branch pair\n"
        "    is the historical ('liquid', 'vapor'))."
    )


def three_feeds_one_tie_line() -> None:
    print("\n3) Three feeds, one tie line (the lever rule)")
    eos = PCSAFTEOS()
    for feed in (FEED, (0.2, 0.8), (0.8, 0.2)):
        result = ct.flash_tp(
            _mixture(feed), temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
        )
        water_rich = result.phases["liquid1"].composition.fractions
        hexane_rich = result.phases["liquid2"].composition.fractions
        beta = result.phase_fractions["liquid2"]
        total = sum(feed)
        recombined = [
            (1.0 - beta) * water_rich[index] + beta * hexane_rich[index] for index in range(2)
        ]
        error = max(abs(value / total - found) for value, found in zip(feed, recombined))
        print(
            f"    z = ({feed[0]:.1f}, {feed[1]:.1f}): liquid1 fraction "
            f"{result.phase_fractions['liquid1']:.8f}, liquid2 fraction {beta:.8f}, "
            f"x_hexane(liquid1) = {water_rich[1]:.8e}, x_water(liquid2) = "
            f"{hexane_rich[0]:.8e}, |z - recombined| = {error:.1e}"
        )
    print(
        "    The two compositions do not move with the feed, only the amounts do - which is\n"
        "    what makes them a tie line. Because ADR-0019 orders the two names by the first\n"
        "    component's mole fraction, 'liquid1' is the water-rich phase at every feed."
    )


def the_same_tie_line_at_one_megapascal() -> None:
    print("\n4) The same tie line at 1 MPa, where the vapour root is gone")
    eos = PCSAFTEOS()
    roots = eos.density_roots(
        mixture=_mixture(),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=HIGH_PRESSURE_PA,
        composition=list(FEED),
    )
    print(f"    density roots at the feed: {len(roots)} -> {roots[0]:,.2f} mol/m^3")
    result = ct.flash_tp(
        _mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=HIGH_PRESSURE_PA, eos=eos
    )
    _print_phases(result, eos, HIGH_PRESSURE_PA)
    _residuals(result)
    print(
        "    Before ADR-0019 this pair came back named 'liquid' / 'vapor' with\n"
        "    vapor_fraction = 0.503164 - which was really the hexane-rich LIQUID's\n"
        "    fraction - because both phases measured as liquids and the two-phase naming\n"
        "    rule fell through to its Wilson-ranking fallback. A liquid tie line barely\n"
        "    moves between 1 atm and 1 MPa, and this one does not: compare section 2."
    )


def main() -> None:
    print("Liquid-liquid equilibrium from PC-SAFT (ADR-0019, validation Case P-8)")
    print("=" * 78)
    the_state_and_the_old_failure()
    the_liquid_liquid_split()
    three_feeds_one_tie_line()
    the_same_tie_line_at_one_megapascal()
    print("\n" + "=" * 78)
    print(
        "k_ij = 0 throughout. PC-SAFT without a fitted binary parameter is known to be\n"
        "poor for water / hydrocarbon mutual solubilities: this model gives about\n"
        "1.7e-05 mole fraction hexane in water and 6.3e-03 water in hexane, against\n"
        "commonly tabulated experimental figures near 2e-06 and 5e-04 at 298 K (which\n"
        "were NOT verified against a primary source here, and which nothing above\n"
        "asserts). This demo shows the solver works, not that the model is accurate;\n"
        "examples/validation/17_pcsaft_lle_vs_feos.py checks the code against FeOs."
    )


if __name__ == "__main__":
    main()
