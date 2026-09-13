"""Three phases from an equation of state, discovered rather than assumed.

What this shows
---------------
ADR-0011 gave ``flash_tp`` a phase addition / removal search - solve a phase
set, add the incipient phase a failed post-split stability test found, remove
a phase whose amount converges to zero or below - and wired it to the
low-pressure ``modified-raoult`` path only. ADR-0020 wires it to the
equation-of-state (phi-phi) path as well, with each phase pinned to its own
density root (ADR-0019), so a vapour and two liquids can coexist in one
answer.

Nothing here asks for three phases. ``flash_tp`` is called exactly as it always
is::

    flash_tp(mixture, temperature_K=..., pressure_Pa=101325.0, eos=PCSAFTEOS())

and the phase count comes back in ``result.phases``.

The three sections
------------------
1. **Three liquid phases from Peng-Robinson.** Water / ethanol / n-hexane at
   280 K and 1 atm with ``k_ij = 0``: ``liquid1`` / ``liquid2`` / ``liquid3``,
   found by the search in about a tenth of a second. ADR-0019 recorded "no
   pure Peng-Robinson three-phase case found in the databank"; this is one.
2. **The water / n-hexane window at 1 atm.** Just *below* the three-phase
   temperature ``T3`` the deepest tangent-plane minimum from a 50/50 feed is a
   **vapour**, so the search starts ``V -> LV``; that pair is unstable towards
   a second liquid, ``LV -> LLV``; and the three-phase solve then drives the
   vapour amount negative, ``LLV -> LL``. The two conjugate liquids are
   reached by *removing* a phase that addition had to add first. Before this
   slice those temperatures raised ``ConvergenceError``. Just *above* ``T3``
   the ordinary vapour-liquid answer is returned and the search is never
   entered.
3. **A vapour-liquid-liquid tie triangle** (``--full``): water / ethanol /
   n-hexane with PC-SAFT at 333 K and 1 atm returns ``liquid1`` / ``liquid2``
   / ``vapor`` with ``phase_regime = "VLLE"``. It is behind the flag only
   because one such flash costs about half a minute.

A binary at a fixed pressure has no three-phase *region*: Gibbs' phase rule
gives one degree of freedom, so three phases meet at a single temperature, and
there the three phase amounts are not fixed by the mass balance either. That is
why section 2 shows the two-phase answers on either side of ``T3`` rather than
a three-phase answer at it, and why section 3 uses a ternary.

Run::

    python examples/basic/flash_tp_pcsaft_vlle_demo.py
    python examples/basic/flash_tp_pcsaft_vlle_demo.py --full

Model caveat, printed and never asserted: ``k_ij = 0`` between water and a
hydrocarbon is not a serious parameterization, so the numbers below are this
model's, not measurements. ``examples/validation/18_pcsaft_vlle_water_hexane.py``
checks the code against FeOs.
"""

from __future__ import annotations

import argparse

import numpy as np

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

ATMOSPHERE_PA = 101325.0
BINARY = ("Water", "n-Hexane")
TERNARY = ("Water", "Ethanol", "n-Hexane")
PR_T_K = 280.0
PCSAFT_TERNARY_T_K = 333.0
OFFSET_K = 0.05


def _print_phases(result: ct.FlashResult) -> None:
    for name in sorted(result.phases):
        fractions = result.phases[name].composition.fractions
        amounts = "  ".join(f"{value:>12.8f}" for value in fractions)
        print(f"    {name:<8} beta = {result.phase_fractions[name]:.8f}   x = {amounts}")


def _print_search(result: ct.FlashResult) -> None:
    diagnostics = result.diagnostics
    history = diagnostics.get("phase_set_history", "(search not entered)")
    print(f"    phase_regime          {diagnostics['phase_regime']}")
    print(f"    phase_set_history     {history}")
    if "phases_added" in diagnostics:
        print(
            f"    phases added/removed  {diagnostics['phases_added']} / "
            f"{diagnostics['phases_removed']}"
        )
    residual_key = (
        "equilibrium_residual" if "equilibrium_residual" in diagnostics else "fugacity_residual"
    )
    print(f"    {residual_key:<21} {float(diagnostics[residual_key]):.3e}")
    print(f"    mass_balance_residual {float(diagnostics['mass_balance_residual']):.3e}")
    print(f"    delta_g_split_rt      {float(diagnostics['delta_g_split_rt']):+.6f}")
    if "delta_g_vs_two_phase_rt" in diagnostics:
        print(
            "    delta_g_vs_two_phase  "
            f"{float(diagnostics['delta_g_vs_two_phase_rt']):+.6e}  (< 0: better than the "
            "two-phase answer the search started from)"
        )
    print(f"    post_split_status     {diagnostics['post_split_status']}")
    print(f"    vapor_fraction        {result.vapor_fraction}")


def three_liquids_from_a_cubic() -> None:
    print("1) Three liquid phases from Peng-Robinson (k_ij = 0)")
    print("-" * 78)
    print(f"   water / ethanol / n-hexane, z = (0.2, 0.4, 0.4), {PR_T_K} K, 1 atm")
    result = ct.flash_tp(
        ct.Mixture.from_database(list(TERNARY), [0.2, 0.4, 0.4], normalize=True),
        temperature_K=PR_T_K,
        pressure_Pa=ATMOSPHERE_PA,
        eos=ct.PengRobinsonEOS(),
    )
    print(f"   phases: {sorted(result.phases)}")
    _print_phases(result)
    _print_search(result)
    print(
        "\n   liquid1 / liquid2 / liquid3 are ordered by the first component's mole\n"
        "   fraction (ADR-0019), so the names mean the same thing at every feed inside\n"
        "   the triangle - unlike the gamma-gamma path, where they are seed roles."
    )


def _three_phase_temperature() -> float:
    """The binary three-phase temperature from a 4-equation Newton.

    Two liquids on the model's liquid branch and a vapour on its vapour branch,
    all at 1 atm: four equal-fugacity equations in ``(x^I, x^II, y, T)``. The
    temperature is an unknown because Gibbs' phase rule makes it one.
    """
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True)

    def ln_f(x1: float, temperature: float, branch: str) -> np.ndarray:
        x = np.array([x1, 1.0 - x1])
        phi = np.asarray(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=ATMOSPHERE_PA,
                composition=x.tolist(),
                phase=branch,
            ),
            dtype=float,
        )
        return np.log(x) + np.log(phi)

    def residual(u: np.ndarray) -> np.ndarray:
        first = ln_f(u[0], u[3], "liquid")
        second = ln_f(u[1], u[3], "liquid")
        vapor = ln_f(u[2], u[3], "vapor")
        return np.concatenate([first - second, first - vapor])

    u = np.array([0.9999, 0.02, 0.20, 334.5])
    steps = (1e-8, 1e-8, 1e-8, 3.345e-5)
    for _iteration in range(40):
        f = residual(u)
        worst = float(np.max(np.abs(f)))
        if worst < 1e-11:
            break
        jacobian = np.zeros((4, 4))
        for column in range(4):
            shifted = u.copy()
            shifted[column] += steps[column]
            jacobian[:, column] = (residual(shifted) - f) / steps[column]
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while scale > 1e-10:
            candidate = u + scale * direction
            inside = bool(np.all(candidate[:3] > 0.0) and np.all(candidate[:3] < 1.0))
            if inside and float(np.max(np.abs(residual(candidate)))) < worst:
                break
            scale *= 0.5
        u = u + scale * direction

    print(f"   three-phase point (4-equation Newton, residual {worst:.2e}):")
    print(f"     T3          = {u[3]:.9f} K  ({u[3] - 273.15:.4f} C)")
    print(f"     x_water(I)  = {u[0]:.12f}   (water-rich liquid)")
    print(f"     x_water(II) = {u[1]:.12f}   (hexane-rich liquid)")
    print(f"     y_water     = {u[2]:.12f}   (vapour)")
    for label, value in (("I", u[0]), ("II", u[1]), ("V", u[2])):
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=float(u[3]),
            pressure_Pa=ATMOSPHERE_PA,
            composition=[float(value), 1.0 - float(value)],
        )
        pretty = ", ".join(f"{root:,.2f}" for root in roots)
        print(f"     density roots of {label:<2}: {pretty} mol/m^3")
    print(
        "\n   Every one of the three compositions has BOTH a vapour and a liquid root,\n"
        "   which is why each phase has to be pinned to its own (ADR-0019): one\n"
        "   liquid/vapour branch assignment cannot describe two liquids and a vapour."
    )
    return float(u[3])


def the_binary_window() -> float:
    print("\n2) The water / n-hexane window at 1 atm, z = 0.5 / 0.5")
    print("-" * 78)
    t3 = _three_phase_temperature()

    for offset, comment in (
        (-OFFSET_K, "below T3: two liquids, reached by adding a phase and removing one"),
        (+OFFSET_K, "above T3: the ordinary vapour-liquid answer; no search"),
    ):
        temperature = t3 + offset
        print(f"\n   T = T3 {offset:+.2f} K = {temperature:.6f} K  ({comment})")
        result = ct.flash_tp(
            ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
            temperature_K=temperature,
            pressure_Pa=ATMOSPHERE_PA,
            eos=PCSAFTEOS(),
        )
        print(f"   phases: {sorted(result.phases)}")
        _print_phases(result)
        _print_search(result)

    print(
        "\n   Before this slice the first of those two raised ConvergenceError: the\n"
        "   two-phase set was provably not the answer and there was no way to say so."
    )
    return t3


def the_ternary_tie_triangle() -> None:
    print("\n3) A vapour-liquid-liquid tie triangle from PC-SAFT")
    print("-" * 78)
    print(f"   water / ethanol / n-hexane, {PCSAFT_TERNARY_T_K} K, 1 atm (2B water and ethanol)")
    vertices = None
    for feed in ((0.4, 0.3, 0.3), (0.5, 0.2, 0.3)):
        result = ct.flash_tp(
            ct.Mixture.from_database(list(TERNARY), list(feed), normalize=True),
            temperature_K=PCSAFT_TERNARY_T_K,
            pressure_Pa=ATMOSPHERE_PA,
            eos=PCSAFTEOS(),
        )
        print(f"\n   z = {feed} -> {sorted(result.phases)}")
        _print_phases(result)
        _print_search(result)
        current = {
            name: tuple(result.phases[name].composition.fractions) for name in sorted(result.phases)
        }
        if vertices is None:
            vertices = current
        else:
            moved = max(
                abs(a - b) for name in current for a, b in zip(current[name], vertices[name])
            )
            print(f"\n   same triangle as the first feed to {moved:.2e} in mole fraction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the PC-SAFT ternary vapour-liquid-liquid tie triangle (~1 min)",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("Three phases from an equation of state (ADR-0020, Cases P-9 and P-10)")
    print("=" * 78)
    three_liquids_from_a_cubic()
    the_binary_window()
    if args.full:
        the_ternary_tie_triangle()
    else:
        print("\n  (pass --full for the PC-SAFT vapour-liquid-liquid tie triangle)")

    print("\n" + "=" * 78)
    print(
        "k_ij = 0 throughout, and nothing above is compared against measurement.\n"
        "The phase count is an output of the tangent-plane stability test, so it is\n"
        "never better than that test: a three-phase region the deterministic trial\n"
        "set misses comes back as two phases with no error. See ADR-0020 and\n"
        "validation Cases P-9 and P-10 for what is and is not covered."
    )


if __name__ == "__main__":
    main()
