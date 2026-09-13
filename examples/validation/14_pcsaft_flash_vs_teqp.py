"""PC-SAFT flash against teqp's own 300 K isotherm, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0014 validated chemthermo's PC-SAFT *properties* against ``teqp`` at fixed
``(T, rho, x)``. This script validates the *equilibrium*: the density roots
(ADR-0015), the tangent-plane stability verdict and the phi-phi flash, for the
methane / n-hexane binary at 300 K.

The two sides share the published model and its 42 universal constants and
nothing else:

  * teqp obtains every derivative by automatic differentiation of one
    hand-written ``alphar``, and gets its tie lines by numerically continuing
    along the isotherm (``trace_VLE_isotherm_binary``) and then polishing with
    a Newton solve on its own residual (``mix_VLE_Tp``, ``mix_VLE_Tx``);
  * chemthermo writes analytic derivatives, solves ``P_model(T, rho, x) = P``
    for the mechanically stable density roots, and gets its tie lines from
    Michelsen's tangent-plane test followed by a Rachford-Rice /
    successive-substitution split.

What is checked here
--------------------
1. Pure n-hexane at 300 K and 400 K: ``density_roots`` at teqp's own saturation
   pressure returns teqp's two saturation densities, and reports that a third
   root was bracketed and discarded as mechanically unstable.
2. Seven pressures between 0.5 and 8.5 MPa: the flash's ``x1`` and ``y1``
   against teqp's polished tie line, the phase densities against teqp's, and
   the three residuals every split must carry.
3. **The decisive check**: teqp's own ``get_fugacity_coefficients``, evaluated
   at chemthermo's converged compositions and densities, must give the two
   phases equal fugacities. This uses no reference tie line at all - it asks
   teqp directly whether chemthermo's answer is an equilibrium.
4. Bubble pressures at three liquid compositions, found by **bisecting
   ``stability_tp``'s verdict** (unstable below, stable above), against teqp's
   ``mix_VLE_Tx``. chemthermo solves no bubble-point equation here.
5. Methane / n-decane at 350 K with ``kij = 0.03``. That value is
   **illustrative**, chosen to be nonzero; it is *not* a literature-validated
   binary parameter for this pair. What is validated is that the two codes
   agree once both are given it.

Requires the optional ``teqp`` dependency::

    pip install -e ".[validation]"

The script prints a message and exits 0 when it is missing.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

try:
    import teqp
except ImportError:  # pragma: no cover - exercised only without the extra
    teqp = None  # type: ignore[assignment]

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos._pcsaft_density import solve_density_roots
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

# Gross & Sadowski (2001) Table 1; written out here so the reference model is
# built from this file rather than from the package under test.
PARAMETERS: dict[str, tuple[float, float, float]] = {
    "Methane": (1.0000, 3.7039, 150.03),
    "n-Hexane": (3.0576, 3.7983, 236.77),
    "n-Decane": (4.6627, 3.8384, 243.87),
}

TEMPERATURE_K = 300.0
BINARY = ("Methane", "n-Hexane")
PRESSURES = (5.0e5, 1.0e6, 2.0e6, 3.0e6, 5.0e6, 7.0e6, 8.5e6)
BUBBLE_COMPOSITIONS = (0.10, 0.20, 0.30)

TIE_LINE_TOL = 1e-6
DENSITY_RTOL = 1e-5
FUGACITY_RTOL = 1e-8

failures: list[str] = []


def record(label: str, passed: bool) -> None:
    status = "PASS" if passed else "FAIL"
    print(f"    [{status}] {label}")
    if not passed:
        failures.append(label)


def build_reference(components: Sequence[str], kij_matrix: np.ndarray):
    coefficients = [
        {
            "name": name,
            "m": PARAMETERS[name][0],
            "sigma_Angstrom": PARAMETERS[name][1],
            "epsilon_over_k": PARAMETERS[name][2],
            "BibTeXKey": "Gross-IECR-2001",
        }
        for name in components
    ]
    return teqp.make_model(
        {"kind": "PCSAFT", "model": {"coeffs": coefficients, "kmat": kij_matrix.tolist()}}
    )


def teqp_fugacity_mismatch(
    model, temperature: float, pressure: float, phases: Sequence[tuple[np.ndarray, float]]
) -> float:
    """Worst relative difference between two phases' fugacities *in teqp's model*."""
    first, second = (
        np.asarray(model.get_fugacity_coefficients(temperature, density * x)) * x * pressure
        for x, density in phases
    )
    return float(np.max(np.abs(first / second - 1.0)))


def bisect_verdict(unstable: Callable[[float], bool], low: float, high: float, steps: int) -> float:
    for _ in range(steps):
        mid = 0.5 * (low + high)
        if unstable(mid):
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def check_pure_saturation_roots(eos: PCSAFTEOS) -> None:
    print("\n1) Pure n-hexane density roots at teqp's saturation pressure")
    print("-" * 78)
    model = build_reference(("n-Hexane",), np.zeros((1, 1)))
    for temperature, guesses in ((300.0, (7700.0, 10.0)), (400.0, (6800.0, 150.0))):
        rho_l_ref, rho_v_ref = model.pure_VLE_T(temperature, guesses[0], guesses[1], 200)
        p_ref = float(
            rho_l_ref
            * R_J_PER_MOL_K
            * temperature
            * (1.0 + model.get_Ar01(temperature, rho_l_ref, np.array([1.0])))
        )
        isotherm = eos._isotherm(names=("n-Hexane",), temperature_K=temperature, composition=[1.0])
        detail = solve_density_roots(isotherm, p_ref)
        rho_v, rho_l = detail.densities[0], detail.densities[-1]
        print(f"  T = {temperature:.0f} K, Psat = {p_ref:,.6f} Pa (teqp)")
        print(
            f"    rho_V : {rho_v:>16.8f}  teqp {float(rho_v_ref):>16.8f}  "
            f"rel {abs(rho_v / rho_v_ref - 1.0):.2e}"
        )
        print(
            f"    rho_L : {rho_l:>16.8f}  teqp {float(rho_l_ref):>16.8f}  "
            f"rel {abs(rho_l / rho_l_ref - 1.0):.2e}"
        )
        print(
            f"    brackets found {detail.bracket_count}, roots returned "
            f"{len(detail.densities)} (the spinodal branch is discarded)"
        )
        record(
            f"pure n-hexane saturation densities at {temperature:.0f} K to 1e-8 relative",
            len(detail.densities) == 2
            and detail.bracket_count == 3
            and abs(rho_v / rho_v_ref - 1.0) < 1e-8
            and abs(rho_l / rho_l_ref - 1.0) < 1e-8,
        )


def check_tie_lines(eos: PCSAFTEOS, model, trace, pressures: np.ndarray) -> None:
    print("\n2) Tie lines at 300 K: flash_tp against teqp's traced + polished isotherm")
    print("-" * 78)
    print(
        f"  {'P / MPa':>8} {'x1 (teqp)':>11} {'dx1':>10} {'y1 (teqp)':>11} {'dy1':>10} "
        f"{'d rhoL':>9} {'d rhoV':>9} {'teqp f':>9} {'dG/RT':>8}"
    )
    rising = int(np.argmax(pressures))
    for pressure in PRESSURES:
        index = int(np.argmin(np.abs(pressures[: rising + 1] - pressure)))
        solved = model.mix_VLE_Tp(
            TEMPERATURE_K,
            pressure,
            np.array(trace[index]["rhoL / mol/m^3"]),
            np.array(trace[index]["rhoV / mol/m^3"]),
        )
        rho_liquid_ref = np.asarray(solved.rhovecL, dtype=float)
        rho_vapor_ref = np.asarray(solved.rhovecV, dtype=float)
        x1_ref = float(rho_liquid_ref[0] / rho_liquid_ref.sum())
        y1_ref = float(rho_vapor_ref[0] / rho_vapor_ref.sum())

        z1 = 0.5 * (x1_ref + y1_ref)
        mixture = ct.Mixture.from_database(list(BINARY), [z1, 1.0 - z1])
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=eos)
        x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
        y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
        rho_liquid = eos.density_roots(
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure,
            composition=x.tolist(),
            mixture=mixture,
        )[-1]
        rho_vapor = eos.density_roots(
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure,
            composition=y.tolist(),
            mixture=mixture,
        )[0]
        mismatch = teqp_fugacity_mismatch(
            model, TEMPERATURE_K, pressure, ((x, rho_liquid), (y, rho_vapor))
        )
        delta_g = float(result.diagnostics["delta_g_split_rt"])
        print(
            f"  {pressure / 1e6:>8.2f} {x1_ref:>11.6f} {float(x[0]) - x1_ref:>+10.2e} "
            f"{y1_ref:>11.6f} {float(y[0]) - y1_ref:>+10.2e} "
            f"{rho_liquid / rho_liquid_ref.sum() - 1.0:>+9.1e} "
            f"{rho_vapor / rho_vapor_ref.sum() - 1.0:>+9.1e} {mismatch:>9.1e} {delta_g:>+8.4f}"
        )
        label = f"P = {pressure / 1e6:.2f} MPa"
        record(
            f"{label}: tie line within {TIE_LINE_TOL:g} mole fraction of teqp",
            abs(float(x[0]) - x1_ref) <= TIE_LINE_TOL and abs(float(y[0]) - y1_ref) <= TIE_LINE_TOL,
        )
        record(
            f"{label}: phase densities within {DENSITY_RTOL:g} relative of teqp",
            abs(rho_liquid / rho_liquid_ref.sum() - 1.0) <= DENSITY_RTOL
            and abs(rho_vapor / rho_vapor_ref.sum() - 1.0) <= DENSITY_RTOL,
        )
        record(
            f"{label}: teqp's own fugacity coefficients make the split an equilibrium",
            mismatch < FUGACITY_RTOL,
        )
        record(
            f"{label}: split verified (mass balance, equal fugacity, dG < 0, post-split stable)",
            float(result.diagnostics["mass_balance_residual"]) < 1e-12
            and float(result.diagnostics["fugacity_residual"]) < 1e-8
            and delta_g < 0.0
            and result.diagnostics["post_split_stable"] is True,
        )


def check_bubble_pressures(eos: PCSAFTEOS, model, trace, pressures: np.ndarray) -> None:
    print("\n3) Bubble pressures from the stability verdict alone, against teqp mix_VLE_Tx")
    print("-" * 78)
    rising = int(np.argmax(pressures))
    liquid_fractions = np.array(
        [point["xL_0 / mole frac."] for point in trace[: rising + 1]], dtype=float
    )
    for x1 in BUBBLE_COMPOSITIONS:
        seed = int(np.argmin(np.abs(liquid_fractions - x1)))
        _code, rho_liquid_ref, _rho_vapor_ref = model.mix_VLE_Tx(
            TEMPERATURE_K,
            np.array(trace[seed]["rhoL / mol/m^3"]),
            np.array(trace[seed]["rhoV / mol/m^3"]),
            np.array([x1, 1.0 - x1]),
            1e-10,
            1e-10,
            1e-10,
            1e-10,
            100,
        )
        rho_liquid_ref = np.asarray(rho_liquid_ref, dtype=float)
        total = float(rho_liquid_ref.sum())
        p_ref = float(
            total
            * R_J_PER_MOL_K
            * TEMPERATURE_K
            * (1.0 + model.get_Ar01(TEMPERATURE_K, total, rho_liquid_ref / total))
        )

        mixture = ct.Mixture.from_database(list(BINARY), [x1, 1.0 - x1])

        def unstable(pressure: float) -> bool:
            return (
                ct.stability_tp(
                    mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=eos
                ).status
                == "unstable"
            )

        p_bubble = bisect_verdict(unstable, 0.5 * p_ref, 2.0 * p_ref, 40)
        relative = abs(p_bubble / p_ref - 1.0)
        print(
            f"  x1 = {x1:.2f}: verdict flip {p_bubble:>16,.6f} Pa   "
            f"teqp {p_ref:>16,.6f} Pa   rel {relative:.2e}"
        )
        record(f"bubble pressure at x1 = {x1:.2f} within 1e-6 relative of teqp", relative < 1e-6)


def check_kij_binary() -> None:
    print("\n4) Methane / n-decane at 350 K with an illustrative kij = 0.03")
    print("-" * 78)
    components = ("Methane", "n-Decane")
    kij = 0.03
    eos = PCSAFTEOS(kij={components: kij})
    model = build_reference(components, np.array([[0.0, kij], [kij, 0.0]]))
    temperature, pressure = 350.0, 5.0e6

    mixture = ct.Mixture.from_database(list(components), [0.4, 0.6])
    result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
    rho_liquid = eos.density_roots(
        temperature_K=temperature,
        pressure_Pa=pressure,
        composition=x.tolist(),
        mixture=mixture,
    )[-1]
    rho_vapor = eos.density_roots(
        temperature_K=temperature,
        pressure_Pa=pressure,
        composition=y.tolist(),
        mixture=mixture,
    )[0]
    mismatch = teqp_fugacity_mismatch(
        model, temperature, pressure, ((x, rho_liquid), (y, rho_vapor))
    )
    print(f"  vapor fraction {result.vapor_fraction: .8f}")
    print(f"  liquid x(methane) = {float(x[0]):.8f}   rho = {rho_liquid:,.4f} mol/m^3")
    print(f"  vapor  y(methane) = {float(y[0]):.8f}   rho = {rho_vapor:,.4f} mol/m^3")
    print(f"  teqp equal-fugacity mismatch: {mismatch:.2e} relative")
    record("kij = 0.03 split is an equilibrium in teqp's model", mismatch < FUGACITY_RTOL)
    print(
        "  kij = 0.03 is illustrative, not a literature-validated parameter for this\n"
        "  pair. The check is that two independent codes agree once both are given it."
    )


def main() -> None:
    if teqp is None:
        print("teqp is not installed; skipping. Install with: pip install -e '.[validation]'")
        return

    print("=" * 78)
    print("PC-SAFT equilibrium against teqp (validation Cases P-3 and P-5)")
    print("=" * 78)

    eos = PCSAFTEOS()
    check_pure_saturation_roots(eos)

    model = build_reference(BINARY, np.zeros((2, 2)))
    hexane = build_reference(("n-Hexane",), np.zeros((1, 1)))
    rho_liquid, rho_vapor = hexane.pure_VLE_T(TEMPERATURE_K, 7700.0, 10.0, 200)
    options = teqp.TVLEOptions()
    options.polish = True
    options.max_steps = 10000
    options.integration_order = 5
    trace = model.trace_VLE_isotherm_binary(
        TEMPERATURE_K,
        np.array([0.0, rho_liquid]),
        np.array([0.0, rho_vapor]),
        options,
    )
    pressures = np.array([point["pL / Pa"] for point in trace])

    check_tie_lines(eos, model, trace, pressures)
    check_bubble_pressures(eos, model, trace, pressures)
    check_kij_binary()

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
