"""PC-SAFT density roots, stability and flash against teqp (Cases P-3 and P-5).

What makes this an independent check
------------------------------------
``teqp`` (NIST, MIT licence) implements the same published model - Gross &
Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244, non-associating - but shares no
code with chemthermo: teqp obtains every derivative by automatic
differentiation of one hand-written ``alphar``, chemthermo writes them
analytically. On top of that the *equilibrium* machinery is entirely separate:

- teqp gets its tie line from ``trace_VLE_isotherm_binary`` (a numerical
  continuation along the isotherm) polished by ``mix_VLE_Tp`` (a Newton solve
  on its own equal-chemical-potential residual);
- chemthermo gets its tie line from Michelsen's tangent-plane stability test
  followed by a Rachford-Rice / successive-substitution split, with the
  densities coming from the root solver of ADR-0015.

So an agreement here is an agreement between two different equilibrium
formulations of one model, not a restatement of one of them.

The decisive check is the last one in each case: teqp's **own**
``get_fugacity_coefficients`` is evaluated at chemthermo's converged phase
compositions and densities, and the two phases must have equal fugacities *in
teqp's model*. That needs no reference tie line at all.

Skipped when ``teqp`` is not installed (``pip install -e ".[validation]"``).
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos.pcsaft import PCSAFTEOS, R_J_PER_MOL_K

teqp = pytest.importorskip("teqp")

#: Gross & Sadowski (2001) Table 1, repeated here so the reference model is
#: built from values written in this file, not read out of the code under test.
_PARAMETERS: dict[str, tuple[float, float, float]] = {
    "Methane": (1.0000, 3.7039, 150.03),
    "n-Hexane": (3.0576, 3.7983, 236.77),
    "n-Decane": (4.6627, 3.8384, 243.87),
}

TEMPERATURE_K = 300.0
BINARY = ("Methane", "n-Hexane")

#: Tie-line agreement required of ``x1`` and ``y1`` (mole fraction).
_TIE_LINE_TOL = 1e-6
#: Relative agreement required of the phase densities.
_DENSITY_RTOL = 1e-5
#: Relative agreement required of teqp's own equal-fugacity check.
_FUGACITY_RTOL = 1e-8


def _reference_model(components: Sequence[str], kij_matrix: np.ndarray):
    coefficients = [
        {
            "name": name,
            "m": _PARAMETERS[name][0],
            "sigma_Angstrom": _PARAMETERS[name][1],
            "epsilon_over_k": _PARAMETERS[name][2],
            "BibTeXKey": "Gross-IECR-2001",
        }
        for name in components
    ]
    return teqp.make_model(
        {
            "kind": "PCSAFT",
            "model": {"coeffs": coefficients, "kmat": kij_matrix.tolist()},
        }
    )


@pytest.fixture(scope="module")
def reference():
    """teqp's methane / n-hexane model plus its traced 300 K isotherm."""
    model = _reference_model(BINARY, np.zeros((2, 2)))
    hexane = _reference_model(("n-Hexane",), np.zeros((1, 1)))
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
    return model, trace, pressures, float(rho_liquid), float(rho_vapor)


def _teqp_tie_line(reference, pressure: float) -> tuple[np.ndarray, np.ndarray]:
    """teqp's own tie line at ``(300 K, pressure)``, from the trace then polished.

    The seed is the traced point nearest in pressure on the rising branch (the
    trace continues past the mixture critical region and comes back down, so
    only the stretch up to the pressure maximum is single valued); teqp's
    ``mix_VLE_Tp`` then solves its own residual at exactly this pressure.
    """
    model, trace, pressures, _rho_l, _rho_v = reference
    rising = pressures[: int(np.argmax(pressures)) + 1]
    index = int(np.argmin(np.abs(rising - pressure)))
    solved = model.mix_VLE_Tp(
        TEMPERATURE_K,
        pressure,
        np.array(trace[index]["rhoL / mol/m^3"]),
        np.array(trace[index]["rhoV / mol/m^3"]),
    )
    return np.array(solved.rhovecL, dtype=float), np.array(solved.rhovecV, dtype=float)


def _teqp_fugacity_mismatch(
    model, temperature: float, pressure: float, phases: Sequence[tuple[np.ndarray, float]]
) -> float:
    """Worst relative difference between two phases' fugacities *in teqp's model*."""
    fugacities = [
        np.asarray(model.get_fugacity_coefficients(temperature, density * x)) * x * pressure
        for x, density in phases
    ]
    first, second = fugacities
    return float(np.max(np.abs(first / second - 1.0)))


def _bisect(predicate: Callable[[float], bool], low: float, high: float, iterations: int) -> float:
    """Return the ``P`` where ``predicate`` flips from True at ``low`` to False at ``high``."""
    assert predicate(low) and not predicate(high)
    for _ in range(iterations):
        mid = 0.5 * (low + high)
        if predicate(mid):
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


# ---------------------------------------------------------------------------
# Case P-3: the density roots themselves
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("temperature", [300.0, 400.0])
def test_pure_hexane_saturation_roots_match_teqp(temperature: float) -> None:
    """``density_roots`` at teqp's saturation pressure returns teqp's two densities.

    Case P-2 established the saturation state with a test-local bisection root
    finder; this restates it through the shipped solver, and additionally
    checks that the mechanically unstable middle root is bracketed and then
    discarded rather than never found.
    """
    eos = PCSAFTEOS(components=("n-Hexane",))
    model = _reference_model(("n-Hexane",), np.zeros((1, 1)))
    guesses = (7700.0, 10.0) if temperature < 350.0 else (6800.0, 150.0)
    rho_l_ref, rho_v_ref = model.pure_VLE_T(temperature, guesses[0], guesses[1], 200)
    p_ref = float(
        rho_l_ref
        * R_J_PER_MOL_K
        * temperature
        * (1.0 + model.get_Ar01(temperature, rho_l_ref, np.array([1.0])))
    )

    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=p_ref, composition=[1.0])
    assert len(roots) == 2
    assert roots[0] == pytest.approx(rho_v_ref, rel=1e-8)
    assert roots[1] == pytest.approx(rho_l_ref, rel=1e-8)

    isotherm = eos._isotherm(names=("n-Hexane",), temperature_K=temperature, composition=[1.0])
    from chemthermo.eos._pcsaft_density import solve_density_roots

    detail = solve_density_roots(isotherm, p_ref)
    assert detail.bracket_count == 3, "the spinodal root must be found and then rejected"
    assert len(detail.densities) == 2


def test_the_reference_is_not_vacuous() -> None:
    """Negative control: a 1 % change in one sigma must break the agreement."""
    eos = PCSAFTEOS(components=("n-Hexane",))
    perturbed = teqp.make_model(
        {
            "kind": "PCSAFT",
            "model": {
                "coeffs": [
                    {
                        "name": "n-Hexane",
                        "m": _PARAMETERS["n-Hexane"][0],
                        "sigma_Angstrom": _PARAMETERS["n-Hexane"][1] * 1.01,
                        "epsilon_over_k": _PARAMETERS["n-Hexane"][2],
                        "BibTeXKey": "perturbed",
                    }
                ],
                "kmat": [[0.0]],
            },
        }
    )
    rho_l, _rho_v = perturbed.pure_VLE_T(300.0, 7700.0, 10.0, 200)
    p_ref = float(
        rho_l * R_J_PER_MOL_K * 300.0 * (1.0 + perturbed.get_Ar01(300.0, rho_l, np.array([1.0])))
    )
    roots = eos.density_roots(temperature_K=300.0, pressure_Pa=p_ref, composition=[1.0])
    assert abs(roots[-1] / rho_l - 1.0) > 1e-3


# ---------------------------------------------------------------------------
# Case P-5: the tie lines
# ---------------------------------------------------------------------------

_PRESSURES = [5.0e5, 1.0e6, 2.0e6, 3.0e6, 5.0e6, 7.0e6, 8.5e6]


@pytest.mark.parametrize("pressure", _PRESSURES, ids=[f"{p:.3g} Pa" for p in _PRESSURES])
def test_flash_tie_line_matches_teqp_at_300_K(reference, pressure: float) -> None:
    model = reference[0]
    rho_liquid_ref, rho_vapor_ref = _teqp_tie_line(reference, pressure)
    x1_ref = float(rho_liquid_ref[0] / rho_liquid_ref.sum())
    y1_ref = float(rho_vapor_ref[0] / rho_vapor_ref.sum())
    assert 0.0 < x1_ref < y1_ref < 1.0

    # The reference itself must satisfy teqp's equal-fugacity condition, or the
    # comparison below would be against an arbitrary pair of states.
    assert (
        _teqp_fugacity_mismatch(
            model,
            TEMPERATURE_K,
            pressure,
            (
                (rho_liquid_ref / rho_liquid_ref.sum(), float(rho_liquid_ref.sum())),
                (rho_vapor_ref / rho_vapor_ref.sum(), float(rho_vapor_ref.sum())),
            ),
        )
        < _FUGACITY_RTOL
    )

    z1 = 0.5 * (x1_ref + y1_ref)
    mixture = ct.Mixture.from_database(list(BINARY), [z1, 1.0 - z1])
    eos = PCSAFTEOS()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=eos)

    assert set(result.phases) == {"liquid", "vapor"}
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["fugacity_residual"]) < 1e-8
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_stable"] is True

    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
    assert abs(float(x[0]) - x1_ref) <= _TIE_LINE_TOL
    assert abs(float(y[0]) - y1_ref) <= _TIE_LINE_TOL

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
    assert rho_liquid == pytest.approx(float(rho_liquid_ref.sum()), rel=_DENSITY_RTOL)
    assert rho_vapor == pytest.approx(float(rho_vapor_ref.sum()), rel=_DENSITY_RTOL)

    # The decisive check: chemthermo's two phases are in equilibrium according
    # to teqp's own fugacity coefficients, evaluated at chemthermo's densities.
    assert (
        _teqp_fugacity_mismatch(model, TEMPERATURE_K, pressure, ((x, rho_liquid), (y, rho_vapor)))
        < _FUGACITY_RTOL
    )


_BUBBLE_COMPOSITIONS = [0.10, 0.20, 0.30]


@pytest.mark.parametrize("x1", _BUBBLE_COMPOSITIONS, ids=[f"x1={v}" for v in _BUBBLE_COMPOSITIONS])
def test_bubble_pressure_from_stability_verdicts_matches_teqp(reference, x1: float) -> None:
    """Where ``stability_tp`` flips verdict at fixed feed is the bubble pressure.

    chemthermo never solves a bubble-point equation here: the pressure is found
    by bisecting the *verdict* of the tangent-plane test on a feed of
    composition ``x1``, which is unstable below the bubble pressure and stable
    above it. teqp solves its own ``mix_VLE_Tx`` at the same liquid
    composition.
    """
    model, trace, pressures, _rho_l, _rho_v = reference
    rising = pressures[: int(np.argmax(pressures)) + 1]
    liquid_fractions = np.array(
        [point["xL_0 / mole frac."] for point in trace[: rising.size]], dtype=float
    )
    seed = int(np.argmin(np.abs(liquid_fractions - x1)))
    code, rho_liquid_ref, rho_vapor_ref = model.mix_VLE_Tx(
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
    assert code != teqp.VLE_return_code.notfinite_step
    rho_liquid_ref = np.asarray(rho_liquid_ref, dtype=float)
    assert float(rho_liquid_ref[0] / rho_liquid_ref.sum()) == pytest.approx(x1, abs=1e-10)
    p_ref = float(
        rho_liquid_ref.sum()
        * R_J_PER_MOL_K
        * TEMPERATURE_K
        * (
            1.0
            + model.get_Ar01(
                TEMPERATURE_K, rho_liquid_ref.sum(), rho_liquid_ref / rho_liquid_ref.sum()
            )
        )
    )

    mixture = ct.Mixture.from_database(list(BINARY), [x1, 1.0 - x1])
    eos = PCSAFTEOS()

    def unstable(pressure: float) -> bool:
        return (
            ct.stability_tp(
                mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=eos
            ).status
            == "unstable"
        )

    # 40 bisections on a bracket of half a decade take the interval well below
    # 1e-6 relative.
    p_bubble = _bisect(unstable, 0.5 * p_ref, 2.0 * p_ref, 40)
    assert p_bubble == pytest.approx(p_ref, rel=1e-6)


def test_a_kij_binary_reaches_the_same_equilibrium_as_teqp() -> None:
    """Methane / n-decane at 350 K with ``kij = 0.03``.

    That value is **illustrative**: it is the same number used in the ADR-0014
    cross-check states, chosen to be nonzero, and is *not* a
    literature-validated binary parameter for this pair. What is being
    validated is that chemthermo and teqp agree once both are given it.
    """
    components = ("Methane", "n-Decane")
    kij = 0.03
    eos = PCSAFTEOS(kij={components: kij})
    model = _reference_model(components, np.array([[0.0, kij], [kij, 0.0]]))

    temperature, pressure = 350.0, 5.0e6
    mixture = ct.Mixture.from_database(list(components), [0.4, 0.6])
    result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)

    assert set(result.phases) == {"liquid", "vapor"}
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_stable"] is True

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
    assert (
        _teqp_fugacity_mismatch(model, temperature, pressure, ((x, rho_liquid), (y, rho_vapor)))
        < _FUGACITY_RTOL
    )

    # teqp's own solver, seeded from chemthermo's answer, must not move it.
    solved = model.mix_VLE_Tp(temperature, pressure, rho_liquid * x, rho_vapor * y)
    rho_liquid_ref = np.asarray(solved.rhovecL, dtype=float)
    rho_vapor_ref = np.asarray(solved.rhovecV, dtype=float)
    assert float(rho_liquid_ref[0] / rho_liquid_ref.sum()) == pytest.approx(
        float(x[0]), abs=_TIE_LINE_TOL
    )
    assert float(rho_vapor_ref[0] / rho_vapor_ref.sum()) == pytest.approx(
        float(y[0]), abs=_TIE_LINE_TOL
    )
