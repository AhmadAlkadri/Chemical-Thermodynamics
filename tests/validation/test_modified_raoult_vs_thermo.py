"""Modified-Raoult VLE against `thermo` 0.6.0 (validation Case R-4).

External route: `thermo.GibbsExcessLiquid` (an activity-coefficient liquid with
pure-component vapor pressures as the reference fugacity, `use_Poynting=False`
and `use_phis_sat=False`, which is exactly modified Raoult) over `thermo.NRTL`
with the same taus and alphas, an `thermo.IdealGas` vapor, and `thermo.FlashVL`.

The comparison is only meaningful if both packages evaluate the *same*
`Psat_i(T)`, so that is established first, to machine precision: `thermo`'s
`VaporPressure` is given the chemthermo databank's own Antoine coefficients,
converted from the base-e / bar form
(`ln(P/bar) = A - B/(T + C)`) to the base-e / Pa form `thermo` and
`chemicals.vapor_pressure.Antoine` expect by adding `ln(1e5)` to `A`.

System: 1-propanol(1) / water(2) at 101325 Pa, NRTL pair 1-3 of Tessier,
Brennecke and Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1. Those
parameters were fitted to *liquid-liquid* data and are temperature independent;
this module compares two implementations of one model, not either against
experiment.

Recorded limitation: `thermo`'s vapor-fraction-specified flash
(`flash(P=..., VF=0)`) does not work on this phase pair - it returns
temperatures of 1.8e5 K and 9.2e3 K. Bubble and dew temperatures from `thermo`
are therefore obtained by bisecting `thermo`'s own T,P flashes, which work
normally. See `test_thermo_vapor_fraction_specified_flash_is_unusable_here`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.models._antoine import antoine_saturation_pressures

PRESSURE_PA = 101325.0
NAMES = ("1-Propanol", "Water")
THERMO_IDS = ("1-propanol", "water")
_PROPANOL, _WATER = 0, 2

#: Flash settings used whenever a phase *boundary* is located by bisecting
#: chemthermo verdicts: the default `tpd_tol = 1e-8` displaces the boundary by
#: a 1e-8 relative offset in `sum_i x_i gamma_i Psat_i / P`.
SHARP = ct.FlashSettings(stability_settings=ct.StabilitySettings(tpd_tol=1e-14))

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "nrtl" / "tessier2000_problem1.json"


def _parameters() -> tuple[list[list[float]], list[list[float]]]:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    tau, alpha = payload["tau"], payload["alpha"]
    return (
        [[0.0, float(tau[_PROPANOL][_WATER])], [float(tau[_WATER][_PROPANOL]), 0.0]],
        [[0.0, float(alpha[_PROPANOL][_WATER])], [float(alpha[_WATER][_PROPANOL]), 0.0]],
    )


def _mixture(z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(NAMES), list(z), normalize=True)


def _model() -> ct.NRTL:
    tau, alpha = _parameters()
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(NAMES[0], NAMES[1], tau[0][1], tau[1][0], alpha[0][1], alpha[1][0])]
        )
    )


def _flash(z: Sequence[float], temperature_K: float, settings=None) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(z),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=_model(),
        flash_mode="modified-raoult",
        settings=settings,
    )


def _thermo_vapor_pressures(thermo: Any) -> list[Any]:
    """`thermo.VaporPressure` objects carrying the chemthermo Antoine records.

    Unit/base conversion: chemthermo stores `ln(P/bar) = A - B/(T + C)`;
    `chemicals.vapor_pressure.Antoine(T, A, B, C, base)` computes
    `base**(A - B/(T + C))` in **Pa**. With `base = e` the only change needed is
    `A -> A + ln(1e5)`.
    """
    mixture = _mixture((0.5, 0.5))
    pressures = []
    for component in mixture.components:
        antoine = component.antoine
        assert antoine is not None
        vapor_pressure = thermo.VaporPressure(
            Tb=None, Tc=component.tc_k, Pc=component.pc_pa, omega=component.omega, CASRN=None
        )
        vapor_pressure.add_correlation(
            name="chemthermo",
            model="Antoine",
            Tmin=antoine.Tmin_K,
            Tmax=antoine.Tmax_K,
            A=antoine.A + math.log(1.0e5),
            B=antoine.B,
            C=antoine.C,
            base=math.e,
        )
        pressures.append(vapor_pressure)
    return pressures


def _thermo_flasher(thermo: Any, temperature_K: float) -> Any:
    tau, alpha = _parameters()
    constants, correlations = thermo.ChemicalConstantsPackage.from_IDs(list(THERMO_IDS))
    excess = thermo.NRTL(
        T=temperature_K,
        xs=[0.5, 0.5],
        tau_coeffs=[[[tau[i][j], 0, 0, 0, 0, 0] for j in range(2)] for i in range(2)],
        alpha_coeffs=[[[alpha[i][j], 0] for j in range(2)] for i in range(2)],
    )
    liquid = thermo.GibbsExcessLiquid(
        VaporPressures=_thermo_vapor_pressures(thermo),
        GibbsExcessModel=excess,
        HeatCapacityGases=correlations.HeatCapacityGases,
        T=temperature_K,
        P=PRESSURE_PA,
        zs=[0.5, 0.5],
    )
    gas = thermo.IdealGas(
        HeatCapacityGases=correlations.HeatCapacityGases,
        T=temperature_K,
        P=PRESSURE_PA,
        zs=[0.5, 0.5],
    )
    return thermo.FlashVL(constants, correlations, liquid=liquid, gas=gas)


def _bisect(low: float, high: float, predicate: Callable[[float], bool]) -> float:
    """Temperature where `predicate` flips between `low` and `high`."""
    low_value = predicate(low)
    assert low_value != predicate(high), "bracket does not straddle the boundary"
    for _ in range(60):
        middle = 0.5 * (low + high)
        if predicate(middle) == low_value:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _modified_raoult_residual(x: np.ndarray, y: np.ndarray, temperature_K: float) -> float:
    """`max_i |ln(x_i gamma_i Psat_i / P) - ln y_i|`, the equilibrium condition."""
    mixture = _mixture((0.5, 0.5))
    gamma = np.asarray(
        _model().activity_coefficients(
            mixture=mixture, temperature_K=temperature_K, composition=x.tolist()
        ),
        dtype=float,
    )
    psat = antoine_saturation_pressures(mixture, temperature_K)
    return float(np.max(np.abs(np.log(x * gamma * psat / PRESSURE_PA) - np.log(y))))


# --------------------------------------------------------------------------


@pytest.mark.parametrize("temperature_K", [300.0, 330.0, 361.0, 380.0, 399.0])
def test_saturation_pressures_match_thermo_to_machine_precision(temperature_K: float) -> None:
    """Establish the shared reference fugacity before comparing any equilibrium.

    Achieved: worst relative difference 2.3e-15 over the five temperatures.
    """
    thermo = pytest.importorskip("thermo")
    ours = antoine_saturation_pressures(_mixture((0.5, 0.5)), temperature_K)
    theirs = np.array(
        [vp.T_dependent_property(temperature_K) for vp in _thermo_vapor_pressures(thermo)]
    )
    assert float(np.max(np.abs(theirs - ours) / ours)) < 1e-10


def test_vapor_liquid_flash_matches_thermo() -> None:
    """A two-phase state: vapor fraction, x and y.

    Achieved at 361 K, z = (0.5, 0.5): |d beta| = 4.3e-06, worst |dx| = 3.0e-07,
    worst |dy| = 3.4e-13. Both solutions satisfy their own material balance to
    round-off, so the difference is in the equilibrium itself and is adjudicated
    below with the modified-Raoult residual: chemthermo 1.6e-15, `thermo`
    2.0e-07, i.e. the difference is `thermo`'s convergence tolerance.
    """
    thermo = pytest.importorskip("thermo")
    temperature_K = 361.0
    z = np.array([0.5, 0.5])

    ours = _flash((0.5, 0.5), temperature_K)
    x = np.array(ours.phases["liquid"].composition.fractions)
    y = np.array(ours.phases["vapor"].composition.fractions)
    beta = ours.vapor_fraction
    assert beta is not None

    result = _thermo_flasher(thermo, temperature_K).flash(
        T=temperature_K, P=PRESSURE_PA, zs=z.tolist()
    )
    their_x = np.array(result.liquid0.zs)
    their_y = np.array(result.gas.zs)

    assert abs(beta - result.VF) < 1e-5
    assert float(np.max(np.abs(x - their_x))) < 1e-5
    assert float(np.max(np.abs(y - their_y))) < 1e-5

    ours_residual = _modified_raoult_residual(x, y, temperature_K)
    theirs_residual = _modified_raoult_residual(their_x, their_y, temperature_K)
    assert ours_residual < 1e-12
    assert ours_residual < theirs_residual


@pytest.mark.parametrize(
    "z1, bubble_low, bubble_high, dew_low, dew_high",
    [(0.3, 350.0, 361.0, 364.0, 365.0), (0.5, 350.0, 361.0, 361.5, 362.0)],
)
def test_bubble_and_dew_temperatures_match_thermo(
    z1: float, bubble_low: float, bubble_high: float, dew_low: float, dew_high: float
) -> None:
    """Both packages' phase boundaries, each found by bisecting its own verdicts.

    Achieved: |dT| <= 2.8e-08 K on all four boundaries, against an asserted
    1e-05 K.
    """
    thermo = pytest.importorskip("thermo")
    z = [z1, 1.0 - z1]
    flasher = _thermo_flasher(thermo, 361.0)

    def their_vapor_fraction(temperature_K: float) -> float:
        return float(flasher.flash(T=temperature_K, P=PRESSURE_PA, zs=z).VF)

    def ours_two_phase(temperature_K: float) -> bool:
        return len(_flash(z, temperature_K, SHARP).phase_names()) == 2

    our_bubble = _bisect(bubble_low, bubble_high, ours_two_phase)
    their_bubble = _bisect(bubble_low, bubble_high, lambda t: their_vapor_fraction(t) > 0.0)
    assert abs(our_bubble - their_bubble) < 1e-5, (our_bubble, their_bubble)

    our_dew = _bisect(dew_low, dew_high, ours_two_phase)
    their_dew = _bisect(dew_low, dew_high, lambda t: their_vapor_fraction(t) >= 1.0)
    assert abs(our_dew - their_dew) < 1e-5, (our_dew, their_dew)


def test_thermo_vapor_fraction_specified_flash_is_unusable_here() -> None:
    """Honest record of why bubble and dew points are bisected above.

    `thermo.FlashVL.flash(P=..., VF=0)` on this `GibbsExcessLiquid` /
    `IdealGas` pair returns T = 1.8e5 K for the bubble point and T = 9.2e3 K
    for the dew point of an equimolar feed whose true boundaries are 360.99 K
    and 361.61 K. The T,P flashes on the same objects are correct (previous
    test), so this is a solver-path problem in `thermo`, not a model
    disagreement. Asserted so that a future `thermo` which fixes it is noticed.
    """
    thermo = pytest.importorskip("thermo")
    flasher = _thermo_flasher(thermo, 361.0)
    for vapor_fraction in (0.0, 1.0):
        result = flasher.flash(P=PRESSURE_PA, zs=[0.5, 0.5], VF=vapor_fraction)
        assert not (350.0 < result.T < 375.0), (vapor_fraction, result.T)
