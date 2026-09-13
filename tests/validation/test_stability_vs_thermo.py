"""Cross-check chemthermo's tangent-plane stability test against `thermo` 0.6.0.

The reference `thermo` Peng-Robinson mixture is built from EXACTLY the Tc, Pc and
omega values carried by chemthermo's own `Component` objects (not thermo's
databank) with all kijs = 0, so the only remaining differences between the two
codes are implementation details of the EOS itself.

One such difference is known and deliberate: chemthermo uses the rounded
Peng-Robinson constants 0.45724 / 0.07780, while `thermo` uses the exact roots
0.4572355289213822 / 0.0777960739038885. That shifts ln(phi) by up to ~2e-4 at
these states, which is the dominant term in the tolerances below.

Achieved on the states below (Python 3.11, numpy, thermo 0.6.0):
- stable/unstable verdict: identical on all 7 states.
- |tpd_chemthermo - tpd_thermo| <= 5.7e-5.
- max |w_chemthermo - w_thermo| <= 2.5e-5.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.flash._common import wilson_k

thermo = pytest.importorskip("thermo")
CEOSGas = thermo.CEOSGas
CEOSLiquid = thermo.CEOSLiquid
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX
stability_iteration_Michelsen = thermo.flash.flash_utils.stability_iteration_Michelsen

EOS = ct.PengRobinsonEOS()

# (names, z, T [K], P [Pa], expected status)
CASES: list[tuple[list[str], list[float], float, float, str]] = [
    (["Methane", "Ethane", "Propane"], [0.5, 0.3, 0.2], 240.0, 3.0e6, "unstable"),
    (["Methane", "Ethane"], [0.5, 0.5], 450.0, 1.0e5, "stable"),
    (["Methane", "Ethane", "Propane"], [0.5, 0.3, 0.2], 300.0, 5.0e7, "stable"),
    (["Methane", "Ethane", "Propane"], [0.5, 0.3, 0.2], 200.0, 1.0e8, "stable"),
    (["Methane", "Ethane", "Propane"], [0.5, 0.3, 0.2], 220.0, 2.0e6, "unstable"),
    (["Methane", "Ethane"], [0.5, 0.5], 200.0, 2.0e6, "unstable"),
    (["Methane", "Propane"], [0.7, 0.3], 250.0, 5.0e6, "unstable"),
]

TPD_ATOL = 1e-3
COMP_ATOL = 1e-3


def _reference_phases(names: list[str]) -> tuple[Any, Any, Any]:
    """Build thermo PR gas/liquid phases from chemthermo's component constants."""
    components = [ct.Component.from_database(name) for name in names]
    tcs = [component.tc_k for component in components]
    pcs = [component.pc_pa for component in components]
    omegas = [component.omega for component in components]

    base_constants, base_properties = ChemicalConstantsPackage.from_IDs(
        [name.lower() for name in names]
    )
    constants = ChemicalConstantsPackage(
        Tcs=tcs,
        Pcs=pcs,
        omegas=omegas,
        MWs=base_constants.MWs,
        names=base_constants.names,
        CASs=base_constants.CASs,
    )
    size = len(names)
    eos_kwargs = {
        "Tcs": tcs,
        "Pcs": pcs,
        "omegas": omegas,
        "kijs": [[0.0] * size for _ in range(size)],
    }
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=base_properties.HeatCapacityGases)
    liquid = CEOSLiquid(
        PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=base_properties.HeatCapacityGases
    )
    flasher = FlashVL(constants, base_properties, liquid=liquid, gas=gas)
    return flasher, gas, liquid


def _min_and_other_phase(gas: Any, liquid: Any, T: float, P: float, zs: list[float]):
    """Pick the lowest-Gibbs phase at (T, P, z), as thermo's own driver does."""
    gas_state = gas.to(T=T, P=P, zs=zs)
    liquid_state = liquid.to(T=T, P=P, zs=zs)
    if liquid_state.G_dep() < gas_state.G_dep():
        return liquid_state, gas_state
    return gas_state, liquid_state


def _trial_guesses(names: list[str], zs: list[float], T: float, P: float):
    """The same deterministic trial set chemthermo uses."""
    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    k = wilson_k(mixture, T, P)
    z = np.array(mixture.fractions, dtype=float)
    guesses = [
        ("wilson-vapor", (z * k) / float(np.sum(z * k))),
        ("wilson-liquid", (z / k) / float(np.sum(z / k))),
    ]
    size = len(names)
    for index in range(size):
        w = np.full(size, 1e-3 / (size - 1))
        w[index] = 1.0 - 1e-3
        guesses.append((f"pure-{names[index]}", w))
    return guesses


def _reference_stationary_point(names: list[str], zs: list[float], T: float, P: float):
    """Largest sum(W) over thermo's Michelsen iteration from the same trial set."""
    _, gas, liquid = _reference_phases(names)
    min_phase, other_phase = _min_and_other_phase(gas, liquid, T, P, zs)
    fugacities_trial = min_phase.fugacities_lowest_Gibbs()

    def test_phase(composition: list[float]) -> list[float]:
        return other_phase.lnphis_at_zs(composition, most_stable=True)

    best_sum_w = 1.0
    best_w = np.array(zs, dtype=float)
    for _, guess in _trial_guesses(names, zs, T, P):
        try:
            solution = stability_iteration_Michelsen(
                T=T,
                P=P,
                zs_trial=list(zs),
                fugacities_trial=fugacities_trial,
                zs_test=list(guess),
                test_phase=test_phase,
                maxiter=500,
                xtol=1e-14,
            )
        except Exception:  # noqa: BLE001 - thermo raises on non-convergence
            continue
        sum_w = float(solution[0])
        if math.isfinite(sum_w) and sum_w > best_sum_w:
            best_sum_w = sum_w
            best_w = np.array(solution[2], dtype=float)
    return best_sum_w, best_w


@pytest.mark.parametrize(
    ("names", "zs", "T", "P", "expected"),
    CASES,
    ids=[f"{'-'.join(c[0])}@{c[2]:g}K-{c[3]:g}Pa" for c in CASES],
)
def test_stability_verdict_matches_thermo(
    names: list[str], zs: list[float], T: float, P: float, expected: str
) -> None:
    """Verdicts must agree with thermo's own Michelsen stability test."""
    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    result = ct.stability_tp(mixture, temperature_K=T, pressure_Pa=P, eos=EOS)
    assert result.status == expected

    flasher, gas, liquid = _reference_phases(names)
    min_phase, other_phase = _min_and_other_phase(gas, liquid, T, P, zs)
    reference_stable, _ = flasher.stability_test_Michelsen(
        T, P, zs, min_phase=min_phase, other_phase=other_phase
    )
    assert result.stable is bool(reference_stable)


@pytest.mark.parametrize(
    ("names", "zs", "T", "P", "expected"),
    [case for case in CASES if case[4] == "unstable"],
    ids=[f"{'-'.join(c[0])}@{c[2]:g}K-{c[3]:g}Pa" for c in CASES if c[4] == "unstable"],
)
def test_unstable_stationary_point_matches_thermo(
    names: list[str], zs: list[float], T: float, P: float, expected: str
) -> None:
    """Compare the converged stationary point (sum W, tpd and w) with thermo."""
    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    result = ct.stability_tp(mixture, temperature_K=T, pressure_Pa=P, eos=EOS)
    assert result.trial_composition is not None

    reference_sum_w, reference_w = _reference_stationary_point(names, zs, T, P)
    assert reference_sum_w > 1.0

    reference_tpd = -math.log(reference_sum_w)
    assert result.tpd_min == pytest.approx(reference_tpd, abs=TPD_ATOL)
    assert np.allclose(np.array(result.trial_composition), reference_w, atol=COMP_ATOL)

    sum_w = float(result.diagnostics["sum_W"])
    assert sum_w == pytest.approx(reference_sum_w, abs=TPD_ATOL)
