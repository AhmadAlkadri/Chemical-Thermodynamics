"""Deterministic cross-check of tangent-plane stability against `thermo`.

The reference Peng-Robinson mixture in `thermo` is built from EXACTLY the Tc, Pc
and omega values carried by chemthermo's own `Component` objects, with all
kijs = 0, so the comparison isolates the stability algorithm rather than the
component databank.

Known residual difference: chemthermo uses the rounded Peng-Robinson constants
0.45724 / 0.07780 while `thermo` uses the exact roots 0.4572355289213822 /
0.0777960739038885. That shifts ln(phi) by up to ~2e-4 here and dominates the
tolerances below.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

import chemthermo as ct
from chemthermo.flash._common import wilson_k

try:
    import thermo
except ImportError:
    print("Error: 'thermo' package is required. Install with: pip install -e '.[validation]'")
    raise SystemExit(1)


CEOSGas = thermo.CEOSGas
CEOSLiquid = thermo.CEOSLiquid
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX
stability_iteration_Michelsen = thermo.flash.flash_utils.stability_iteration_Michelsen

TPD_TOL = 1e-3
COMP_TOL = 1e-3

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


def _reference_phases(names: Sequence[str]) -> tuple[Any, Any, Any]:
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


def _min_and_other_phase(gas: Any, liquid: Any, T: float, P: float, zs: Sequence[float]):
    gas_state = gas.to(T=T, P=P, zs=list(zs))
    liquid_state = liquid.to(T=T, P=P, zs=list(zs))
    if liquid_state.G_dep() < gas_state.G_dep():
        return liquid_state, gas_state
    return gas_state, liquid_state


def _trial_guesses(names: Sequence[str], zs: Sequence[float], T: float, P: float):
    mixture = ct.Mixture.from_database(list(names), list(zs), normalize=True)
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


def _reference_stationary_point(names: Sequence[str], zs: Sequence[float], T: float, P: float):
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


def main() -> int:
    print("Stability validation (chemthermo vs thermo, Michelsen tangent-plane test)")
    print(f"Tolerances: |d tpd| <= {TPD_TOL}, max |d w| <= {COMP_TOL}")
    failures = 0

    for names, zs, T, P, expected in CASES:
        mixture = ct.Mixture.from_database(names, zs, normalize=True)
        result = ct.stability_tp(mixture, temperature_K=T, pressure_Pa=P, eos=ct.PengRobinsonEOS())

        flasher, gas, liquid = _reference_phases(names)
        min_phase, other_phase = _min_and_other_phase(gas, liquid, T, P, zs)
        reference_stable, _ = flasher.stability_test_Michelsen(
            T, P, list(zs), min_phase=min_phase, other_phase=other_phase
        )

        label = f"{'/'.join(names)} z={zs} T={T:.1f} K P={P:.3e} Pa"
        print("-" * 78)
        print(label)
        print(
            f"  chemthermo: status={result.status:<10} tpd_min={result.tpd_min: .8e}  "
            f"thermo stable={bool(reference_stable)}"
        )

        ok = result.status == expected and result.stable is bool(reference_stable)
        if not ok:
            print(f"  FAIL: expected {expected!r}, verdicts disagree or mismatch.")
            failures += 1
            continue

        if expected == "unstable":
            reference_sum_w, reference_w = _reference_stationary_point(names, zs, T, P)
            reference_tpd = -math.log(reference_sum_w)
            tpd_diff = abs(result.tpd_min - reference_tpd)
            assert result.trial_composition is not None
            comp_diff = float(np.max(np.abs(np.array(result.trial_composition) - reference_w)))
            print(
                f"  sum(W): chemthermo={float(result.diagnostics['sum_W']):.10f} "
                f"thermo={reference_sum_w:.10f}"
            )
            print(f"  |d tpd| = {tpd_diff:.3e}   max |d w| = {comp_diff:.3e}")
            if tpd_diff > TPD_TOL or comp_diff > COMP_TOL:
                print("  FAIL: stationary point exceeded tolerance.")
                failures += 1
                continue

        print("  PASS")

    print("=" * 78)
    if failures:
        print(f"FAIL: {failures} case(s) outside tolerance.")
        return 1
    print(f"PASS: all {len(CASES)} stability cases agree with thermo.")
    return 0


if __name__ == "__main__":
    exit_code = main()
    if exit_code != 0:
        raise SystemExit(exit_code)
