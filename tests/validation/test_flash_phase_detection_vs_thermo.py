"""`flash_tp` phase detection against `thermo`'s `FlashVL` (validation Cases F-1, F-2).

Both sides run Peng-Robinson with **the same** `Tc`, `Pc` and `omega` (read from
the chemthermo databank and handed to `thermo`'s `ChemicalConstantsPackage`) and
`kij = 0`, so a verdict difference is a solver difference, not a parameter
difference. The `thermo` ideal-gas heat capacities are pulled from its own
databank; they do not enter an isothermal-isobaric VLE flash verdict.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import chemthermo as ct

thermo = pytest.importorskip("thermo")
CEOSGas = thermo.CEOSGas
CEOSLiquid = thermo.CEOSLiquid
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX

EOS = ct.PengRobinsonEOS()
LEGACY = ct.FlashSettings(phase_detection="wilson-heuristic")

#: (names, CAS numbers, feed). CAS numbers only select `thermo`'s heat-capacity
#: correlations; every EOS parameter comes from the chemthermo databank below.
SYSTEMS: tuple[tuple[tuple[str, ...], tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), ("74-82-8", "74-84-0"), (0.5, 0.5)),
    (("Methane", "n-Pentane"), ("74-82-8", "109-66-0"), (0.6, 0.4)),
    (("Ethane", "n-Heptane"), ("74-84-0", "142-82-5"), (0.7, 0.3)),
    (
        ("Methane", "Ethane", "Propane"),
        ("74-82-8", "74-84-0", "74-98-6"),
        (0.5, 0.3, 0.2),
    ),
    (
        ("Propane", "n-Butane", "n-Pentane"),
        ("74-98-6", "106-97-8", "109-66-0"),
        (0.4, 0.3, 0.3),
    ),
)
GRID_T_K = (170.0, 175.0, 200.0, 240.0, 280.0, 320.0, 360.0)
GRID_P_PA = (2.0e5, 1.0e6, 1.778e6, 3.0e6, 8.0e6)


def _flasher(names: tuple[str, ...], cas: tuple[str, ...]) -> Any:
    """`thermo` FlashVL using chemthermo's own critical constants."""
    components = [ct.Component.from_database(name) for name in names]
    constants = ChemicalConstantsPackage(
        Tcs=[c.tc_k for c in components],
        Pcs=[c.pc_pa for c in components],
        omegas=[c.omega for c in components],
        MWs=[c.mw_kg_per_mol * 1000.0 for c in components],
        names=list(names),
        CASs=list(cas),
    )
    _, properties = ChemicalConstantsPackage.from_IDs(list(cas))
    size = len(names)
    eos_kwargs = {
        "Tcs": constants.Tcs,
        "Pcs": constants.Pcs,
        "omegas": constants.omegas,
        "kijs": [[0.0] * size for _ in range(size)],
    }
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(
        PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases
    )
    return FlashVL(constants, properties, liquid=liquid, gas=gas)


def _phase_count(result: ct.FlashResult) -> int:
    beta = result.vapor_fraction
    return 2 if (beta is not None and 0.0 < beta < 1.0) else 1


def test_the_disagreement_state_is_two_phase_in_thermo_and_lowers_the_gibbs_energy() -> None:
    """Case F-2: Methane/n-Pentane at 175 K, 1.778 MPa.

    The legacy heuristic returns a single liquid (`rr_no_root`). The
    tangent-plane path returns a two-phase split. Two independent witnesses say
    the split is the right answer:

    1. `thermo`'s `FlashVL` with the same constants reports VF = 0.0791276.
    2. The Gibbs energy of the split, evaluated at *thermo's* compositions with
       chemthermo's own fugacity coefficients on the minimum-Gibbs root, is
       below the single-phase feed by dG/RT = -1.767e-3.
    """
    from chemthermo.stability.tp import _ln_phi_min_gibbs

    names = ("Methane", "n-Pentane")
    z = (0.6, 0.4)
    temperature_K, pressure_Pa = 175.0, 1.778e6
    mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)

    legacy = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=EOS,
        settings=LEGACY,
    )
    assert legacy.phase_names() == ["liquid"]

    result = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS)
    assert set(result.phase_names()) == {"liquid", "vapor"}
    beta = result.vapor_fraction
    assert beta is not None

    reference = _flasher(names, ("74-82-8", "109-66-0")).flash(
        T=temperature_K, P=pressure_Pa, zs=list(z)
    )
    reference_beta = float(reference.VF)
    assert 0.0 < reference_beta < 1.0
    # Achieved: |d beta| = 1.08e-4, max |dx| = 6.8e-5, max |dy| = 1.9e-7.
    assert beta == pytest.approx(reference_beta, abs=5e-4)
    assert np.allclose(
        result.phases["liquid"].composition.fractions, reference.liquid0.zs, atol=1e-3
    )
    assert np.allclose(result.phases["vapor"].composition.fractions, reference.gas.zs, atol=1e-3)

    # Independent Gibbs-energy check at thermo's split, using chemthermo's EOS.
    def reduced_g(fractions: np.ndarray) -> float:
        ln_phi, _branch = _ln_phi_min_gibbs(
            EOS,
            mixture=mixture,
            temperature=temperature_K,
            pressure=pressure_Pa,
            composition=fractions,
        )
        mask = fractions > 0.0
        return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_phi[mask])))

    feed = np.array(mixture.fractions, dtype=float)
    reference_x = np.array(reference.liquid0.zs, dtype=float)
    reference_y = np.array(reference.gas.zs, dtype=float)
    delta_g = (
        reference_beta * reduced_g(reference_y)
        + (1.0 - reference_beta) * reduced_g(reference_x)
        - reduced_g(feed)
    )
    assert delta_g == pytest.approx(-1.7666e-3, rel=1e-3)
    assert delta_g < 0.0

    # And the solver's own reported reduction agrees in sign and magnitude.
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert float(result.diagnostics["delta_g_split_rt"]) == pytest.approx(delta_g, rel=5e-2)


def test_phase_verdicts_over_a_grid_match_thermo_more_often_than_the_heuristic() -> None:
    """Case F-1: verdict and vapor-fraction sweep against `FlashVL`.

    Counts recorded over the 175 states of this grid (5 systems x 7
    temperatures x 5 pressures); the assertions below enforce the load-bearing
    ones:

    - tangent-plane path: 175/175 verdicts match `thermo`; 0 mismatches.
    - legacy heuristic: 166/175. Its 9 misses are 8 states where it raises
      `ConvergenceError` (thermo says single phase on all 8, and the
      tangent-plane path now returns a single phase there) and 1 state where it
      returns a single phase and thermo says two (the Case F-2 state).
    - 56 states are two-phase on both sides; worst |beta - beta_thermo| there is
      1.601e-4.
    """
    total = 0
    tangent_matches = 0
    legacy_matches = 0
    legacy_nonconvergence_now_answered = 0
    worst_beta_deviation = 0.0
    two_phase_states = 0

    for names, cas, z in SYSTEMS:
        mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)
        flasher = _flasher(names, cas)
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                try:
                    reference = flasher.flash(T=temperature_K, P=pressure_Pa, zs=list(z))
                except Exception:  # pragma: no cover - reference solver failure
                    continue
                reference_beta = float(reference.VF)
                reference_count = 2 if 0.0 < reference_beta < 1.0 else 1
                total += 1

                try:
                    result = ct.flash_tp(
                        mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS
                    )
                    tangent_count = _phase_count(result)
                except ct.ConvergenceError:
                    result = None
                    tangent_count = -1

                try:
                    legacy = ct.flash_tp(
                        mixture,
                        temperature_K=temperature_K,
                        pressure_Pa=pressure_Pa,
                        eos=EOS,
                        settings=LEGACY,
                    )
                    legacy_count = _phase_count(legacy)
                except ct.ConvergenceError:
                    legacy_count = -1

                if tangent_count == reference_count:
                    tangent_matches += 1
                if legacy_count == reference_count:
                    legacy_matches += 1
                if legacy_count == -1 and tangent_count == 1:
                    legacy_nonconvergence_now_answered += 1
                    assert reference_count == 1, (names, temperature_K, pressure_Pa)

                if result is not None and tangent_count == 2 and reference_count == 2:
                    two_phase_states += 1
                    beta = result.vapor_fraction
                    assert beta is not None
                    worst_beta_deviation = max(worst_beta_deviation, abs(beta - reference_beta))

    assert total >= 150
    assert two_phase_states >= 30
    assert tangent_matches >= legacy_matches
    assert tangent_matches == total
    assert legacy_nonconvergence_now_answered >= 5
    assert tangent_matches > legacy_matches
    assert worst_beta_deviation < 1e-3
