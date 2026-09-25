"""Adjudicate the 5 unresolved Case F-1 phase-count disagreements with `thermo`.

Case F-1 (`tests/validation/test_flash_phase_detection_vs_thermo.py`) found 5
states of Methane/Propane/n-Decane (0.7, 0.2, 0.1) where `flash_tp` splits into
two phases (`stability_tp` reports `status == "unstable"`) while `thermo`'s
`FlashVL` reports `VF = 0.0`, a single phase. That left the direction of the
disagreement open.

This module adjudicates it *inside thermo's own model*: build a `thermo`
`PRMIX` `CEOSLiquid` with chemthermo's own `Tc`, `Pc` and `omega` (`kij = 0`,
the same constants `flash_tp` used), take chemthermo's `stability_tp`
minimizing trial composition `w`, and evaluate the reduced tangent-plane
distance

    tpd(w) = sum_i w_i * [ln w_i + ln phi_i(w) - ln z_i - ln phi_i(z)]

using *thermo's* fugacity coefficients (`lnphis_at_zs(..., most_stable=True)`,
the minimum-Gibbs root) at both `w` and the feed `z`. If thermo's own
fugacity-coefficient model puts a negative tpd at chemthermo's stationary
point, then `FlashVL`'s `VF = 0` answer is a miss by thermo's stability
*search*, not evidence that the feed is stable under thermo's *model*.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import chemthermo as ct

thermo = pytest.importorskip("thermo")
CEOSLiquid = thermo.CEOSLiquid
CEOSGas = thermo.CEOSGas
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX

NAMES = ("Methane", "Propane", "n-Decane")
CAS = ("74-82-8", "74-98-6", "124-18-5")
FEED = (0.7, 0.2, 0.1)

#: The 5 Case F-1 states where chemthermo splits and `thermo.FlashVL` reports
#: VF = 0. Orchestrator's recomputed numbers (this module re-derives them, see
#: `test_all_five_states_are_unstable_under_thermos_own_fugacities` below):
#: T=270K,P=1.22e7Pa tpd_ct=-5.2781e-2 tpd_thermo=-5.2788e-2;
#: T=300K,P=1.5e7Pa tpd_ct=-3.9096e-2 tpd_thermo=-3.9107e-2;
#: T=330K,P=1.98e7Pa tpd_ct=-6.9083e-3 tpd_thermo=-6.9142e-3;
#: T=360K,P=1.5e7Pa tpd_ct=-4.7123e-2 tpd_thermo=-4.7135e-2;
#: T=450K,P=1.5e7Pa tpd_ct=-3.5904e-2 tpd_thermo=-3.5939e-2. `thermo` VF = 0.0
#: at all 5.
STATES: tuple[tuple[float, float], ...] = (
    (270.0, 1.22e7),
    (300.0, 1.5e7),
    (330.0, 1.98e7),
    (360.0, 1.5e7),
    (450.0, 1.5e7),
)

EOS = ct.PengRobinsonEOS()


def _thermo_liquid_and_flasher() -> tuple[Any, Any]:
    """A thermo `CEOSLiquid` and `FlashVL` built from chemthermo's own EOS constants."""
    components = [ct.Component.from_database(name) for name in NAMES]
    size = len(NAMES)
    eos_kwargs = {
        "Tcs": [c.tc_k for c in components],
        "Pcs": [c.pc_pa for c in components],
        "omegas": [c.omega for c in components],
        "kijs": [[0.0] * size for _ in range(size)],
    }
    constants, properties = ChemicalConstantsPackage.from_IDs(list(CAS))
    liquid = CEOSLiquid(
        PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases
    )
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
    return liquid, flasher


def _thermo_tpd_at(liquid: Any, temperature_K: float, pressure_Pa: float, w: np.ndarray) -> float:
    """Reduced tangent-plane distance at trial composition `w`, using thermo's own
    fugacity coefficients (minimum-Gibbs root) for both `w` and the feed."""
    feed = np.array(FEED, dtype=float)
    liquid_state = liquid.to_TP_zs(T=temperature_K, P=pressure_Pa, zs=list(FEED))
    ln_phi_feed = np.array(liquid_state.lnphis_at_zs(list(FEED), most_stable=True))
    ln_phi_w = np.array(liquid_state.lnphis_at_zs(w.tolist(), most_stable=True))
    ln_d_feed = np.log(feed) + ln_phi_feed
    return float(np.sum(w * (np.log(w) + ln_phi_w - ln_d_feed)))


def test_all_five_states_are_unstable_under_thermos_own_fugacities() -> None:
    """Case F-1 adjudication: all 5 disagreement states are unstable in thermo's model.

    For every state: (a) chemthermo's `stability_tp` reports `"unstable"`; (b)
    thermo's own fugacity coefficients, evaluated at chemthermo's minimizing
    trial composition, give a negative reduced tpd within 5e-5 (absolute) of
    chemthermo's `tpd_min`; (c) chemthermo's `flash_tp` two-phase split lowers
    the Gibbs energy (`delta_g_split_rt < 0`). thermo's `VF` is recorded, not
    asserted: today it is 0.0 (a single-phase miss) at all 5 states, but the
    adjudication above does not depend on that, so the test stays meaningful
    if a future thermo release finds the split too.
    """
    mixture = ct.Mixture.from_database(list(NAMES), list(FEED), normalize=True)
    liquid, flasher = _thermo_liquid_and_flasher()

    adjudicated = 0
    already_resolved = 0
    for temperature_K, pressure_Pa in STATES:
        stability = ct.stability_tp(
            mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS
        )
        assert stability.status == "unstable", (temperature_K, pressure_Pa)
        assert stability.trial_composition is not None
        w = np.array(stability.trial_composition, dtype=float)

        split = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS)
        assert set(split.phase_names()) == {"liquid", "vapor"}, (temperature_K, pressure_Pa)
        delta_g_split_rt = float(split.diagnostics["delta_g_split_rt"])
        assert delta_g_split_rt < 0.0, (temperature_K, pressure_Pa)

        thermo_tpd_at_w = _thermo_tpd_at(liquid, temperature_K, pressure_Pa, w)
        assert thermo_tpd_at_w < 0.0, (temperature_K, pressure_Pa)
        assert thermo_tpd_at_w == pytest.approx(stability.tpd_min, abs=5e-5), (
            temperature_K,
            pressure_Pa,
        )

        reference = flasher.flash(T=temperature_K, P=pressure_Pa, zs=list(FEED))
        thermo_vf = float(reference.VF)  # recorded, not asserted
        if thermo_vf == 0.0:
            adjudicated += 1
        else:
            already_resolved += 1
            print(
                f"NOTE: thermo now finds a two-phase split (VF={thermo_vf:.6g}) at "
                f"T={temperature_K} K, P={pressure_Pa:.4g} Pa; this disagreement is "
                "no longer open on thermo's side."
            )

    assert adjudicated + already_resolved == len(STATES)
    assert adjudicated + already_resolved >= 4
