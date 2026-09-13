"""Phase addition and removal in `flash_tp` (ADR-0011, validation Case V-3).

This module covers the *mechanism*: the `FlashSettings.max_phases` contract,
the add-then-remove route that resolves the refusal window of Case R-3, the
result contract for a three-phase answer, and the regression that none of it
reaches the paths it was not wired to.

The tie-triangle itself, with its independent Newton solve and its Gibbs
ordering, is validation Cases V-1 and V-2 in
`tests/validation/test_vlle_water_propanol_butanol.py`.

System for the refusal window
-----------------------------
n-Butanol(1) / Water(2) at P = 101325 Pa, `flash_mode="modified-raoult"`, NRTL
pair 2-3 of Table 1 of Tessier, Brennecke and Stadtherr, Chem. Eng. Sci. 55
(2000) 1785-1796 (tau_12 = 0.90047, tau_21 = 3.51307, alpha = 0.48 implied by
the printed `G`), Antoine from the packaged databank. The three-phase
temperature of this model is T3 = 366.213774 K (Case R-3), and the binodal is
temperature independent here because the fitted tau are, so the two-liquid
answer below T3 is the *same* tie-line at every temperature:
x1 = 0.019998419467 / 0.359999661508.

Why the window needs *removal*
------------------------------
Just below T3 the deepest tangent-plane minimum from the feed is a vapor, so
the first two-phase iterate is a vapor-liquid pair that is not the equilibrium.
Adding the second liquid the post-split test finds gives a three-phase set, and
a binary cannot have three phases at any temperature other than T3 (Gibbs'
phase rule: F = 2 - 3 + 2 = 1), so the three-phase Rachford-Rice has no finite
solution and the vapor's amount runs negative. Removing it leaves the two
liquids. The history is `"V -> LV -> LLV -> LL"`.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct

PRESSURE_PA = 101325.0

NAMES = ("n-Butanol", "Water")
TAU_12 = 0.90047
TAU_21 = 3.51307
ALPHA = 0.48
FEED = (0.20, 0.80)

#: Three-phase temperature of this model, from validation Case R-3.
T3_K = 366.213774
#: Temperature-independent binodal of this NRTL pair (Cases L-3, R-1, R-3).
BINODAL_X1 = (0.019998419467, 0.359999661508)

# Mirrors the phi-phi grid of tests/test_flash_phase_detection.py, kept in
# sync by hand so this module stays self-contained.
GRID_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
GRID_T_K = (170.0, 200.0, 240.0, 280.0, 320.0, 360.0)
GRID_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)


@pytest.fixture(scope="module")
def butanol_water() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)]
        )
    )


def _flash(
    model: ct.NRTL,
    temperature_K: float,
    z: Sequence[float] = FEED,
    settings: ct.FlashSettings | None = None,
) -> ct.FlashResult:
    return ct.flash_tp(
        ct.Mixture.from_database(list(NAMES), list(z), normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
        settings=settings,
    )


def _ln_gamma(model: ct.NRTL):
    mixture = ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True)

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        return np.log(
            np.array(
                model.activity_coefficients(
                    mixture=mixture,
                    temperature_K=298.15,
                    composition=[float(value) for value in values],
                ),
                dtype=float,
            )
        )

    return ln_gamma


def _psat(temperature_K: float) -> np.ndarray:
    values = []
    for name in NAMES:
        antoine = ct.Component.from_database(name).antoine
        assert antoine is not None, name
        values.append(np.exp(antoine.A - antoine.B / (temperature_K + antoine.C)) * 1.0e5)
    return np.array(values, dtype=float)


# --------------------------------------------------------------------------
# FlashSettings.max_phases
# --------------------------------------------------------------------------


def test_max_phases_defaults_to_three_and_is_validated() -> None:
    assert ct.FlashSettings().max_phases == 3
    assert ct.FlashSettings(max_phases=5).max_phases == 5
    for invalid in (0, -1):
        with pytest.raises(ct.InputRangeError, match="max_phases"):
            ct.FlashSettings(max_phases=invalid)


def test_max_phases_two_reproduces_the_pre_adr_0011_refusal(butanol_water) -> None:
    """The documented behavior of Case R-3 is still reachable, unchanged."""
    with pytest.raises(ct.ConvergenceError, match="third phase is required"):
        _flash(butanol_water, T3_K - 0.05, settings=ct.FlashSettings(max_phases=2))
    # max_phases=1 cannot un-split a feed the stability test proved unstable,
    # so it behaves exactly like max_phases=2 here.
    with pytest.raises(ct.ConvergenceError, match="third phase is required"):
        _flash(butanol_water, T3_K - 0.05, settings=ct.FlashSettings(max_phases=1))


def test_post_split_stability_false_still_returns_the_two_phase_pair(butanol_water) -> None:
    """The escape hatch is not overridden by the search.

    ``post_split_stability=False`` means "do not police the phase set", so it
    returns the converged two-phase answer with the failure in diagnostics and
    never enters the addition/removal loop. This is what Case R-3 inspects.
    """
    result = _flash(
        butanol_water, T3_K - 0.05, settings=ct.FlashSettings(post_split_stability=False)
    )
    assert sorted(result.phase_names()) == ["liquid", "vapor"]
    assert result.diagnostics["post_split_status"] == "unstable"
    assert "phase_set_history" not in result.diagnostics


# --------------------------------------------------------------------------
# Case V-3: the refusal window is resolved by removal
# --------------------------------------------------------------------------


@pytest.mark.parametrize("offset_K", (-0.05, -0.10))
def test_below_t3_the_window_resolves_to_the_two_liquids(butanol_water, offset_K: float) -> None:
    """Case V-3. Achieved: tie-line to 1.4e-16 of the independent binodal.

    Phase fractions 0.470585520652 / 0.529414479348 at both temperatures - the
    tie-line of this model is temperature independent, so the lever rule gives
    the same answer at each.
    """
    result = _flash(butanol_water, T3_K + offset_K)

    assert sorted(result.phase_names()) == ["liquid1", "liquid2"]
    assert result.diagnostics["phase_regime"] == "LLE"
    assert result.diagnostics["phase_count"] == 2
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_set_history"] == "V -> LV -> LLV -> LL"
    assert result.diagnostics["phases_added"] == 1
    assert result.diagnostics["phases_removed"] == 1

    tie_line = sorted(
        float(result.phases[name].composition.fractions[0]) for name in result.phase_names()
    )
    assert tie_line[0] == pytest.approx(BINODAL_X1[0], abs=1e-8)
    assert tie_line[1] == pytest.approx(BINODAL_X1[1], abs=1e-8)

    lever = (FEED[0] - tie_line[0]) / (tie_line[1] - tie_line[0])
    heavier = max(result.phase_names(), key=lambda n: result.phases[n].composition.fractions[0])
    assert result.phase_fractions[heavier] == pytest.approx(lever, abs=1e-10)

    assert float(result.diagnostics["equilibrium_residual"]) < 1e-10
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_stable"] is True


def test_just_above_t3_the_answer_is_a_vapor_liquid_pair(butanol_water) -> None:
    """Case V-3, the control on the other side of T3.

    Case R-3 recorded a *single vapor* two kelvin above T3; 0.05 K above it the
    feed is still between its bubble and dew points, so the model gives a
    vapor-liquid pair. That is checked independently here against modified
    Raoult's law and against the liquid being exactly at its bubble point, both
    written out in this module.
    """
    result = _flash(butanol_water, T3_K + 0.05)

    assert sorted(result.phase_names()) == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert "phase_set_history" not in result.diagnostics
    assert result.diagnostics["post_split_stable"] is True

    ln_gamma = _ln_gamma(butanol_water)
    psat = _psat(T3_K + 0.05)
    x = np.array(result.phases["liquid"].composition.fractions)
    y = np.array(result.phases["vapor"].composition.fractions)

    bubble = float(np.sum(x * np.exp(ln_gamma(x)) * psat / PRESSURE_PA))
    assert bubble == pytest.approx(1.0, abs=1e-10)
    raoult = np.max(np.abs(y * PRESSURE_PA - x * np.exp(ln_gamma(x)) * psat)) / PRESSURE_PA
    assert raoult < 1e-10, raoult

    # And the single liquid is inside the miscibility gap no longer: it is
    # stable against a second liquid.
    liquid = ct.Mixture.from_database(list(NAMES), list(x), normalize=True)
    assert (
        ct.stability_tp(
            liquid,
            temperature_K=T3_K + 0.05,
            pressure_Pa=PRESSURE_PA,
            activity_model=butanol_water,
        ).status
        == "stable"
    )


def test_the_window_answer_is_the_same_as_a_direct_liquid_liquid_flash(butanol_water) -> None:
    """The route does not change the answer.

    Below T3 the two-liquid state is also reachable without any vapor at all,
    by the `gamma-gamma` path. The tie-line the addition/removal search returns
    must be that tie-line, not merely something near it: achieved worst
    composition deviation 1.4e-16.
    """
    through_the_window = _flash(butanol_water, T3_K - 0.05)
    direct = ct.flash_tp(
        ct.Mixture.from_database(list(NAMES), list(FEED), normalize=True),
        temperature_K=T3_K - 0.05,
        pressure_Pa=PRESSURE_PA,
        activity_model=butanol_water,
    )

    assert sorted(direct.phase_names()) == ["liquid1", "liquid2"]
    windowed = sorted(
        tuple(through_the_window.phases[name].composition.fractions)
        for name in through_the_window.phase_names()
    )
    reference = sorted(
        tuple(direct.phases[name].composition.fractions) for name in direct.phase_names()
    )
    assert np.allclose(np.array(windowed), np.array(reference), rtol=0.0, atol=1e-9)


def test_the_window_result_is_deterministic_and_permutation_invariant(butanol_water) -> None:
    temperature_K = T3_K - 0.05
    base = _flash(butanol_water, temperature_K)
    assert dict(_flash(butanol_water, temperature_K).diagnostics) == dict(base.diagnostics)

    swapped_model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)]
        )
    )
    swapped = ct.flash_tp(
        ct.Mixture.from_database([NAMES[1], NAMES[0]], [FEED[1], FEED[0]], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=swapped_model,
        flash_mode="modified-raoult",
    )
    assert sorted(swapped.phase_names()) == ["liquid1", "liquid2"]
    restored = sorted(
        tuple(reversed(swapped.phases[name].composition.fractions))
        for name in swapped.phase_names()
    )
    reference = sorted(
        tuple(base.phases[name].composition.fractions) for name in base.phase_names()
    )
    assert np.allclose(np.array(restored), np.array(reference), rtol=0.0, atol=1e-9)
    assert swapped.diagnostics["phase_set_history"] == base.diagnostics["phase_set_history"]


# --------------------------------------------------------------------------
# The result contract
# --------------------------------------------------------------------------


def test_a_three_phase_result_satisfies_the_flash_result_invariants(
    tessier2000_names, tessier2000_model
) -> None:
    """Case V-1 invariants at the API level, for the ternary tie-triangle."""
    result = ct.flash_tp(
        ct.Mixture.from_database(
            tessier2000_names, [0.13418838, 0.08427618, 0.78153544], normalize=True
        ),
        temperature_K=364.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
        flash_mode="modified-raoult",
    )

    assert sorted(result.phase_names()) == ["liquid1", "liquid2", "vapor"]
    assert set(result.phase_fractions) == set(result.phases)
    total = sum(result.phase_fractions.values())
    assert total == pytest.approx(1.0, abs=1e-15)
    assert all(0.0 <= value <= 1.0 for value in result.phase_fractions.values())
    assert result.vapor_fraction == result.phase_fractions["vapor"]
    for name, phase in result.phases.items():
        assert phase.name == name
        assert sum(phase.composition.fractions) == pytest.approx(1.0, abs=1e-12)

    diagnostics = result.diagnostics
    for key in (
        "phase_count",
        "phase_regime",
        "phase_state",
        "phase_set_history",
        "phases_added",
        "phases_removed",
        "post_split_stable",
        "post_split_status",
        "post_split_tpd_min",
        "equilibrium_residual",
        "mass_balance_residual",
        "delta_g_split_rt",
        "delta_g_vs_two_phase_rt",
        "ssi_iterations",
        "second_order_iterations",
        "rachford_rice_iterations",
        "converged_stage",
    ):
        assert key in diagnostics, key
    for name in result.phase_names():
        assert f"phase_stability_{name}" in diagnostics
        assert f"phase_stability_tpd_min_{name}" in diagnostics
    assert diagnostics["flash_mode"] == "modified-raoult"
    assert diagnostics["phase_detection"] == "tangent-plane"


# --------------------------------------------------------------------------
# Regression: nothing else moved
# --------------------------------------------------------------------------


def test_the_phi_phi_grid_never_reaches_a_third_phase() -> None:
    """Requirement of ADR-0011: `max_phases = 3` changes no phi-phi state.

    No state on the in-repo Peng-Robinson grid needs a third phase (Case L-4),
    so none of them enters the search: `phase_set_history` is absent from every
    one and no state has more than two phases. The stronger statement - that
    every number is bit-identical - is
    `tests/test_flash_refactor_bit_identity.py`, whose fixture predates this
    slice and is unchanged by it.
    """
    eos = ct.PengRobinsonEOS()
    counts: dict[int, int] = {}
    for names, z in GRID_MIXTURES:
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                result = ct.flash_tp(
                    ct.Mixture.from_database(list(names), list(z), normalize=True),
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=eos,
                )
                count = int(result.diagnostics["phase_count"])
                counts[count] = counts.get(count, 0) + 1
                assert count <= 2, (names, temperature_K, pressure_Pa, count)
                assert "phase_set_history" not in result.diagnostics
                assert "phases_added" not in result.diagnostics
    assert sum(counts.values()) == 144
    assert counts == {1: 97, 2: 47}, counts


def test_the_other_paths_still_stop_at_two_phases_whatever_max_phases_says(
    butanol_water,
) -> None:
    """ADR-0011 wires the search to the modified-Raoult path only.

    No in-repo phi-phi or gamma-gamma state needs a third phase (Case L-4), so
    wiring the loop there would ship a path nothing exercises. Both failures
    below are manufactured by loosening the split tolerance, exactly as in
    `tests/test_flash_lle.py::test_post_split_failure_raises_and_post_split_stability_false_returns`;
    what this test pins is that `max_phases` does not change either of them.
    """
    phi_phi = ct.Mixture.from_database(["Ethane", "n-Heptane"], [0.7, 0.3], normalize=True)
    with pytest.raises(ct.ConvergenceError, match="not a stable phase set"):
        ct.flash_tp(
            phi_phi,
            temperature_K=360.0,
            pressure_Pa=1.0e6,
            eos=ct.PengRobinsonEOS(),
            settings=ct.FlashSettings(tol=1e-3, max_phases=4),
        )

    gamma_gamma = ct.Mixture.from_database(list(NAMES), list(FEED), normalize=True)
    with pytest.raises(ct.ConvergenceError, match="not a stable phase set"):
        ct.flash_tp(
            gamma_gamma,
            temperature_K=298.15,
            pressure_Pa=PRESSURE_PA,
            activity_model=butanol_water,
            settings=ct.FlashSettings(tol=1e-3, second_order=False, max_phases=4),
        )
