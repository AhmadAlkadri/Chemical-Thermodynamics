"""Validation Case F-4: phi-phi split robustness with PC-SAFT (ADR-0016).

Two independent things are checked.

**The grid.** 188 PC-SAFT states - carbon dioxide / n-decane and methane /
n-hexane, over feed composition, temperature and pressure - are flashed on the
tangent-plane path. None may raise, every two-phase answer must carry the four
invariants (mass balance, equal fugacities, a negative Gibbs-energy change and
a stable phase set) and every single-phase answer must be a stability verdict.
The four states that raised before ADR-0016 are named explicitly so the count
is a measurement and not a moving target.

**The cross-check.** At the reference state (240 K, 1.0 MPa) ``teqp``'s *own*
``get_fugacity_coefficients``, evaluated at chemthermo's converged compositions
and densities, must give the two phases equal fugacities. That uses no
reference tie line: it asks a different implementation of the same published
model, with no shared derivative code, whether chemthermo's answer is an
equilibrium. Skipped when ``teqp`` is not installed
(``pip install -e ".[validation]"``).
"""

from __future__ import annotations

import numpy as np
import pytest

import chemthermo as ct

#: Gross & Sadowski (2001) Table 1, written out here so the reference model is
#: built from values in this file rather than from the package under test.
_PARAMETERS: dict[str, tuple[float, float, float]] = {
    "Methane": (1.0000, 3.7039, 150.03),
    "n-Hexane": (3.0576, 3.7983, 236.77),
    "n-Decane": (4.6627, 3.8384, 243.87),
    "Carbon dioxide": (2.0729, 2.7852, 169.21),
}

#: The Case F-4 grid: (component pair, feed fractions of the first component,
#: temperatures in K, pressures in Pa).
GRID: tuple[
    tuple[tuple[str, str], tuple[float, ...], tuple[float, ...], tuple[float, ...]], ...
] = (
    (
        ("Carbon dioxide", "n-Decane"),
        (0.6, 0.8, 0.9),
        (230.0, 240.0, 250.0, 260.0),
        (1.0e6, 1.5e6, 2.0e6, 2.5e6),
    ),
    (
        ("Methane", "n-Hexane"),
        (0.5, 0.8, 0.9, 0.95),
        (170.0, 180.0, 190.0, 195.0, 200.0),
        (0.5e6, 1.0e6, 1.5e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6),
    ),
)

#: States of the grid that raised ``ConvergenceError`` before ADR-0016, with
#: the "Rachford-Rice failed to bracket a vapor fraction" message. Measured on
#: this repository at commit cf846fe; all four are still failures on the legacy
#: ``phase_detection="wilson-heuristic"`` path, which ADR-0016 left alone, and
#: that is what the last test below re-measures rather than trusting this list.
PREVIOUSLY_FAILING = (
    (("Carbon dioxide", "n-Decane"), 0.8, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 260.0, 1.5e6),
)

REFERENCE_COMPONENTS = ("Carbon dioxide", "n-Decane")
REFERENCE_FEED = (0.9, 0.1)
REFERENCE_T_K = 240.0
REFERENCE_P_PA = 1.0e6

#: Relative agreement required of teqp's own equal-fugacity check, matching
#: `tests/validation/test_pcsaft_flash_vs_teqp.py`.
_FUGACITY_RTOL = 1e-8


def _grid_states():
    for components, feeds, temperatures, pressures in GRID:
        for z1 in feeds:
            for temperature_K in temperatures:
                for pressure_Pa in pressures:
                    yield components, z1, temperature_K, pressure_Pa


def _flash(
    components: tuple[str, str],
    z1: float,
    temperature_K: float,
    pressure_Pa: float,
    settings: ct.FlashSettings | None = None,
) -> ct.FlashResult:
    mixture = ct.Mixture.from_database(list(components), [z1, 1.0 - z1], normalize=True)
    return ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PCSAFTEOS(),
        settings=settings,
    )


# ---------------------------------------------------------------------------
# The grid.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_the_whole_grid_answers_and_every_answer_is_verified() -> None:
    """Case F-4: 0 `ConvergenceError`s, and every answer carries its invariants.

    The full 188-state grid; deselected by default (`pyproject.toml`
    `addopts = "-m 'not slow'"`), run explicitly with `pytest -q -m slow`. A
    24-state representative subset runs in CI's default `pytest -q` instead -
    see `tests/validation/test_flash_split_robustness_pcsaft_subset.py`.
    """
    failures: list[str] = []
    two_phase = 0
    single_phase = 0
    rescued = 0
    worst_mass_balance = 0.0
    worst_fugacity = 0.0
    worst_delta_g = -np.inf

    for components, z1, temperature_K, pressure_Pa in _grid_states():
        label = f"{'/'.join(components)} z1={z1} T={temperature_K} P={pressure_Pa}"
        try:
            result = _flash(components, z1, temperature_K, pressure_Pa)
        except ct.ConvergenceError as error:
            failures.append(f"{label}: {error}")
            continue

        diagnostics = result.diagnostics
        if diagnostics["phase_count"] == 1:
            single_phase += 1
            # A single phase on this path is a stability verdict, never a
            # heuristic and never a fallback.
            assert diagnostics["stability_status"] == "stable", label
            assert diagnostics["termination_reason"] == "feed_stable_tangent_plane", label
            continue

        two_phase += 1
        beta = result.vapor_fraction
        assert beta is not None and 0.0 < beta < 1.0, label
        mass_balance = float(diagnostics["mass_balance_residual"])
        fugacity = float(diagnostics["fugacity_residual"])
        delta_g = float(diagnostics["delta_g_split_rt"])
        assert mass_balance < 1e-12, label
        assert fugacity < 1e-6, label
        assert delta_g < 0.0, label
        assert diagnostics["post_split_status"] == "stable", label
        worst_mass_balance = max(worst_mass_balance, mass_balance)
        worst_fugacity = max(worst_fugacity, fugacity)
        worst_delta_g = max(worst_delta_g, delta_g)
        if diagnostics.get("converged_stage") == "second-order":
            rescued += 1

    assert not failures, "\n".join(failures)
    assert two_phase + single_phase == 188, (two_phase, single_phase)
    # Measured on this repository: 123 two-phase, 65 single-phase, and the four
    # PREVIOUSLY_FAILING states are exactly the ones that need the stage.
    assert two_phase == 123, two_phase
    assert single_phase == 65, single_phase
    assert rescued == len(PREVIOUSLY_FAILING), rescued
    assert worst_mass_balance < 1e-12
    assert worst_fugacity < 1e-6
    assert worst_delta_g < 0.0


@pytest.mark.parametrize(("components", "z1", "temperature_K", "pressure_Pa"), PREVIOUSLY_FAILING)
@pytest.mark.slow  # ADR-0020 runtime trim: the same four states run in tests/test_flash_phi_phi_second_order.py
def test_each_previously_failing_state_needs_and_gets_the_second_order_stage(
    components: tuple[str, str], z1: float, temperature_K: float, pressure_Pa: float
) -> None:
    """Every one of the four is rescued, and fails again without the stage."""
    result = _flash(components, z1, temperature_K, pressure_Pa)
    assert result.diagnostics["converged_stage"] == "second-order"
    assert int(result.diagnostics["negative_flash_steps"]) > 0
    assert result.diagnostics["post_split_status"] == "stable"

    with pytest.raises(ct.ConvergenceError):
        _flash(
            components,
            z1,
            temperature_K,
            pressure_Pa,
            settings=ct.FlashSettings(second_order=False),
        )
    with pytest.raises(ct.ConvergenceError, match="failed to bracket"):
        _flash(
            components,
            z1,
            temperature_K,
            pressure_Pa,
            settings=ct.FlashSettings(phase_detection="wilson-heuristic"),
        )


# ---------------------------------------------------------------------------
# The teqp cross-check at the reference state.
# ---------------------------------------------------------------------------


def _reference_model(components, kij_matrix):
    teqp = pytest.importorskip("teqp")
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
        {"kind": "PCSAFT", "model": {"coeffs": coefficients, "kmat": kij_matrix.tolist()}}
    )


def test_teqp_calls_the_rescued_split_an_equilibrium() -> None:
    """teqp's own fugacity coefficients at chemthermo's phases and densities.

    No reference tie line is used. The two phases' fugacities must be equal
    *in teqp's model*, which is what makes this a check of the answer rather
    than of the route to it.
    """
    pytest.importorskip("teqp")
    mixture = ct.Mixture.from_database(
        list(REFERENCE_COMPONENTS), list(REFERENCE_FEED), normalize=True
    )
    eos = ct.PCSAFTEOS()
    result = ct.flash_tp(
        mixture,
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        eos=ct.PCSAFTEOS(),
    )
    assert result.diagnostics["converged_stage"] == "second-order"

    model = _reference_model(REFERENCE_COMPONENTS, np.zeros((2, 2)))
    fugacities = []
    for name, root in (("liquid", -1), ("vapor", 0)):
        composition = np.array(result.phases[name].composition.fractions, dtype=float)
        roots = eos.density_roots(
            temperature_K=REFERENCE_T_K,
            pressure_Pa=REFERENCE_P_PA,
            composition=composition.tolist(),
            mixture=mixture,
        )
        density = float(roots[root])
        phi = np.asarray(model.get_fugacity_coefficients(REFERENCE_T_K, density * composition))
        fugacities.append(phi * composition * REFERENCE_P_PA)

    liquid, vapor = fugacities
    mismatch = float(np.max(np.abs(liquid / vapor - 1.0)))
    assert mismatch < _FUGACITY_RTOL, mismatch


def test_the_teqp_cross_check_is_not_vacuous() -> None:
    """Negative control: a 1 % change in one sigma must break the agreement."""
    pytest.importorskip("teqp")
    mixture = ct.Mixture.from_database(
        list(REFERENCE_COMPONENTS), list(REFERENCE_FEED), normalize=True
    )
    eos = ct.PCSAFTEOS()
    result = ct.flash_tp(
        mixture,
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        eos=ct.PCSAFTEOS(),
    )

    teqp = pytest.importorskip("teqp")
    perturbed = teqp.make_model(
        {
            "kind": "PCSAFT",
            "model": {
                "coeffs": [
                    {
                        "name": name,
                        "m": _PARAMETERS[name][0],
                        "sigma_Angstrom": _PARAMETERS[name][1] * (1.01 if index == 0 else 1.0),
                        "epsilon_over_k": _PARAMETERS[name][2],
                        "BibTeXKey": "perturbed",
                    }
                    for index, name in enumerate(REFERENCE_COMPONENTS)
                ],
                "kmat": np.zeros((2, 2)).tolist(),
            },
        }
    )
    fugacities = []
    for name, root in (("liquid", -1), ("vapor", 0)):
        composition = np.array(result.phases[name].composition.fractions, dtype=float)
        roots = eos.density_roots(
            temperature_K=REFERENCE_T_K,
            pressure_Pa=REFERENCE_P_PA,
            composition=composition.tolist(),
            mixture=mixture,
        )
        phi = np.asarray(
            perturbed.get_fugacity_coefficients(REFERENCE_T_K, float(roots[root]) * composition)
        )
        fugacities.append(phi * composition * REFERENCE_P_PA)
    liquid, vapor = fugacities
    assert float(np.max(np.abs(liquid / vapor - 1.0))) > 1e-3
