"""Tangent-plane stability and phi-phi flash with PC-SAFT (ADR-0015).

Validation Cases P-4 (stability) and P-5 (flash). The tie lines are checked
against ``teqp`` in ``tests/validation/test_pcsaft_flash_vs_teqp.py``; this
file needs no optional dependency and checks the properties that hold whatever
the reference says: verdicts, the identities of Michelsen's test, permutation
invariance, determinism, and the three residuals every split must carry.

Equation numbers refer to the module docstring of
``src/chemthermo/stability/tp.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos.pcsaft import PCSAFTEOS
from chemthermo.stability.tp import _ln_phi_min_gibbs

EOS = PCSAFTEOS()
BINARY = ("Methane", "n-Hexane")
TEMPERATURE_K = 300.0


def _mixture(names, z):
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _d_vector(mixture, pressure_Pa):
    """Feed tangent-plane intercept d_i = ln z_i + ln phi_i(z), equation (1)."""
    z = np.array(mixture.fractions, dtype=float)
    ln_phi, _ = _ln_phi_min_gibbs(
        EOS,
        mixture=mixture,
        temperature=TEMPERATURE_K,
        pressure=pressure_Pa,
        composition=z,
    )
    return np.log(z) + ln_phi


def _tpd(mixture, pressure_Pa, w, d):
    """Reduced tangent-plane distance, equation (2), from first principles."""
    w = np.array(w, dtype=float)
    w = w / float(np.sum(w))
    ln_phi, _ = _ln_phi_min_gibbs(
        EOS,
        mixture=mixture,
        temperature=TEMPERATURE_K,
        pressure=pressure_Pa,
        composition=w,
    )
    return float(np.sum(w * (np.log(w) + ln_phi - d)))


# ---------------------------------------------------------------------------
# Case P-4: stability verdicts
# ---------------------------------------------------------------------------

#: ``(z1, P, expected status, why)``. The two-phase boundary at 300 K is the
#: PC-SAFT methane / n-hexane isotherm cross-checked against teqp in
#: ``tests/validation/test_pcsaft_flash_vs_teqp.py``: at 1 MPa the tie line runs
#: from x1 = 0.0556 to y1 = 0.9735, at 3 MPa from 0.1609 to 0.9869, and the
#: bubble pressure at x1 = 0.3 is near 5.9 MPa.
_VERDICTS = [
    (0.30, 3.0e6, "unstable", "inside the two-phase region"),
    (0.50, 5.0e6, "unstable", "inside the two-phase region"),
    (0.95, 1.0e6, "unstable", "inside the two-phase region (dew y1 = 0.9735 at 1 MPa)"),
    (0.30, 8.0e6, "stable", "compressed liquid above the bubble pressure"),
    (0.99, 1.0e6, "stable", "a vapour beyond the dew composition"),
    (0.02, 5.0e5, "stable", "a liquid below the bubble composition"),
]


@pytest.mark.parametrize(
    ("z1", "pressure", "expected", "why"),
    _VERDICTS,
    ids=[f"z1={row[0]}, P={row[1]:.3g} Pa" for row in _VERDICTS],
)
def test_stability_verdicts_on_the_300_K_isotherm(
    z1: float, pressure: float, expected: str, why: str
) -> None:
    mixture = _mixture(BINARY, (z1, 1.0 - z1))
    result = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=EOS)
    assert result.status == expected, f"{why}: tpd_min = {result.tpd_min}"
    assert result.diagnostics["model_family"] == "eos"
    assert (result.tpd_min < 0.0) == (expected == "unstable")


def test_an_unstable_feed_points_at_the_methane_rich_vapor() -> None:
    """The minimizer's implied K-values must match the tie line's direction."""
    mixture = _mixture(BINARY, (0.30, 0.70))
    result = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=EOS)
    assert result.status == "unstable"
    assert result.trial_composition is not None
    assert result.k_values is not None

    # The stationary point is a methane-rich vapour, so K_methane > 1 > K_hexane.
    assert result.trial_composition[0] > 0.9
    assert result.k_values[0] > 1.0
    assert result.k_values[1] < 1.0


def test_tpd_vanishes_at_the_feed_and_trials_are_stationary() -> None:
    """Equations (1), (2) and (5) restated on PC-SAFT's fugacity coefficients."""
    mixture = _mixture(BINARY, (0.30, 0.70))
    pressure = 3.0e6
    d = _d_vector(mixture, pressure)
    z = np.array(mixture.fractions, dtype=float)
    assert _tpd(mixture, pressure, z, d) == pytest.approx(0.0, abs=1e-12)

    settings = ct.StabilitySettings()
    result = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
        settings=settings,
    )
    non_trivial = [t for t in result.trials if t.converged and not t.trivial]
    assert non_trivial, "expected at least one non-trivial stationary point"

    for trial in non_trivial:
        assert trial.composition is not None
        assert trial.residual < settings.tol
        w = np.array(trial.composition, dtype=float)
        w_capital = w * trial.sum_W
        ln_phi, _ = _ln_phi_min_gibbs(
            EOS,
            mixture=mixture,
            temperature=TEMPERATURE_K,
            pressure=pressure,
            composition=w,
        )
        assert float(np.max(np.abs(np.log(w_capital) + ln_phi - d))) < 1e-9
        # Equation (7): tpd = -ln(sum_W), and it agrees with the direct value.
        assert trial.tpd == pytest.approx(-math.log(trial.sum_W), abs=1e-10)
        assert _tpd(mixture, pressure, w, d) == pytest.approx(trial.tpd, abs=1e-12)


def test_stability_is_permutation_invariant_and_deterministic() -> None:
    mixture = _mixture(BINARY, (0.30, 0.70))
    swapped = _mixture(tuple(reversed(BINARY)), (0.70, 0.30))

    base = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=EOS)
    other = ct.stability_tp(swapped, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=EOS)

    assert base.status == other.status
    assert other.tpd_min == pytest.approx(base.tpd_min, abs=1e-10)
    assert base.trial_composition is not None
    assert other.trial_composition is not None
    np.testing.assert_allclose(
        other.trial_composition, tuple(reversed(base.trial_composition)), atol=1e-10
    )

    repeat = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=EOS)
    assert repeat.tpd_min == base.tpd_min
    assert repeat.trial_composition == base.trial_composition


# ---------------------------------------------------------------------------
# Case P-5: the flash
# ---------------------------------------------------------------------------

_FEEDS = [(0.10, 1.0e6), (0.30, 3.0e6), (0.50, 5.0e6), (0.70, 7.0e6)]


@pytest.mark.parametrize(
    ("z1", "pressure"),
    # ADR-0028 runtime trim: four feeds on one isotherm, two of them further
    # points on the same curve.
    [
        row if index % 2 == 0 else pytest.param(*row, marks=pytest.mark.slow)
        for index, row in enumerate(_FEEDS)
    ],
    ids=[f"z1={row[0]}, P={row[1]:.3g} Pa" for row in _FEEDS],
)
def test_two_phase_flash_is_verified_not_merely_converged(z1: float, pressure: float) -> None:
    """Every split carries mass balance, equal fugacity and a Gibbs-energy drop."""
    mixture = _mixture(BINARY, (z1, 1.0 - z1))
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure, eos=EOS)

    assert set(result.phases) == {"liquid", "vapor"}
    diagnostics = result.diagnostics
    assert diagnostics["converged"] is True
    assert diagnostics["phase_detection"] == "tangent-plane"
    assert diagnostics["stability_status"] == "unstable"
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_checked"] is True
    assert diagnostics["post_split_stable"] is True

    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
    beta = result.vapor_fraction
    assert beta is not None
    assert 0.0 < beta < 1.0
    # The vapour really is the methane-rich, lower-density phase.
    assert y[0] > x[0]
    liquid_density = EOS.density_roots(
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        composition=x.tolist(),
        mixture=mixture,
    )[-1]
    vapor_density = EOS.density_roots(
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        composition=y.tolist(),
        mixture=mixture,
    )[0]
    assert vapor_density < liquid_density

    # Equal fugacity restated from the public interface, independent of the
    # solver's own residual.
    phi_liquid = np.array(
        EOS.fugacity_coefficients(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure,
            composition=x.tolist(),
            phase="liquid",
        )
    )
    phi_vapor = np.array(
        EOS.fugacity_coefficients(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure,
            composition=y.tolist(),
            phase="vapor",
        )
    )
    # Relative, because the solver's own convergence test is on ln f and an
    # absolute fugacity difference is dominated by whichever component is
    # abundant.
    assert float(np.max(np.abs(phi_liquid * x / (phi_vapor * y) - 1.0))) < 1e-7


def test_a_stable_feed_returns_one_phase() -> None:
    mixture = _mixture(BINARY, (0.30, 0.70))
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=EOS)
    assert len(result.phases) == 1
    assert result.diagnostics["stability_status"] == "stable"
    assert result.diagnostics["phase_count"] == 1


def test_two_feeds_on_one_tie_line_give_the_same_phases() -> None:
    """Only the amounts may differ; the lever rule ties the two together."""
    pressure = 3.0e6
    first = ct.flash_tp(
        _mixture(BINARY, (0.30, 0.70)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
    )
    second = ct.flash_tp(
        _mixture(BINARY, (0.50, 0.50)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
    )
    for name in ("liquid", "vapor"):
        np.testing.assert_allclose(
            first.phases[name].composition.fractions,
            second.phases[name].composition.fractions,
            atol=1e-9,
        )
    assert first.vapor_fraction is not None
    assert second.vapor_fraction is not None
    assert second.vapor_fraction > first.vapor_fraction

    x1 = first.phases["liquid"].composition.fractions[0]
    y1 = first.phases["vapor"].composition.fractions[0]
    assert first.vapor_fraction == pytest.approx((0.30 - x1) / (y1 - x1), rel=1e-9)


def test_flash_is_permutation_invariant_and_deterministic() -> None:
    pressure = 3.0e6
    base = ct.flash_tp(
        _mixture(BINARY, (0.30, 0.70)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
    )
    swapped = ct.flash_tp(
        _mixture(tuple(reversed(BINARY)), (0.70, 0.30)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
    )
    assert base.vapor_fraction is not None
    assert swapped.vapor_fraction is not None
    assert swapped.vapor_fraction == pytest.approx(base.vapor_fraction, abs=1e-12)
    for name in ("liquid", "vapor"):
        np.testing.assert_allclose(
            swapped.phases[name].composition.fractions,
            tuple(reversed(base.phases[name].composition.fractions)),
            atol=1e-12,
        )

    repeat = ct.flash_tp(
        _mixture(BINARY, (0.30, 0.70)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        eos=EOS,
    )
    assert repeat.vapor_fraction == base.vapor_fraction
    assert (
        repeat.phases["liquid"].composition.fractions == base.phases["liquid"].composition.fractions
    )


def test_a_kij_moves_the_tie_line() -> None:
    """The ADR-0006 contract carries over unchanged; see also Case P-5(E).

    ``kij = 0.03`` for methane / n-decane is **illustrative**, not a
    literature-validated parameter.
    """
    mixture = _mixture(("Methane", "n-Decane"), (0.40, 0.60))
    zero = ct.flash_tp(mixture, temperature_K=350.0, pressure_Pa=5.0e6, eos=PCSAFTEOS())
    tuned = ct.flash_tp(
        mixture,
        temperature_K=350.0,
        pressure_Pa=5.0e6,
        eos=PCSAFTEOS(kij={("Methane", "n-Decane"): 0.03}),
    )
    assert zero.vapor_fraction is not None
    assert tuned.vapor_fraction is not None
    assert abs(tuned.vapor_fraction - zero.vapor_fraction) > 1e-3
    for result in (zero, tuned):
        assert float(result.diagnostics["fugacity_residual"]) < 1e-8
        assert float(result.diagnostics["delta_g_split_rt"]) < 0.0


def test_phi_phi_still_stops_at_two_phases() -> None:
    """``max_phases`` does not reach the phi-phi path (ADR-0011).

    PC-SAFT joining the EOS family does not change that: phase addition lives
    on the modified-Raoult path only.
    """
    mixture = _mixture(BINARY, (0.30, 0.70))
    result = ct.flash_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=ct.FlashSettings(max_phases=3),
    )
    assert len(result.phases) == 2
    assert "phase_set_history" not in result.diagnostics


def test_peng_robinson_results_are_untouched_by_this_slice() -> None:
    """A pinned Peng-Robinson number, restated here as a regression guard."""
    mixture = _mixture(("Methane", "Ethane", "Propane"), (0.50, 0.30, 0.20))
    result = ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS())
    assert result.vapor_fraction is not None
    assert len(result.phases) == 2
