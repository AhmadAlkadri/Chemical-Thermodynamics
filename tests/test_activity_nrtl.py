"""Unit tests for the NRTL activity coefficient model.

The Tessier (2000) Problem 1 parameters are deliberately used for the
thermodynamic-consistency checks: they are strongly asymmetric
(tau_ij != tau_ji, alpha_23 = 0.48 != alpha_12 = 0.3), which is exactly the
regime in which a row-sum/column-sum mix-up in the NRTL equation is visible.
A symmetric binary cannot detect it.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct

PACKAGED_NRTL_JSON = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "chemthermo"
    / "parameters"
    / "data"
    / "activity"
    / "nrtl.json"
)

# Directions along the composition simplex (each sums to zero, so a step keeps
# sum(x) = 1) and interior compositions used for the Gibbs-Duhem check.
GIBBS_DUHEM_DIRECTIONS = (
    (1.0, -1.0, 0.0),
    (1.0, 0.0, -1.0),
    (0.0, 1.0, -1.0),
    (1.0, 1.0, -2.0),
)
GIBBS_DUHEM_COMPOSITIONS = (
    (0.12, 0.08, 0.80),
    (0.50, 0.30, 0.20),
    (0.20, 0.20, 0.60),
    (0.05, 0.05, 0.90),
)


def _ln_gamma_callable(
    model: ct.NRTL, names: Sequence[str], temperature_K: float = 298.15
) -> Callable[[np.ndarray], np.ndarray]:
    mixture = ct.Mixture.from_database(list(names), [1.0 / len(names)] * len(names), normalize=True)

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        gamma = model.activity_coefficients(
            mixture=mixture,
            temperature_K=temperature_K,
            composition=[float(value) for value in x],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma


def test_nrtl_gamma_unity_with_zero_tau() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.4, 0.6), basis="mole", normalize=False),
    )

    params = ct.NRTLParameters.from_pairs([("Methane", "Ethane", 0.0, 0.0, 0.3, 0.3)])
    model = ct.NRTL(parameters=params)

    gamma = model.activity_coefficients(
        mixture=mixture, temperature_K=300.0, composition=mixture.fractions
    )

    assert gamma == pytest.approx((1.0, 1.0))


def test_nrtl_symmetric_pair_yields_equal_gamma() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
    )

    params = ct.NRTLParameters.from_pairs([("Methane", "Ethane", 0.2, 0.2, 0.3, 0.3)])
    model = ct.NRTL(parameters=params)

    gamma = model.activity_coefficients(
        mixture=mixture, temperature_K=300.0, composition=mixture.fractions
    )

    assert gamma[0] == pytest.approx(gamma[1])
    assert gamma[0] > 0.0


def test_nrtl_fully_symmetric_ternary_is_permutation_symmetric() -> None:
    """tau_ij = tau_ji and alpha symmetric at an equimolar ternary.

    With every ordered pair carrying the same tau and alpha, no component is
    distinguishable from any other at x = (1/3, 1/3, 1/3), so all three
    activity coefficients must coincide.
    """
    names = ["Methane", "Ethane", "Propane"]
    params = ct.NRTLParameters.from_pairs(
        [(a, b, 0.35, 0.35, 0.3, 0.3) for a, b in itertools.combinations(names, 2)]
    )
    mixture = ct.Mixture.from_database(names, [1 / 3, 1 / 3, 1 / 3], normalize=True)

    gamma = ct.NRTL(parameters=params).activity_coefficients(
        mixture=mixture, temperature_K=300.0, composition=mixture.fractions
    )

    assert gamma[0] == pytest.approx(gamma[1], rel=1e-14)
    assert gamma[0] == pytest.approx(gamma[2], rel=1e-14)


def test_nrtl_parameter_loader_matches_database_pair() -> None:
    params = ct.NRTLParameters.load()
    tau, alpha = params.for_components(["Methane", "Ethane"])

    assert tau[0, 1] == pytest.approx(0.2)
    assert tau[1, 0] == pytest.approx(0.1)
    assert alpha[0, 1] == pytest.approx(0.3)
    assert alpha[1, 0] == pytest.approx(0.3)


def test_packaged_nrtl_pairs_declare_synthetic_provenance() -> None:
    """The shipped default pairs are placeholders and must say so.

    They were added ad hoc to make the gamma-phi demo/CLI path runnable; they
    are not fitted to data. The loader must also tolerate the extra
    documentation keys without a schema bump.
    """
    with PACKAGED_NRTL_JSON.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    assert payload["schema_version"] == 1
    provenance = payload["provenance"]
    assert provenance["status"] == "synthetic-demo"
    assert provenance["fitted_to_data"] is False
    for entry in payload["pairs"]:
        assert entry["source"] == "synthetic-demo"

    # Extra keys must not break the loader (no schema_version change).
    params = ct.NRTLParameters.load()
    tau, _ = params.for_components(["Benzene", "Water"])
    assert tau[0, 1] == pytest.approx(3.0)


def test_nrtl_binary_reduces_to_textbook_two_component_formulas() -> None:
    """Independent two-component NRTL expressions (Renon & Prausnitz 1968).

    ln g1 = x2^2 [ tau21 (G21/(x1 + x2 G21))^2 + tau12 G12/(x2 + x1 G12)^2 ]
    ln g2 = x1^2 [ tau12 (G12/(x2 + x1 G12))^2 + tau21 G21/(x1 + x2 G21)^2 ]

    These are written out here from the binary form, not derived from the
    multicomponent code under test.
    """
    tau_12, tau_21 = 0.6, -0.35
    alpha_12 = 0.3
    g_12 = float(np.exp(-alpha_12 * tau_12))
    g_21 = float(np.exp(-alpha_12 * tau_21))

    params = ct.NRTLParameters.from_pairs(
        [("Methane", "Ethane", tau_12, tau_21, alpha_12, alpha_12)]
    )
    model = ct.NRTL(parameters=params)
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)

    for x1 in (0.01, 0.1, 0.3, 0.5, 0.75, 0.99):
        x2 = 1.0 - x1
        gamma = model.activity_coefficients(
            mixture=mixture, temperature_K=310.0, composition=[x1, x2]
        )
        expected_1 = x2**2 * (
            tau_21 * (g_21 / (x1 + x2 * g_21)) ** 2 + tau_12 * g_12 / (x2 + x1 * g_12) ** 2
        )
        expected_2 = x1**2 * (
            tau_12 * (g_12 / (x2 + x1 * g_12)) ** 2 + tau_21 * g_21 / (x1 + x2 * g_21) ** 2
        )
        assert float(np.log(gamma[0])) == pytest.approx(expected_1, rel=0.0, abs=1e-14)
        assert float(np.log(gamma[1])) == pytest.approx(expected_2, rel=0.0, abs=1e-14)


def test_nrtl_satisfies_gibbs_duhem_for_asymmetric_parameters(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """sum_i x_i d(ln gamma_i)/ds = 0 along directions in the simplex.

    This is the isothermal-isobaric Gibbs-Duhem relation. It holds for any
    ln gamma that is the composition derivative of a single g^E function and
    fails otherwise, so it is the sharpest cheap check that the NRTL equation
    is written correctly. The pre-fix implementation returned residuals of
    order 1e-1 here.
    """
    ln_gamma = _ln_gamma_callable(tessier2000_model, tessier2000_names)
    step = 1e-6
    worst = 0.0

    for composition in GIBBS_DUHEM_COMPOSITIONS:
        x = np.asarray(composition, dtype=float)
        for direction in GIBBS_DUHEM_DIRECTIONS:
            e = np.asarray(direction, dtype=float)
            derivative = (ln_gamma(x + step * e) - ln_gamma(x - step * e)) / (2.0 * step)
            residual = float(np.sum(x * derivative))
            worst = max(worst, abs(residual))
            assert abs(residual) < 1e-7, (
                f"Gibbs-Duhem residual {residual:.3e} at x={composition} along {direction}"
            )

    # Achieved (double precision, step 1e-6): worst residual ~2.1e-10.
    assert worst < 1e-8


def test_nrtl_is_permutation_invariant(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """Reordering components permutes ln gamma and nothing else."""
    base_x = (0.12, 0.08, 0.80)
    ln_gamma = _ln_gamma_callable(tessier2000_model, tessier2000_names)
    reference = ln_gamma(np.asarray(base_x, dtype=float))

    for order in itertools.permutations(range(3)):
        names = [tessier2000_names[i] for i in order]
        permuted_ln_gamma = _ln_gamma_callable(tessier2000_model, names)
        value = permuted_ln_gamma(np.asarray([base_x[i] for i in order], dtype=float))
        expected = np.asarray([reference[i] for i in order], dtype=float)
        assert np.allclose(value, expected, rtol=1e-12, atol=0.0)


def test_nrtl_regression_guard_against_row_sum_bug(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """Hard-coded ln gamma from the standard Renon-Prausnitz equation.

    Values produced by an independent implementation of

        ln gamma_i = C_i/S_i + sum_j x_j G_ij/S_j (tau_ij - C_j/S_j),
        S_j = sum_k G_kj x_k,  C_j = sum_k tau_kj G_kj x_k,

    and cross-checked against `thermo.NRTL` (agreement 8.9e-16; see
    `tests/validation/test_nrtl_tessier2000.py`). The pre-fix implementation
    returned (0.857845, 1.045984, 0.174959) here, i.e. off by up to 0.113.
    """
    ln_gamma = _ln_gamma_callable(tessier2000_model, tessier2000_names)
    value = ln_gamma(np.asarray([0.12, 0.08, 0.80], dtype=float))
    expected = np.asarray(
        [0.9120183964066216, 1.1593117704861213, 0.18497649540992964], dtype=float
    )
    assert np.allclose(value, expected, rtol=0.0, atol=1e-9)
