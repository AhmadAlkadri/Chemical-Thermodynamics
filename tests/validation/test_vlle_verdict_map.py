"""Validation Case V-5: the verdict map of the ternary VLLE region.

System
------
1-Propanol(1) / n-Butanol(2) / Water(3) at P = 101325 Pa, modified Raoult: an
NRTL liquid with Antoine pure-liquid reference fugacities against an ideal
vapor (ADR-0010). Parameters are Table 1 of S. R. Tessier, J. F. Brennecke and
M. A. Stadtherr, Chem. Eng. Sci. 55 (2000) 1785-1796, loaded from
`tests/fixtures/nrtl/tessier2000_problem1.json`; Antoine coefficients from the
packaged databank (Koretsky 2012). **These parameters were fitted to
liquid-liquid data and no experimental ternary VLLE data is used here.**
Nothing in this module is a comparison against measurement: what is validated
is that `flash_tp` returns the phase state *this model* has, over a whole grid
of feeds rather than at a handful of hand-picked ones.

Why a map and not more points
-----------------------------
Case V-2 recorded a feed whose three-phase answer the deterministic stability
trial set missed at 363 K. The fix (ADR-0012: one fixed phase-candidate surface
per stability trial) changes *which* stationary points are reachable, so the
honest check is not "the one feed that used to fail now passes" but "the whole
verdict map is right". A grid answered one feed at a time by an independent
route, compared against `flash_tp` feed by feed, is what this module is.

The independent route, written here and sharing no code with chemthermo.flash
-----------------------------------------------------------------------------
For each temperature:

1. the **tie-triangle** from the six-equation damped Newton solve of
   `tests/validation/test_vlle_water_propanol_butanol.py` - three equal
   activities, two normalizations, and the bubble condition on *one* liquid;
2. a **vapor-liquid** state from a four-equation Newton solve of
   ``z_i - (1 - b) x_i - b K_i(x) x_i = 0`` with ``K_i = gamma_i(x) Psat_i / P``
   and ``sum_i x_i = 1``;
3. a **liquid-liquid** state from a seven-equation Newton solve of the three
   equal activities, the two normalizations and two mass balances;
4. the **single-phase** state, the feed itself on whichever candidate has the
   lower Gibbs energy.

Every candidate that exists is scored with the reduced molar Gibbs energy

    G/RT = sum_j beta_j sum_i x_i^j [ ln x_i^j + t_i^j ]

(``t_i = ln gamma_i + ln(Psat_i/P)`` for a liquid, ``0`` for the vapor) and the
**lowest-Gibbs feasible candidate is the expected answer**. That is the whole
adjudication rule: equal fugacities alone are satisfied by more than one state,
and the equilibrium state is the one of least Gibbs energy.

The two geometric statements the slice claims are then consequences that this
module checks rather than assumes: a feed with all-positive barycentric weights
in the tie-triangle comes out three-phase, and a feed outside it does not.

The grid
--------
Per temperature, ``21`` feeds strictly inside the triangle (a barycentric
lattice with all weights >= 1/8) and every point of a ``1/12`` mole-fraction
lattice that is not inside or within ``1e-3`` of the triangle, which is 54-55
feeds outside it depending on the temperature. Deterministic, and 75-76 feeds
per temperature in total.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct

PRESSURE_PA = 101325.0
TEMPERATURES = (363.0, 364.0, 365.0)
FIXTURE_PATH = (
    Path(__file__).resolve().parents[1] / "fixtures" / "nrtl" / "tessier2000_problem1.json"
)

#: Newton tolerance of every independent solve in this module.
NEWTON_TOL = 1e-13
#: Tie-triangle tolerance (achieved <= 1.2e-15).
TRIANGLE_TOL = 1e-14
#: Newton starting points for the tie-triangle, one per temperature. Shared with
#: `tests/validation/test_vlle_water_propanol_butanol.py`.
TRIANGLE_SEEDS: dict[float, tuple[float, ...]] = {
    365.0: (0.023616, 0.024407, 0.951977, 0.098737, 0.201955, 0.699308),
    364.0: (0.052144, 0.029172, 0.918684, 0.139639, 0.123478, 0.736883),
    363.0: (0.102828, 0.035390, 0.861782, 0.156303, 0.064227, 0.779470),
}

#: Barycentric lattice denominator for the feeds inside the triangle.
INSIDE_LATTICE = 8
#: Mole-fraction lattice denominator for the feeds outside it.
OUTSIDE_LATTICE = 12
#: How far outside the triangle a lattice point must be to count as "outside".
OUTSIDE_MARGIN = 1e-3

LnGamma = Callable[[np.ndarray], np.ndarray]
Psat = Callable[[float], np.ndarray]


# --------------------------------------------------------------------------
# The model, written here from the fixture, not taken from chemthermo.models
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def system() -> tuple[list[str], ct.NRTL, LnGamma, Psat]:
    with FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]
    tau = np.array(payload["tau"], dtype=float)
    alpha = np.array(payload["alpha"], dtype=float)
    capital_g = np.exp(-alpha * tau)

    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [
                (
                    names[i],
                    names[j],
                    float(tau[i][j]),
                    float(tau[j][i]),
                    float(alpha[i][j]),
                    float(alpha[j][i]),
                )
                for i in range(len(names))
                for j in range(i + 1, len(names))
            ]
        )
    )

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        """Renon-Prausnitz NRTL with column sums, written out here."""
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        s = capital_g.T @ values
        c = (tau * capital_g).T @ values
        return np.array(
            [
                c[i] / s[i] + float(np.sum(values * capital_g[i, :] / s * (tau[i, :] - c / s)))
                for i in range(len(names))
            ]
        )

    antoine = []
    for name in names:
        record = ct.Component.from_database(name).antoine
        assert record is not None, name
        antoine.append(record)

    def psat(temperature_K: float) -> np.ndarray:
        """``ln(P^sat / bar) = A - B / (T + C)`` from the packaged databank."""
        return np.array(
            [np.exp(a.A - a.B / (temperature_K + a.C)) * 1.0e5 for a in antoine],
            dtype=float,
        )

    return names, model, ln_gamma, psat


# --------------------------------------------------------------------------
# Independent numerics
# --------------------------------------------------------------------------


def _newton(
    residual: Callable[[np.ndarray], np.ndarray],
    start: np.ndarray,
    *,
    positive: Sequence[int],
    tol: float = NEWTON_TOL,
    max_iter: int = 200,
) -> tuple[np.ndarray, float]:
    """Damped Newton with a forward-difference Jacobian and a line search.

    ``positive`` lists the indices that must stay strictly positive (mole
    fractions, and the phase fraction where one is being solved for); the line
    search halves until they do and until the max-norm residual decreases.
    """
    u = np.array(start, dtype=float)
    step = 1e-8
    for _ in range(max_iter):
        f = residual(u)
        norm = float(np.max(np.abs(f)))
        if not np.isfinite(norm):
            return u, float("inf")
        if norm < tol:
            return u, norm
        jacobian = np.zeros((u.size, u.size), dtype=float)
        for column in range(u.size):
            forward = u.copy()
            forward[column] += step
            jacobian[:, column] = (residual(forward) - f) / step
        try:
            direction = np.linalg.solve(jacobian, -f)
        except np.linalg.LinAlgError:
            return u, norm
        if not np.all(np.isfinite(direction)):
            return u, norm
        scale = 1.0
        for _ in range(60):
            candidate = u + scale * direction
            if np.any(candidate[list(positive)] <= 0.0):
                scale *= 0.5
                continue
            trial = float(np.max(np.abs(residual(candidate))))
            if np.isfinite(trial) and trial < norm:
                u = candidate
                break
            scale *= 0.5
        else:
            return u, norm
    return u, float(np.max(np.abs(residual(u))))


def _liquid_terms(x: np.ndarray, temperature_K: float, ln_gamma: LnGamma, psat: Psat) -> np.ndarray:
    return ln_gamma(x) + np.log(psat(temperature_K) / PRESSURE_PA)


def _reduced_g(
    x: np.ndarray, is_vapor: bool, temperature_K: float, ln_gamma: LnGamma, psat: Psat
) -> float:
    """``sum_i x_i (ln x_i + t_i)``, the composition-dependent part of ``G/RT``."""
    values = np.asarray(x, dtype=float)
    terms = (
        np.zeros_like(values) if is_vapor else _liquid_terms(values, temperature_K, ln_gamma, psat)
    )
    mask = values > 0.0
    return float(np.sum(values[mask] * (np.log(values[mask]) + terms[mask])))


def _tie_triangle(
    temperature_K: float, ln_gamma: LnGamma, psat: Psat
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """The six-equation solve of Case V-1, repeated here."""

    def residual(u: np.ndarray) -> np.ndarray:
        x_i, x_ii = u[:3], u[3:]
        activity_i = np.log(np.abs(x_i)) + ln_gamma(x_i)
        activity_ii = np.log(np.abs(x_ii)) + ln_gamma(x_ii)
        y = np.exp(activity_i) * psat(temperature_K) / PRESSURE_PA
        return np.concatenate(
            [
                activity_i - activity_ii,
                [float(np.sum(x_i)) - 1.0, float(np.sum(x_ii)) - 1.0, float(np.sum(y)) - 1.0],
            ]
        )

    u, norm = _newton(
        residual,
        np.array(TRIANGLE_SEEDS[temperature_K], dtype=float),
        positive=range(6),
        tol=TRIANGLE_TOL,
    )
    x_i, x_ii = u[:3], u[3:]
    y = np.exp(np.log(x_i) + ln_gamma(x_i)) * psat(temperature_K) / PRESSURE_PA
    return x_i, x_ii, y, norm


def _vapor_liquid(
    z: np.ndarray, temperature_K: float, seeds: Sequence[np.ndarray], ln_gamma: LnGamma, psat: Psat
) -> tuple[np.ndarray, np.ndarray, float, float] | None:
    """``z_i = (1 - b) x_i + b K_i(x) x_i`` with ``sum x = 1``, in ``(x, b)``."""
    saturation = psat(temperature_K)

    def residual(u: np.ndarray) -> np.ndarray:
        x = u[:3]
        beta = u[3]
        k = np.exp(ln_gamma(x)) * saturation / PRESSURE_PA
        return np.concatenate([z - (1.0 - beta) * x - beta * k * x, [float(np.sum(x)) - 1.0]])

    for seed in seeds:
        u, norm = _newton(residual, np.concatenate([seed, [0.5]]), positive=range(4))
        x = u[:3]
        beta = float(u[3])
        if norm >= 1e-12 or not (1e-8 < beta < 1.0 - 1e-8) or not np.all(x > 0.0):
            continue
        y = np.exp(ln_gamma(x)) * saturation / PRESSURE_PA * x
        y = y / float(np.sum(y))
        if float(np.max(np.abs(y - x))) > 1e-7:
            return x, y, beta, norm
    return None


def _liquid_liquid(
    z: np.ndarray,
    temperature_K: float,
    seeds: Sequence[tuple[np.ndarray, np.ndarray]],
    ln_gamma: LnGamma,
    psat: Psat,
) -> tuple[np.ndarray, np.ndarray, float, float] | None:
    """Three equal activities, two normalizations and two mass balances."""

    def residual(u: np.ndarray) -> np.ndarray:
        x_i, x_ii, beta = u[:3], u[3:6], u[6]
        activity_i = np.log(x_i) + _liquid_terms(x_i, temperature_K, ln_gamma, psat)
        activity_ii = np.log(x_ii) + _liquid_terms(x_ii, temperature_K, ln_gamma, psat)
        balance = z - (1.0 - beta) * x_i - beta * x_ii
        return np.concatenate(
            [
                activity_i - activity_ii,
                [float(np.sum(x_i)) - 1.0, float(np.sum(x_ii)) - 1.0],
                balance[:2],
            ]
        )

    for first, second in seeds:
        u, norm = _newton(residual, np.concatenate([first, second, [0.5]]), positive=range(7))
        x_i, x_ii, beta = u[:3], u[3:6], float(u[6])
        if norm >= 1e-12 or not (1e-8 < beta < 1.0 - 1e-8):
            continue
        if not (np.all(x_i > 0.0) and np.all(x_ii > 0.0)):
            continue
        if float(np.max(np.abs(x_i - x_ii))) <= 1e-6:
            continue
        balance = float(np.max(np.abs(z - (1.0 - beta) * x_i - beta * x_ii)))
        if balance < 1e-10:
            return x_i, x_ii, beta, max(norm, balance)
    return None


# --------------------------------------------------------------------------
# The independent classifier
# --------------------------------------------------------------------------


class _Candidate:
    """One admissible phase state of a feed, with its reduced Gibbs energy."""

    def __init__(self, kind: str, phases: int, g_rt: float, detail: object) -> None:
        self.kind = kind
        self.phases = phases
        self.g_rt = g_rt
        self.detail = detail


def _classify(
    z: np.ndarray,
    temperature_K: float,
    triangle: tuple[np.ndarray, np.ndarray, np.ndarray],
    ln_gamma: LnGamma,
    psat: Psat,
) -> tuple[_Candidate, list[_Candidate], np.ndarray]:
    """Return the lowest-Gibbs admissible state of ``z``, and all of them."""
    x_i, x_ii, y = triangle
    vertices = np.column_stack([x_i, x_ii, y])
    weights = np.linalg.solve(vertices, z)

    candidates: list[_Candidate] = []
    if np.all(weights > 0.0):
        g_three = float(
            weights[0] * _reduced_g(x_i, False, temperature_K, ln_gamma, psat)
            + weights[1] * _reduced_g(x_ii, False, temperature_K, ln_gamma, psat)
            + weights[2] * _reduced_g(y, True, temperature_K, ln_gamma, psat)
        )
        candidates.append(_Candidate("VLL", 3, g_three, weights))

    g_liquid = _reduced_g(z, False, temperature_K, ln_gamma, psat)
    g_vapor = _reduced_g(z, True, temperature_K, ln_gamma, psat)
    candidates.append(
        _Candidate("V" if g_vapor < g_liquid else "L", 1, min(g_liquid, g_vapor), None)
    )

    ideal_x = z * PRESSURE_PA / psat(temperature_K)
    ideal_x = ideal_x / float(np.sum(ideal_x))
    vapor_liquid = _vapor_liquid(z, temperature_K, (z, ideal_x, x_i, x_ii, y), ln_gamma, psat)
    if vapor_liquid is not None:
        x, y_vl, beta, _residual = vapor_liquid
        candidates.append(
            _Candidate(
                "VL",
                2,
                (1.0 - beta) * _reduced_g(x, False, temperature_K, ln_gamma, psat)
                + beta * _reduced_g(y_vl, True, temperature_K, ln_gamma, psat),
                vapor_liquid,
            )
        )

    liquid_liquid = _liquid_liquid(z, temperature_K, ((x_i, x_ii), (x_ii, x_i)), ln_gamma, psat)
    if liquid_liquid is not None:
        first, second, beta, _residual = liquid_liquid
        candidates.append(
            _Candidate(
                "LL",
                2,
                (1.0 - beta) * _reduced_g(first, False, temperature_K, ln_gamma, psat)
                + beta * _reduced_g(second, False, temperature_K, ln_gamma, psat),
                liquid_liquid,
            )
        )

    candidates.sort(key=lambda candidate: candidate.g_rt)
    return candidates[0], candidates, weights


def _grid(triangle: tuple[np.ndarray, np.ndarray, np.ndarray]) -> list[tuple[str, np.ndarray]]:
    """A deterministic feed grid: a barycentric lattice inside, a mole-fraction one outside."""
    vertices = np.column_stack(list(triangle))
    feeds: list[tuple[str, np.ndarray]] = []

    n = INSIDE_LATTICE
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            if k < 1:
                continue
            z = vertices @ (np.array([i, j, k], dtype=float) / n)
            feeds.append(("inside", z / float(np.sum(z))))

    m = OUTSIDE_LATTICE
    for i in range(1, m):
        for j in range(1, m - i):
            k = m - i - j
            if k < 1:
                continue
            z = np.array([i, j, k], dtype=float) / m
            if np.all(np.linalg.solve(vertices, z) > -OUTSIDE_MARGIN):
                continue
            feeds.append(("outside", z))
    return feeds


# --------------------------------------------------------------------------
# The exam
# --------------------------------------------------------------------------


class _Row:
    """One feed of the map: what was expected, what `flash_tp` returned."""

    def __init__(
        self,
        region: str,
        z: np.ndarray,
        expected: _Candidate,
        weights: np.ndarray,
        result: ct.FlashResult,
    ) -> None:
        self.region = region
        self.z = z
        self.expected = expected
        self.weights = weights
        self.result = result
        self.obtained = len(result.phase_names())


_CACHE: dict[float, tuple[tuple[np.ndarray, np.ndarray, np.ndarray], float, list[_Row]]] = {}


def _exam(
    temperature_K: float, system: tuple[list[str], ct.NRTL, LnGamma, Psat]
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], float, list[_Row]]:
    if temperature_K in _CACHE:
        return _CACHE[temperature_K]
    names, model, ln_gamma, psat = system
    x_i, x_ii, y, norm = _tie_triangle(temperature_K, ln_gamma, psat)
    triangle = (x_i, x_ii, y)

    rows: list[_Row] = []
    for region, z in _grid(triangle):
        expected, _all, weights = _classify(z, temperature_K, triangle, ln_gamma, psat)
        result = ct.flash_tp(
            ct.Mixture.from_database(names, [float(value) for value in z], normalize=True),
            temperature_K=temperature_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
            flash_mode="modified-raoult",
        )
        rows.append(_Row(region, z, expected, weights, result))

    _CACHE[temperature_K] = (triangle, norm, rows)
    return _CACHE[temperature_K]


def _confusion(rows: Sequence[_Row]) -> dict[tuple[int, int], int]:
    counts: dict[tuple[int, int], int] = {}
    for row in rows:
        key = (row.expected.phases, row.obtained)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_the_independent_tie_triangle_converges(temperature_K: float, system) -> None:
    _triangle, norm, _rows = _exam(temperature_K, system)
    assert norm < TRIANGLE_TOL, norm


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_the_grid_is_large_and_spans_both_regions(temperature_K: float, system) -> None:
    _triangle, _norm, rows = _exam(temperature_K, system)
    assert len(rows) >= 60, len(rows)
    inside = [row for row in rows if row.region == "inside"]
    outside = [row for row in rows if row.region == "outside"]
    assert len(inside) == 21, len(inside)
    assert len(outside) >= 39, len(outside)
    # "Inside" is a statement about the reference triangle, checked here.
    for row in inside:
        assert np.all(row.weights > 0.0), (row.z, row.weights)
    for row in outside:
        assert np.any(row.weights <= -OUTSIDE_MARGIN), (row.z, row.weights)


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_the_verdict_map_agrees_with_the_independent_classifier(
    temperature_K: float, system
) -> None:
    """Zero disagreements between the lowest-Gibbs candidate and `flash_tp`.

    The confusion matrix is keyed ``(expected phase count, obtained phase
    count)``; only diagonal entries are admissible. Measured (see Case V-5):

    ===========  =========  =========  =========
    T / K        (1, 1)     (2, 2)     (3, 3)
    ===========  =========  =========  =========
    363.0        48         7          21
    364.0        44         10         21
    365.0        41         14         21
    ===========  =========  =========  =========
    """
    _triangle, _norm, rows = _exam(temperature_K, system)
    counts = _confusion(rows)
    disagreements = [
        (tuple(np.round(row.z, 8)), row.expected.kind, row.expected.phases, row.obtained)
        for row in rows
        if row.expected.phases != row.obtained
    ]
    assert not disagreements, (counts, disagreements)
    assert set(counts) <= {(1, 1), (2, 2), (3, 3)}, counts


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_every_three_phase_answer_is_the_same_tie_triangle(temperature_K: float, system) -> None:
    """The three compositions do not depend on the feed; only the amounts do."""
    triangle, _norm, rows = _exam(temperature_K, system)
    reference = list(triangle)
    three_phase = [row for row in rows if row.obtained == 3]
    assert len(three_phase) == 21, len(three_phase)

    for row in three_phase:
        names = row.result.phase_names()
        assert sorted(names) == ["liquid1", "liquid2", "vapor"], names
        compositions = [
            np.array(row.result.phases[name].composition.fractions, dtype=float) for name in names
        ]
        best_error = float("inf")
        best_order: tuple[int, ...] = (0, 1, 2)
        for order in itertools.permutations(range(3)):
            error = max(
                float(np.max(np.abs(compositions[index] - reference[order[index]])))
                for index in range(3)
            )
            if error < best_error:
                best_error, best_order = error, order
        assert best_error < 1e-8, (row.z, best_error)

        fractions = [row.result.phase_fractions[name] for name in names]
        for index in range(3):
            assert abs(fractions[index] - row.weights[best_order[index]]) < 1e-8, (
                row.z,
                fractions,
                row.weights,
            )
        assert abs(sum(fractions) - 1.0) < 1e-12


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_every_two_phase_answer_satisfies_equal_fugacity_and_mass_balance(
    temperature_K: float, system
) -> None:
    """``ln x_i + t_i`` is equal in both phases, and the amounts add back to ``z``."""
    _names, _model, ln_gamma, psat = system
    _triangle, _norm, rows = _exam(temperature_K, system)
    two_phase = [row for row in rows if row.obtained == 2]
    assert two_phase, "the grid must contain two-phase feeds"

    for row in two_phase:
        names = row.result.phase_names()
        compositions = [
            np.array(row.result.phases[name].composition.fractions, dtype=float) for name in names
        ]
        terms = [
            np.zeros(3)
            if name.startswith("vapor")
            else _liquid_terms(composition, temperature_K, ln_gamma, psat)
            for name, composition in zip(names, compositions)
        ]
        equilibrium = float(
            np.max(np.abs(np.log(compositions[0]) + terms[0] - np.log(compositions[1]) - terms[1]))
        )
        assert equilibrium < 1e-8, (row.z, names, equilibrium)

        fractions = [row.result.phase_fractions[name] for name in names]
        balance = float(
            np.max(np.abs(fractions[0] * compositions[0] + fractions[1] * compositions[1] - row.z))
        )
        assert balance < 1e-10, (row.z, balance)


@pytest.mark.parametrize("temperature_K", TEMPERATURES)
def test_every_single_phase_answer_is_the_feed_on_its_own_candidate(
    temperature_K: float, system
) -> None:
    """One phase means the feed itself, on the candidate of lower Gibbs energy."""
    _names, _model, ln_gamma, psat = system
    _triangle, _norm, rows = _exam(temperature_K, system)
    single = [row for row in rows if row.obtained == 1]
    assert single, "the grid must contain single-phase feeds"

    for row in single:
        name = row.result.phase_names()[0]
        composition = np.array(row.result.phases[name].composition.fractions, dtype=float)
        assert float(np.max(np.abs(composition - row.z))) < 1e-12, (row.z, composition)
        g_liquid = _reduced_g(row.z, False, temperature_K, ln_gamma, psat)
        g_vapor = _reduced_g(row.z, True, temperature_K, ln_gamma, psat)
        expected_name = "vapor" if g_vapor < g_liquid else "liquid"
        assert name == expected_name, (row.z, name, g_liquid, g_vapor)
        assert row.result.diagnostics["phase_regime"] == "single-phase"
