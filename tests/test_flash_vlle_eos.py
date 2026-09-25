"""Phase addition and removal on the phi-phi path (ADR-0020, Cases P-9, P-10).

ADR-0011 gave `flash_tp` a phase addition/removal search and wired it to the
`modified-raoult` path only, because no state in this repository exercised a
third phase on an equation of state. ADR-0018 (association) and ADR-0019
(per-phase density roots) supplied the states; this file is the evidence that
the search now serves an EOS as well.

Three kinds of state are checked, and they are deliberately different shapes:

1. **A binary refusal window.** PC-SAFT (2B water, `k_ij = 0`), water /
   n-hexane at 1 atm. Just *below* the three-phase temperature `T3` the
   deepest tangent-plane minimum from the feed is a vapour, so the search
   starts `V -> LV`, the pair is unstable towards a second liquid,
   `LV -> LLV`, and the three-phase solve drives the vapour fraction negative,
   `LLV -> LL`. The correct two-liquid answer is reached by *removing* a phase
   that addition had to add first - the ADR-0011 Case R-3 mechanism, now on an
   equation of state. Before this slice these temperatures raised
   `ConvergenceError`.

   A binary at a fixed pressure has **no three-phase region**: Gibbs' phase
   rule gives `F = 2 - 3 + 2 = 1`, so three phases coexist at exactly one
   temperature, and there the three phase *amounts* are not determined by the
   mass balance either (three unknowns, two independent equations). What is
   checked at `T3` is therefore the geometry - the two-phase answers on either
   side meet there - and not a three-phase `FlashResult`, which would be
   arithmetic on an underdetermined system.

2. **A ternary tie-triangle with a vapour**, which does have a finite
   three-phase region: PC-SAFT water / ethanol / n-hexane at 1 atm and 333 K
   returns `liquid1` / `liquid2` / `vapor` with `phase_regime = "VLLE"`. That
   is the headline capability and it is `slow`-marked here only because one
   such flash costs ~35 s; the cheap three-phase check is (3).

3. **Three liquids from a cubic.** Peng-Robinson with `k_ij = 0` on the same
   ternary at 280 K returns `liquid1` / `liquid2` / `liquid3` in ~0.1 s. ADR-0019
   and the first draft of ADR-0020 both recorded "no pure Peng-Robinson
   three-phase case found in the databank"; this is one. Whether the *model* is
   right there is not asserted - what is asserted is that the returned set is a
   verified equilibrium **of that model**, checked against an independent
   Newton written here.

Every reference number in this file is computed here, by Newton solves on the
public `EquationOfState.fugacity_coefficients` interface, and never taken from
`flash_tp`'s own output. The FeOs cross-check lives in
`tests/validation/test_pcsaft_vlle_water_hexane.py`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

BINARY = ("Water", "n-Hexane")
TERNARY = ("Water", "Ethanol", "n-Hexane")
PRESSURE_PA = 101325.0

#: Offset from `T3` used for the two-sided window checks. 0.05 K is ~200,000
#: times the width of the band in which the post-split tangent-plane distance
#: of the third phase is inside `tpd_tol`, so both sides are unambiguous.
OFFSET_K = 0.05

#: Peng-Robinson three-liquid state (3): 280 K, 1 atm, water / ethanol /
#: n-hexane. Cheap enough to run several feeds in the default suite.
PR_TERNARY_T_K = 280.0
PR_TERNARY_FEEDS = ((0.2, 0.4, 0.4), (0.3, 0.3, 0.4), (0.2, 0.5, 0.3), (0.1, 0.5, 0.4))

#: PC-SAFT ternary VLLE state (2).
VLLE_TERNARY_T_K = 333.0
VLLE_TERNARY_FEEDS = ((0.4, 0.3, 0.3), (0.3, 0.3, 0.4), (0.5, 0.2, 0.3), (0.8, 0.1, 0.1))

#: Newton tolerances used by the reference solves written in this file.
_NEWTON_TOL = 1e-11


# ---------------------------------------------------------------------------
# Independent references: Newton solves on the public EOS interface
# ---------------------------------------------------------------------------


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _ln_f(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    composition: np.ndarray,
    branch: str,
) -> np.ndarray:
    """``ln(x_i phi_i)`` on one named density/compressibility branch."""
    values = np.asarray(composition, dtype=float)
    values = values / float(np.sum(values))
    phi = np.asarray(
        eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=PRESSURE_PA,
            composition=values.tolist(),
            phase=branch,
        ),
        dtype=float,
    )
    return np.log(values) + np.log(phi)


def _reduced_g(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    composition: np.ndarray,
    branch: str,
) -> float:
    """``sum_i x_i ln(x_i phi_i)``: the composition-dependent part of ``G/RT``."""
    values = np.asarray(composition, dtype=float)
    return float(np.sum(values * _ln_f(eos, mixture, temperature_K, values, branch)))


def _damped_newton(
    residual: Callable[[np.ndarray], np.ndarray],
    start: np.ndarray,
    *,
    steps: tuple[float, ...],
    admissible: Callable[[np.ndarray], bool] = lambda _u: True,
    tol: float = _NEWTON_TOL,
    max_iter: int = 60,
) -> tuple[np.ndarray, float, int]:
    """Damped Newton with a finite-difference Jacobian. Returns ``(u, residual, iterations)``."""
    u = np.array(start, dtype=float)
    worst = float(np.max(np.abs(residual(u))))
    iteration = 0
    for iteration in range(1, max_iter + 1):
        f = residual(u)
        worst = float(np.max(np.abs(f)))
        if worst < tol:
            break
        jacobian = np.zeros((f.size, u.size))
        for column in range(u.size):
            shifted = u.copy()
            shifted[column] += steps[column]
            jacobian[:, column] = (residual(shifted) - f) / steps[column]
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while scale > 1e-10:
            candidate = u + scale * direction
            if admissible(candidate) and float(np.max(np.abs(residual(candidate)))) < worst:
                break
            scale *= 0.5
        u = u + scale * direction
    return u, float(np.max(np.abs(residual(u)))), iteration


def _three_phase_point(
    eos: ct.EquationOfState, mixture: ct.Mixture
) -> tuple[float, float, float, float]:
    """The binary three-phase point from a 4-equation Newton in ``(x^I, x^II, y, T)``.

    Two liquids on the model's liquid branch and a vapour on its vapour branch,
    all at 1 atm:

        ln(x_i^I phi_i^L(x^I)) = ln(x_i^II phi_i^L(x^II)) = ln(y_i phi_i^V(y))

    Four equations (two components, two independent equalities) in four
    unknowns, the fourth being the temperature - which is what Gibbs' phase
    rule says is determined once the pressure is fixed.

    Returns ``(T3, x_water^I, x_water^II, y_water)``.
    """

    def residual(u: np.ndarray) -> np.ndarray:
        first = _ln_f(eos, mixture, u[3], np.array([u[0], 1.0 - u[0]]), "liquid")
        second = _ln_f(eos, mixture, u[3], np.array([u[1], 1.0 - u[1]]), "liquid")
        vapor = _ln_f(eos, mixture, u[3], np.array([u[2], 1.0 - u[2]]), "vapor")
        return np.concatenate([first - second, first - vapor])

    def admissible(u: np.ndarray) -> bool:
        return bool(np.all(u[:3] > 0.0) and np.all(u[:3] < 1.0) and 200.0 < u[3] < 600.0)

    # A coarse bracket, not the answer: 334.5 K with rough binodal and vapour
    # estimates. The solve below moves it by 0.3 K and the compositions by
    # 2e-2 / 1e-2.
    start = np.array([0.9999, 0.02, 0.20, 334.5])
    u, worst, _iterations = _damped_newton(
        residual, start, steps=(1e-8, 1e-8, 1e-8, 1e-7 * 334.5), admissible=admissible
    )
    assert worst < 1e-10, worst
    return float(u[3]), float(u[0]), float(u[1]), float(u[2])


def _binary_pair(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    branches: tuple[str, str],
    start: tuple[float, float],
) -> tuple[float, float]:
    """Equal fugacity for one binary pair on two named branches; returns ``(x1_a, x1_b)``."""

    def residual(u: np.ndarray) -> np.ndarray:
        a = _ln_f(eos, mixture, temperature_K, np.array([u[0], 1.0 - u[0]]), branches[0])
        b = _ln_f(eos, mixture, temperature_K, np.array([u[1], 1.0 - u[1]]), branches[1])
        return a - b

    u, worst, _iterations = _damped_newton(
        residual,
        np.array(start, dtype=float),
        steps=(1e-8, 1e-8),
        admissible=lambda v: bool(np.all(v > 0.0) and np.all(v < 1.0)),
    )
    assert worst < 1e-10, worst
    return float(u[0]), float(u[1])


def _multi_liquid_newton(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    compositions: Sequence[np.ndarray],
    branches: Sequence[str],
) -> tuple[list[np.ndarray], float]:
    """Equal fugacity across N ternary phases, each on its own named branch.

    Parameterized by ``ln(x_k / x_last)`` per phase so that every iterate is a
    positive composition - one converged phase here has a mole fraction of
    5e-14 and a Newton on the mole fractions themselves steps straight out of
    the simplex.
    """
    count = len(compositions)
    width = len(mixture.components) - 1

    def unpack(u: np.ndarray) -> list[np.ndarray]:
        phases = []
        for index in range(count):
            logs = np.concatenate([u[index * width : (index + 1) * width], [0.0]])
            weights = np.exp(logs - float(np.max(logs)))
            phases.append(weights / float(np.sum(weights)))
        return phases

    def residual(u: np.ndarray) -> np.ndarray:
        phases = unpack(u)
        terms = [
            _ln_f(eos, mixture, temperature_K, x, branch) for x, branch in zip(phases, branches)
        ]
        return np.concatenate([terms[0] - terms[index] for index in range(1, count)])

    start = np.concatenate(
        [
            np.log(np.maximum(np.asarray(x, dtype=float)[:width], 1e-300) / x[width])
            for x in compositions
        ]
    )
    u, worst, _iterations = _damped_newton(
        residual, start, steps=tuple(1e-7 * max(1.0, abs(value)) for value in start), tol=1e-13
    )
    return unpack(u), worst


def _fractions_from_mass_balance(
    compositions: Sequence[np.ndarray], z: Sequence[float]
) -> np.ndarray:
    """Phase amounts from ``z = sum_j beta_j x^j`` and ``sum_j beta_j = 1``."""
    matrix = np.vstack([np.column_stack(compositions), np.ones(len(compositions))])
    target = np.concatenate([np.asarray(z, dtype=float), [1.0]])
    beta, *_rest = np.linalg.lstsq(matrix, target, rcond=None)
    return beta


# ---------------------------------------------------------------------------
# Session references (one PC-SAFT Newton, reused)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def binary_reference() -> dict[str, object]:
    eos = PCSAFTEOS()
    mixture = _mixture(BINARY, (0.5, 0.5))
    t3, x_first, x_second, y = _three_phase_point(eos, mixture)
    return {"eos": eos, "mixture": mixture, "T3": t3, "xI": x_first, "xII": x_second, "y": y}


@pytest.fixture(scope="module")
def below_t3(binary_reference) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(BINARY, (0.5, 0.5)),
        temperature_K=float(binary_reference["T3"]) - OFFSET_K,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )


@pytest.fixture(scope="module")
def above_t3(binary_reference) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(BINARY, (0.5, 0.5)),
        temperature_K=float(binary_reference["T3"]) + OFFSET_K,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )


# ---------------------------------------------------------------------------
# (1) The binary three-phase point and the refusal window around it
# ---------------------------------------------------------------------------


def test_the_three_phase_point_is_a_real_stationary_state(binary_reference) -> None:
    """The 4-equation Newton, and why per-phase density roots are needed at all.

    Every one of the three phase compositions has **two** mechanically stable
    density roots at this state, so "liquid" and "vapor" are two different
    numbers for each of them and a split that pinned one phase to each branch
    could not describe two liquids and a vapour at once (ADR-0019).
    """
    eos, mixture = binary_reference["eos"], binary_reference["mixture"]
    t3 = float(binary_reference["T3"])

    # Weak external sanity check: the water / n-hexane heteroazeotrope at 1 atm
    # is commonly tabulated near 61.6 C with y_water ~ 0.21. Not a reference
    # value (no primary source verified), so the band is wide.
    assert 330.0 < t3 < 340.0, t3
    assert 0.15 < float(binary_reference["y"]) < 0.28

    for composition in (
        [binary_reference["xI"], 1.0 - float(binary_reference["xI"])],
        [binary_reference["xII"], 1.0 - float(binary_reference["xII"])],
        [binary_reference["y"], 1.0 - float(binary_reference["y"])],
    ):
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=t3,
            pressure_Pa=PRESSURE_PA,
            composition=[float(value) for value in composition],
        )
        assert len(roots) == 2, roots
        assert roots[0] < 100.0 < 1000.0 < roots[-1]


def test_below_the_three_phase_temperature_the_search_returns_the_conjugate_liquids(
    binary_reference, below_t3
) -> None:
    """The refusal window, resolved by adding a phase and then removing one.

    Case P-9 (i). Before ADR-0020 this state raised `ConvergenceError` (the
    two-phase set is provably not the answer); now the two liquids come back,
    and they are the tie line an independent two-liquid Newton finds.
    """
    eos, mixture = binary_reference["eos"], binary_reference["mixture"]
    temperature = float(binary_reference["T3"]) - OFFSET_K

    assert sorted(below_t3.phases) == ["liquid1", "liquid2"]
    assert below_t3.vapor_fraction is None
    assert below_t3.diagnostics["phase_regime"] == "LLE"
    assert below_t3.diagnostics["phase_set_history"] == "V -> LV -> LLV -> LL"
    assert below_t3.diagnostics["phases_added"] == 1
    assert below_t3.diagnostics["phases_removed"] == 1
    assert below_t3.diagnostics["phase_label_method"] == "compressibility"
    assert below_t3.diagnostics["post_split_status"] == "stable"
    assert float(below_t3.diagnostics["equilibrium_residual"]) < 1e-9
    assert float(below_t3.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(below_t3.diagnostics["delta_g_split_rt"]) < 0.0

    reference = _binary_pair(eos, mixture, temperature, ("liquid", "liquid"), (0.9999, 0.022))
    ours = (
        float(below_t3.phases["liquid1"].composition.fractions[0]),
        float(below_t3.phases["liquid2"].composition.fractions[0]),
    )
    assert ours[0] == pytest.approx(reference[0], abs=1e-9)
    assert ours[1] == pytest.approx(reference[1], abs=1e-9)

    # Amounts from the lever rule on the independent tie line.
    beta = (0.5 - reference[0]) / (reference[1] - reference[0])
    assert below_t3.phase_fractions["liquid2"] == pytest.approx(beta, abs=1e-9)

    # Both phases stable, including against the vapour candidate: `stability_tp`
    # selects the lowest-Gibbs of the two density roots at every iterate, so a
    # vapour stationary point is inside the trial set it searches.
    for name, phase in below_t3.phases.items():
        verdict = ct.stability_tp(
            _mixture(BINARY, phase.composition.fractions),
            temperature_K=temperature,
            pressure_Pa=PRESSURE_PA,
            eos=PCSAFTEOS(),
        )
        assert verdict.status == "stable", (name, verdict.tpd_min)


def test_below_t3_the_two_liquids_beat_the_vapor_liquid_pair_in_gibbs_energy(
    binary_reference, below_t3
) -> None:
    """`delta_g_vs_two_phase_rt` is the two-phase answer the search started from.

    Both pairs exist below `T3` - the vapour-liquid one is the pair the search
    converged first and then rejected - so which is the equilibrium is a Gibbs
    energy comparison, and it is done here from independent Newton solves.
    """
    eos, mixture = binary_reference["eos"], binary_reference["mixture"]
    temperature = float(binary_reference["T3"]) - OFFSET_K

    liquids = _binary_pair(eos, mixture, temperature, ("liquid", "liquid"), (0.9999, 0.022))
    vapor_liquid = _binary_pair(eos, mixture, temperature, ("liquid", "vapor"), (0.9999, 0.213))

    def pair_energy(pair: tuple[float, float], branches: tuple[str, str]) -> float:
        beta = (0.5 - pair[0]) / (pair[1] - pair[0])
        first = np.array([pair[0], 1.0 - pair[0]])
        second = np.array([pair[1], 1.0 - pair[1]])
        return (1.0 - beta) * _reduced_g(
            eos, mixture, temperature, first, branches[0]
        ) + beta * _reduced_g(eos, mixture, temperature, second, branches[1])

    g_liquids = pair_energy(liquids, ("liquid", "liquid"))
    g_vapor_liquid = pair_energy(vapor_liquid, ("liquid", "vapor"))
    g_feed = min(
        _reduced_g(eos, mixture, temperature, np.array([0.5, 0.5]), branch)
        for branch in ("liquid", "vapor")
    )

    assert g_liquids < g_vapor_liquid < g_feed
    assert float(below_t3.diagnostics["delta_g_vs_two_phase_rt"]) == pytest.approx(
        g_liquids - g_vapor_liquid, abs=1e-9
    )
    assert float(below_t3.diagnostics["delta_g_split_rt"]) == pytest.approx(
        g_liquids - g_feed, abs=1e-9
    )


def test_just_above_the_three_phase_temperature_the_answer_is_vapor_liquid(
    binary_reference, above_t3
) -> None:
    """Case P-9 (ii): the search is not entered, and Gibbs says it should not be."""
    eos, mixture = binary_reference["eos"], binary_reference["mixture"]
    temperature = float(binary_reference["T3"]) + OFFSET_K

    assert sorted(above_t3.phases) == ["liquid", "vapor"]
    assert above_t3.diagnostics["phase_regime"] == "VLE"
    assert "phase_set_history" not in above_t3.diagnostics
    assert above_t3.diagnostics["post_split_status"] == "stable"

    vapor_liquid = _binary_pair(eos, mixture, temperature, ("liquid", "vapor"), (0.9999, 0.213))
    assert float(above_t3.phases["liquid"].composition.fractions[0]) == pytest.approx(
        vapor_liquid[0], abs=1e-9
    )
    assert float(above_t3.phases["vapor"].composition.fractions[0]) == pytest.approx(
        vapor_liquid[1], abs=1e-9
    )

    liquids = _binary_pair(eos, mixture, temperature, ("liquid", "liquid"), (0.9999, 0.0226))

    def pair_energy(pair: tuple[float, float], branches: tuple[str, str]) -> float:
        beta = (0.5 - pair[0]) / (pair[1] - pair[0])
        return (1.0 - beta) * _reduced_g(
            eos, mixture, temperature, np.array([pair[0], 1.0 - pair[0]]), branches[0]
        ) + beta * _reduced_g(
            eos, mixture, temperature, np.array([pair[1], 1.0 - pair[1]]), branches[1]
        )

    assert pair_energy(vapor_liquid, ("liquid", "vapor")) < pair_energy(
        liquids, ("liquid", "liquid")
    )


def test_at_the_three_phase_temperature_the_two_answers_meet(binary_reference) -> None:
    """Case P-9 (iii): the geometry at `T3`, and what is *not* claimed there.

    A binary at a fixed pressure has one degree of freedom, so three phases
    coexist at exactly one temperature - and there the three phase *amounts*
    solve an underdetermined system (three unknowns, two independent balances).
    `flash_tp` therefore returns a two-phase answer at `T3`: the vapour-liquid
    edge of the tie triangle, whose two compositions are two of the three
    vertices the 4-equation Newton found. The third vertex is the incipient
    phase whose tangent-plane distance is zero there, which is why the
    post-split test reports the pair stable.
    """
    t3 = float(binary_reference["T3"])
    result = ct.flash_tp(
        _mixture(BINARY, (0.5, 0.5)),
        temperature_K=t3,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert float(result.phases["liquid"].composition.fractions[0]) == pytest.approx(
        float(binary_reference["xI"]), abs=1e-6
    )
    assert float(result.phases["vapor"].composition.fractions[0]) == pytest.approx(
        float(binary_reference["y"]), abs=1e-6
    )

    # The missing vertex is a zero of the tangent-plane distance from either
    # returned phase, not a negative one - which is the statement "three phases
    # coexist here" written in the stability test's own terms.
    second_liquid = np.array([float(binary_reference["xII"]), 1.0 - float(binary_reference["xII"])])
    eos, mixture = binary_reference["eos"], binary_reference["mixture"]
    liquid = np.asarray(result.phases["liquid"].composition.fractions, dtype=float)
    plane = _ln_f(eos, mixture, t3, liquid, "liquid")
    incipient = _ln_f(eos, mixture, t3, second_liquid, "liquid")
    tpd = float(np.sum(second_liquid * (incipient - plane)))
    assert abs(tpd) < 1e-9, tpd


def test_max_phases_two_still_refuses_the_window(binary_reference) -> None:
    """Case P-9: `max_phases = 2` reproduces the pre-ADR-0020 refusal."""
    with pytest.raises(ct.ConvergenceError, match="not a stable phase set"):
        ct.flash_tp(
            _mixture(BINARY, (0.5, 0.5)),
            temperature_K=float(binary_reference["T3"]) - OFFSET_K,
            pressure_Pa=PRESSURE_PA,
            eos=PCSAFTEOS(),
            settings=ct.FlashSettings(max_phases=2),
        )


# ---------------------------------------------------------------------------
# (3) Three liquids from a cubic - the cheap three-phase evidence
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pr_three_liquid() -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(TERNARY, PR_TERNARY_FEEDS[0]),
        temperature_K=PR_TERNARY_T_K,
        pressure_Pa=PRESSURE_PA,
        eos=ct.PengRobinsonEOS(),
    )


def test_peng_robinson_returns_three_liquid_phases(pr_three_liquid) -> None:
    """A pure cubic three-phase state, which ADR-0019 recorded as not found.

    Water / ethanol / n-hexane at 280 K and 1 atm with `k_ij = 0`. Nothing here
    asserts that the model is *right* - `k_ij = 0` between water and a
    hydrocarbon is not a serious parameterization - only that what comes back
    is a verified equilibrium of that model, with three phases discovered
    rather than assumed.
    """
    result = pr_three_liquid
    assert sorted(result.phases) == ["liquid1", "liquid2", "liquid3"]
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_count"] == 3
    assert result.diagnostics["phase_state"] == "three_phase"
    assert result.diagnostics["phase_regime"] == "LLE"
    assert result.diagnostics["phase_set_history"] == "L -> LL -> LLL"
    assert result.diagnostics["phases_added"] == 1
    assert result.diagnostics["phases_removed"] == 0
    assert result.diagnostics["phase_label_method"] == "compressibility"
    assert result.diagnostics["post_split_status"] == "stable"
    assert float(result.diagnostics["equilibrium_residual"]) < 1e-8
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert float(result.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0

    # liquid1 / liquid2 / liquid3 are ordered by the first component's mole
    # fraction (ADR-0019 decision 3 generalized), so the names are comparable
    # between feeds rather than being solve-order roles.
    water = [float(result.phases[name].composition.fractions[0]) for name in sorted(result.phases)]
    assert water == sorted(water, reverse=True)


def test_the_peng_robinson_triangle_matches_an_independent_newton(pr_three_liquid) -> None:
    """Six equal-fugacity equations, solved here, against the search's answer."""
    eos = ct.PengRobinsonEOS()
    mixture = _mixture(TERNARY, (1 / 3, 1 / 3, 1 / 3))
    names = ("liquid1", "liquid2", "liquid3")
    ours = [
        np.asarray(pr_three_liquid.phases[name].composition.fractions, dtype=float)
        for name in names
    ]
    solved, worst = _multi_liquid_newton(
        eos, mixture, PR_TERNARY_T_K, ours, ("liquid", "liquid", "liquid")
    )
    assert worst < 1e-12, worst
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(solved, ours)) < 1e-7

    # Every feed inside the triangle returns the same three vertices, with the
    # amounts the mass balance alone fixes.
    for feed in PR_TERNARY_FEEDS:
        result = ct.flash_tp(
            _mixture(TERNARY, feed),
            temperature_K=PR_TERNARY_T_K,
            pressure_Pa=PRESSURE_PA,
            eos=eos,
        )
        assert sorted(result.phases) == list(names)
        phases = [
            np.asarray(result.phases[name].composition.fractions, dtype=float) for name in names
        ]
        assert max(float(np.max(np.abs(a - b))) for a, b in zip(phases, solved)) < 1e-7
        beta = _fractions_from_mass_balance(solved, feed)
        got = np.array([result.phase_fractions[name] for name in names])
        assert float(np.max(np.abs(beta - got))) < 1e-7

        # G3 < G2 < G1 for every one of them.
        assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
        assert float(result.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def test_a_superheated_feed_is_a_single_vapor(binary_reference) -> None:
    result = ct.flash_tp(
        _mixture(BINARY, (0.5, 0.5)),
        temperature_K=float(binary_reference["T3"]) + 20.0,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert sorted(result.phases) == ["vapor"]
    assert result.diagnostics["phase_regime"] == "single-phase"
    assert "phase_set_history" not in result.diagnostics


def test_a_water_rich_subcooled_feed_is_a_single_liquid() -> None:
    """`z_water = 0.99999` at 300 K is inside the water-rich binodal.

    `z_water = 0.999` - the obvious choice - is **not**: this model puts only
    1.67e-05 mole fraction of n-hexane in the water-rich liquid at 298.15 K
    (validation Case P-8), so a feed with 1e-03 of it is genuinely two liquids.
    That half is `test_a_feed_inside_the_binodal_is_two_liquids` below, kept
    `slow` only for its cost.
    """
    single = ct.flash_tp(
        _mixture(BINARY, (0.99999, 0.00001)),
        temperature_K=300.0,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert sorted(single.phases) == ["liquid"]
    assert single.diagnostics["phase_regime"] == "single-phase"
    assert "phase_set_history" not in single.diagnostics


@pytest.mark.slow
def test_a_feed_inside_the_binodal_is_two_liquids() -> None:
    """The other half of the control above: `z_water = 0.999` is *inside* the gap."""
    inside = ct.flash_tp(
        _mixture(BINARY, (0.999, 0.001)),
        temperature_K=300.0,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert sorted(inside.phases) == ["liquid1", "liquid2"]
    assert "phase_set_history" not in inside.diagnostics


def test_no_bit_identity_fixture_state_carries_a_phase_search_key() -> None:
    """ADR-0020 decision 6: the search is never entered on the pinned states.

    `tests/test_flash_refactor_bit_identity.py` compares every float of all 155
    states with `==`, so a state that started entering the search would fail
    there. This asserts the complementary half directly on the fixture: none of
    those states carries a diagnostics key the search produces, so the fixture
    itself is evidence that the new code path is not reached by any of them.
    """
    path = Path(__file__).resolve().parent / "fixtures" / "flash" / "refactor_bit_identity_v3.json"
    with path.open(encoding="utf-8") as handle:
        fixture = json.load(handle)
    assert len(fixture) == 155
    search_keys = {
        "phase_set_history",
        "phases_added",
        "phases_removed",
        "rachford_rice_iterations",
        "delta_g_vs_two_phase_rt",
    }
    for label, state in fixture.items():
        assert not (search_keys & set(state["diagnostics"])), label


# ---------------------------------------------------------------------------
# Slow: the full window scan, the ternary VLLE triangle, permutation invariance
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("z_water", [0.3, 0.7])
def test_the_window_scan_never_raises(binary_reference, z_water: float) -> None:
    """Case P-9 (iv): 41 temperatures per feed across `T3`, no `ConvergenceError`.

    Both feeds switch `LLE -> VLE` exactly once, at `T3`. At `z = 0.7` that is
    what ADR-0021 repaired: until the equation-of-state trials were pinned to a
    density root, the vapour stationary point reachable from the hexane-rich
    liquid was missed and every temperature above `T3` came back as two liquids
    that are metastable by up to 6.5e-03 RT (validation Cases P-9 (iv) and
    P-11 carry the measured Gibbs-energy gaps). The two feeds are checked
    together, and the switch is located against the independently computed `T3`.
    """
    t3 = float(binary_reference["T3"])
    verdicts = []
    for temperature in np.linspace(t3 - 1.0, t3 + 1.0, 41):
        result = ct.flash_tp(
            _mixture(BINARY, (z_water, 1.0 - z_water)),
            temperature_K=float(temperature),
            pressure_Pa=PRESSURE_PA,
            eos=PCSAFTEOS(),
        )
        verdicts.append(str(result.diagnostics["phase_regime"]))

    assert len(verdicts) == 41
    assert set(verdicts) == {"LLE", "VLE"}
    assert verdicts[0] == "LLE" and verdicts[-1] == "VLE"
    switches = [i for i in range(1, 41) if verdicts[i] != verdicts[i - 1]]
    assert len(switches) == 1
    # The one switch straddles `T3`: the last `LLE` grid point is at or below
    # it and the first `VLE` one at or above it. Asserted as a bracket rather
    # than as `|boundary - T3| <= 0.05`, because the grid step *is* 0.05 K and
    # that comparison fails on the last bit whenever the switch lands on the
    # point above `T3` (it does at `z = 0.7`). The sub-microkelvin statement is
    # the bisection test below, and Case P-11.
    grid = np.linspace(t3 - 1.0, t3 + 1.0, 41)
    last_lle = float(grid[switches[0] - 1])
    first_vle = float(grid[switches[0]])
    assert last_lle <= t3 + 1e-9
    assert first_vle >= t3 - 1e-9
    assert first_vle - last_lle == pytest.approx(2.0 / 40.0)


@pytest.mark.slow
def test_the_verdict_boundary_locates_the_three_phase_temperature(binary_reference) -> None:
    """Bisect the LLE/VLE verdict; it must land on the independent `T3`."""
    t3 = float(binary_reference["T3"])
    low, high = t3 - 1.0e-3, t3 + 1.0e-3
    for _step in range(10):
        middle = 0.5 * (low + high)
        result = ct.flash_tp(
            _mixture(BINARY, (0.5, 0.5)),
            temperature_K=middle,
            pressure_Pa=PRESSURE_PA,
            eos=PCSAFTEOS(),
        )
        if result.diagnostics["phase_regime"] == "LLE":
            low = middle
        else:
            high = middle
    assert abs(0.5 * (low + high) - t3) < 1e-5


@pytest.mark.slow
def test_the_pcsaft_ternary_returns_a_vapor_liquid_liquid_tie_triangle() -> None:
    """Case P-10: a finite three-phase region with a vapour in it.

    Water / ethanol / n-hexane, PC-SAFT with 2B water and 2B ethanol,
    `k_ij = 0`, 1 atm, 333 K (below the binary `T3`). Every feed inside the
    triangle returns the same three vertices; the reference is a 6-equation
    Newton written here on the same public interface.
    """
    eos = PCSAFTEOS()
    mixture = _mixture(TERNARY, (1 / 3, 1 / 3, 1 / 3))
    names = ("liquid1", "liquid2", "vapor")
    first = ct.flash_tp(
        _mixture(TERNARY, VLLE_TERNARY_FEEDS[0]),
        temperature_K=VLLE_TERNARY_T_K,
        pressure_Pa=PRESSURE_PA,
        eos=eos,
    )
    assert sorted(first.phases) == ["liquid1", "liquid2", "vapor"]
    assert first.diagnostics["phase_regime"] == "VLLE"
    assert first.diagnostics["phase_count"] == 3
    assert first.diagnostics["phase_set_history"] == "L -> LL -> LLV"
    assert first.diagnostics["phase_label_method"] == "compressibility"
    assert first.diagnostics["post_split_status"] == "stable"
    assert first.vapor_fraction == pytest.approx(first.phase_fractions["vapor"])
    assert float(first.diagnostics["equilibrium_residual"]) < 1e-9
    assert float(first.diagnostics["delta_g_split_rt"]) < 0.0
    assert float(first.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0

    ours = [np.asarray(first.phases[name].composition.fractions, dtype=float) for name in names]
    solved, worst = _multi_liquid_newton(
        eos, mixture, VLLE_TERNARY_T_K, ours, ("liquid", "liquid", "vapor")
    )
    assert worst < 1e-11, worst
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(solved, ours)) < 1e-8

    for feed in VLLE_TERNARY_FEEDS[1:]:
        result = ct.flash_tp(
            _mixture(TERNARY, feed),
            temperature_K=VLLE_TERNARY_T_K,
            pressure_Pa=PRESSURE_PA,
            eos=eos,
        )
        assert sorted(result.phases) == ["liquid1", "liquid2", "vapor"]
        phases = [
            np.asarray(result.phases[name].composition.fractions, dtype=float) for name in names
        ]
        assert max(float(np.max(np.abs(a - b))) for a, b in zip(phases, solved)) < 1e-8
        beta = _fractions_from_mass_balance(solved, feed)
        got = np.array([result.phase_fractions[name] for name in names])
        assert float(np.max(np.abs(beta - got))) < 1e-8


@pytest.mark.slow
def test_the_window_answer_is_permutation_invariant(binary_reference) -> None:
    """Reordering the components reorders the names, not the phase set."""
    temperature = float(binary_reference["T3"]) - OFFSET_K
    forward = ct.flash_tp(
        _mixture(BINARY, (0.5, 0.5)),
        temperature_K=temperature,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    reversed_result = ct.flash_tp(
        _mixture(tuple(reversed(BINARY)), (0.5, 0.5)),
        temperature_K=temperature,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert sorted(reversed_result.phases) == ["liquid1", "liquid2"]
    forward_set = sorted(
        tuple(round(value, 12) for value in phase.composition.fractions)
        for phase in forward.phases.values()
    )
    reversed_set = sorted(
        tuple(round(value, 12) for value in reversed(phase.composition.fractions))
        for phase in reversed_result.phases.values()
    )
    assert forward_set == reversed_set
