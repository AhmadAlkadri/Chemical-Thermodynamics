"""The nine `multiphase-solver-failure` states of the map, repaired (ADR-0029).

The robustness map at `74820b8` (ADR-0027 amendment, ledger Case R-MAP-2) found
fourteen refusals in 2505 states. Nine of them were one class,
`multiphase-solver-failure`, in one family, `eos-three-phase`, and they had two
distinct causes. This file is the evidence that both are fixed and that each
repaired state is a *verified* answer rather than merely a returned one.

**(i) An irreversible removal** - 5 states, refusal stage `collapsed`. PC-SAFT
(2B water, `k_ij = 0`) water / n-hexane at 1 atm and `z_water = 0.05`, at
`T3 + 0.01` through `T3 + 1.0 K`, plus the PC-SAFT ternary feed
`(0.1, 0.1, 0.8)` at 333 K. The feed is unstable, the two-liquid split is
post-split unstable on both phases, a vapour is added, and the three-phase set
has no Rachford-Rice solution at all (Gibbs' phase rule: a binary at fixed
pressure has three phases at one temperature only). The recession direction
then names the **hexane-rich** liquid as the phase leaving, the water-rich pair
that remains negative-flashes, and the search - which could not take a removal
back - raised. ADR-0029 makes removal reversible: the pair the *other* removal
leaves is the vapour-liquid answer, and it is post-split stable.

**(ii) A second-order stage that cannot form its Hessian** - 4 states, refusal
stage `split`. Peng-Robinson (`k_ij = 0`) water / ethanol / n-hexane at 1 atm,
feeds `(0.2, 0.6, 0.2)` at 280 K and 300 K and `(0.4, 0.4, 0.2)` /
`(0.5, 0.3, 0.2)` at 300 K. These are genuine three-liquid states whose
water-rich phase holds n-hexane at `x ~ 1e-12`; the linear second-order stage
perturbs mole numbers by `1e-7`, its first Hessian column leaves the box, and
it aborts at iteration 1 leaving successive substitution's residual (1.8e-08 to
1.5e-04) to be reported as a failure. ADR-0029 adds the multiphase log-space
stage, which re-chooses the reference phase and carries `u = ln n`.

Every reference number here is computed in this file by Newton solves on the
public `EquationOfState.fugacity_coefficients` interface and never read out of
`flash_tp`'s own output. The FeOs cross-check for the PC-SAFT states lives in
`tests/validation/test_pcsaft_vlle_water_hexane.py`.
"""

from __future__ import annotations

import itertools
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.flash._multiphase import _ranked_removals, _RemovalChoice, _take_next_removal
from chemthermo.flash._multiphase_log_space import reference_phase_index

BINARY = ("Water", "n-Hexane")
TERNARY = ("Water", "Ethanol", "n-Hexane")
PRESSURE_PA = 101325.0

#: The water / n-hexane three-phase temperature, the constant the map sweeps
#: around (`chemthermo.bench.robustness.EOS3P_T3_K`, itself from the
#: independent 4-equation Newton of Case P-9).
T3_K = 334.807826336

#: The four water-lean offsets above `T3` that refused at `74820b8`. Only the
#: widest runs by default; see the `slow` marks below.
WATER_LEAN_OFFSETS_K = (0.01, 0.1, 0.5, 1.0)
WATER_LEAN_FEED = (0.05, 0.95)

#: The four Peng-Robinson three-liquid states that refused at `74820b8`.
PR_THREE_LIQUID_STATES: tuple[tuple[tuple[float, float, float], float], ...] = (
    ((0.2, 0.6, 0.2), 280.0),
    ((0.2, 0.6, 0.2), 300.0),
    ((0.4, 0.4, 0.2), 300.0),
    ((0.5, 0.3, 0.2), 300.0),
)

#: Tolerances asserted below (ledger Case P-18 records the measured values).
NEWTON_TOL = 1e-10
COMPOSITION_TOL = 1e-8
MASS_BALANCE_TOL = 1e-12
EQUILIBRIUM_TOL = 1e-8


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
    steps: Sequence[float],
    admissible: Callable[[np.ndarray], bool] = lambda _u: True,
    tol: float = 1e-13,
    max_iter: int = 120,
) -> tuple[np.ndarray, float]:
    """Damped Newton with a finite-difference Jacobian. Returns ``(u, residual)``."""
    u = np.array(start, dtype=float)
    for _iteration in range(max_iter):
        f = residual(u)
        worst = float(np.max(np.abs(f)))
        if worst < tol:
            break
        jacobian = np.zeros((f.size, u.size))
        for column in range(u.size):
            shifted = u.copy()
            shifted[column] += steps[column]
            jacobian[:, column] = (residual(shifted) - f) / steps[column]
        try:
            direction = np.linalg.solve(jacobian, -f)
        except np.linalg.LinAlgError:  # pragma: no cover - a singular Jacobian
            direction = -f
        scale = 1.0
        while scale > 1e-12:
            candidate = u + scale * direction
            if admissible(candidate) and float(np.max(np.abs(residual(candidate)))) < worst:
                break
            scale *= 0.5
        u = u + scale * direction
    return u, float(np.max(np.abs(residual(u))))


def _equilibrium_newton(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    z: np.ndarray,
    compositions: Sequence[np.ndarray],
    branches: Sequence[str],
    *,
    perturb: float = 1e-3,
) -> tuple[list[np.ndarray], np.ndarray, float]:
    """The full ``NC - 1`` equilibrium system for ``N`` phases and ``C`` components.

    ``(N - 1) C`` equal-fugacity equations plus ``C - 1`` independent mass
    balances, in ``N (C - 1)`` composition degrees of freedom plus ``N - 1``
    phase fractions - 3 equations for a two-phase binary, 5 for a two-phase
    ternary, 8 for a three-phase ternary. Compositions are carried as
    ``ln(x_k / x_last)`` so that every iterate is a positive composition: one
    converged phase here holds a mole fraction of 5e-14 and a Newton on the
    mole fractions themselves steps straight out of the simplex.

    ``perturb`` moves the start off `flash_tp`'s answer so that convergence
    back onto it is a result and not an identity.
    """
    count = len(compositions)
    width = len(mixture.components) - 1

    def unpack(u: np.ndarray) -> tuple[list[np.ndarray], np.ndarray]:
        phases = []
        for index in range(count):
            logs = np.concatenate([u[index * width : (index + 1) * width], [0.0]])
            weights = np.exp(logs - float(np.max(logs)))
            phases.append(weights / float(np.sum(weights)))
        tail = u[count * width :]
        return phases, np.concatenate([tail, [1.0 - float(np.sum(tail))]])

    def residual(u: np.ndarray) -> np.ndarray:
        phases, beta = unpack(u)
        terms = [
            _ln_f(eos, mixture, temperature_K, x, branch) for x, branch in zip(phases, branches)
        ]
        equal = [terms[0] - terms[index] for index in range(1, count)]
        recombined = sum(weight * x for weight, x in zip(beta, phases))
        return np.concatenate(equal + [(recombined - z)[:width]])

    start = np.concatenate(
        [
            np.log(np.maximum(np.asarray(x, dtype=float)[:width], 1e-300) / x[width])
            for x in compositions
        ]
        + [np.full(count - 1, 1.0 / count)]
    )
    start = start + perturb * np.array([(-1.0) ** k for k in range(start.size)])

    def admissible(u: np.ndarray) -> bool:
        _phases, beta = unpack(u)
        return bool(np.all(beta > 0.0) and np.all(beta < 1.0))

    u, worst = _damped_newton(
        residual,
        start,
        steps=[1e-7 * max(1.0, abs(value)) for value in start],
        admissible=admissible,
    )
    phases, beta = unpack(u)
    return phases, beta, worst


def _binary_pair(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    branches: tuple[str, str],
    start: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Equal fugacity for one binary pair on two named branches.

    Two equations in two unknowns - the mass balance is *not* imposed, so the
    tie line is found first and the lever rule then says whether the feed lies
    on it (:func:`_lever_rule`). That separation is the point: the pair the
    pre-ADR-0029 search was left holding is a perfectly good tie line that
    simply does not contain this feed.
    """

    def residual(u: np.ndarray) -> np.ndarray:
        first = _ln_f(eos, mixture, temperature_K, np.array([u[0], 1.0 - u[0]]), branches[0])
        second = _ln_f(eos, mixture, temperature_K, np.array([u[1], 1.0 - u[1]]), branches[1])
        return first - second

    u, worst = _damped_newton(
        residual,
        np.array(start, dtype=float),
        steps=(1e-8, 1e-8),
        admissible=lambda v: bool(np.all(v > 0.0) and np.all(v < 1.0)),
        # 1e-11 rather than 1e-13: the last two digits cost 75 more iterations
        # of PC-SAFT density solves and move nothing this file asserts.
        tol=1e-11,
    )
    assert worst < NEWTON_TOL, worst
    assert abs(u[0] - u[1]) > 1e-6, u
    return np.array([u[0], 1.0 - u[0]]), np.array([u[1], 1.0 - u[1]])


def _lever_rule(pair: tuple[np.ndarray, np.ndarray], z: np.ndarray) -> float:
    """Fraction of the second phase that puts the feed on this tie line."""
    first, second = pair
    return float((z[0] - first[0]) / (second[0] - first[0]))


def _pair_reduced_g(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    pair: tuple[np.ndarray, np.ndarray],
    branches: tuple[str, str],
    beta: float,
) -> float:
    first, second = pair
    return (1.0 - beta) * _reduced_g(
        eos, mixture, temperature_K, first, branches[0]
    ) + beta * _reduced_g(eos, mixture, temperature_K, second, branches[1])


def _phases_of(result: ct.FlashResult) -> tuple[list[str], list[np.ndarray], np.ndarray, list[str]]:
    """``(names, compositions, fractions, branches)`` of a `FlashResult`.

    The branch of a phase is its *name*: on the EOS path a phase named
    ``"vapor"`` converged on the vapour root and everything else on a liquid
    root (ADR-0017 measures the name from the root, so reading it back this way
    is the inverse of that measurement, not an assumption about it).
    """
    names = list(result.phases)
    compositions = [
        np.asarray(result.phases[name].composition.fractions, dtype=float) for name in names
    ]
    fractions = np.array([float(result.phase_fractions[name]) for name in names])
    branches = ["vapor" if name == "vapor" else "liquid" for name in names]
    return names, compositions, fractions, branches


def _flash(
    names: Sequence[str],
    z: Sequence[float],
    temperature_K: float,
    eos: ct.EquationOfState,
) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(names, z), temperature_K=temperature_K, pressure_Pa=PRESSURE_PA, eos=eos
    )


# ---------------------------------------------------------------------------
# (i) The water-lean band above T3: an irreversible removal, made reversible
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "offset_K",
    # `T3 + 1.0 K` is the state the robustness map's quick subset pins, so it
    # runs by default; the three narrower offsets are repetitions of the same
    # scan (`.agents/dev-contract.md`, "The `slow` marker") at ~2 s each.
    [
        pytest.param(offset, marks=pytest.mark.slow)
        for offset in WATER_LEAN_OFFSETS_K
        if offset != 1.0
    ]
    + [1.0],
)
def test_the_water_lean_band_above_t3_returns_a_verified_vapor_liquid_state(
    offset_K: float,
) -> None:
    """z_water = 0.05 just above T3: vapour plus the *hexane-rich* liquid.

    Before ADR-0029 the search removed the hexane-rich liquid from the LLV set
    and refused on what was left. The answer is the other removal, and the
    checks here are that it is an equilibrium (an independent Newton), that it
    is *the* equilibrium (Gibbs against every competing two-phase set and
    against one phase), and that it closes the mass balance.
    """
    eos = PCSAFTEOS()
    temperature = T3_K + offset_K
    mixture = _mixture(BINARY, WATER_LEAN_FEED)
    result = _flash(BINARY, WATER_LEAN_FEED, temperature, eos)

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert result.diagnostics["post_split_stable"] is True
    assert float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL
    assert float(result.diagnostics["equilibrium_residual"]) < EQUILIBRIUM_TOL
    # The search reached it by adding a vapour to the two liquids and then
    # undoing the removal the recession direction named first.
    assert str(result.diagnostics["phase_set_history"]).startswith("L -> LL -> LLV -> LV")
    assert int(result.diagnostics["phases_removed"]) >= 2

    names, compositions, fractions, branches = _phases_of(result)
    # It is the hexane-rich liquid, not the water-rich one.
    liquid = compositions[names.index("liquid")]
    assert liquid[0] < 0.05, liquid

    z = np.asarray(_mixture(BINARY, WATER_LEAN_FEED).composition.fractions, dtype=float)
    solved, beta, worst = _equilibrium_newton(eos, mixture, temperature, z, compositions, branches)
    assert worst < NEWTON_TOL, worst
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(solved, compositions)) < COMPOSITION_TOL
    assert float(np.max(np.abs(beta - fractions))) < COMPOSITION_TOL

    returned_g = float(
        sum(
            weight * _reduced_g(eos, mixture, temperature, x, branch)
            for weight, x, branch in zip(fractions, compositions, branches)
        )
    )
    for branch in ("liquid", "vapor"):
        assert returned_g < _reduced_g(eos, mixture, temperature, z, branch)

    # The two competing two-phase sets, each solved here on its own.
    # `VL water-rich` is the pair the pre-ADR-0029 search was left holding:
    # its lever rule puts the feed *outside* the tie line, which is the
    # negative flash the old code then raised on. `LL` is a real pair and it
    # is the higher Gibbs energy.
    liquid_liquid = _binary_pair(eos, mixture, temperature, ("liquid", "liquid"), (0.9999, 0.0226))
    water_rich_vl = _binary_pair(eos, mixture, temperature, ("liquid", "vapor"), (0.9999, 0.213))
    hexane_rich_vl = _binary_pair(eos, mixture, temperature, ("liquid", "vapor"), (0.0226, 0.213))

    inadmissible = _lever_rule(water_rich_vl, z)
    assert not 0.0 < inadmissible < 1.0, inadmissible

    two_liquid_beta = _lever_rule(liquid_liquid, z)
    assert 0.0 < two_liquid_beta < 1.0, two_liquid_beta
    two_liquid_g = _pair_reduced_g(
        eos, mixture, temperature, liquid_liquid, ("liquid", "liquid"), two_liquid_beta
    )
    assert returned_g < two_liquid_g

    # ... and the pair the repaired search returns is the one it found.
    assert abs(hexane_rich_vl[0][0] - liquid[0]) < COMPOSITION_TOL


@pytest.mark.slow  # ~8 s; the binary band above covers the same repair by default
def test_the_hexane_rich_ternary_corner_returns_a_verified_vapor_liquid_state() -> None:
    """PC-SAFT (0.1, 0.1, 0.8) at 333 K: Case P-10 (i)'s pre-existing failure.

    The same shape as the binary band - the LLV set has no Rachford-Rice
    solution, the first removal is the wrong one - on a ternary feed at the
    tie-triangle's hexane-rich edge, where the equilibrium is two phases rather
    than three.
    """
    eos = PCSAFTEOS()
    feed = (0.1, 0.1, 0.8)
    mixture = _mixture(TERNARY, feed)
    result = _flash(TERNARY, feed, 333.0, eos)

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["post_split_stable"] is True
    assert float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL
    assert float(result.diagnostics["equilibrium_residual"]) < EQUILIBRIUM_TOL

    _names, compositions, fractions, branches = _phases_of(result)
    z = np.asarray(mixture.composition.fractions, dtype=float)
    solved, beta, worst = _equilibrium_newton(eos, mixture, 333.0, z, compositions, branches)
    assert worst < NEWTON_TOL, worst
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(solved, compositions)) < COMPOSITION_TOL
    assert float(np.max(np.abs(beta - fractions))) < COMPOSITION_TOL

    returned_g = float(
        sum(
            weight * _reduced_g(eos, mixture, 333.0, x, branch)
            for weight, x, branch in zip(fractions, compositions, branches)
        )
    )
    for branch in ("liquid", "vapor"):
        assert returned_g < _reduced_g(eos, mixture, 333.0, z, branch)


# ---------------------------------------------------------------------------
# (ii) The Peng-Robinson three-liquid states: the multiphase log-space stage
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("feed", "temperature_K"), PR_THREE_LIQUID_STATES)
def test_the_peng_robinson_three_liquid_states_converge_and_are_verified(
    feed: tuple[float, float, float], temperature_K: float
) -> None:
    """Three liquids from a cubic, each verified by an 8-equation Newton.

    Whether the *model* is right at these conditions is not asserted - what is
    asserted is that the returned set is an equilibrium of that model, that it
    has the lowest reduced Gibbs energy of every admissible alternative, and
    that it came out of the stage ADR-0029 added.
    """
    eos = ct.PengRobinsonEOS()
    mixture = _mixture(TERNARY, feed)
    result = _flash(TERNARY, feed, temperature_K, eos)

    assert sorted(result.phases) == ["liquid1", "liquid2", "liquid3"]
    assert result.diagnostics["post_split_stable"] is True
    assert float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL
    assert float(result.diagnostics["equilibrium_residual"]) < EQUILIBRIUM_TOL
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert float(result.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0
    # The repair is the log-space stage, and it is what met the tolerance.
    assert result.diagnostics["converged_stage"] == "second-order-log"
    assert int(result.diagnostics["log_space_iterations"]) >= 1

    names, compositions, fractions, branches = _phases_of(result)
    z = np.asarray(mixture.composition.fractions, dtype=float)
    solved, beta, worst = _equilibrium_newton(
        eos, mixture, temperature_K, z, compositions, branches
    )
    assert worst < NEWTON_TOL, worst
    assert max(float(np.max(np.abs(a - b))) for a, b in zip(solved, compositions)) < COMPOSITION_TOL
    assert float(np.max(np.abs(beta - fractions))) < COMPOSITION_TOL

    returned_g = float(
        sum(
            weight * _reduced_g(eos, mixture, temperature_K, x, branch)
            for weight, x, branch in zip(fractions, compositions, branches)
        )
    )
    assert returned_g < _reduced_g(eos, mixture, temperature_K, z, "liquid")
    assert returned_g < _reduced_g(eos, mixture, temperature_K, z, "vapor")

    # Every two-phase subset of the three, re-solved here on its own. A subset
    # whose lever rule puts the feed outside the tie line is not a candidate at
    # all; the ones that are all sit above the three-phase answer.
    admissible = 0
    for first, second in itertools.combinations(range(3), 2):
        pair, pair_beta, pair_residual = _equilibrium_newton(
            eos,
            mixture,
            temperature_K,
            z,
            [compositions[first], compositions[second]],
            [branches[first], branches[second]],
            perturb=0.0,
        )
        if pair_residual > NEWTON_TOL or not np.all((pair_beta > 1e-9) & (pair_beta < 1.0 - 1e-9)):
            continue
        if float(np.max(np.abs(pair[0] - pair[1]))) < 1e-6:
            continue
        admissible += 1
        pair_g = float(
            sum(
                weight * _reduced_g(eos, mixture, temperature_K, x, branch)
                for weight, x, branch in zip(pair_beta, pair, [branches[first], branches[second]])
            )
        )
        assert returned_g < pair_g, (names[first], names[second], pair_g, returned_g)
    assert admissible >= 1, "no two-phase competitor was admissible, so the ordering is vacuous"


# ---------------------------------------------------------------------------
# The two rules, on their own
# ---------------------------------------------------------------------------


def test_the_removal_ranking_starts_at_the_phase_the_search_always_removed() -> None:
    """`_ranked_removals` must extend the old rule, not replace it."""
    fractions = np.array([0.4, -0.2, 0.9, -0.7])
    ranked = _ranked_removals(fractions)
    assert ranked[0] == int(np.argmin(fractions))
    assert ranked == (3, 1)
    # A set with every fraction positive offers nothing to remove.
    assert _ranked_removals(np.array([0.3, 0.7])) == ()
    # A tie keeps the earlier phase, so the choice is deterministic.
    assert _ranked_removals(np.array([-0.5, -0.5, 0.1]))[0] == 0


def test_a_removal_can_be_taken_back_and_each_candidate_is_tried_once() -> None:
    """The undo stack terminates: a candidate is consumed when it is taken."""
    surfaces: list[Callable[[np.ndarray], np.ndarray]] = [lambda x: x, lambda x: x, lambda x: x]
    undo = [
        _RemovalChoice(
            labels=["liquid", "liquid", "vapor"],
            surfaces=surfaces,
            compositions=[np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.5, 0.5])],
            remaining=[0, 2],
        )
    ]
    first = _take_next_removal(undo)
    assert first is not None
    assert first[0] == ["liquid", "vapor"]
    assert undo[0].remaining == [2]

    second = _take_next_removal(undo)
    assert second is not None
    assert second[0] == ["liquid", "liquid"]
    assert undo[0].remaining == []

    assert _take_next_removal(undo) is None
    assert undo == []


def test_the_reference_phase_is_the_one_that_holds_every_component() -> None:
    """The log-space stage's reference choice, which is what makes it work."""
    active = np.array([True, True, True])
    water_rich = np.array([1.0 - 2e-12, 1e-12, 1e-12])
    balanced = np.array([0.2, 0.5, 0.3])
    other = np.array([0.1, 0.8, 0.1])
    assert reference_phase_index([water_rich, balanced, other], active) == 1
    assert reference_phase_index([balanced, water_rich], active) == 0
    # A component absent from the feed never decides the choice.
    absent = np.array([True, True, False])
    assert reference_phase_index([np.array([0.5, 0.5, 0.0]), water_rich], absent) == 0


def test_a_state_that_already_converged_carries_no_log_space_diagnostics() -> None:
    """Dormancy: the ADR-0029 key appears only where the stage actually ran."""
    eos = ct.PengRobinsonEOS()
    result = _flash(TERNARY, (0.2, 0.4, 0.4), 280.0, eos)
    assert sorted(result.phases) == ["liquid1", "liquid2", "liquid3"]
    assert result.diagnostics["converged_stage"] == "successive-substitution"
    assert "log_space_iterations" not in result.diagnostics
