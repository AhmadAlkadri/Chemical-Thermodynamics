"""Tangent-plane stability in log mole numbers (ADR-0025, validation Case P-15).

Michelsen's stability iteration carries the *unnormalized* mole numbers ``W``,
and the only place they were ever needed as doubles is the normalization
``w = W / sum_j W_j``. Before ADR-0025 that normalization was protected by a
clamp, ``ln W <- clip(ln W, -700, 700)``, which is where ``exp`` stops
existing. The clamp is not a safeguard on a *result*: a stationary point whose
``ln W`` legitimately sits outside that window is simply not reachable, the
iteration parks on the boundary, and the trial is reported as
``second_order_no_progress``. A polyethylene melt against a solvent-vapour feed
needs ``ln W_polymer ~ 1450``.

What is checked here:

- **the gate.** The log-space arithmetic runs where, and only where, the old
  clamp would have engaged; ``StabilityTrial.log_space`` says which route a
  trial took, and every state in this repository that had an answer before
  ADR-0025 still takes the old one. That is asserted over the 144-state
  Peng-Robinson stability grid and the two activity families, not argued.
- **the algebra.** ``tpd = -ln sum_W`` (equation (7)) and the trivial-solution
  test still hold when ``sum_W`` itself has overflowed to ``inf``.
- **the capability, without a polymer.** A synthetic equation of state whose
  ``ln phi`` is large enough to put a stationary point outside the window,
  built so that its stationary point is known in closed form.
- **the capability, with one.** ``Mw = 53 000`` polyethylene in n-pentane at
  453 K and 0.5 / 1 MPa: the melt is found, and ``flash_tp`` returns the split
  it seeds. This retires the pinned miss of validation Case P-14; the deeper
  checks on that split (an independent 1-D solve, FeOs's potentials) live in
  `examples/validation/22_stability_log_space.py`.
- **the hand-over.** ``ln W`` reaches the log-space split stage directly, and
  it agrees with the ``ln w - tpd`` reconstruction it replaces wherever that
  reconstruction is defined at all.

The polymer parameters come from ``tests/fixtures/pcsaft/martini2009_polymers.json``;
read its provenance block before reading any number here.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.flash._detect import _stationary_point_ln_capital_w
from chemthermo.models import EquationOfState
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord
from chemthermo.stability.tp import _LN_W_MAX, _LN_W_MIN, _logsumexp, _normalize

FIXTURE = Path(__file__).parent / "fixtures" / "pcsaft" / "martini2009_polymers.json"

TEMPERATURE_K = 453.0
PENTANE_MW_G_MOL = 72.146
KIJ = -0.006


def _grid_module():
    """The 144-state Peng-Robinson grid, loaded by path.

    By path, not by ``from tests... import``: CI runs the ``pytest`` console
    script, which does not put the working directory on ``sys.path``
    (`.agents/dev-contract.md`).
    """
    path = Path(__file__).parent / "test_stability_eos_surfaces.py"
    spec = importlib.util.spec_from_file_location("_stability_grid", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# The polymer system (identical to `tests/test_flash_log_space_stage.py`)
# ---------------------------------------------------------------------------


def _polymer_record(mw_g_mol: float) -> PCSAFTRecord:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    row = next(p for p in payload["polymers"] if p["name"] == "Polyethylene")
    return PCSAFTRecord(
        name="Polyethylene",
        segments_per_g=row["segments_per_g"],
        MW_g_mol=mw_g_mol,
        sigma_A=row["sigma_A"],
        epsilon_k_K=row["epsilon_k_K"],
        source="Martini et al. 2009 Table 1 citing Gross & Sadowski 2002; see fixture",
    )


def _eos(mw_g_mol: float, kij: float = KIJ) -> PCSAFTEOS:
    parameters = PCSAFTParameters.from_records(
        [
            _polymer_record(mw_g_mol),
            PCSAFTRecord(
                name="n-Pentane",
                m=2.6896,
                sigma_A=3.7729,
                epsilon_k_K=231.20,
                MW_g_mol=PENTANE_MW_G_MOL,
            ),
        ]
    )
    return PCSAFTEOS(parameters=parameters, kij=kij)


def _mixture(mw_g_mol: float, weight_fraction: float = 0.05) -> ct.Mixture:
    polymer = weight_fraction / mw_g_mol
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return ct.Mixture.from_components(
        [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=mw_g_mol / 1000.0,
                formula="(C2H4)n",
                volatile=False,
                source="see tests/fixtures/pcsaft/martini2009_polymers.json",
            ),
            ct.Component.from_database("n-Pentane"),
        ],
        [polymer / total, solvent / total],
        normalize=True,
    )


# ---------------------------------------------------------------------------
# 1. The two normalizations, and the gate between them
# ---------------------------------------------------------------------------


def test_logsumexp_matches_the_direct_sum_where_the_direct_sum_exists() -> None:
    values = np.array([-3.25, 0.0, 12.5, -40.0])
    assert _logsumexp(values) == pytest.approx(math.log(float(np.sum(np.exp(values)))), rel=1e-15)


def test_logsumexp_carries_an_absent_component_and_an_all_absent_input() -> None:
    assert _logsumexp(np.array([-math.inf, 0.0])) == pytest.approx(0.0, abs=1e-15)
    assert _logsumexp(np.array([-math.inf, -math.inf])) == -math.inf


def test_the_gate_runs_the_old_arithmetic_whenever_the_old_clamp_was_dormant() -> None:
    """Bit-identity by construction, stated as an executable claim.

    Inside ``[-700, 700]`` the pre-ADR-0025 ``np.clip`` returned its argument,
    so the expression ``w = exp(ln W) / sum_j exp(ln W_j)`` is what ran then and
    is what runs now - the same operations in the same order on the same
    doubles.
    """
    active = np.array([True, True, True])
    for values in (
        np.array([0.0, -1.0, 2.0]),
        np.array([_LN_W_MIN, 0.0, _LN_W_MAX]),
        np.array([-699.5, 699.5, 12.0]),
    ):
        normalized = _normalize(values, active)
        assert normalized is not None
        w, sum_capital_w, ln_sum_capital_w, log_space = normalized
        assert not log_space
        expected_capital_w = np.exp(values)
        expected_sum = float(np.sum(expected_capital_w))
        assert sum_capital_w == expected_sum
        assert ln_sum_capital_w == math.log(expected_sum)
        assert np.array_equal(w, expected_capital_w / expected_sum)


def test_the_gate_switches_exactly_when_a_clip_would_have_engaged() -> None:
    active = np.array([True, True])
    inside = _normalize(np.array([_LN_W_MAX, 0.0]), active)
    outside = _normalize(np.array([np.nextafter(_LN_W_MAX, math.inf), 0.0]), active)
    assert inside is not None and outside is not None
    assert not inside[3]
    assert outside[3]


def test_log_space_normalization_survives_mole_numbers_that_exp_cannot_hold() -> None:
    """``ln W = (1452.2, 3.6)``: the melt of validation Case P-15.

    ``sum_W`` overflows and ``w`` underflows, and neither is a failure: the
    composition is a genuine ``(1.0, 0.0)`` to the last bit a double has, and
    the quantity the verdict is read from - ``ln sum_W`` - is an ordinary
    number.
    """
    values = np.array([1452.2, 3.6])
    normalized = _normalize(values, np.array([True, True]))
    assert normalized is not None
    w, sum_capital_w, ln_sum_capital_w, log_space = normalized
    assert log_space
    assert ln_sum_capital_w == pytest.approx(1452.2, abs=1e-12)
    assert sum_capital_w == math.inf
    assert tuple(w.tolist()) == (1.0, 0.0)


def test_an_all_absent_or_infinite_sum_is_still_refused() -> None:
    active = np.array([True, True])
    assert _normalize(np.array([-math.inf, -math.inf]), active) is None
    assert _normalize(np.array([math.inf, 0.0]), active) is None


def test_a_sum_below_the_smallest_double_is_normalized_rather_than_lost() -> None:
    """The other end of the same window, where the old clamp also engaged.

    ``sum_W = exp(-799.3)`` is an exact ``0.0``, which the old code could not
    divide by; in logs the composition is still ``(0.5, 0.5)`` to round-off.
    """
    normalized = _normalize(np.array([-800.0, -800.0]), np.array([True, True]))
    assert normalized is not None
    w, sum_capital_w, ln_sum_capital_w, log_space = normalized
    assert log_space
    assert sum_capital_w == 0.0
    assert ln_sum_capital_w == pytest.approx(-800.0 + math.log(2.0), abs=1e-12)
    assert w.tolist() == pytest.approx([0.5, 0.5], abs=1e-12)


# ---------------------------------------------------------------------------
# 2. The capability without a polymer: a synthetic equation of state
# ---------------------------------------------------------------------------


class _TwoBranchLnPhiEOS(EquationOfState):
    """An EOS whose ``ln phi`` is a fixed vector per phase label.

    Not a physical model and not meant to be one. It is the smallest thing that
    puts a tangent-plane stationary point outside the exponential's range, and
    its value is that the answer is known in closed form: with ``ln phi``
    independent of composition, equation (5) is solved by one substitution,
    ``ln W_i = d_i - ln phi_i``, and equation (7) then gives
    ``tpd = -logsumexp(ln W)`` exactly. The test below therefore checks the
    solver against arithmetic written out here, not against another solver, and
    it needs no polymer and no PC-SAFT.

    ``exp(-1500)`` is an exact ``0.0``, so ``fugacity_coefficients`` cannot
    express the liquid branch and ``log_fugacity_coefficients`` (the ADR-0022
    guard) is what carries it - the same division of labour a real model has
    here.
    """

    name = "two-branch-constant-ln-phi"

    def __init__(
        self, ln_phi_vapor: Sequence[float], ln_phi_liquid: Sequence[float] | None = None
    ) -> None:
        self._vapor = tuple(float(value) for value in ln_phi_vapor)
        self._liquid = (
            self._vapor if ln_phi_liquid is None else tuple(float(value) for value in ln_phi_liquid)
        )

    def _terms(self, phase: str) -> tuple[float, ...]:
        return self._liquid if phase == "liquid" else self._vapor

    def fugacity_coefficients(
        self,
        *,
        mixture: ct.Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> tuple[float, ...]:
        del mixture, temperature_K, pressure_Pa, composition
        with np.errstate(under="ignore"):
            return tuple(float(value) for value in np.exp(np.array(self._terms(phase))))

    def log_fugacity_coefficients(
        self,
        *,
        mixture: ct.Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> tuple[float, ...]:
        del mixture, temperature_K, pressure_Pa, composition
        return self._terms(phase)


def _synthetic_stability(
    ln_phi_vapor: Sequence[float],
    ln_phi_liquid: Sequence[float] | None = None,
    z: Sequence[float] = (0.25, 0.75),
) -> ct.StabilityResult:
    return ct.stability_tp(
        ct.Mixture.from_database(["Methane", "Ethane"], list(z), normalize=True),
        temperature_K=200.0,
        pressure_Pa=1.0e6,
        eos=_TwoBranchLnPhiEOS(ln_phi_vapor, ln_phi_liquid),
    )


def test_one_composition_independent_branch_leaves_the_feed_stable() -> None:
    """The control for the test below.

    With a single ``ln phi`` for both labels, ``d_i - ln phi_i(w) = ln z_i`` at
    every ``w``, so every trial lands on ``W = z`` - the trivial solution,
    ``tpd = 0``. Nothing leaves the exponential's range, and nothing claims to,
    even though ``ln phi = -1500`` makes ``phi`` itself an exact zero.
    """
    result = _synthetic_stability((-1500.0, -2.0))
    assert result.status == "stable"
    assert all(trial.trivial for trial in result.trials if trial.converged)
    assert not any(trial.log_space for trial in result.trials)
    assert "log_space_trial_count" not in result.diagnostics


def test_the_log_space_route_reaches_a_stationary_point_at_ln_w_of_1500() -> None:
    """No polymer, no PC-SAFT: two constant branches and one known answer.

    With the liquid branch constant, one substitution solves equation (5), so
    the stationary point is ``ln W_i = ln z_i + ln phi_i^V(z) - ln phi_i^L`` and
    ``tpd = -logsumexp(ln W)`` exactly. Both are checked against arithmetic
    done here, and ``log_space`` records that the trial could not have got
    there under the pre-ADR-0025 clamp.
    """
    z = np.array([0.25, 0.75])
    ln_phi_vapor = np.array([-1.0, -2.0])
    # Chosen so the *feed* is on the vapour branch (its reduced residual Gibbs
    # energy there is -1.75 against 0.0 on the liquid one), while component 1
    # is enormously favoured by the liquid branch.
    ln_phi_liquid = np.array([-1500.0, 500.0])
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], z.tolist(), normalize=True)
    result = ct.stability_tp(
        mixture,
        temperature_K=200.0,
        pressure_Pa=1.0e6,
        eos=_TwoBranchLnPhiEOS(ln_phi_vapor.tolist(), ln_phi_liquid.tolist()),
    )

    d = np.log(z) + ln_phi_vapor
    expected_ln_capital_w = d - ln_phi_liquid
    expected_tpd = -_logsumexp(expected_ln_capital_w)

    assert result.feed_branch == "vapor"
    assert result.status == "unstable"
    assert result.tpd_min == pytest.approx(expected_tpd, abs=1e-9)
    assert result.tpd_min < -1400.0
    assert result.trial_ln_W is not None
    assert np.allclose(np.array(result.trial_ln_W), expected_ln_capital_w, atol=1e-9)
    assert float(result.diagnostics["ln_sum_W"]) == pytest.approx(-expected_tpd, abs=1e-9)
    assert float(result.diagnostics["tpd_from_sum_W"]) == pytest.approx(expected_tpd, abs=1e-9)
    assert result.diagnostics["tm_at_stationary_point"] == -math.inf
    assert result.diagnostics["minimizing_trial_log_space"] is True
    assert int(result.diagnostics["log_space_trial_count"]) >= 1


def test_the_synthetic_case_is_deterministic_and_permutation_invariant() -> None:
    """Swapping the two components swaps the answer and changes nothing else."""
    z = [0.25, 0.75]
    ln_phi_vapor = [-1.0, -2.0]
    ln_phi_liquid = [-1500.0, 500.0]

    def run(order: Sequence[int]) -> ct.StabilityResult:
        names = [["Methane", "Ethane"][i] for i in order]
        return ct.stability_tp(
            ct.Mixture.from_database(names, [z[i] for i in order], normalize=True),
            temperature_K=200.0,
            pressure_Pa=1.0e6,
            eos=_TwoBranchLnPhiEOS(
                [ln_phi_vapor[i] for i in order], [ln_phi_liquid[i] for i in order]
            ),
        )

    forward = run((0, 1))
    reversed_ = run((1, 0))
    assert forward.tpd_min == pytest.approx(reversed_.tpd_min, abs=1e-9)
    assert run((0, 1)).tpd_min == forward.tpd_min  # exactly, not approximately
    assert forward.trial_ln_W is not None and reversed_.trial_ln_W is not None
    assert forward.trial_ln_W[0] == pytest.approx(reversed_.trial_ln_W[1], abs=1e-9)


# ---------------------------------------------------------------------------
# 3. The polymer: validation Case P-14's pinned miss, flipped
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "pressure_Pa",
    [5.0e5, pytest.param(1.0e6, marks=pytest.mark.slow)],
)
def test_the_melt_is_the_minimizer_for_the_longest_chain_below_1_mpa(pressure_Pa: float) -> None:
    """The reversal of `test_the_longest_chain_below_1_mpa_is_still_out_of_reach`.

    ``Mw = 53 000`` is ``m = 1393.9`` segments, and a melt of it against a
    solvent-vapour feed has ``ln phi_polymer = -1550``; the stationary point is
    therefore at ``ln W_polymer ~ 1450``, which the pre-ADR-0025 clamp put 750
    out of reach. The three liquid-surface trials used to spend 51 iterations
    parked on the clamp and end ``second_order_no_progress``; they now converge
    in 3 successive substitutions.

    The 1 MPa case is marked ``slow``: it is the same statement at a second
    pressure, and 0.5 MPa runs by default.
    """
    mixture = _mixture(53000.0)
    eos = _eos(53000.0)
    result = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)

    assert result.status == "unstable"
    assert result.tpd_min < -1000.0
    assert result.feed_branch == "vapor"
    assert result.phase_branch == "liquid"
    assert result.diagnostics["minimizing_trial_surface"] == "liquid"
    assert result.diagnostics["minimizing_trial_log_space"] is True
    # Equation (7) on the stationary point, which `sum_W` can no longer carry.
    assert float(result.diagnostics["tpd_from_sum_W"]) == pytest.approx(result.tpd_min, rel=1e-12)
    assert result.diagnostics["sum_W"] == math.inf
    assert result.diagnostics["tm_at_stationary_point"] == -math.inf

    # The melt is essentially pure polymer *as a normalized composition* - the
    # solvent's share is `exp(-1449)` - and `ln W` is where that magnitude
    # survives.
    assert result.trial_composition == (1.0, 0.0)
    assert result.trial_ln_W is not None
    assert result.trial_ln_W[0] > 1000.0
    assert math.isfinite(result.trial_ln_W[1])

    converged = [trial for trial in result.trials if trial.converged]
    assert len(converged) == 4
    melts = [trial for trial in result.trials if trial.log_space]
    assert sorted(trial.label for trial in melts) == [
        "pure-Polyethylene",
        "pure-n-Pentane",
        "wilson-liquid",
    ]
    assert all(trial.iterations <= 5 for trial in melts)
    assert all(trial.termination_reason == "stationarity_met" for trial in melts)


def test_the_melt_stationary_point_satisfies_its_own_equations() -> None:
    """Equations (5) and (7), re-derived here from the model rather than trusted."""
    mixture = _mixture(53000.0)
    eos = _eos(53000.0)
    result = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=5.0e5, eos=eos)
    z = np.asarray(mixture.fractions, dtype=float)
    assert result.trial_ln_W is not None
    ln_capital_w = np.asarray(result.trial_ln_W, dtype=float)
    w = np.asarray(result.trial_composition, dtype=float)

    def ln_phi(composition: np.ndarray, phase: str) -> np.ndarray:
        return np.asarray(
            eos.log_fugacity_coefficients(
                mixture=mixture,
                temperature_K=TEMPERATURE_K,
                pressure_Pa=5.0e5,
                composition=composition.tolist(),
                phase=phase,
            )
        )

    d = np.log(z) + ln_phi(z, "vapor")
    # Equation (5): `ln W_i + ln phi_i(w) - d_i = 0` on the trial's own surface.
    residual = ln_capital_w + ln_phi(w, "liquid") - d
    assert float(np.max(np.abs(residual))) < 1e-9
    # Equation (7): `tpd = -ln sum_i W_i`.
    assert -_logsumexp(ln_capital_w) == pytest.approx(result.tpd_min, rel=1e-12)


@pytest.mark.parametrize(
    "pressure_Pa",
    [5.0e5, pytest.param(1.0e6, marks=pytest.mark.slow)],
)
def test_flash_tp_returns_the_vapour_liquid_split_the_melt_seeds(pressure_Pa: float) -> None:
    """Case P-15: the split Case P-14 could not reach, now reached and verified.

    The deeper cross-checks - an independently written 1-D equal-fugacity
    solve and FeOs's chemical potentials - are in
    `examples/validation/22_stability_log_space.py`; what is pinned here is
    that ``flash_tp`` answers, that the answer carries its own invariants, and
    that the route it took is the log-space one seeded from the melt.
    """
    mixture = _mixture(53000.0)
    result = ct.flash_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=_eos(53000.0)
    )
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert diagnostics["k_seed"] == "stability-log"
    assert diagnostics["converged_stage"] == "second-order-log"
    assert float(diagnostics["tpd_min"]) < -1000.0
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["log_space_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"
    # The vapour is pure solvent to the last bit a double has; its polymer
    # content is in the diagnostics (ADR-0024 decision 3).
    assert result.phases["vapor"].composition.fractions[0] == 0.0
    assert float(diagnostics["log_space_ln_x_min"]) < -1400.0


def test_the_0_5_mpa_split_is_the_published_melt_composition() -> None:
    """One pinned number, so a silent drift in the answer is visible.

    5.9 wt% solvent in the melt at 0.5 MPa, which is the direction and the
    order of magnitude ADR-0024 reports for the shorter chain.
    """
    mixture = _mixture(53000.0)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=5.0e5, eos=_eos(53000.0))
    melt = result.phases["liquid"].composition.fractions
    assert melt[0] == pytest.approx(2.1187890289570684e-02, rel=1e-9)
    assert result.vapor_fraction == pytest.approx(0.996618853739762, rel=1e-9)
    solvent_mass = melt[1] * PENTANE_MW_G_MOL
    polymer_mass = melt[0] * 53000.0
    assert solvent_mass / (solvent_mass + polymer_mass) == pytest.approx(0.0592, abs=5e-4)


@pytest.mark.slow
# `slow`: a repetition of the 0.5 MPa state above at four further pressures on
# the same map, which the default run already covers.
@pytest.mark.parametrize("pressure_Pa", [4.0e5, 7.5e5, 1.2e6, 1.4e6])
def test_every_state_of_the_vapour_liquid_region_reaches_the_melt(pressure_Pa: float) -> None:
    """All of 0.4-1.4 MPa raised `ConvergenceError` before ADR-0025, not only 0.5 and 1.

    Case P-14 pinned the two states it had measured. The repair is the whole
    band: a 0.3-3.6 MPa sweep at 0.1 MPa steps puts 13 states in the
    "raised, now returns a verified split" column.
    """
    mixture = _mixture(53000.0)
    eos = _eos(53000.0)
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
    )
    assert stability.tpd_min < -1000.0
    assert stability.feed_branch == "vapor"

    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert float(result.diagnostics["fugacity_residual"]) < 1e-8
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"


def test_the_polymer_case_is_deterministic() -> None:
    first = ct.stability_tp(
        _mixture(53000.0), temperature_K=TEMPERATURE_K, pressure_Pa=5.0e5, eos=_eos(53000.0)
    )
    second = ct.stability_tp(
        _mixture(53000.0), temperature_K=TEMPERATURE_K, pressure_Pa=5.0e5, eos=_eos(53000.0)
    )
    assert first.tpd_min == second.tpd_min
    assert first.trial_ln_W == second.trial_ln_W
    assert first.trial_composition == second.trial_composition


def test_the_16400_chain_never_leaves_the_old_arithmetic() -> None:
    """Validation Case P-14's own system is untouched, and that is asserted.

    The shorter chain's melt sits at ``ln W`` of order ``400``, inside the
    window, so every trial of every state ADR-0024 measured runs the
    pre-ADR-0025 expressions character for character.
    """
    mixture = _mixture(16400.0)
    eos = _eos(16400.0)
    for pressure_Pa in (5.0e5, 1.0e6, 2.0e6):
        result = ct.stability_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
        )
        assert not any(trial.log_space for trial in result.trials)
        assert "log_space_trial_count" not in result.diagnostics
        for trial in result.trials:
            if trial.ln_W is None:
                continue
            finite = [value for value in trial.ln_W if math.isfinite(value)]
            assert max(abs(value) for value in finite) <= _LN_W_MAX


# ---------------------------------------------------------------------------
# 4. Dormancy
# ---------------------------------------------------------------------------


def test_the_peng_robinson_stability_grid_never_engages_the_log_space_route() -> None:
    """144 states: same verdicts, same ``tpd_min``, same ``w``, old arithmetic.

    The verdicts and compositions of this grid are pinned elsewhere
    (`tests/test_stability_eos_surfaces.py`, `tests/test_flash_refactor_bit_identity.py`);
    what this adds is the ADR-0025 statement those cannot make - that not one
    of the 144 states left the pre-ADR-0025 expressions, so their agreement is
    by construction and not by luck.
    """
    grid = _grid_module()
    eos = ct.PengRobinsonEOS()
    states = 0
    for names, z in grid.GRID_MIXTURES:
        for temperature_K in grid.GRID_T_K:
            for pressure_Pa in grid.GRID_P_PA:
                result = ct.stability_tp(
                    ct.Mixture.from_database(list(names), list(z), normalize=True),
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=eos,
                )
                states += 1
                assert not any(trial.log_space for trial in result.trials)
                assert "log_space_trial_count" not in result.diagnostics
    assert states == 144


def test_the_activity_families_never_engage_the_log_space_route() -> None:
    """An NRTL liquid-liquid test and a modified-Raoult one, for the same reason."""
    names = ("n-Butanol", "Water")
    # The packaged synthetic pair used by `tests/test_stability_activity.py`.
    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(names[0], names[1], 0.90047, 3.51307, 0.48, 0.48)]
        )
    )
    mixture = ct.Mixture.from_database(list(names), [0.5, 0.5], normalize=True)
    liquid_liquid = ct.stability_tp(
        mixture, temperature_K=298.15, pressure_Pa=101325.0, activity_model=model
    )
    modified_raoult = ct.stability_tp(
        mixture, temperature_K=350.0, pressure_Pa=101325.0, activity_model=model, vapor="ideal"
    )
    for result in (liquid_liquid, modified_raoult):
        assert not any(trial.log_space for trial in result.trials)
        assert "log_space_trial_count" not in result.diagnostics


def test_every_trial_records_ln_w_and_its_sum_consistently() -> None:
    """``ln_sum_W`` is ``log(sum_W)`` wherever ``sum_W`` is a positive double.

    That equality is what makes ``tpd_from_sum_W`` bit-identical to the
    pre-ADR-0025 ``-math.log(sum_W)`` on every state that had one.
    """
    result = ct.stability_tp(
        ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    for trial in result.trials:
        assert trial.ln_W is not None
        assert math.isfinite(trial.sum_W) and trial.sum_W > 0.0
        assert trial.ln_sum_W == math.log(trial.sum_W)
        assert trial.ln_W[0] >= _LN_W_MIN


# ---------------------------------------------------------------------------
# 5. The hand-over to the log-space split stage
# ---------------------------------------------------------------------------


def test_the_split_stage_is_seeded_from_ln_w_only_where_w_cannot_carry_it() -> None:
    """The gate on the seed, which is the gate on the arithmetic again."""
    melt = ct.stability_tp(
        _mixture(53000.0), temperature_K=TEMPERATURE_K, pressure_Pa=5.0e5, eos=_eos(53000.0)
    )
    assert _stationary_point_ln_capital_w(melt) is not None

    ordinary = ct.stability_tp(
        ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True),
        temperature_K=200.0,
        pressure_Pa=1.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    assert _stationary_point_ln_capital_w(ordinary) is None


def test_ln_w_agrees_with_the_reconstruction_it_replaces() -> None:
    """``ln W`` against ``ln w - tpd``, over the PR grid's unstable verdicts.

    The split seed used to be built from the second expression and is built
    from the first only where the second has stopped working (ADR-0025). On an
    **unstable** verdict - the only kind that seeds a split - the two are the
    same quantity, and the measurement below says how nearly: 5.4e-15 over the
    grid.

    On a *stable* verdict they can differ by whole units, and that is not new
    and not a defect: ``tpd_min`` is measured to the lower envelope of the
    phase candidates while ``ln sum_W`` is equation (7) on the surface the
    trial iterated on, and the two part company exactly where a pinned trial
    stops above the other candidate (the note in
    ``chemthermo.stability.tp._summarize``, ADR-0021). A stable verdict seeds
    nothing, so nothing reads either number there. The worst such gap is
    reported here rather than hidden, so that widening the gate would have to
    confront it.
    """
    grid = _grid_module()
    eos = ct.PengRobinsonEOS()
    worst = 0.0
    worst_stable = 0.0
    compared = 0
    for names, z in grid.GRID_MIXTURES:
        for temperature_K in grid.GRID_T_K:
            for pressure_Pa in grid.GRID_P_PA:
                result = ct.stability_tp(
                    ct.Mixture.from_database(list(names), list(z), normalize=True),
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=eos,
                )
                if result.trial_ln_W is None or result.trial_composition is None:
                    continue
                w = np.asarray(result.trial_composition, dtype=float)
                mask = w > 0.0
                rebuilt = np.log(w[mask]) - result.tpd_min
                gap = float(np.max(np.abs(rebuilt - np.asarray(result.trial_ln_W)[mask])))
                if result.status == "unstable":
                    compared += 1
                    worst = max(worst, gap)
                else:
                    worst_stable = max(worst_stable, gap)
    assert compared == 47  # the grid's unstable verdicts (Case P-11)
    # Measured 5.4e-15; the bound is loose because the claim is "the same
    # quantity", not a pinned number.
    assert worst < 1e-10
    # And the documented disagreement on stable verdicts is real, so that a
    # future reader does not mistake the bound above for a global one.
    assert worst_stable > 1.0
