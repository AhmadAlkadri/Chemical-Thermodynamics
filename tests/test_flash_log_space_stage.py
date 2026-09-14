"""The log-space split stage (ADR-0024, validation Case P-14).

``chemthermo.flash._log_space`` carries the two-phase Gibbs minimization in
``u = ln n`` instead of in the mole numbers themselves. What that buys is one
thing: a phase composition whose smallest component is ``exp(-450)`` - the
polymer content of the vapour above a polyethylene melt - which the linear
stage cannot represent, cannot difference and cannot step towards.

What is checked here:

- the seed. Michelsen's ``W`` at an essentially pure polymer melt has
  K-values spanning ``e^460``, and the Rachford-Rice window they define is
  **empty**; the log-space seed is admissible componentwise for any such ``K``
  by construction, and that is asserted rather than asserted-about.
- the stage. The same answer from four unrelated starting points, a converged
  point that really is a local minimum of the Gibbs energy, and a refusal when
  the seed itself is outside the box.
- the two entry points. The stationary-point seed (no successive substitution
  at all) and the hand-over from a failed linear stage, which reach the same
  split to thirteen figures from different directions.
- the reporting convention: a mole fraction outside the exponential's range is
  an exact ``0.0`` in the composition, with its logarithm in the diagnostics.
- that none of it touches an ordinary state.

The polymer parameters come from ``tests/fixtures/pcsaft/martini2009_polymers.json``;
read its provenance block before reading any number here.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.core import Composition
from chemthermo.eos import PCSAFTEOS
from chemthermo.exceptions import ModelError
from chemthermo.flash import _detect
from chemthermo.flash._common import wilson_k
from chemthermo.flash._log_space import (
    _SEED_PHASE_FRACTION,
    TRACE_MOLE_FRACTION,
    _safeguarded_directions,
    has_trace_component,
    lever_rule_phase_fraction,
    log_space_seed,
    log_space_split,
    seed_from_iterate,
    stability_seed_ladder,
)
from chemthermo.flash._split import (
    _PhaseRoot,
    _rachford_rice,
    _rachford_rice_window,
    _solve_k_loop,
)
from chemthermo.flash.settings import FlashSettings
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

FIXTURE = Path(__file__).parent / "fixtures" / "pcsaft" / "martini2009_polymers.json"

TEMPERATURE_K = 453.0
PE_MW_G_MOL = 16400.0
PENTANE_MW_G_MOL = 72.146
KIJ = -0.006


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


def _parameters(mw_g_mol: float = PE_MW_G_MOL) -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
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


def _mixture(weight_fraction: float = 0.05, mw_g_mol: float = PE_MW_G_MOL) -> ct.Mixture:
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


def _eos(kij: float = KIJ, mw_g_mol: float = PE_MW_G_MOL) -> PCSAFTEOS:
    return PCSAFTEOS(parameters=_parameters(mw_g_mol), kij=kij)


def _pinned_roots(pressure_Pa: float, mixture: ct.Mixture, eos: PCSAFTEOS):
    """The two per-phase root holders the flash would build at this state."""
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
    )
    assert stability.status == "unstable"
    assert stability.trial_composition is not None
    z = np.asarray(mixture.fractions, dtype=float)
    w = np.asarray(stability.trial_composition, dtype=float)
    _k, incipient = _detect._stability_k_seed(
        mixture, TEMPERATURE_K, pressure_Pa, z=z, w=w, tpd_min=float(stability.tpd_min)
    )
    roots_x, roots_y = _detect._phi_phi_roots(
        eos,
        mixture,
        TEMPERATURE_K,
        pressure_Pa,
        feed_branch=stability.feed_branch,
        incipient_branch=stability.phase_branch,
        incipient_phase=incipient,
        seed_label="stability",
    )
    return stability, w, incipient, roots_x, roots_y


# ---------------------------------------------------------------------------
# The seed
# ---------------------------------------------------------------------------


def test_the_stationary_points_k_values_bracket_no_vapor_fraction_at_all() -> None:
    """The premise of ADR-0024 decision 2, as numbers.

    The tangent-plane minimizer at 1 MPa is an essentially pure polymer melt,
    its K-values span ``1e+180`` and every one of them is below one - so the
    Rachford-Rice function has no root in ``[0, 1]`` *and* the extended
    Leibovici-Neoschil window is empty. Successive substitution has nowhere to
    start, for a feed the same stability test has just proved unstable.
    """
    mixture = _mixture()
    z = np.asarray(mixture.fractions, dtype=float)
    stability = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=1.0e6, eos=_eos())
    assert stability.status == "unstable"
    assert stability.trial_composition is not None
    w = np.asarray(stability.trial_composition, dtype=float)
    assert w[1] < 1e-100  # measured 2.35e-178: the melt holds no solvent at all here

    k_seed, _incipient = _detect._stability_k_seed(
        mixture, TEMPERATURE_K, 1.0e6, z=z, w=w, tpd_min=float(stability.tpd_min)
    )
    assert float(np.max(k_seed)) < 1.0
    assert _rachford_rice(z, k_seed)[0] is None
    assert _rachford_rice_window(k_seed) is None

    assert has_trace_component(w, z > 0.0)
    assert TRACE_MOLE_FRACTION == 1e-30


def test_the_log_space_seed_is_inside_the_two_phase_box_for_any_k() -> None:
    """``0 < n_i < z_i`` componentwise, including for K-values that are not doubles.

    The seed is ``n_i = beta K_i z_i / ((1 - beta) + beta K_i)``, whose ratio to
    ``z_i`` is a number in ``(0, 1)`` for every positive ``K``; the assertion
    below is that the *implementation* preserves that through the logarithms.
    """
    z = np.array([2.3148042e-04, 0.99976852])
    for exponent in (-4000.0, -460.0, -1.0, 0.0, 1.0, 460.0, 4000.0):
        w = np.array([1.0, math.exp(min(exponent, 0.0))])
        u = log_space_seed(z=z, w=w, tpd_min=-exponent, incipient_vapor=False)
        assert np.all(np.isfinite(u))
        moles = np.exp(u)
        assert np.all(moles >= 0.0)
        assert np.all(moles < z)
    # A stationary point with an exact zero gives ln K = -inf; the clamp keeps
    # the seed finite rather than propagating it.
    u = log_space_seed(z=z, w=np.array([1.0, 0.0]), tpd_min=-10.0, incipient_vapor=False)
    assert np.all(np.isfinite(u))


def test_the_seed_ignores_components_absent_from_the_feed() -> None:
    z = np.array([0.4, 0.6, 0.0])
    u = log_space_seed(z=z, w=np.array([0.5, 0.5, 0.0]), tpd_min=-1.0, incipient_vapor=True)
    assert u[2] == -math.inf
    assert np.all(np.isfinite(u[:2]))


def test_seed_from_iterate_pulls_a_negative_flash_back_into_the_box() -> None:
    z = np.array([0.3, 0.7])
    u = seed_from_iterate(z=z, x_ii=np.array([0.2, 0.8]), beta=-6.3e10)
    moles = np.exp(u)
    assert np.all(moles > 0.0)
    assert np.all(moles < z)


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------

#: The 1 MPa answer, to the digits the stage reproduces from four seeds.
MELT_SOLVENT_1MPA = 0.9714683526646756
LN_Y_POLYMER_1MPA = -450.5307940616726


def test_the_stage_finds_the_same_split_from_four_unrelated_seeds() -> None:
    """Seed-independence, which is what makes the answer a property of the state.

    Three of the four seeds are trivial splits ``n = f z`` that know nothing
    about the model; the fourth is the stationary-point seed the flash actually
    uses. All four must reach the same melt composition and the same
    ``ln y_polymer``, a quantity 450 orders of magnitude below anything the
    linear stage can hold.
    """
    mixture = _mixture()
    eos = _eos()
    z = np.asarray(mixture.fractions, dtype=float)
    stability, w, incipient, roots_x, roots_y = _pinned_roots(1.0e6, mixture, eos)
    settings = FlashSettings()

    seeds = {
        "stability-w": log_space_seed(
            z=z, w=w, tpd_min=float(stability.tpd_min), incipient_vapor=incipient == "vapor"
        ),
        "half": np.log(0.5 * z),
        "nine-tenths": np.log(0.9 * z),
        "ninety-nine-hundredths": np.log(0.99 * z),
    }
    for label, u0 in seeds.items():
        split = log_space_split(
            z=z,
            u0=u0,
            terms_i=roots_x.ln_fugacity_terms,
            terms_ii=roots_y.ln_fugacity_terms,
            settings=settings,
        )
        assert split.residual < 1e-11, label
        assert split.x_i[1] == pytest.approx(MELT_SOLVENT_1MPA, rel=1e-12), label
        assert split.ln_x_ii[0] == pytest.approx(LN_Y_POLYMER_1MPA, rel=1e-12), label
        assert 0.0 < split.beta < 1.0, label


def test_the_converged_point_is_a_minimum_of_the_two_phase_gibbs_energy() -> None:
    """Not merely a stationary point, and in the one direction the energy can see it.

    The gradient in ``u`` is ``n_k r_k`` by the chain rule, so a vanishing
    equal-fugacity residual makes the point stationary. Walking away from it
    along the **solvent** coordinate must raise the energy - that is what
    distinguishes the equilibrium from the trivial solution, which is a
    stationary ridge at the feed's own energy.

    Along the **polymer** coordinate the energy cannot see anything: the
    curvature there is ``1/n`` with ``n ~ 1e-196``, so a step of ``e^0.001``
    changes ``g`` by about ``1e-196`` against a ``g`` of ``0.24``, which is
    nothing in doubles. That is precisely why this stage exists, and it is why
    the minimality of that coordinate is asserted on the *residual* instead -
    the quantity that is still accurate there.
    """
    mixture = _mixture()
    eos = _eos()
    z = np.asarray(mixture.fractions, dtype=float)
    _stability, _w, _incipient, roots_x, roots_y = _pinned_roots(1.0e6, mixture, eos)
    split = log_space_split(
        z=z,
        u0=np.log(0.9 * z),
        terms_i=roots_x.ln_fugacity_terms,
        terms_ii=roots_y.ln_fugacity_terms,
        settings=FlashSettings(),
    )

    def state(moles: np.ndarray) -> tuple[float, np.ndarray]:
        reference = z - moles
        composition_i = reference / float(np.sum(reference))
        ln_composition_ii = np.log(moles) - math.log(float(np.sum(moles)))
        activity_i = np.log(composition_i) + roots_x.ln_fugacity_terms(composition_i)
        activity_ii = ln_composition_ii + roots_y.ln_fugacity_terms(np.exp(ln_composition_ii))
        return float(reference @ activity_i + moles @ activity_ii), activity_ii - activity_i

    converged = split.beta * split.x_ii
    base_energy, base_residual = state(converged)
    assert float(np.max(np.abs(base_residual))) < 1e-11

    solvent = 1
    for step in (+1e-3, -1e-3):
        moved = converged.copy()
        moved[solvent] *= math.exp(step)
        assert state(moved)[0] > base_energy

    polymer = 0
    for step in (+1e-3, -1e-3):
        moved = converged.copy()
        moved[polymer] *= math.exp(step)
        assert state(moved)[0] >= base_energy
        assert abs(state(moved)[1][polymer]) > abs(base_residual[polymer])


def test_the_stage_refuses_a_seed_that_is_not_a_two_phase_split() -> None:
    mixture = _mixture()
    eos = _eos()
    z = np.asarray(mixture.fractions, dtype=float)
    _stability, _w, _incipient, roots_x, roots_y = _pinned_roots(1.0e6, mixture, eos)
    with pytest.raises(ModelError, match="cannot start from this seed"):
        log_space_split(
            z=z,
            u0=np.log(2.0 * z),  # phase II holds more than the feed has
            terms_i=roots_x.ln_fugacity_terms,
            terms_ii=roots_y.ln_fugacity_terms,
            settings=FlashSettings(),
        )


# ---------------------------------------------------------------------------
# The two entry points, through the public API
# ---------------------------------------------------------------------------


def test_flash_tp_solves_the_polymer_vle_from_the_stationary_point_seed() -> None:
    """Entry point 1 (ADR-0024 decision 2): no successive substitution at all."""
    mixture = _mixture()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=1.0e6, eos=_eos())
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert diagnostics["k_seed"] == "stability-log"
    assert diagnostics["converged_stage"] == "second-order-log"
    assert diagnostics["log_space_seed"] == "stability-w"
    assert diagnostics["ssi_iterations"] == 0
    assert diagnostics["log_space_ln_x_min_component"] == "Polyethylene"
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(LN_Y_POLYMER_1MPA, rel=1e-12)
    assert result.phases["liquid"].composition.fractions[1] == pytest.approx(
        MELT_SOLVENT_1MPA, rel=1e-12
    )


def test_a_failed_linear_stage_hands_over_to_log_space(monkeypatch: pytest.MonkeyPatch) -> None:
    """Entry point 2 (ADR-0024 decision 1), and an independent check of the answer.

    Forcing the Wilson K-seed sends the same state down the *other* route: one
    successive-substitution step, a linear second-order stage that spends its
    whole budget without converging, and only then the log-space stage, seeded
    from that failed iterate rather than from the stationary point. It reaches
    the same split - which is two seeds and two code paths agreeing to thirteen
    figures, not one solver repeated.
    """
    mixture = _mixture()
    eos = _eos()

    def no_root(*_args: object, **_kwargs: object) -> tuple[None, float, float]:
        return None, 0.0, 0.0

    # The stability seed is refused a bracket, so `_flash_tp_tangent_plane`
    # takes its documented Wilson fallback; the trace-component gate is
    # side-stepped by leaving the fallback's own bracket intact.
    calls = {"n": 0}
    original = _detect._rachford_rice

    def once(z: np.ndarray, K: np.ndarray):
        calls["n"] += 1
        if calls["n"] == 1:
            return no_root()
        return original(z, K)

    monkeypatch.setattr(_detect, "_rachford_rice", once)
    monkeypatch.setattr(_detect, "has_trace_component", lambda *_a, **_k: False)

    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=1.0e6, eos=eos)
    diagnostics = result.diagnostics
    assert diagnostics["k_seed"] == "wilson"
    assert diagnostics["converged_stage"] == "second-order-log"
    assert diagnostics["log_space_seed"] == "linear-iterate"
    assert int(diagnostics["ssi_iterations"]) >= 1
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.phases["liquid"].composition.fractions[1] == pytest.approx(
        MELT_SOLVENT_1MPA, rel=1e-12
    )
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(LN_Y_POLYMER_1MPA, rel=1e-12)


def test_the_wilson_seeded_split_really_does_fail_without_the_log_space_stage() -> None:
    """The premise of the test above, measured rather than assumed."""
    mixture = _mixture()
    eos = _eos()
    z = np.asarray(mixture.fractions, dtype=float)
    stability, w, incipient, roots_x, roots_y = _pinned_roots(1.0e6, mixture, eos)
    del stability, w
    settings = FlashSettings(second_order=False)
    k_seed = wilson_k(mixture, TEMPERATURE_K, 1.0e6)
    vapor_fraction, _f0, _f1 = _rachford_rice(z, k_seed)
    assert vapor_fraction is not None
    del incipient
    split = _solve_k_loop(
        mixture,
        TEMPERATURE_K,
        1.0e6,
        eos=eos,
        activity_model=None,
        mode="phi-phi",
        settings=settings,
        z=z,
        K=k_seed,
        vapor_fraction=vapor_fraction,
        max_iter=settings.max_iter,
        allow_unconverged=True,
        extended_rachford_rice=True,
        roots_x=roots_x,
        roots_y=roots_y,
    )
    assert not split.converged
    assert split.iterations == 1  # it stops after one step, as Case P-13 (vi) recorded


# ---------------------------------------------------------------------------
# The zero-mole-fraction convention (ADR-0024 decision 3)
# ---------------------------------------------------------------------------


def test_a_mole_fraction_outside_the_exponentials_range_is_reported_as_zero() -> None:
    """Mw = 53 000 at 2 MPa: ``ln y_polymer = -1315``, so ``y_polymer`` is ``0.0``.

    The mole fraction is not lost - its logarithm is in the diagnostics - and
    the mass balance is *exact* rather than approximate, because the melt then
    holds every mole of polymer the feed had.
    """
    mixture = _mixture(0.05, 53000.0)
    result = ct.flash_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=2.0e6, eos=_eos(KIJ, 53000.0)
    )
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.phases["vapor"].composition.fractions[0] == 0.0
    assert int(diagnostics["log_space_zero_fractions"]) == 1
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(-1314.985239549, rel=1e-9)
    assert math.exp(float(diagnostics["log_space_ln_x_min"])) == 0.0
    assert diagnostics["log_space_ln_x_min_component"] == "Polyethylene"
    # `fugacity_residual` is taken over the components present in *both*
    # phases, so the polymer drops out of it; the stage's own residual is the
    # one that says the polymer's condition holds too.
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["log_space_residual"]) < 1e-8
    assert float(diagnostics["mass_balance_residual"]) < 1e-15
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"
    assert result.phases["liquid"].composition.fractions[0] == pytest.approx(
        2.7736965292611865e-03, rel=1e-9
    )


def test_a_composition_accepts_an_exact_zero() -> None:
    """The convention's precondition, stated where a reader will look for it."""
    composition = Composition(fractions=(0.0, 1.0))
    assert composition.fractions == (0.0, 1.0)


@pytest.mark.parametrize("pressure_Pa", [5.0e5, pytest.param(1.0e6, marks=pytest.mark.slow)])
def test_the_longest_chain_below_1_mpa_is_now_in_reach(pressure_Pa: float) -> None:
    """The former limitation, pinned the other way round (ADR-0025, Case P-15).

    Until ADR-0025 this test asserted the opposite: for ``Mw = 53 000`` at 0.5
    and 1 MPa the tangent-plane test did not find the melt at all - its deepest
    stationary point was a shallow vapour-side one (``tpd`` of order 1e-04) -
    so this stage was seeded 1300 orders of magnitude from the answer and spent
    its budget. It was a property of the **stability** iteration, not of the
    split: the melt's stationary point sits at ``ln W_polymer ~ 1450`` and the
    iteration clamped ``ln W`` at ``700``. With the clamp gone the melt is the
    minimizer, and this stage - unchanged - converges from it.

    What this module owns is the *stage*; the stability side of the same
    reversal is in `tests/test_stability_log_space.py`.
    """
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
    )
    assert stability.status == "unstable"
    assert stability.tpd_min < -1000.0  # the melt, not the shallow vapour-side point
    assert stability.phase_branch == "liquid"

    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["k_seed"] == "stability-log"
    assert result.diagnostics["converged_stage"] == "second-order-log"
    assert float(result.diagnostics["fugacity_residual"]) < 1e-8
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"


# ---------------------------------------------------------------------------
# Dormancy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "case",
    [
        "pr-methane-ethane",
        # another model on the same check, which the two unmarked cases cover
        pytest.param("pcsaft-water-hexane", marks=pytest.mark.slow),
        "pcsaft-polymer-lle",
    ],
)
def test_the_log_space_stage_is_dormant_on_ordinary_states(case: str) -> None:
    """No state that converged before ADR-0024 carries a single log-space key."""
    if case == "pr-methane-ethane":
        mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
        result = ct.flash_tp(
            mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS()
        )
    elif case == "pcsaft-water-hexane":
        mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5])
        result = ct.flash_tp(mixture, temperature_K=298.15, pressure_Pa=101325.0, eos=PCSAFTEOS())
    else:
        result = ct.flash_tp(_mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())

    assert not [key for key in result.diagnostics if key.startswith("log_space_")]
    assert result.diagnostics.get("converged_stage") != "second-order-log"
    assert result.diagnostics["k_seed"] != "stability-log"


# ---------------------------------------------------------------------------
# The curvature safeguard (ADR-0026, validation Case P-16)
# ---------------------------------------------------------------------------

#: ``Mw = 53 000`` at 0.3 MPa: the one state where the ADR-0024 iteration
#: parks next to the *trivial* solution instead of converging. The numbers are
#: measured at HEAD 584c508 and reproduced by
#: ``examples/validation/22_stability_log_space.py``.
STALL_PRESSURE_PA = 3.0e5
STALL_RESIDUAL = 3.900817e-05
#: The answer the safeguarded stage reaches, confirmed by a one-dimensional
#: equal-fugacity solve that shares no code with the flash (the melt's solvent
#: mole fraction agrees to 1.4e-14).
STALL_MELT_SOLVENT = 0.9631916542803276


_STALL_CACHE: list[tuple[np.ndarray, np.ndarray, _PhaseRoot, _PhaseRoot]] = []


def _stall_setup() -> tuple[np.ndarray, np.ndarray, _PhaseRoot, _PhaseRoot]:
    """The seed and the two root holders of the 0.3 MPa state.

    Cached: building it runs a stability test, and the two tests below want the
    *same* seed, which is the point they are making.
    """
    if _STALL_CACHE:
        return _STALL_CACHE[0]
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    z = np.asarray(mixture.fractions, dtype=float)
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=STALL_PRESSURE_PA, eos=eos
    )
    assert stability.trial_composition is not None
    w = np.asarray(stability.trial_composition, dtype=float)
    ln_capital_w = _detect._stationary_point_ln_capital_w(stability)
    _k, incipient = _detect._stability_k_seed(
        mixture,
        TEMPERATURE_K,
        STALL_PRESSURE_PA,
        z=z,
        w=w,
        tpd_min=float(stability.tpd_min),
        ln_capital_w=ln_capital_w,
    )
    u0 = log_space_seed(
        z=z,
        w=w,
        tpd_min=float(stability.tpd_min),
        incipient_vapor=incipient == "vapor",
        ln_capital_w=ln_capital_w,
    )
    roots_x, roots_y = _detect._phi_phi_roots(
        eos,
        mixture,
        TEMPERATURE_K,
        STALL_PRESSURE_PA,
        feed_branch=stability.feed_branch,
        incipient_branch=stability.phase_branch,
        incipient_phase=incipient,
        seed_label="stability-log",
    )
    _STALL_CACHE.append((z, u0, roots_x, roots_y))
    return _STALL_CACHE[0]


def test_the_unsafeguarded_stage_parks_next_to_the_trivial_solution() -> None:
    """The defect, measured rather than described.

    ``-r`` is a descent direction for ``g`` but its *length* is ``|r|``: once
    the residual is small the step is small, whatever distance is left. Here
    the Newton direction is rejected at every iterate (the Gibbs Hessian is
    indefinite next to the trivial solution) and the stage spends its whole
    budget with the two phases at the same composition.
    """
    z, u0, roots_x, roots_y = _stall_setup()
    split = log_space_split(
        z=z,
        u0=u0,
        terms_i=roots_x.ln_fugacity_terms,
        terms_ii=roots_y.ln_fugacity_terms,
        settings=FlashSettings(),
    )
    assert split.iterations == FlashSettings().second_order_max_iter
    assert split.residual == pytest.approx(STALL_RESIDUAL, rel=1e-5)
    assert split.residual > FlashSettings().tol
    # Both phases are still within a few percent of the feed in the solvent,
    # which is what "next to the trivial solution" means here: the polymer-rich
    # phase holds 7.5e-05 polymer where the answer is 3.7e-02, five hundred
    # times more, and the phase fraction is 0.042 where the answer is 0.998.
    assert split.x_i[1] == pytest.approx(float(z[1]), rel=1e-3)
    assert split.x_i[0] < 1e-4 < 1e-2 < 1.0 - STALL_MELT_SOLVENT
    assert split.beta == pytest.approx(0.0417, rel=1e-2)


def test_the_safeguard_converges_the_same_seed_from_the_same_place() -> None:
    """Same seed, same model, same tolerance - only the direction rule changes."""
    z, u0, roots_x, roots_y = _stall_setup()
    split = log_space_split(
        z=z,
        u0=u0,
        terms_i=roots_x.ln_fugacity_terms,
        terms_ii=roots_y.ln_fugacity_terms,
        settings=FlashSettings(),
        curvature_safeguard=True,
    )
    assert split.residual < 1e-12
    assert split.iterations <= 20
    assert split.x_i[1] == pytest.approx(STALL_MELT_SOLVENT, rel=1e-12)
    assert 0.0 < split.beta < 1.0
    # The melt really is polymer-rich and the incipient phase really is not.
    assert split.x_i[0] > 100.0 * float(z[0])
    assert split.ln_x_ii[0] < -1000.0


def test_the_safeguard_leaves_a_converging_state_on_the_same_answer() -> None:
    """Continuity: where the default rule converges, the safeguard agrees.

    It is a different sequence of iterates - it is allowed to be - so this is
    an agreement to the stage's own tolerance, not a bit-identity claim. The
    bit-identity claim is that the flag is *off* unless the default rule has
    already failed, which `test_the_safeguard_is_reported_and_is_off_by_default`
    pins on the flash.
    """
    mixture = _mixture()
    eos = _eos()
    z = np.asarray(mixture.fractions, dtype=float)
    stability, w, incipient, roots_x, roots_y = _pinned_roots(1.0e6, mixture, eos)
    u0 = log_space_seed(
        z=z, w=w, tpd_min=float(stability.tpd_min), incipient_vapor=incipient == "vapor"
    )
    split = log_space_split(
        z=z,
        u0=u0,
        terms_i=roots_x.ln_fugacity_terms,
        terms_ii=roots_y.ln_fugacity_terms,
        settings=FlashSettings(),
        curvature_safeguard=True,
    )
    assert split.residual < 1e-11
    assert split.x_i[1] == pytest.approx(MELT_SOLVENT_1MPA, rel=1e-12)
    assert split.ln_x_ii[0] == pytest.approx(LN_Y_POLYMER_1MPA, rel=1e-12)


def test_every_safeguarded_direction_descends_and_the_list_ends_with_ssi() -> None:
    """The contract of :func:`_safeguarded_directions`, on a Jacobian by hand.

    ``J`` below is the log-space Jacobian of a state whose Gibbs Hessian has
    one negative eigenvalue - the shape measured at the 0.3 MPa stall - so the
    Newton direction is *not* a descent direction and is dropped, while the
    modified-Newton direction of the same Hessian is kept.
    """
    jacobian = np.array([[0.8058676, -0.7265849], [-3.549692e-06, -2.360102e-06]])
    residual = np.array([-8.602786e-05, -3.898835e-05])
    moles = np.array([2.031340e-07, 4.158566e-02])
    gradient = moles * residual

    hessian = 0.5 * (
        (moles[:, None] * jacobian + np.diag(gradient))
        + (moles[:, None] * jacobian + np.diag(gradient)).T
    )
    assert float(np.min(np.linalg.eigvalsh(hessian))) < 0.0  # indefinite

    newton = np.linalg.solve(jacobian, -residual)
    assert float(gradient @ newton) > 0.0  # not a descent direction

    directions = _safeguarded_directions(jacobian=jacobian, residual=residual, moles=moles)
    assert len(directions) >= 2
    for direction, slope in directions:
        assert slope < 0.0
        assert slope == pytest.approx(float(gradient @ direction), rel=1e-12)
        assert not np.allclose(direction, newton)
    # The last one is always the successive-substitution step.
    assert np.array_equal(directions[-1][0], -residual)

    # Without a Jacobian there is only that step.
    assert len(_safeguarded_directions(jacobian=None, residual=residual, moles=moles)) == 1


def test_the_safeguard_is_off_on_a_state_that_never_needed_it() -> None:
    """`flash_tp` says which rule produced the answer, and here it is the old one.

    The other half of this statement - the flag reported `True`, on a state that
    raised before ADR-0026 - is pinned with that state's answer in
    `tests/test_pcsaft_polymer.py::test_the_0_3_mpa_vapour_liquid_split_the_stage_used_to_stall_on`.
    """
    result = ct.flash_tp(
        _mixture(0.05, 53000.0),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=2.0e6,
        eos=_eos(KIJ, 53000.0),
    )
    assert result.diagnostics["converged_stage"] == "second-order-log"
    assert result.diagnostics["log_space_curvature_safeguard"] is False


# ---------------------------------------------------------------------------
# The seed ladder (ADR-0028)
# ---------------------------------------------------------------------------


def test_the_lever_rule_fraction_is_the_bound_a_non_negative_complement_allows() -> None:
    """``lambda <= min_i z_i / w_i``, and the complement at it is non-negative.

    This is the whole content of :func:`lever_rule_phase_fraction`: if the
    ``w``-phase holds a fraction ``lambda`` of one mole of feed, the other
    phase holds ``(z - lambda w) / (1 - lambda)``, and that has to be a
    composition. Nothing about the model enters.
    """
    z = np.array([1.37497633e-05, 0.999986250])
    w = np.array([1.0, 1.7e-197])  # an essentially pure polymer melt

    # The convention follows `log_space_seed`: when `w` is phase II the
    # fraction returned *is* ``lambda``, and when it is phase I the complement's
    # fraction is returned instead.
    lam = lever_rule_phase_fraction(z=z, w=w, incipient_vapor=True)
    assert lam == float(np.min(z / w))
    assert lever_rule_phase_fraction(z=z, w=w, incipient_vapor=False) == 1.0 - lam

    complement = (z - lam * w) / (1.0 - lam)
    assert np.all(complement >= 0.0)
    # A hair past the bound and the complement goes negative, which is what
    # makes this a bound and not a guess.
    assert np.any((z - 1.01 * lam * w) < 0.0)

    # A stationary point at the feed gives the neutral answer back: the bound
    # is 1, so the complement's fraction is 0, and the clamp keeps it inside.
    edge = lever_rule_phase_fraction(z=z, w=z, incipient_vapor=True)
    assert 0.0 < edge < 1.0


def test_the_lever_rule_seed_is_inside_the_two_phase_box_like_the_neutral_one() -> None:
    """The property that makes any ``beta`` in ``(0, 1)`` admissible, at this one."""
    z = np.array([2.3148042e-04, 0.99976852])
    w = np.array([1.0, 1e-200])
    beta = lever_rule_phase_fraction(z=z, w=w, incipient_vapor=False)
    u = log_space_seed(z=z, w=w, tpd_min=-455.0, incipient_vapor=False, phase_fraction=beta)
    assert np.all(np.isfinite(u))
    moles = np.exp(u)
    assert np.all(moles > 0.0)
    assert np.all(moles < z)


def test_the_ladders_first_entry_is_the_seed_the_route_already_formed() -> None:
    """Entry 1 is bit-identical to the pre-ADR-0028 seed, which is why it is skippable.

    ``_phi_phi_log_space`` runs entries 1 and 2 itself and walks the ladder from
    entry 3; that is only sound if entry 1 is the very array it was given.
    """
    z = np.array([1.37497633e-05, 0.999986250])
    w = np.array([1.0, 2.3e-178])
    ladder = stability_seed_ladder(z=z, w=w, tpd_min=-1492.0, incipient_vapor=False)
    neutral = log_space_seed(z=z, w=w, tpd_min=-1492.0, incipient_vapor=False)

    assert len(ladder) == 3
    assert np.array_equal(ladder[0][0], neutral)
    assert np.array_equal(ladder[1][0], neutral)
    assert [safeguard for _seed, safeguard in ladder] == [False, True, True]
    # Entry 3 is a different seed, and it is the lever-rule one.
    assert not np.array_equal(ladder[2][0], neutral)
    assert np.array_equal(
        ladder[2][0],
        log_space_seed(
            z=z,
            w=w,
            tpd_min=-1492.0,
            incipient_vapor=False,
            phase_fraction=lever_rule_phase_fraction(z=z, w=w, incipient_vapor=False),
        ),
    )
    assert _SEED_PHASE_FRACTION == 0.5
