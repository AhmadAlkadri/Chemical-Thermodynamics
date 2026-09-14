"""Liquid-liquid equilibrium from an equation of state (ADR-0019).

The defect this pins: before this slice ``chemthermo.flash._split`` evaluated
one phase of a phi-phi split with ``phase="liquid"`` and the other with
``phase="vapor"``, so the two phases could never share a density branch. At
298.15 K and 1 atm water / n-hexane has both a liquid and a vapour root at most
compositions, the true answer is two *liquids* (the two pure vapour pressures
sum to about 23 kPa, far below 1 atm), and the split therefore converged on a
water-rich liquid against a hexane-rich vapour whose Gibbs energy is **above**
the feed's - correctly refused by the post-split stability test, but refused
rather than solved. Validation Case P-7(iii) recorded it as a limitation.

What replaces it: each phase is pinned to the branch the tangent-plane
stability test found *that phase* on, so the pair may be
``("liquid", "liquid")``, and the two converged phases are then named from
``EquationOfState.phase_identity`` (ADR-0017) - ``liquid1`` / ``liquid2`` with
``vapor_fraction = None`` when both measure as liquids.

Nothing here is compared against a number copied out of the implementation.
The tie line is checked against the lever rule and against a second and third
feed on the same tie line; the phase identities are re-derived from a finite
difference of the public ``PCSAFTEOS.pressure_Pa``; the Gibbs decrease and the
stability of each converged phase are recomputed from the public API. The
comparison against an independent implementation (FeOs) lives in
``tests/validation/test_pcsaft_lle_vs_feos.py``.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD

NAMES = ("Water", "n-Hexane")
TEMPERATURE_K = 298.15
ATMOSPHERE_PA = 101325.0
HIGH_PRESSURE_PA = 1.0e6

#: The phi-phi grid of ``tests/test_flash_phase_detection.py``, which is also
#: the 144 Peng-Robinson states of the bit-identity fixture. Duplicated rather
#: than imported so this module stays self-contained.
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


def _mixture(names: Sequence[str] = NAMES, z: Sequence[float] = (0.5, 0.5)) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


@lru_cache(maxsize=None)
def _lle_result(pressure_Pa: float, z: tuple[float, ...]) -> ct.FlashResult:
    """One PC-SAFT water / n-hexane flash, cached for the module.

    `flash_tp` is deterministic for fixed inputs, models and settings (an
    invariant this file also tests directly), and one of these costs about four
    seconds, so the tests that only need to *read* a converged result share
    one. ``test_the_liquid_liquid_answer_is_deterministic`` deliberately does
    not use this.
    """
    return ct.flash_tp(
        _mixture(z=z),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PCSAFTEOS(),
    )


def _pcsaft_kappa(
    eos: ct.PCSAFTEOS,
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float] | np.ndarray,
) -> float:
    """``kappa = P / (rho dP/drho)`` at the phase's own root, independently.

    Built from the **public** ``density_roots`` and ``pressure_Pa`` with a
    central difference, never from ``phase_identity``'s own analytic
    derivative, so this is a cross-check of the label and not a restatement of
    it (the same construction ``tests/test_phase_identity.py`` uses).
    """
    values = np.asarray(composition, dtype=float).tolist()
    roots = eos.density_roots(
        mixture=mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=values,
    )
    # The liquid-like root is the largest; a single root is both.
    density = roots[-1]
    step = 1e-4 * density
    bound = ct.PCSAFTEOS(components=tuple(mixture.component_names))
    plus = bound.pressure_Pa(
        temperature_K=temperature_K,
        density_mol_m3=density + step,
        composition=values,
    )
    minus = bound.pressure_Pa(
        temperature_K=temperature_K,
        density_mol_m3=density - step,
        composition=values,
    )
    slope = (plus - minus) / (2.0 * step)
    return pressure_Pa / (density * slope)


def _reduced_g(fractions: np.ndarray, ln_phi: np.ndarray) -> float:
    mask = fractions > 0.0
    return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_phi[mask])))


def _ln_phi(
    eos: ct.EquationOfState,
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float] | np.ndarray,
    phase: str,
) -> np.ndarray:
    return np.log(
        np.array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=np.asarray(composition, dtype=float).tolist(),
                phase=phase,
            ),
            dtype=float,
        )
    )


# --------------------------------------------------------------------------
# A. Water / n-hexane at 1 atm: the state the old split could not express.
# --------------------------------------------------------------------------


def test_water_hexane_at_one_atm_returns_two_liquids() -> None:
    """Case P-8(i): the defect of Case P-7(iii), resolved."""
    eos = ct.PCSAFTEOS()
    mixture = _mixture()

    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
    )
    assert stability.status == "unstable"
    # The premise of the defect: at 1 atm there really are two density roots,
    # so the old fixed liquid/vapour pairing had a vapour root to land on.
    roots = eos.density_roots(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        composition=[0.5, 0.5],
    )
    assert len(roots) == 2

    result = _lle_result(ATMOSPHERE_PA, (0.5, 0.5))
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert diagnostics["phase_regime"] == "LLE"
    assert diagnostics["phase_label_method"] == "compressibility"
    assert diagnostics["phase_i_branch"] == "liquid"
    assert diagnostics["phase_ii_branch"] == "liquid"

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-9
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"

    water_rich = list(result.phases["liquid1"].composition.fractions)
    hexane_rich = list(result.phases["liquid2"].composition.fractions)
    assert water_rich[0] > 0.99 and hexane_rich[1] > 0.99

    # Both phases are liquids by a kappa recomputed here from the public
    # pressure routine, not by asking `phase_identity` again.
    for composition in (water_rich, hexane_rich):
        kappa = _pcsaft_kappa(eos, mixture, TEMPERATURE_K, ATMOSPHERE_PA, composition)
        assert 0.0 < kappa < KAPPA_LIQUID_THRESHOLD, kappa
        # ... and each sits on a liquid-like density, not the 43 mol/m^3
        # vapour root the old split paired it with.
        density = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=composition,
        )[-1]
        assert density > 5000.0

    # The Gibbs decrease, recomputed from the public API on the roots the two
    # phases actually converged on.
    z = np.array(mixture.composition.fractions, dtype=float)
    beta = result.phase_fractions["liquid2"]
    x = np.array(water_rich, dtype=float)
    y = np.array(hexane_rich, dtype=float)
    feed_liquid = _ln_phi(eos, mixture, TEMPERATURE_K, ATMOSPHERE_PA, z, "liquid")
    feed_vapor = _ln_phi(eos, mixture, TEMPERATURE_K, ATMOSPHERE_PA, z, "vapor")
    feed_g = min(_reduced_g(z, feed_liquid), _reduced_g(z, feed_vapor))
    delta_g = (
        beta * _reduced_g(y, _ln_phi(eos, mixture, TEMPERATURE_K, ATMOSPHERE_PA, y, "liquid"))
        + (1.0 - beta)
        * _reduced_g(x, _ln_phi(eos, mixture, TEMPERATURE_K, ATMOSPHERE_PA, x, "liquid"))
        - feed_g
    )
    assert delta_g < 0.0
    assert delta_g == pytest.approx(float(diagnostics["delta_g_split_rt"]), abs=1e-12)


@pytest.mark.slow  # ADR-0020 runtime trim: the 1-atm tie line is covered by test_water_hexane_at_one_atm_returns_two_liquids
def test_the_same_tie_line_comes_back_from_three_feeds_and_obeys_the_lever_rule() -> None:
    """Case P-8(ii): a tie line is a property of the state, not of the feed."""
    reference: list[np.ndarray] | None = None

    for feed in ((0.5, 0.5), (0.2, 0.8), (0.8, 0.2)):
        result = _lle_result(ATMOSPHERE_PA, feed)
        assert sorted(result.phases) == ["liquid1", "liquid2"]
        phases = [
            np.array(result.phases[name].composition.fractions, dtype=float)
            for name in ("liquid1", "liquid2")
        ]
        if reference is None:
            reference = phases
        else:
            for found, expected in zip(phases, reference):
                assert found == pytest.approx(expected, abs=1e-10)

        # Lever rule, restated on the reported fractions.
        z = np.array(feed, dtype=float) / sum(feed)
        beta = result.phase_fractions["liquid2"]
        assert np.max(np.abs(z - ((1.0 - beta) * phases[0] + beta * phases[1]))) < 1e-12


@pytest.mark.slow  # ADR-0020 runtime trim: naming is also asserted by the two fast tie-line tests
def test_liquid1_is_the_phase_richer_in_the_first_component() -> None:
    """ADR-0019 decision 3: the LLE order is composition-based, so it is stable."""
    eos = ct.PCSAFTEOS()
    result = ct.flash_tp(
        _mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
    )
    first = result.phases["liquid1"].composition.fractions[0]
    second = result.phases["liquid2"].composition.fractions[0]
    assert first > second

    # Reversing the component order reverses which phase carries the label,
    # and nothing else: the phase *set* is the same pair of compositions.
    reversed_result = ct.flash_tp(
        _mixture(names=NAMES[::-1], z=(0.5, 0.5)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        eos=eos,
    )
    assert sorted(reversed_result.phases) == ["liquid1", "liquid2"]
    assert (
        reversed_result.phases["liquid1"].composition.fractions[0]
        > (reversed_result.phases["liquid2"].composition.fractions[0])
    )
    forward = {
        tuple(round(value, 12) for value in result.phases[name].composition.fractions)
        for name in result.phases
    }
    backward = {
        tuple(round(value, 12) for value in reversed(phase.composition.fractions))
        for phase in reversed_result.phases.values()
    }
    assert forward == backward


# --------------------------------------------------------------------------
# B. The 1 MPa tie line, previously mislabelled by the Wilson fallback.
# --------------------------------------------------------------------------


def test_at_one_megapascal_the_pair_is_now_named_liquid1_liquid2() -> None:
    """Case P-8(iii): the same tie line, without the ``wilson-ranking`` fallback.

    Above about 0.6 MPa the isotherm has a single density root, so ``"vapor"``
    and ``"liquid"`` name the same root and the pre-slice machinery already
    found this split - but both phases measured as liquids, the two-phase
    naming rule fell through to the Wilson ranking, and ``vapor_fraction`` was
    the hexane-rich **liquid**'s fraction. That mislabel is what ADR-0019
    removes.
    """
    eos = ct.PCSAFTEOS()
    mixture = _mixture()
    for composition in ([0.5, 0.5], [0.999, 0.001], [0.001, 0.999]):
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=HIGH_PRESSURE_PA,
            composition=composition,
        )
        assert len(roots) == 1

    result = ct.flash_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=HIGH_PRESSURE_PA, eos=eos
    )
    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_regime"] == "LLE"
    assert result.diagnostics["phase_label_method"] == "compressibility"
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"

    # The tie line barely moves between 1 atm and 1 MPa, as a liquid one should.
    at_one_atm = _lle_result(ATMOSPHERE_PA, (0.5, 0.5))
    for name in ("liquid1", "liquid2"):
        here = np.array(result.phases[name].composition.fractions, dtype=float)
        there = np.array(at_one_atm.phases[name].composition.fractions, dtype=float)
        assert np.max(np.abs(here - there)) < 1e-5


# --------------------------------------------------------------------------
# C. Peng-Robinson, the same binary.
# --------------------------------------------------------------------------


def test_peng_robinson_water_hexane_is_a_verified_liquid_liquid_split() -> None:
    """Case P-8(iv): what Peng-Robinson with ``kij = 0`` does at the same state.

    Reported, not forced: nothing here asserts a solubility, only that whatever
    the model returns is a verified equilibrium of *that* model - a Gibbs
    decrease, equal fugacities on the roots the phases sit on, and every phase
    stable when re-tested.
    """
    eos = ct.PengRobinsonEOS()
    mixture = _mixture()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos)

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_regime"] == "LLE"
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert float(result.diagnostics["fugacity_residual"]) < 1e-9
    assert result.diagnostics["post_split_status"] == "stable"

    # Each converged phase, re-tested on its own: stable or marginal, never
    # unstable. (`_post_split_status` says the same; this asks again from the
    # public API so the check is not a restatement of a diagnostics key.)
    for phase in result.phases.values():
        single = ct.Mixture.from_database(
            list(NAMES), list(phase.composition.fractions), normalize=True
        )
        verdict = ct.stability_tp(
            single, temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
        )
        assert verdict.tpd_min > -1e-6, verdict.tpd_min
        assert (
            eos.phase_identity(
                mixture=mixture,
                temperature_K=TEMPERATURE_K,
                pressure_Pa=ATMOSPHERE_PA,
                composition=list(phase.composition.fractions),
                phase="liquid",
            )
            == "liquid"
        )


# --------------------------------------------------------------------------
# D. The audit: vapour-liquid behaviour is untouched.
# --------------------------------------------------------------------------


def test_the_pinned_root_is_the_historical_branch_on_the_whole_phi_phi_grid(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """ADR-0019's bit-identity evidence, stated as a property rather than a fixture.

    ``tests/test_flash_refactor_bit_identity.py`` proves the *outputs* of the
    155 pinned states are unchanged. This proves *why*: on every one of the 144
    Peng-Robinson phi-phi states in that fixture, at **every** iterate of the
    split, the fugacity coefficients each phase was evaluated with are exactly
    (``==``, not ``approx``) the ones the pre-ADR-0019 fixed assignment -
    ``phase="liquid"`` for phase I, ``phase="vapor"`` for phase II - would have
    produced. If any state ever stops satisfying this, the fixture has to be
    audited state by state and regenerated, not relaxed.
    """
    from chemthermo.flash import _detect
    from chemthermo.flash._common import EosBranchTerms
    from chemthermo.flash._split import _PhaseRoot

    created: list["_Recording"] = []

    class _Recording(_PhaseRoot):
        def __init__(self, *args: object, **kwargs: object) -> None:
            super().__init__(*args, **kwargs)  # type: ignore[arg-type]
            self.seen: list[tuple[np.ndarray, np.ndarray]] = []
            created.append(self)

        def branch_terms(self, composition: np.ndarray) -> EosBranchTerms:
            terms = super().branch_terms(composition)
            # ADR-0022's log-space guard must stay dormant here: Peng-Robinson
            # fugacity coefficients on this grid are all representable, so the
            # split takes the same `phi_l / phi_v` it always took.
            assert terms.phi is not None
            self.seen.append((np.asarray(composition, dtype=float).copy(), terms.phi.copy()))
            return terms

    monkeypatch.setattr(_detect, "_PhaseRoot", _Recording)

    eos = ct.PengRobinsonEOS()
    two_phase_states = 0
    evaluations = 0
    for names, z in GRID_MIXTURES:
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                created.clear()
                mixture = _mixture(names, z)
                result = ct.flash_tp(
                    mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=eos
                )
                if len(result.phases) == 1:
                    assert not created
                    continue
                two_phase_states += 1
                assert len(created) == 2
                for holder, historical in zip(created, ("liquid", "vapor")):
                    assert holder.seen
                    for composition, phi in holder.seen:
                        evaluations += 1
                        expected = np.array(
                            eos.fugacity_coefficients(
                                mixture=mixture,
                                temperature_K=temperature_K,
                                pressure_Pa=pressure_Pa,
                                composition=composition.tolist(),
                                phase=historical,
                            ),
                            dtype=float,
                        )
                        assert np.array_equal(phi, expected), (
                            names,
                            temperature_K,
                            pressure_Pa,
                            historical,
                        )

    assert two_phase_states == 47, two_phase_states
    assert evaluations > 100, evaluations


def test_a_vapor_liquid_result_carries_none_of_the_new_diagnostics_keys() -> None:
    """ADR-0019 decision 4: the branch keys are conditional, deliberately."""
    result = ct.flash_tp(
        _mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert result.vapor_fraction is not None
    for key in ("phase_i_branch", "phase_ii_branch"):
        assert key not in result.diagnostics


def test_a_model_without_phase_identity_keeps_the_wilson_ranking_fallback() -> None:
    """The documented last resort is still reachable and still recorded."""

    class _Unmeasured(ct.PengRobinsonEOS):
        """Peng-Robinson with ADR-0017's measurement switched off."""

        def phase_identity(self, **kwargs: object) -> str | None:
            return None

    result = ct.flash_tp(
        _mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=_Unmeasured(),
    )
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["phase_label_method"] == "wilson-ranking"
    assert result.diagnostics["phase_regime"] == "VLE"


@pytest.mark.slow  # ADR-0020 runtime trim: determinism of the phi-phi path is covered on the PR grid
def test_the_liquid_liquid_answer_is_deterministic() -> None:
    """Same inputs, same models, same settings - same doubles, every time."""
    eos = ct.PCSAFTEOS()
    runs = [
        ct.flash_tp(_mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos)
        for _ in range(2)
    ]
    for result in runs[1:]:
        assert result.phase_names() == runs[0].phase_names()
        assert dict(result.phase_fractions) == dict(runs[0].phase_fractions)
        assert dict(result.diagnostics) == dict(runs[0].diagnostics)
        for name, phase in result.phases.items():
            assert phase.composition.fractions == runs[0].phases[name].composition.fractions
