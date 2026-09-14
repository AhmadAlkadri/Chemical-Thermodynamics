"""PC-SAFT for a polymer/solvent mixture (ADR-0022, validation Cases P-12, P-13).

The system is polyethylene / n-pentane at 453 K, with the parameters the cited
fixture ``tests/fixtures/pcsaft/martini2009_polymers.json`` carries. **Read
that file's provenance block before reading any number here**: the polymer
parameters come from one open secondary source citing a paywalled primary
table, they are not verified against that table, and chemthermo packages no
polymer parameters at all. Everything below is therefore a statement about
*this model with these inputs*, never about polyethylene.

What is checked here (the external cross-checks are in
``tests/validation/test_pcsaft_polymer_vs_feos.py``):

- the segments-per-mass record arithmetic, ``m = (m/M) * Mw``;
- the pure melt's density root and its compressibility identity;
- the ADR-0022 log-space guard: dormant on ordinary states, active and
  *correct* for a chain long enough that ``exp(ln phi)`` underflows;
- the liquid-liquid split itself - verified residuals, both phases liquid, a
  negative Gibbs change, post-split stable;
- the qualitative behaviour the source describes: an LCST-type switch with
  temperature, a cloud-point pressure rising with ``k_ij``;
- permutation invariance and determinism under a mass ratio of 1e5 : 1;
- and two states where the split does **not** converge, pinned as limitations
  rather than worked around.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.exceptions import ConvergenceError
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

FIXTURE = Path(__file__).parent / "fixtures" / "pcsaft" / "martini2009_polymers.json"

#: The exam state (validation Case P-13): 453 K, the k_ij the fixture's Table 3
#: gives for the Mw = 16400 sample, and a feed of 5 wt% polymer.
TEMPERATURE_K = 453.0
PE_MW_G_MOL = 16400.0
KIJ = -0.006
PENTANE_MW_G_MOL = 72.146


def _fixture() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _polymer_record(mw_g_mol: float) -> PCSAFTRecord:
    """The fixture's polyethylene row, as a record, with no ``m`` written out."""
    row = next(p for p in _fixture()["polymers"] if p["name"] == "Polyethylene")
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


def _polymer_component(mw_g_mol: float = PE_MW_G_MOL) -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=mw_g_mol / 1000.0,
        formula="(C2H4)n",
        volatile=False,
        source="see tests/fixtures/pcsaft/martini2009_polymers.json",
    )


def _weight_to_mole_fractions(weight_fraction: float, mw_g_mol: float) -> list[float]:
    """``(x_polymer, x_solvent)`` from a polymer **mass** fraction."""
    polymer = weight_fraction / mw_g_mol
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return [polymer / total, solvent / total]


def _mixture(
    weight_fraction: float = 0.05,
    mw_g_mol: float = PE_MW_G_MOL,
    *,
    reversed_order: bool = False,
) -> ct.Mixture:
    components = [_polymer_component(mw_g_mol), ct.Component.from_database("n-Pentane")]
    z = _weight_to_mole_fractions(weight_fraction, mw_g_mol)
    if reversed_order:
        components = components[::-1]
        z = z[::-1]
    return ct.Mixture.from_components(components, z, normalize=True)


def _eos(kij: float = KIJ, mw_g_mol: float = PE_MW_G_MOL) -> PCSAFTEOS:
    return PCSAFTEOS(parameters=_parameters(mw_g_mol), kij=kij)


def _phase_by_polymer(result: ct.FlashResult, polymer_index: int = 0) -> tuple[str, str]:
    """``(polymer-rich name, solvent-rich name)`` - labels are roles, not identities."""
    order = sorted(
        result.phases, key=lambda name: result.phases[name].composition.fractions[polymer_index]
    )
    return order[-1], order[0]


# ---------------------------------------------------------------------------
# The record arithmetic
# ---------------------------------------------------------------------------


def test_the_segment_number_is_the_mass_based_parameter_times_the_molar_mass() -> None:
    for mw, expected in ((16400.0, 431.32), (53000.0, 1393.9)):
        record = _parameters(mw).record("Polyethylene")
        assert record.m == pytest.approx(0.0263 * mw, rel=1e-15)
        assert record.m == pytest.approx(expected, rel=1e-12)
        assert record.MW_g_mol == mw


def test_the_fixture_states_its_provenance_caveat() -> None:
    """The parameters must never be used without the caveat travelling with them."""
    citation = _fixture()["citation"]
    assert citation["primary_source_read"] is False
    assert "NOT VERIFIED AGAINST THE PRIMARY TABLE" in citation["provenance_caveat"]
    assert "NOT packaged runtime data" in citation["redistribution_note"]


def test_no_polymer_is_packaged_as_runtime_data() -> None:
    from chemthermo.parameters import get_pcsaft_parameters

    packaged = get_pcsaft_parameters().names()
    assert not any("poly" in name for name in packaged)


# ---------------------------------------------------------------------------
# The pure melt (validation Case P-12 (ii))
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mw_g_mol", [16400.0, 53000.0])
def test_the_pure_melt_has_one_density_root_of_a_plausible_magnitude(mw_g_mol: float) -> None:
    """One liquid-like root over 1-30 MPa, and no overflow in the eta scan.

    The ~0.757-0.793 g/cm^3 window below is a **sanity bracket**, not an
    assertion against data: commonly tabulated polyethylene melt densities near
    450 K are around 0.77-0.80 g/cm^3, and that remark is unverified here.
    """
    eos = PCSAFTEOS(components=("Polyethylene",), parameters=_parameters(mw_g_mol))
    previous = 0.0
    for pressure_Pa in (1.0e6, 5.0e6, 1.0e7, 1.5e7, 2.0e7, 3.0e7):
        roots = eos.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=[1.0]
        )
        assert len(roots) == 1
        density_g_cm3 = roots[0] * mw_g_mol / 1.0e6
        assert 0.75 < density_g_cm3 < 0.80
        assert density_g_cm3 > previous  # a melt is compressed by pressure
        previous = density_g_cm3

        compressibility = eos.compressibility_factor(
            temperature_K=TEMPERATURE_K, density_mol_m3=roots[0], composition=[1.0]
        )
        assert math.isfinite(compressibility) and compressibility > 0.0


def test_the_melt_root_measures_as_a_liquid() -> None:
    mixture = ct.Mixture.from_components([_polymer_component()], [1.0])
    eos = _eos()
    bound = PCSAFTEOS(components=("Polyethylene",), parameters=_parameters())
    for pressure_Pa in (1.0e6, 1.0e7, 3.0e7):
        assert (
            eos.phase_identity(
                mixture=mixture,
                temperature_K=TEMPERATURE_K,
                pressure_Pa=pressure_Pa,
                composition=[1.0],
                phase="liquid",
            )
            == "liquid"
        )
        # ... and independently, from the public pressure routine.
        density = eos.density_roots(
            mixture=mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=[1.0]
        )[0]
        step = 1e-4 * density
        slope = (
            bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density + step, composition=[1.0]
            )
            - bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density - step, composition=[1.0]
            )
        ) / (2.0 * step)
        assert 0.0 < pressure_Pa / (density * slope) < KAPPA_LIQUID_THRESHOLD


# ---------------------------------------------------------------------------
# The ADR-0022 log-space guard
# ---------------------------------------------------------------------------


def _count_branch_evaluations(monkeypatch: pytest.MonkeyPatch, call: "object") -> tuple[int, int]:
    """Run ``call()`` and report ``(branch evaluations, of which in log space)``."""
    from chemthermo.flash import _common, _split
    from chemthermo.stability import _evaluator

    original = _common.eos_branch_terms
    counters = [0, 0]

    def counting(*args: object, **kwargs: object) -> _common.EosBranchTerms:
        terms = original(*args, **kwargs)  # type: ignore[arg-type]
        counters[0] += 1
        if terms.phi is None:
            counters[1] += 1
        return terms

    monkeypatch.setattr(_evaluator, "eos_branch_terms", counting)
    monkeypatch.setattr(_split, "eos_branch_terms", counting)
    call()  # type: ignore[operator]
    return counters[0], counters[1]


def test_a_long_chains_fugacity_coefficient_underflows_and_its_logarithm_does_not() -> None:
    """The measured failure the guard exists for, stated as a number."""
    parameters = _parameters(53000.0)
    eos = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=parameters, kij=KIJ)
    x = _weight_to_mole_fractions(0.05, 53000.0)
    density = eos.density_roots(temperature_K=TEMPERATURE_K, pressure_Pa=1.0e7, composition=x)[-1]
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
    )
    assert ln_phi[0] < -1600.0  # measured: -1690.6
    # The reason: the smallest positive double is exp(-744.44), so exp of
    # anything below that is an exact zero.
    assert math.log(5e-324) > -745.0
    assert np.exp(np.asarray(ln_phi))[0] == 0.0

    mixture = _mixture(0.05, 53000.0)
    direct = PCSAFTEOS(parameters=parameters, kij=KIJ).log_fugacity_coefficients(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=1.0e7,
        composition=x,
        phase="liquid",
    )
    assert direct == ln_phi


def test_the_guard_is_dormant_on_an_ordinary_state(monkeypatch: pytest.MonkeyPatch) -> None:
    """Nothing that converged before ADR-0022 reaches the logarithmic route."""
    mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5])
    evaluations, log_space = _count_branch_evaluations(
        monkeypatch,
        lambda: ct.flash_tp(mixture, temperature_K=298.15, pressure_Pa=101325.0, eos=PCSAFTEOS()),
    )
    assert evaluations > 100
    assert log_space == 0


def test_the_guard_is_dormant_for_the_16400_polymer(monkeypatch: pytest.MonkeyPatch) -> None:
    """m = 431.32 is large but ``exp(ln phi)`` still exists; nothing is guarded."""
    mixture = _mixture()
    evaluations, log_space = _count_branch_evaluations(
        monkeypatch,
        lambda: ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos()),
    )
    assert evaluations > 100
    assert log_space == 0


def test_the_guard_carries_the_53000_polymer_through_a_whole_flash(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """m = 1393.9: every branch evaluation is in log space, and the answer is verified.

    Before ADR-0022 this state raised ``ModelError`` ("No usable
    fugacity-coefficient branch for stability analysis ... non-finite or
    non-positive fugacity coefficients") from ``stability_tp``, because
    ``exp(-1690.6)`` is ``0.0``.
    """
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    holder: list[ct.FlashResult] = []
    evaluations, log_space = _count_branch_evaluations(
        monkeypatch,
        lambda: holder.append(
            ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=eos)
        ),
    )
    assert evaluations > 100
    assert log_space == evaluations

    result = holder[0]
    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["fugacity_residual"]) < 1e-8
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"
    polymer_rich, solvent_rich = _phase_by_polymer(result)
    assert result.phases[polymer_rich].composition.fractions[0] == pytest.approx(
        3.729971e-04, rel=1e-6
    )
    assert result.phases[solvent_rich].composition.fractions[0] == pytest.approx(
        6.269597e-10, rel=1e-6
    )


# ---------------------------------------------------------------------------
# The equilibrium exam, Case P-13
# ---------------------------------------------------------------------------


def _verdict(pressure_Pa: float, *, kij: float = KIJ, temperature_K: float = TEMPERATURE_K) -> str:
    return ct.stability_tp(
        _mixture(),
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=_eos(kij),
    ).status


def _bisect(low: float, high: float, unstable_side: str, key, tolerance: float) -> float:
    """Bisect a verdict switch; ``low`` is the ``unstable_side`` end."""
    assert key(low) == unstable_side
    assert key(high) != unstable_side
    while abs(high - low) > tolerance:
        middle = 0.5 * (low + high)
        if key(middle) == unstable_side:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def test_the_feed_is_two_phase_below_a_cloud_point_pressure_and_one_phase_above() -> None:
    """Case P-13 (i). 5 wt% polymer, 453 K, k_ij = -0.006.

    The verdict switches once between 8 and 10 MPa and never switches back up
    to 30 MPa; ``test_the_cloud_point_pressure_is_bisected_from_the_verdict``
    refines the switch, and is the same check at a finer step.
    """
    verdicts = {P: _verdict(P) for P in (5.0e6, 8.0e6, 1.0e7, 3.0e7)}
    assert verdicts[5.0e6] == "unstable"
    assert verdicts[8.0e6] == "unstable"
    assert verdicts[1.0e7] == "stable"
    assert verdicts[3.0e7] == "stable"


@pytest.mark.slow  # a finer step on the switch the default run already brackets
def test_the_cloud_point_pressure_is_bisected_from_the_verdict() -> None:
    cloud_point = _bisect(5.0e6, 3.0e7, "unstable", _verdict, 1.0)
    assert cloud_point == pytest.approx(9.7518e6, rel=1e-4)


def test_the_split_inside_the_two_phase_region_is_two_verified_liquids() -> None:
    """Case P-13 (ii) at 8 MPa: the headline answer, with its three residuals."""
    mixture = _mixture()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert diagnostics["phase_regime"] == "LLE"
    assert diagnostics["phase_label_method"] == "compressibility"
    assert (diagnostics["phase_i_branch"], diagnostics["phase_ii_branch"]) == ("liquid", "liquid")

    # The polymer's K-value is ~1e-2 and its mole fraction ~1e-5, and the mass
    # balance still closes to the last bit: the asymmetry costs nothing here.
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"
    assert all(0.0 < value < 1.0 for value in result.phase_fractions.values())
    assert sum(result.phase_fractions.values()) == pytest.approx(1.0, abs=1e-12)

    polymer_rich, solvent_rich = _phase_by_polymer(result)
    assert result.phases[polymer_rich].composition.fractions[0] == pytest.approx(
        9.751377e-04, rel=1e-6
    )
    assert result.phases[solvent_rich].composition.fractions[0] == pytest.approx(
        1.1478711e-05, rel=1e-6
    )
    assert result.phase_fractions[polymer_rich] == pytest.approx(0.2282983, rel=1e-6)

    # Both phases are liquids by a kappa computed here from the public pressure
    # routine, not from `phase_identity`'s own derivative (measured: 0.0680 and
    # 0.1225 against the 0.5 threshold).
    bound = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=_parameters(), kij=KIJ)
    eos = _eos()
    for phase in result.phases.values():
        x = list(phase.composition.fractions)
        density = eos.density_roots(
            mixture=mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, composition=x
        )[-1]
        step = 1e-4 * density
        slope = (
            bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density + step, composition=x
            )
            - bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density - step, composition=x
            )
        ) / (2.0 * step)
        assert 0.0 < 8.0e6 / (density * slope) < KAPPA_LIQUID_THRESHOLD


@pytest.mark.slow  # another feed on the tie line the default run already solves
def test_the_tie_line_is_a_property_of_the_state_and_not_of_the_feed() -> None:
    """Two feeds inside the same two-phase region give the same two phases."""
    reference = ct.flash_tp(
        _mixture(0.05), temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos()
    )
    other = ct.flash_tp(_mixture(0.10), temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())
    for result in (reference, other):
        assert sorted(result.phases) == ["liquid1", "liquid2"]
    for source, target in zip(_phase_by_polymer(reference), _phase_by_polymer(other)):
        assert reference.phases[source].composition.fractions == pytest.approx(
            other.phases[target].composition.fractions, rel=1e-9
        )


def test_the_split_matches_an_independently_written_newton_solve() -> None:
    """Case P-13 (ii): the same tie line from a solver written here, not reused.

    Two unknowns - the polymer's mole fraction in each phase, carried as
    logarithms so the 1e-05 one keeps its digits - and the two equal-fugacity
    equations ``ln x_i + ln phi_i(x) = ln y_i + ln phi_i(y)``, solved by a
    finite-difference Newton started a thousandth away from ``flash_tp``'s
    answer so that it has to do real work. Nothing of the flash's Rachford-Rice,
    successive substitution, second-order stage or phase-count logic is used;
    only ``PCSAFTEOS``'s own ``(T, rho, x)`` interface and ``density_roots``.
    """
    pressure_Pa = 8.0e6
    result = ct.flash_tp(
        _mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=_eos()
    )
    bound = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=_parameters(), kij=KIJ)

    def ln_activity(polymer_fraction: float) -> np.ndarray:
        x = [polymer_fraction, 1.0 - polymer_fraction]
        density = bound.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
        )[-1]
        return np.log(np.asarray(x)) + np.asarray(
            bound.ln_fugacity_coefficients(
                temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
            )
        )

    def residual(u: np.ndarray) -> np.ndarray:
        return ln_activity(math.exp(u[0])) - ln_activity(math.exp(u[1]))

    names = list(result.phases)
    u = np.array(
        [math.log(result.phases[name].composition.fractions[0]) for name in names]
    ) + np.array([1e-3, -1e-3])
    for _ in range(60):
        value = residual(u)
        jacobian = np.zeros((2, 2))
        for column in range(2):
            step = 1e-7
            plus = u.copy()
            minus = u.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2.0 * step)
        correction = np.linalg.solve(jacobian, -value)
        u = u + correction
        if float(np.max(np.abs(correction))) < 1e-14:
            break

    assert float(np.max(np.abs(residual(u)))) < 1e-10
    ours = sorted(result.phases[name].composition.fractions[0] for name in names)
    theirs = sorted(math.exp(value) for value in u)
    # Measured: 1.7e-15 absolute on both branches.
    assert max(abs(a - b) for a, b in zip(ours, theirs)) < 1e-12


def test_the_model_switches_to_two_phases_as_the_temperature_rises() -> None:
    """Case P-13 (iii): the LCST-type direction, at a fixed 10 MPa.

    Qualitative only, and deliberately so. The cited manuscript's Figure 1
    shows the region *above* the cloud-point curve to be single phase, i.e.
    demixing appears on heating at fixed pressure. That is the direction
    checked here; the figure was not digitized and no pressure or temperature
    is asserted against it.
    """
    verdicts = {T: _verdict(1.0e7, temperature_K=T) for T in (400.0, 450.0, 460.0)}
    assert verdicts[400.0] == "stable"
    assert verdicts[450.0] == "stable"
    assert verdicts[460.0] == "unstable"


@pytest.mark.slow  # a finer step on the switch the default run already brackets
def test_the_cloud_point_temperature_is_bisected_from_the_verdict() -> None:
    switch = _bisect(460.0, 400.0, "unstable", lambda T: _verdict(1.0e7, temperature_K=T), 1e-3)
    assert switch == pytest.approx(454.56, abs=0.01)


def test_the_cloud_point_pressure_rises_with_kij() -> None:
    """Case P-13 (iv), qualitative: the manuscript's own statement, reproduced.

    At 15 MPa the feed is already single phase for ``k_ij = -0.006`` and for
    ``k_ij = 0``, and still two phases for ``k_ij = +0.02``: raising ``k_ij``
    pushes the cloud point up in pressure.
    ``test_the_cloud_point_pressures_are_bisected_for_three_kij`` puts numbers
    on the same three curves.
    """
    assert _verdict(1.5e7, kij=-0.006) == "stable"
    assert _verdict(1.5e7, kij=0.0) == "stable"
    assert _verdict(1.5e7, kij=0.02) == "unstable"


@pytest.mark.slow  # numbers for the ordering the default run already establishes
def test_the_cloud_point_pressures_are_bisected_for_three_kij() -> None:
    points = []
    for kij in (-0.006, 0.0, 0.02):
        high = 3.0e7
        while _verdict(high, kij=kij) == "unstable":
            high *= 1.5
        points.append(_bisect(5.0e6, high, "unstable", lambda P, k=kij: _verdict(P, kij=k), 1e3))

    assert points[0] < points[1] < points[2]
    assert points[0] == pytest.approx(9.75e6, rel=1e-3)
    assert points[1] == pytest.approx(1.077e7, rel=1e-3)
    assert points[2] == pytest.approx(2.151e7, rel=1e-3)


# ---------------------------------------------------------------------------
# Asymmetry stress (validation Case P-13 (v))
# ---------------------------------------------------------------------------


def test_the_answer_does_not_depend_on_the_component_order() -> None:
    forward = ct.flash_tp(_mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())
    backward = ct.flash_tp(
        _mixture(reversed_order=True),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=8.0e6,
        eos=_eos(),
    )
    for source, target in zip(_phase_by_polymer(forward), _phase_by_polymer(backward, 1)):
        ours = np.asarray(forward.phases[source].composition.fractions)
        theirs = np.asarray(backward.phases[target].composition.fractions)[::-1]
        assert float(np.max(np.abs(ours - theirs))) < 1e-13
        assert abs(forward.phase_fractions[source] - backward.phase_fractions[target]) < 1e-11


def test_the_answer_is_deterministic() -> None:
    results = [
        ct.flash_tp(_mixture(), temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())
        for _ in range(2)
    ]
    assert results[0].phase_fractions == results[1].phase_fractions
    for name in results[0].phases:
        assert (
            results[0].phases[name].composition.fractions
            == results[1].phases[name].composition.fractions
        )


@pytest.mark.parametrize("weight_fraction", [0.001, 0.40])
def test_the_dilute_and_the_concentrated_feed_both_converge(weight_fraction: float) -> None:
    """Both are outside the two-phase region at 8 MPa and return one phase.

    0.1 wt% is ``x_polymer = 4.4e-06`` and 40 wt% is ``x_polymer = 2.9e-03``;
    the binodal at this state runs from ``1.1e-05`` to ``9.8e-04``, so both
    feeds are single phase - the dilute one below the solvent-rich branch, the
    concentrated one above the polymer-rich branch. Neither raises, which is
    what this pins.
    """
    mixture = _mixture(weight_fraction)
    assert (
        ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos()).status
        == "stable"
    )
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=_eos())
    assert list(result.phases) == ["liquid"]
    assert result.phases["liquid"].composition.fractions == pytest.approx(
        tuple(mixture.fractions), rel=1e-12
    )


def test_a_polymer_in_two_solvents_runs() -> None:
    """Three components, a mass ratio of 1e5 : 1, and a verified split."""
    parameters = PCSAFTParameters.from_records(
        [
            _polymer_record(PE_MW_G_MOL),
            PCSAFTRecord(name="n-Pentane", m=2.6896, sigma_A=3.7729, epsilon_k_K=231.20),
            PCSAFTRecord(name="n-Hexane", m=3.0576, sigma_A=3.7983, epsilon_k_K=236.77),
        ]
    )
    amounts = [0.05 / PE_MW_G_MOL, 0.475 / PENTANE_MW_G_MOL, 0.475 / 86.177]
    total = sum(amounts)
    mixture = ct.Mixture.from_components(
        [
            _polymer_component(),
            ct.Component.from_database("n-Pentane"),
            ct.Component.from_database("n-Hexane"),
        ],
        [value / total for value in amounts],
        normalize=True,
    )
    eos = PCSAFTEOS(parameters=parameters, kij=KIJ)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=5.0e6, eos=eos)

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["fugacity_residual"]) < 1e-6
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"


# ---------------------------------------------------------------------------
# Limitations, pinned rather than worked around
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "pressure_Pa",
    [1.0e6, pytest.param(2.0e6, marks=pytest.mark.slow)],  # the second is the same gap
)
def test_the_split_does_not_converge_below_the_solvents_saturation_pressure(
    pressure_Pa: float,
) -> None:
    """A vapour root exists there, and the split fails - recorded, not patched.

    n-pentane is subcritical at 453 K, so below roughly 3 MPa the mixture has a
    vapour density root as well as a liquid one and the equilibrium in question
    is vapour-liquid, not liquid-liquid. ``stability_tp`` still reports the feed
    unstable; the split then stops after one successive-substitution step with
    an equal-fugacity residual of order 1e+02 and ``flash_tp`` raises. This is a
    genuine gap in the phi-phi split for this system, not a property of the
    polymer support added by ADR-0022, and it is pinned here so that any future
    change to it is deliberate. Polymer/solvent *vapour*-liquid equilibrium is
    named as a later slice candidate, not delivered here.
    """
    mixture = _mixture()
    eos = _eos()
    assert (
        ct.stability_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
        ).status
        == "unstable"
    )
    assert (
        len(
            eos.density_roots(
                mixture=mixture,
                temperature_K=TEMPERATURE_K,
                pressure_Pa=pressure_Pa,
                composition=list(mixture.fractions),
            )
        )
        == 2
    )
    with pytest.raises(ConvergenceError):
        ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)


def test_the_ternary_split_runs_away_at_3_mpa() -> None:
    """The same gap, in the ternary: pinned with the ``beta`` the split reached."""
    parameters = PCSAFTParameters.from_records(
        [
            _polymer_record(PE_MW_G_MOL),
            PCSAFTRecord(name="n-Pentane", m=2.6896, sigma_A=3.7729, epsilon_k_K=231.20),
            PCSAFTRecord(name="n-Hexane", m=3.0576, sigma_A=3.7983, epsilon_k_K=236.77),
        ]
    )
    amounts = [0.05 / PE_MW_G_MOL, 0.475 / PENTANE_MW_G_MOL, 0.475 / 86.177]
    total = sum(amounts)
    mixture = ct.Mixture.from_components(
        [
            _polymer_component(),
            ct.Component.from_database("n-Pentane"),
            ct.Component.from_database("n-Hexane"),
        ],
        [value / total for value in amounts],
        normalize=True,
    )
    eos = PCSAFTEOS(parameters=parameters, kij=KIJ)
    with pytest.raises(ConvergenceError, match="outside"):
        ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=eos)


def test_a_record_without_a_molar_mass_cannot_use_the_mass_based_parameter() -> None:
    from chemthermo.parameters import PCSAFTParameterError

    with pytest.raises(PCSAFTParameterError, match="MW_g_mol"):
        PCSAFTRecord(name="Polyethylene", segments_per_g=0.0263, sigma_A=4.0217, epsilon_k_K=247.5)
