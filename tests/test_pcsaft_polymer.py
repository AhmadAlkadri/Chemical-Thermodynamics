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
- and the **vapour**-liquid split below the solvent's saturation pressure
  (validation Case P-14, ADR-0024), which was pinned here as a defect until
  the split learned to run in log mole numbers. The stage that does that is
  tested in ``tests/test_flash_log_space_stage.py``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from _capture_identity import on_capture_platform

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.flash import _detect
from chemthermo.flash._split import _rachford_rice
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
# Vapour-liquid below the solvent's saturation pressure (Case P-14)
#
# Both states below were pinned as *defects* by ADR-0022 (Case P-13 (vi)) and
# are resolved by the log-space split stage of ADR-0024. The stage itself is
# tested in ``tests/test_flash_log_space_stage.py``; what is pinned here is the
# answer for this system.
# ---------------------------------------------------------------------------


#: ``(melt composition, melt phase fraction, ln y_polymer)`` per pressure, at
#: 5 wt% polymer, 453 K and ``k_ij = -0.006``. The melt's solvent content is
#: reproduced to 1e-10 by an independent one-dimensional equal-fugacity solve
#: in ``examples/validation/21_pcsaft_polymer_vle.py``, which shares no code
#: with the flash.
VLE_STATES = {
    1.0e6: ((0.028531647335324484, 0.9714683526646756), 0.008113111018756114, -450.5307940616726),
    2.0e6: ((0.008772681701263427, 0.9912273182987366), 0.026386506459723290, -408.9441806099687),
}


@pytest.mark.parametrize(
    "pressure_Pa",
    [1.0e6, pytest.param(2.0e6, marks=pytest.mark.slow)],  # the second is the same state
)
def test_the_split_below_the_solvents_saturation_pressure_is_a_vapour_and_a_melt(
    pressure_Pa: float,
) -> None:
    """Case P-14 (i): a solvent vapour over a solvent-swollen melt.

    n-pentane is subcritical at 453 K, so below roughly 2.6 MPa the mixture has
    a vapour density root as well as a liquid one and the equilibrium in
    question is vapour-liquid, not liquid-liquid. Before ADR-0024 the split
    stopped after one successive-substitution step with an equal-fugacity
    residual of order 1e+02 and ``flash_tp`` raised: the tangent-plane
    minimizer here is an essentially pure polymer melt whose K-values bracket
    no vapour fraction at all, and the equilibrium vapour's polymer content is
    ``exp(-450)``. Both are now carried in log mole numbers.
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

    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    diagnostics = result.diagnostics
    melt, melt_fraction, ln_y_polymer = VLE_STATES[pressure_Pa]

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert diagnostics["phase_regime"] == "VLE"
    assert diagnostics["phase_label_method"] == "compressibility"
    assert diagnostics["k_seed"] == "stability-log"
    assert diagnostics["converged_stage"] == "second-order-log"

    assert result.phases["liquid"].composition.fractions == pytest.approx(melt, rel=1e-12)
    assert result.phase_fractions["liquid"] == pytest.approx(melt_fraction, rel=1e-12)
    assert result.phases["vapor"].composition.fractions[1] == 1.0
    assert math.log(result.phases["vapor"].composition.fractions[0]) == pytest.approx(
        ln_y_polymer, rel=1e-12
    )
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(ln_y_polymer, rel=1e-12)
    assert result.vapor_fraction == pytest.approx(1.0 - melt_fraction, rel=1e-12)

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["log_space_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"

    # The compressibility identity of ADR-0017, computed here from the public
    # pressure routine: the vapour is one side of the 0.5 threshold and the
    # melt is far on the other (measured 1.16 and 0.0024 at 1 MPa).
    bound = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=_parameters(), kij=KIJ)
    kappa = {}
    for name, index in (("vapor", 0), ("liquid", -1)):
        x = list(result.phases[name].composition.fractions)
        density = eos.density_roots(
            mixture=mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
        )[index]
        step = 1e-4 * density
        slope = (
            bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density + step, composition=x
            )
            - bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density - step, composition=x
            )
        ) / (2.0 * step)
        kappa[name] = pressure_Pa / (density * slope)
    assert kappa["vapor"] > KAPPA_LIQUID_THRESHOLD
    assert 0.0 < kappa["liquid"] < KAPPA_LIQUID_THRESHOLD


def test_the_ternary_converges_at_3_mpa() -> None:
    """Case P-14 (ii): the state that used to run away to ``beta = -6.3e+10``.

    Successive substitution converges here - on the **trivial** solution, whose
    vapour fraction is outside ``[0, 1]``; before ADR-0024 that was refused and
    ``flash_tp`` raised. It is now handed to the second-order stage, which
    finds the real split. The stage is the *linear* one: no composition here is
    outside machine range, and this state is in the ledger as the one that
    shows the two ADR-0024 entry points are separate.
    """
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
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=eos)
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert diagnostics["converged_stage"] == "second-order"
    assert not [key for key in diagnostics if key.startswith("log_space_")]
    polymer_rich, solvent_rich = _phase_by_polymer(result)
    assert result.phases[polymer_rich].composition.fractions == pytest.approx(
        (0.0013272336238799715, 0.5402484470355627, 0.45842431934055733), rel=1e-9
    )
    assert result.phases[solvent_rich].composition.fractions == pytest.approx(
        (1.863504765793684e-06, 0.5450873602321722, 0.45491077626306187), rel=1e-9
    )
    assert result.phase_fractions[polymer_rich] == pytest.approx(0.18872173701842598, rel=1e-9)
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"


# ---------------------------------------------------------------------------
# The last six refusals of the Mw = 53 000 chain (validation Case P-16,
# ADR-0026)
#
# A 0.3-3.6 MPa sweep at 0.1 MPa steps over both molar masses - 68 states - had
# exactly six states raising `ConvergenceError` at HEAD 584c508, all of them on
# the Mw = 53 000 chain, and for two different reasons:
#
#   0.3 MPa        the log-space stage parked next to the trivial solution and
#   2.8, 2.9 MPa   spent its budget (residual 3.9e-05 / 1.1e-08 / 6.0e+00);
#                  repaired by the ADR-0026 curvature safeguard.
#   3.0-3.2 MPa    "neither the stability-seeded nor the Wilson K-values
#                  bracket a Rachford-Rice root", because `K_polymer` is under
#                  the spacing of doubles at one and `1 + (K - 1)` cancels to
#                  an exact zero at `beta = 1`; repaired by the ADR-0026
#                  convex denominator.
#
# Every number below is reproduced by a solve that shares no code with
# `flash_tp`: the 0.3 MPa melt by the one-dimensional equal-fugacity solve of
# Case P-14, the five liquid-liquid tie lines by a two-equation Newton, both in
# `examples/validation/22_stability_log_space.py`.
# ---------------------------------------------------------------------------


#: ``pressure -> (polymer-rich x_polymer, polymer-lean ln x_polymer, phase
#: fraction of the polymer-rich phase)`` for the five liquid-liquid states.
P16_LLE_STATES = {
    2.8e6: (9.231250339583601e-04, -93.81711157873, 0.07760525755403569),
    2.9e6: (9.084501098780754e-04, -90.20777508681, 0.07885887759376387),
    3.0e6: (8.940738921079231e-04, -86.82803360271, 0.08012688509001775),
    3.1e6: (8.799780999223999e-04, -83.65315387491, 0.08141038512348653),
    3.2e6: (8.661464766005232e-04, -80.66231922913, 0.08271043980469528),
}


def test_the_0_3_mpa_vapour_liquid_split_the_stage_used_to_stall_on() -> None:
    """Case P-16 (i): the shallowest state of the vapour-liquid region.

    Before ADR-0026 the log-space stage spent all 100 iterations here moving
    ``ln n`` of the solvent by 4e-03 while the answer was 3.2 away, and
    ``flash_tp`` raised with a residual of 3.9e-05. The melt's solvent content
    below is reproduced to 1.4e-14 by an independent one-dimensional
    equal-fugacity solve.
    """
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e5, eos=eos)
    diagnostics = result.diagnostics

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert diagnostics["phase_regime"] == "VLE"
    assert diagnostics["k_seed"] == "stability-log"
    assert diagnostics["converged_stage"] == "second-order-log"
    assert diagnostics["log_space_curvature_safeguard"] is True

    assert result.phases["liquid"].composition.fractions == pytest.approx(
        (0.036808345719672377, 0.9631916542803276), rel=1e-12
    )
    assert result.phases["vapor"].composition.fractions[0] == 0.0
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(-1528.39130415, rel=1e-9)
    assert result.vapor_fraction == pytest.approx(0.9980537197579997, rel=1e-12)

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["log_space_residual"]) < 1e-12
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"


@pytest.mark.parametrize(
    "pressure_Pa",
    [
        2.8e6,
        # `slow`: the same statement at four further pressures on the same tie
        # line, one of which (3.0 MPa) runs by default just below.
        pytest.param(2.9e6, marks=pytest.mark.slow),
        3.0e6,
        pytest.param(3.1e6, marks=pytest.mark.slow),
        pytest.param(3.2e6, marks=pytest.mark.slow),
    ],
)
def test_the_liquid_liquid_states_between_2_8_and_3_2_mpa(pressure_Pa: float) -> None:
    """Case P-16 (ii) and (iii): the two remaining refusals, both repaired.

    2.8 and 2.9 MPa reach the answer through the curvature safeguard; 3.0 to
    3.2 MPa reach it through the successive-substitution loop, because the
    stability seed's K-values bracket a Rachford-Rice root once the
    denominator stops cancelling. The two routes meet in the middle: the
    polymer-rich composition is smooth in pressure across the boundary between
    them, and so is the polymer-lean one, which runs from ``exp(-93.8)`` to
    ``exp(-80.7)`` over the five states.
    """
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    diagnostics = result.diagnostics
    rich_x, lean_ln_x, rich_fraction = P16_LLE_STATES[pressure_Pa]

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert diagnostics["phase_regime"] == "LLE"
    assert result.vapor_fraction is None

    polymer_rich, solvent_rich = _phase_by_polymer(result)
    assert result.phases[polymer_rich].composition.fractions[0] == pytest.approx(rich_x, rel=1e-11)
    assert result.phase_fractions[polymer_rich] == pytest.approx(rich_fraction, rel=1e-11)
    lean = result.phases[solvent_rich].composition.fractions[0]
    assert math.log(lean) == pytest.approx(lean_ln_x, rel=1e-9)

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"

    if pressure_Pa < 3.0e6:
        assert diagnostics["k_seed"] == "stability-log"
        assert diagnostics["converged_stage"] == "second-order-log"
        assert diagnostics["log_space_curvature_safeguard"] is True
    else:
        assert diagnostics["k_seed"] == "stability"
        assert diagnostics["converged_stage"] == "second-order"
        assert diagnostics["rachford_rice_convex_denominators"] is True


def test_the_underflowed_k_is_the_stability_seed_at_3_mpa() -> None:
    """Where ``tests/test_rachford_rice_extended.py``'s ``TRACE_K`` comes from.

    That module pins the equation with plain numbers so it needs no equation of
    state; this is the state those numbers are a snapshot of, recomputed from
    the model. The polymer's K is 1.1e-18, which is under the spacing of
    doubles at one - and that, not the physics, is what made the split refuse.
    """
    mixture = _mixture(0.05, 53000.0)
    eos = _eos(KIJ, 53000.0)
    z = np.asarray(mixture.fractions, dtype=float)
    stability = ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=eos)
    assert stability.status == "unstable"
    assert stability.trial_composition is not None

    k_seed, incipient = _detect._stability_k_seed(
        mixture,
        TEMPERATURE_K,
        3.0e6,
        z=z,
        w=np.asarray(stability.trial_composition, dtype=float),
        tpd_min=float(stability.tpd_min),
        ln_capital_w=_detect._stationary_point_ln_capital_w(stability),
    )
    assert incipient == "vapor"
    assert k_seed[0] == pytest.approx(1.11871119e-18, rel=1e-6)
    assert k_seed[1] == pytest.approx(1.00132622e00, rel=1e-8)
    assert k_seed[0] - 1.0 == -1.0  # the cancellation, in one line
    assert _rachford_rice(z, k_seed)[0] is None
    assert _rachford_rice(z, k_seed, convex_denominators=True)[0] is not None


def test_a_record_without_a_molar_mass_cannot_use_the_mass_based_parameter() -> None:
    from chemthermo.parameters import PCSAFTParameterError

    with pytest.raises(PCSAFTParameterError, match="MW_g_mol"):
        PCSAFTRecord(name="Polyethylene", segments_per_g=0.0263, sigma_A=4.0217, epsilon_k_K=247.5)


# ---------------------------------------------------------------------------
# What the robustness map left, and ADR-0028 retires (validation Case P-17)
#
# The 2110-state map at `87f0820` refused 36 states and every one of them was
# this system (ledger Case R-MAP-1). Three causes, all of them about where a
# stage *starts* rather than how it steps:
#
#   1 wt%, 0.3-1.2 MPa    the log-space stage converged on the **trivial**
#   (both molar masses)   solution - equal compositions, residual 8e-13, phase
#                         fraction collapsing to 0 - because its seed puts the
#                         two phases half and half while the melt here holds
#                         4e-04 of the feed. Repaired by the third ladder
#                         entry, whose phase fraction is the lever rule's.
#   3.6-8.7 MPa           the linear second-order stage and its log-space
#   (1, 5 and 15 wt%)     retry both continue from the K-loop's last iterate,
#                         and the K-loop had diverged (`max_delta_k` to
#                         1e+128). Repaired by running the log-space stage from
#                         the stationary point instead - ladder entries 1 and 2.
#   15 wt%, 10.5/10.8 MPa every stability trial stalled in the Newton stage.
#                         Repaired by the ADR-0028 substitution-budget retry;
#                         the verdict is `stable`, which is what the pressures
#                         on both sides of it say too.
#
# Every tie line below is reproduced by a solver that shares no code with
# `flash_tp` (the two-equation Newton above, the one-dimensional solve of Case
# P-14), and the three that are checked against FeOs are in
# `tests/validation/test_pcsaft_polymer_vs_feos.py`.
# ---------------------------------------------------------------------------


#: ``pressure -> (melt composition, melt phase fraction, ln y_polymer)`` at
#: 1 wt% polymer, ``Mw = 16400``. The melt fraction at 0.3 MPa, 4.0742e-04, is
#: the lever rule's on the tie line the 5 wt% feed converges on - which is how
#: the map's diagnosis knew what these states were refusing to find.
P17_DILUTE_VLE = {
    3.0e5: ((0.10906250886173996, 0.89093749113826), 0.00040741633069130145, -472.87659070916015),
    6.0e5: ((0.053053918285085024, 0.946946081714915), 0.0008375224415597682, -463.8731421615793),
    9.0e5: ((0.032720262948131754, 0.9672797370518682), 0.001357991751070986, -454.00264952079226),
    1.2e6: ((0.022154085500088714, 0.9778459144999112), 0.00200567282166797, -443.2920259578849),
}


@pytest.mark.parametrize(
    "pressure_Pa",
    [
        3.0e5,
        # `slow`: the same statement at three further pressures on the same
        # 1 wt% isopleth, all through the same ladder entry.
        pytest.param(6.0e5, marks=pytest.mark.slow),
        pytest.param(9.0e5, marks=pytest.mark.slow),
        pytest.param(1.2e6, marks=pytest.mark.slow),
    ],
)
def test_the_dilute_feed_converges_on_the_lever_rule_seed(pressure_Pa: float) -> None:
    """Case P-17 (i): the four ``beta``-outside-window refusals of the map.

    Before ADR-0028 the log-space stage reached a residual of 8e-13 here and
    `flash_tp` still raised, because what it had converged on was the trivial
    solution: both phases at the feed composition, where every equal-fugacity
    residual is zero by construction and the phase fraction is an exact ``0``.
    The seed's phase fraction, not the stage, is what put it there.
    """
    mixture = _mixture(0.01)
    eos = _eos()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    diagnostics = result.diagnostics
    melt, melt_fraction, ln_y_polymer = P17_DILUTE_VLE[pressure_Pa]

    assert sorted(result.phases) == ["liquid", "vapor"]
    assert diagnostics["phase_regime"] == "VLE"
    assert diagnostics["k_seed"] == "stability-log"
    assert diagnostics["converged_stage"] == "second-order-log"
    assert diagnostics["log_space_seed"] == "stability-w-lever-rule"
    assert diagnostics["log_space_curvature_safeguard"] is True

    assert result.phases["liquid"].composition.fractions == pytest.approx(melt, rel=1e-11)
    assert result.phase_fractions["liquid"] == pytest.approx(melt_fraction, rel=1e-11)
    assert math.log(result.phases["vapor"].composition.fractions[0]) == pytest.approx(
        ln_y_polymer, rel=1e-9
    )
    assert result.vapor_fraction == pytest.approx(1.0 - melt_fraction, rel=1e-11)

    # The lever rule on the converged tie line, which is the number the map's
    # diagnosis predicted before any of this converged.
    z = np.asarray(mixture.fractions, dtype=float)
    lean = result.phases["vapor"].composition.fractions[0]
    assert result.phase_fractions["liquid"] == pytest.approx(
        (z[0] - lean) / (melt[0] - lean), rel=1e-12
    )

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["log_space_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"


def test_the_53000_dilute_state_lands_on_the_tie_line_case_p_16_already_pinned() -> None:
    """Case P-17 (ii): the three ``log-space`` refusals, checked against a fixed point.

    The ``Mw = 53 000`` chain at 1 wt% and 0.3 MPa refused with a residual of
    8.0e-02 before ADR-0028. The answer it now returns is not checked against a
    number invented for it: 0.3 MPa at **5 wt%** is a state ADR-0026 already
    pinned, verified there against a one-dimensional equal-fugacity solve, and
    a tie line is a property of the state and not of the feed. The two melts
    agree to 1.1e-12 relative, from two different ladder entries.
    """
    eos = _eos(KIJ, 53000.0)
    dilute = ct.flash_tp(
        _mixture(0.01, 53000.0), temperature_K=TEMPERATURE_K, pressure_Pa=3.0e5, eos=eos
    )
    diagnostics = dilute.diagnostics

    assert sorted(dilute.phases) == ["liquid", "vapor"]
    assert diagnostics["log_space_seed"] == "stability-w-lever-rule"
    assert dilute.phases["liquid"].composition.fractions == pytest.approx(
        (0.036808345719631735, 0.9631916542803682), rel=1e-11
    )
    # The 5 wt% melt of Case P-16, from the block above, on the same tie line.
    assert dilute.phases["liquid"].composition.fractions[0] == pytest.approx(
        0.036808345719672377, rel=1e-11
    )
    assert dilute.phase_fractions["liquid"] == pytest.approx(0.00037355015625717414, rel=1e-11)
    assert dilute.phases["vapor"].composition.fractions[0] == 0.0
    assert float(diagnostics["log_space_ln_x_min"]) == pytest.approx(-1528.39130415, rel=1e-9)

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["log_space_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"


#: ``(Mw, weight fraction, pressure) -> (polymer-rich composition, polymer-lean
#: ln x_polymer, polymer-rich phase fraction, curvature safeguard)`` for the
#: band the K-loop used to hand over in ruins.
P17_LADDER_LLE = {
    (53000.0, 0.15, 8.1e6): (
        (0.0003647955115724242, 0.9996352044884275),
        -20.715500022759343,
        0.6583461479606705,
        False,
    ),
    (16400.0, 0.05, 7.5e6): (
        (0.001129717115748221, 0.9988702828842517),
        -12.185046049951254,
        0.20129111537244893,
        True,
    ),
    (53000.0, 0.01, 3.6e6): (
        (0.0008131909284340129, 0.999186809071566),
        -70.220361234825,
        0.016908407133382708,
        True,
    ),
    (53000.0, 0.01, 8.1e6): (
        (0.0003647955115718861, 0.9996352044884281),
        -20.715500022776144,
        0.03768904234241921,
        True,
    ),
}

#: Ladder states whose *route* (not answer) is decided by last-bit noise.
#: Over the 17 pressures within +-8 ULP of 8.1 MPa, one Linux x86_64 host took
#: the `stability-w` rung 12 times, the `linear-iterate` rung 3 times, the
#: linear stage once, and refused once; every converged answer was this tie
#: line to <= 4.7e-12. The 16400 g/mol state took `stability-w` 17 times of 17
#: and stays pinned everywhere. ADR-0032; ledger Case P-17 "cross-platform".
P17_ROUTE_BY_NOISE = {(53000.0, 0.15, 8.1e6)}


@pytest.mark.parametrize(
    ("mw_g_mol", "weight_fraction", "pressure_Pa"),
    [
        # One state per ladder entry runs by default: the 15 wt% state needs
        # only the stationary-point seed, the 16400 state needs the curvature
        # safeguard with it.
        (53000.0, 0.15, 8.1e6),
        (16400.0, 0.05, 7.5e6),
        # `slow`: two further pressures of the same 26-state band on the same
        # chain, both through the entry the 16400 state above already covers.
        pytest.param(53000.0, 0.01, 3.6e6, marks=pytest.mark.slow),
        pytest.param(53000.0, 0.01, 8.1e6, marks=pytest.mark.slow),
    ],
)
def test_the_band_the_diverged_k_loop_used_to_end(
    mw_g_mol: float, weight_fraction: float, pressure_Pa: float
) -> None:
    """Case P-17 (iii): the 27 ``phi-phi`` refusals, the map's largest class.

    Successive substitution runs away here - ``max_delta_k`` reaches 1e+89 to
    1e+128 - and both stages that follow it were started from what it left.
    They are now started from the tangent-plane stationary point instead, which
    is a phase the stability test actually measured.
    """
    mixture = _mixture(weight_fraction, mw_g_mol)
    eos = _eos(KIJ, mw_g_mol)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    diagnostics = result.diagnostics
    rich_x, lean_ln_x, rich_fraction, safeguard = P17_LADDER_LLE[
        (mw_g_mol, weight_fraction, pressure_Pa)
    ]

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert diagnostics["phase_regime"] == "LLE"
    assert result.vapor_fraction is None
    # The seed is the stationary point's, and the successive-substitution loop
    # is what failed - so `k_seed` is still `"stability"`, not `"stability-log"`.
    assert diagnostics["k_seed"] == "stability"
    assert diagnostics["converged_stage"] == "second-order-log"
    if on_capture_platform() or (mw_g_mol, weight_fraction, pressure_Pa) not in P17_ROUTE_BY_NOISE:
        assert diagnostics["log_space_seed"] == "stability-w"
    else:
        # Which rung lands first is decided by the wreckage of a diverged
        # K-loop, so it is last-bit noise here (ADR-0032, ledger Case P-17
        # "cross-platform"); the tie line asserted below is not.
        assert diagnostics["log_space_seed"] in ("stability-w", "linear-iterate")
    assert diagnostics["log_space_curvature_safeguard"] is safeguard

    polymer_rich, solvent_rich = _phase_by_polymer(result)
    assert result.phases[polymer_rich].composition.fractions == pytest.approx(rich_x, rel=1e-11)
    assert result.phase_fractions[polymer_rich] == pytest.approx(rich_fraction, rel=1e-11)
    assert math.log(result.phases[solvent_rich].composition.fractions[0]) == pytest.approx(
        lean_ln_x, rel=1e-9
    )

    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-8
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"


@pytest.mark.slow  # another feed on a tie line the default run already solves
def test_the_recovered_tie_line_is_the_same_from_a_1_and_a_15_wt_percent_feed() -> None:
    """Case P-17 (iv): two feeds, two ladder entries, one tie line.

    8.1 MPa is the pressure the quick subset of the robustness map pins, and
    the map refuses it at both 1 wt% and 15 wt%. They converge through
    different ladder entries - the 15 wt% feed needs no curvature safeguard and
    the 1 wt% feed does - so agreeing to 1.5e-12 is a statement about the tie
    line rather than about a shared code path.
    """
    eos = _eos(KIJ, 53000.0)
    lines = []
    for weight_fraction in (0.01, 0.15):
        result = ct.flash_tp(
            _mixture(weight_fraction, 53000.0),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=8.1e6,
            eos=eos,
        )
        polymer_rich, solvent_rich = _phase_by_polymer(result)
        lines.append(
            (
                result.phases[polymer_rich].composition.fractions[0],
                result.phases[solvent_rich].composition.fractions[0],
                result.diagnostics["log_space_curvature_safeguard"],
            )
        )
    assert lines[0][2] is not lines[1][2]
    assert lines[0][0] == pytest.approx(lines[1][0], rel=1e-11)
    assert lines[0][1] == pytest.approx(lines[1][1], rel=1e-8)


@pytest.mark.parametrize(
    "pressure_Pa",
    [1.05e7, pytest.param(1.08e7, marks=pytest.mark.slow)],  # the second is the same statement
)
def test_the_stalled_stability_pair_is_resolved_as_stable(pressure_Pa: float) -> None:
    """Case P-17 (v): the two ``stability-inconclusive`` refusals of the map.

    All four trials ended ``second_order_no_progress`` at a residual of 0.685
    and 1.746 - the Newton line search finding no admissible step at all - from
    a 50-substitution iterate that was still travelling. With the substitutions
    given their full ``max_iter`` budget the same unchanged Newton stage
    converges, on the trivial solution, and the verdict is ``stable``.

    That verdict is not taken on the solver's word here: 9.9 and 10.2 MPa
    below it are ``stable`` with a positive tangent-plane minimum, 11.1 and
    11.4 MPa above it are ``stable`` on the trivial solution, and the scan in
    ``test_no_trial_free_scan_finds_a_negative_tangent_plane_distance`` finds
    no negative distance anywhere.

    ``flash_tp``'s end of it - a single liquid returned rather than a refusal -
    is the test just below; it runs the analysis a second time, so it is
    ``slow`` and this one is what runs by default.
    """
    stability = ct.stability_tp(
        _mixture(0.15, 53000.0),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        eos=_eos(KIJ, 53000.0),
    )
    assert stability.status == "stable"
    assert stability.diagnostics["substitution_budget_retry"] is True
    assert stability.diagnostics["substitution_budget_retry_from"] == 50
    assert stability.diagnostics["ssi_iterations_budget"] == 300
    converged = [trial for trial in stability.trials if trial.converged]
    assert converged, "the retry is only kept when it produced a verdict"
    assert all(trial.trivial for trial in converged)


@pytest.mark.slow  # the same two states through `flash_tp`, which re-runs the analysis
@pytest.mark.parametrize("pressure_Pa", [1.05e7, 1.08e7])
def test_the_resolved_pair_returns_a_single_liquid_rather_than_a_refusal(
    pressure_Pa: float,
) -> None:
    """Case P-17 (v), end to end: what the user gets where the map recorded a refusal."""
    result = ct.flash_tp(
        _mixture(0.15, 53000.0),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        eos=_eos(KIJ, 53000.0),
    )
    assert list(result.phases) == ["liquid"]
    assert result.diagnostics["termination_reason"] == "feed_stable_tangent_plane"
    assert result.diagnostics["phase_regime"] == "single-phase"
    assert result.diagnostics["stability_status"] == "stable"


@pytest.mark.parametrize("pressure_Pa", [9.9e6, 1.11e7])
def test_the_substitution_budget_retry_is_dormant_on_the_pressures_either_side(
    pressure_Pa: float,
) -> None:
    """Both neighbours reach a verdict on the first pass, so the retry never runs.

    9.9 MPa converges to a non-trivial stationary point with ``tpd > 0`` and
    11.1 MPa to the trivial one; the key's *absence* is the assertion that
    nothing about either was re-run.
    """
    stability = ct.stability_tp(
        _mixture(0.15, 53000.0),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        eos=_eos(KIJ, 53000.0),
    )
    assert stability.status == "stable"
    assert "substitution_budget_retry" not in stability.diagnostics


@pytest.mark.slow  # a 220-point scan of the same verdict the default run asserts
def test_no_trial_free_scan_finds_a_negative_tangent_plane_distance() -> None:
    """Case P-17 (v), independently: ``tpd(w) >= 0`` on a grid, with no solver at all.

    Michelsen's test is a local stationary-point search, so "stable" is only
    ever "no negative distance was found from this trial set". This walks
    equation (2) directly over a grid of trial compositions spanning twelve
    decades of polymer content, using nothing but ``density_roots`` and
    ``ln_fugacity_coefficients``. Measured minimum: +8.0e-09 at 10.5 MPa and
    +8.8e-09 at 10.8 MPa, both at the feed composition itself.
    """
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=_parameters(53000.0), kij=KIJ
    )

    def min_gibbs_ln_phi(pressure_Pa: float, x: list[float]) -> np.ndarray:
        best: tuple[float, np.ndarray] | None = None
        for density in bound.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
        ):
            terms = np.asarray(
                bound.ln_fugacity_coefficients(
                    temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
                )
            )
            energy = float(np.dot(x, terms))
            if best is None or energy < best[0]:
                best = (energy, terms)
        assert best is not None
        return best[1]

    z = np.asarray(_mixture(0.15, 53000.0).fractions, dtype=float)
    grid = np.concatenate(
        [np.logspace(-12.0, math.log10(200.0 * z[0]), 200), np.linspace(0.5, 1.0 - 1e-9, 20)]
    )
    for pressure_Pa in (1.05e7, 1.08e7):
        plane = np.log(z) + min_gibbs_ln_phi(pressure_Pa, z.tolist())
        worst = math.inf
        for polymer in grid:
            w = np.array([polymer, 1.0 - polymer])
            terms = min_gibbs_ln_phi(pressure_Pa, w.tolist())
            worst = min(worst, float(np.sum(w * (np.log(w) + terms - plane))))
        assert worst >= 0.0, (pressure_Pa, worst)
