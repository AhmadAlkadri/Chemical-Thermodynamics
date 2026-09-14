"""The one-solve branch capability returns the same doubles as two model calls.

ADR-0023 adds ``EquationOfState.ln_fugacity_branches``: one density-root solve
answering for every branch, where the minimum-Gibbs rule of ADR-0005 and the
lowest-Gibbs fallback of ADR-0019 used to ask the model twice at the same
``(T, P, x)``. That is a **performance** change and its whole gate is
bit-identity, so this module compares the two routes with ``==``, never with a
tolerance:

1. ``exp(ln_fugacity_branches(...)[phase])`` must equal
   ``fugacity_coefficients(..., phase=phase)`` exactly, for both equations of
   state, over the states the library's own validation grids use;
2. :func:`chemthermo.flash._common.eos_branch_terms_all` - the caller-side
   wrapper that re-applies the ADR-0022 guard - must return exactly what one
   :func:`~chemthermo.flash._common.eos_branch_terms` call per branch returns,
   including ``phi is None`` where ``exp(ln phi)`` has underflowed;
3. a model that does not implement the capability must go on being served by
   the per-branch route, unchanged.

The pinned fixture ``tests/test_flash_refactor_bit_identity.py`` is the other
half of the evidence: it proves the *flash* answers did not move. This file
proves the model calls underneath them did not either.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.flash._common import eos_branch_terms, eos_branch_terms_all
from chemthermo.models import PengRobinsonEOS
from chemthermo.models.base import EquationOfState

BRANCHES = ("vapor", "liquid")

#: Verbatim from ``tests/validation/test_flash_split_robustness_pcsaft_subset.py``:
#: the default-run subset of the Case F-4 grid, including the four states that
#: need the ADR-0016 second-order stage.
F4_SUBSET: tuple[tuple[tuple[str, str], float, float, float], ...] = (
    (("Carbon dioxide", "n-Decane"), 0.8, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 260.0, 1.5e6),
    (("Carbon dioxide", "n-Decane"), 0.6, 230.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.6, 260.0, 2.5e6),
    (("Carbon dioxide", "n-Decane"), 0.8, 230.0, 1.5e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 230.0, 2.0e6),
    (("Methane", "n-Hexane"), 0.5, 170.0, 0.5e6),
    (("Methane", "n-Hexane"), 0.5, 200.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.8, 195.0, 2.5e6),
    (("Methane", "n-Hexane"), 0.8, 200.0, 3.0e6),
    (("Methane", "n-Hexane"), 0.9, 190.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.95, 170.0, 1.0e6),
    (("Methane", "n-Hexane"), 0.95, 195.0, 2.5e6),
    (("Methane", "n-Hexane"), 0.95, 200.0, 3.5e6),
)

#: The ADR-0019 water / n-hexane liquid-liquid states, at the feed and on both
#: sides of its tie line, plus the pressure where the vapour root disappears.
WATER_HEXANE_STATES: tuple[tuple[float, float, tuple[float, float]], ...] = (
    (298.15, 101325.0, (0.5, 0.5)),
    (298.15, 101325.0, (0.2, 0.8)),
    (298.15, 101325.0, (0.8, 0.2)),
    (298.15, 101325.0, (0.999983, 1.7e-05)),
    (298.15, 1.0e6, (0.5, 0.5)),
    (335.0, 101325.0, (0.7, 0.3)),
)

#: A slice of the ADR-0017 Peng-Robinson grid: both binaries and both ternaries,
#: at the pressure and temperatures where the cubic has one root and three.
PR_STATES: tuple[tuple[tuple[str, ...], tuple[float, ...], float, float], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5), 170.0, 2.0e5),
    (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6),
    (("Methane", "Ethane"), (0.5, 0.5), 360.0, 8.0e6),
    (("Ethane", "n-Heptane"), (0.7, 0.3), 280.0, 1.0e6),
    (("Methane", "n-Pentane"), (0.6, 0.4), 200.0, 8.0e6),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3), 320.0, 2.0e5),
)

POLYMER_FIXTURE = (
    Path(__file__).resolve().parent / "fixtures" / "pcsaft" / "martini2009_polymers.json"
)


def _assert_branches_match(
    eos: EquationOfState,
    *,
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float],
    label: str,
) -> int:
    """Compare both routes at one state; return how many branches were compared."""
    branches = eos.ln_fugacity_branches(
        mixture=mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=list(composition),
    )
    assert branches is not None, f"{label}: the model must offer the capability"

    compared = 0
    for phase in BRANCHES:
        expected = np.asarray(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=list(composition),
                phase=phase,
            ),
            dtype=float,
        )
        assert phase in branches, f"{label}: branch {phase!r} missing"
        found = np.exp(np.asarray(branches[phase], dtype=float))
        assert np.array_equal(found, expected), (
            f"{label} / {phase}: one-solve branch differs from the per-branch call\n"
            f"  per-branch: {expected!r}\n  one-solve:  {found!r}"
        )
        compared += 1
    return compared


def _assert_terms_match(
    eos: EquationOfState,
    *,
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float],
    label: str,
) -> None:
    """``eos_branch_terms_all`` must equal one ``eos_branch_terms`` per branch."""
    w = np.array(composition, dtype=float)
    together = eos_branch_terms_all(
        eos,
        mixture=mixture,
        temperature=temperature_K,
        pressure=pressure_Pa,
        composition=w,
        phases=BRANCHES,
    )
    for phase in BRANCHES:
        try:
            alone = eos_branch_terms(
                eos,
                mixture=mixture,
                temperature=temperature_K,
                pressure=pressure_Pa,
                composition=w,
                phase=phase,
            )
        except ct.ModelError as exc:
            assert phase not in together.terms, f"{label} / {phase}: should have been refused"
            assert together.failures[phase] == str(exc)
            continue

        assert phase in together.terms, f"{label} / {phase}: refused, but the single call worked"
        found = together.terms[phase]
        assert np.array_equal(found.ln_phi, alone.ln_phi), f"{label} / {phase}: ln phi moved"
        if alone.phi is None:
            assert found.phi is None, f"{label} / {phase}: phi appeared where exp underflows"
        else:
            assert found.phi is not None
            assert np.array_equal(found.phi, alone.phi), f"{label} / {phase}: phi moved"


def test_pcsaft_branches_are_bit_identical_over_the_case_f4_subset() -> None:
    eos = PCSAFTEOS()
    compared = 0
    for components, z1, temperature_K, pressure_Pa in F4_SUBSET:
        z = (z1, 1.0 - z1)
        mixture = ct.Mixture.from_database(list(components), list(z), normalize=True)
        label = f"{components} z={z1} T={temperature_K} P={pressure_Pa}"
        compared += _assert_branches_match(
            eos,
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=mixture.fractions,
            label=label,
        )
        _assert_terms_match(
            eos,
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=mixture.fractions,
            label=label,
        )
    assert compared == 2 * len(F4_SUBSET)


def test_pcsaft_branches_are_bit_identical_on_the_water_hexane_states() -> None:
    eos = PCSAFTEOS()
    for temperature_K, pressure_Pa, z in WATER_HEXANE_STATES:
        mixture = ct.Mixture.from_database(["Water", "n-Hexane"], list(z), normalize=True)
        label = f"water/n-hexane T={temperature_K} P={pressure_Pa} z={z}"
        _assert_branches_match(
            eos,
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=mixture.fractions,
            label=label,
        )
        _assert_terms_match(
            eos,
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=mixture.fractions,
            label=label,
        )


def test_peng_robinson_branches_are_bit_identical() -> None:
    for kij in (0.0, 0.0411):
        eos = PengRobinsonEOS(kij=kij)
        for components, z, temperature_K, pressure_Pa in PR_STATES:
            mixture = ct.Mixture.from_database(list(components), list(z), normalize=True)
            _assert_branches_match(
                eos,
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=mixture.fractions,
                label=f"PR kij={kij} {components} T={temperature_K} P={pressure_Pa}",
            )
            _assert_terms_match(
                eos,
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=mixture.fractions,
                label=f"PR kij={kij} {components} T={temperature_K} P={pressure_Pa}",
            )


def test_the_single_root_case_answers_both_labels_with_the_same_root() -> None:
    """Where the model has one admissible root both labels carry it (ADR-0005).

    The one-solve route evaluates that root once and hands the same numbers to
    both labels; the per-branch route evaluates it twice. They must agree
    exactly, or ``_select_density_root_surface``'s degeneracy test - which
    compares the two branches with ``np.array_equal`` - would stop recognising
    a single-root region.
    """
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5], normalize=True)
    roots = eos.density_roots(
        mixture=mixture, temperature_K=298.15, pressure_Pa=1.0e6, composition=[0.5, 0.5]
    )
    assert len(roots) == 1, "this state is chosen for having exactly one root"

    branches = eos.ln_fugacity_branches(
        mixture=mixture, temperature_K=298.15, pressure_Pa=1.0e6, composition=[0.5, 0.5]
    )
    assert branches is not None
    assert np.array_equal(
        np.asarray(branches["vapor"], dtype=float),
        np.asarray(branches["liquid"], dtype=float),
    )


def test_the_polymer_state_keeps_taking_the_logarithmic_route() -> None:
    """Where ``exp(ln phi)`` underflows, both routes must still refuse ``phi``.

    The ADR-0022 guard lives in the *caller*, so the one-solve wrapper has to
    re-apply it and reach the same conclusion: ``phi is None`` and ``ln phi``
    the model's own, not ``log(0.0)``.
    """
    payload = json.loads(POLYMER_FIXTURE.read_text(encoding="utf-8"))
    row = next(entry for entry in payload["polymers"] if entry["name"] == "Polyethylene")
    mw_g_mol = 53000.0
    parameters = ct.PCSAFTParameters.from_records(
        [
            ct.PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=row["segments_per_g"],
                MW_g_mol=mw_g_mol,
                sigma_A=row["sigma_A"],
                epsilon_k_K=row["epsilon_k_K"],
                source="Martini et al. 2009 Table 1 citing Gross & Sadowski 2002",
            ),
            ct.PCSAFTRecord(
                name="n-Pentane",
                m=2.6896,
                sigma_A=3.7729,
                epsilon_k_K=231.20,
                MW_g_mol=72.146,
                source="Gross & Sadowski 2001 Table 1 (the packaged record)",
            ),
        ]
    )
    eos = PCSAFTEOS(parameters=parameters, kij=-0.006)
    moles_polymer = 0.05 / mw_g_mol
    moles_solvent = 0.95 / 72.146
    total = moles_polymer + moles_solvent
    z = [moles_polymer / total, moles_solvent / total]
    mixture = ct.Mixture.from_components(
        [
            ct.Component.custom("Polyethylene", mw_kg_per_mol=mw_g_mol / 1000.0, volatile=False),
            ct.Component.from_database("n-Pentane"),
        ],
        z,
        normalize=True,
    )

    w = np.array(mixture.fractions, dtype=float)
    alone = eos_branch_terms(
        eos,
        mixture=mixture,
        temperature=453.0,
        pressure=8.0e6,
        composition=w,
        phase="liquid",
    )
    assert alone.phi is None, "this state is chosen because exp(ln phi) underflows"

    together = eos_branch_terms_all(
        eos,
        mixture=mixture,
        temperature=453.0,
        pressure=8.0e6,
        composition=w,
        phases=BRANCHES,
    )
    assert together.terms["liquid"].phi is None
    assert np.array_equal(together.terms["liquid"].ln_phi, alone.ln_phi)


class _NoCapabilityEOS(EquationOfState):
    """A model that implements only the pre-ADR-0023 interface."""

    name = "no-capability"

    def __init__(self) -> None:
        self.calls = 0

    def fugacity_coefficients(
        self,
        *,
        mixture: ct.Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> Sequence[float]:
        self.calls += 1
        offset = 1.0 if phase == "vapor" else 2.0
        return [offset + index for index in range(len(composition))]


def test_a_model_without_the_capability_is_served_by_the_per_branch_route() -> None:
    eos = _NoCapabilityEOS()
    assert (
        eos.ln_fugacity_branches(
            mixture=ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5]),
            temperature_K=200.0,
            pressure_Pa=1.0e6,
            composition=[0.5, 0.5],
        )
        is None
    )

    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    branches = eos_branch_terms_all(
        eos,
        mixture=mixture,
        temperature=200.0,
        pressure=1.0e6,
        composition=np.array([0.5, 0.5]),
        phases=BRANCHES,
    )
    assert eos.calls == 2, "one per-branch call per requested label"
    assert set(branches.terms) == set(BRANCHES)
    assert branches.failures == {}
    assert branches.terms["vapor"].phi is not None
    assert branches.terms["vapor"].phi[0] == pytest.approx(1.0)
    assert branches.terms["liquid"].phi is not None
    assert branches.terms["liquid"].phi[0] == pytest.approx(2.0)
