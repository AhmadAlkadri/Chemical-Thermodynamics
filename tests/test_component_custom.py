"""``Component.custom``: a component the databank does not carry (ADR-0022).

Three things are pinned here, and they are the whole of the new contract:

1. a custom component works in a ``Mixture`` with no databank lookup, and its
   molar mass is in **kg/mol** like every other molar mass in the package;
2. the critical constants are optional, and asking for one that was not given
   raises ``PropertyNotFoundError`` instead of returning a placeholder - which
   is the point, because a polymer has no critical point;
3. ``volatile=False`` changes exactly one number, the Wilson K-value
   *estimate*, and changes nothing for any databank component.

The databank schema is checked to be untouched: ``ComponentData`` still
requires ``Tc`` / ``Pc`` / ``omega``, and the relaxed model is a different
class that ``Database`` cannot hold.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from pydantic import ValidationError

import chemthermo as ct
from chemthermo.exceptions import InputRangeError, PropertyNotFoundError
from chemthermo.flash._common import NON_VOLATILE_WILSON_K, wilson_k
from chemthermo.schemas import ComponentData, CustomComponentData, Database, Parameter


def _polymer() -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=16.4,
        formula="(C2H4)n",
        volatile=False,
        source="test fixture",
    )


# ---------------------------------------------------------------------------
# The constructor
# ---------------------------------------------------------------------------


def test_a_custom_component_needs_no_databank_entry() -> None:
    with pytest.raises(Exception):
        ct.Component.from_database("Polyethylene")

    polymer = _polymer()
    assert polymer.name == "Polyethylene"
    assert polymer.formula == "(C2H4)n"
    assert polymer.mw_kg_per_mol == 16.4
    assert polymer.volatile is False
    assert polymer.antoine is None


def test_a_custom_component_enters_a_mixture_beside_a_databank_one() -> None:
    mixture = ct.Mixture.from_components(
        [_polymer(), ct.Component.from_database("n-Pentane")],
        [0.001, 0.999],
    )
    assert mixture.component_names == ["Polyethylene", "n-Pentane"]
    assert mixture.fractions == (0.001, 0.999)


def test_the_critical_constants_are_optional_and_their_absence_is_an_error_not_a_default() -> None:
    polymer = _polymer()
    for name in ("tc_k", "pc_pa", "omega"):
        with pytest.raises(PropertyNotFoundError, match="Component.custom"):
            getattr(polymer, name)


def test_the_critical_constants_are_returned_when_they_were_given() -> None:
    oligomer = ct.Component.custom(
        "Squalane",
        mw_kg_per_mol=0.42282,
        tc_k=863.0,
        pc_pa=8.7e5,
        omega=1.0,
    )
    assert oligomer.tc_k == 863.0
    assert oligomer.pc_pa == 8.7e5
    assert oligomer.omega == 1.0
    assert oligomer.volatile is True


@pytest.mark.parametrize(
    "kwargs",
    [
        {"mw_kg_per_mol": 0.0},
        {"mw_kg_per_mol": -1.0},
        {"mw_kg_per_mol": math.nan},
        {"mw_kg_per_mol": 1.0, "tc_k": -5.0},
        {"mw_kg_per_mol": 1.0, "pc_pa": 0.0},
        {"mw_kg_per_mol": 1.0, "omega": math.inf},
    ],
)
def test_unphysical_inputs_are_refused(kwargs: dict[str, float]) -> None:
    with pytest.raises(InputRangeError):
        ct.Component.custom(
            "X",
            mw_kg_per_mol=kwargs["mw_kg_per_mol"],
            tc_k=kwargs.get("tc_k"),
            pc_pa=kwargs.get("pc_pa"),
            omega=kwargs.get("omega"),
        )


# ---------------------------------------------------------------------------
# The Wilson estimate
# ---------------------------------------------------------------------------


def test_a_non_volatile_component_gets_the_fixed_wilson_estimate() -> None:
    mixture = ct.Mixture.from_components(
        [_polymer(), ct.Component.from_database("n-Pentane")],
        [0.001, 0.999],
    )
    k = wilson_k(mixture, 453.0, 1.0e7)
    assert k[0] == NON_VOLATILE_WILSON_K
    assert k[1] > 0.0

    # The vapour-like start puts the polymer at ~1e-13 and the liquid-like one
    # makes it dominant: the two estimates a polymer/solvent feed needs.
    z = np.asarray(mixture.fractions)
    vapor_like = (z * k) / float(np.sum(z * k))
    liquid_like = (z / k) / float(np.sum(z / k))
    assert vapor_like[0] < 1e-10
    assert liquid_like[0] > 0.99


def test_the_wilson_estimate_is_unchanged_for_databank_components() -> None:
    """Every packaged component is volatile, so nothing pre-ADR-0022 moved."""
    mixture = ct.Mixture.from_database(["Methane", "n-Decane"], [0.9, 0.1])
    k = wilson_k(mixture, 240.0, 1.0e6)
    expected = [
        math.exp(
            math.log(component.pc_pa / 1.0e6)
            + 5.373 * (1.0 + component.omega) * (1.0 - component.tc_k / 240.0)
        )
        for component in mixture.components
    ]
    assert k.tolist() == expected
    assert all(component.volatile for component in mixture.components)


def test_a_volatile_custom_component_still_uses_the_correlation() -> None:
    volatile = ct.Component.custom("Pseudo", mw_kg_per_mol=0.1, tc_k=500.0, pc_pa=3.0e6, omega=0.3)
    mixture = ct.Mixture.from_components(
        [volatile, ct.Component.from_database("n-Pentane")], [0.5, 0.5]
    )
    k = wilson_k(mixture, 400.0, 1.0e6)
    assert k[0] == pytest.approx(
        math.exp(math.log(3.0e6 / 1.0e6) + 5.373 * 1.3 * (1.0 - 500.0 / 400.0))
    )


# ---------------------------------------------------------------------------
# The databank schema is untouched
# ---------------------------------------------------------------------------


def test_the_databank_record_still_requires_the_critical_constants() -> None:
    payload = {
        "name": "X",
        "formula": "X",
        "MW": {"value": 0.1, "units": "kg/mol"},
    }
    with pytest.raises(ValidationError):
        ComponentData(**payload)
    # The relaxed model accepts the same payload; that is the only difference.
    assert CustomComponentData(**payload).Tc is None


def test_a_custom_record_cannot_enter_the_databank_container() -> None:
    relaxed = CustomComponentData(name="X", formula="X", MW=Parameter(value=0.1, units="kg/mol"))
    with pytest.raises(ValidationError):
        Database(schema_version=1, components=[relaxed])  # type: ignore[list-item]


def test_the_packaged_databank_still_loads_at_schema_version_one() -> None:
    from chemthermo.data import load_component_database

    records = load_component_database()
    assert records
    for record in records.values():
        model = ComponentData(**record)
        assert model.Tc.value > 0.0
        assert model.Pc.value > 0.0
