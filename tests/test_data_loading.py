import pytest

from chemthermo.data import get_component_record, list_component_names, load_component_database
from chemthermo.exceptions import PropertyNotFoundError


def test_load_component_database_contains_methane() -> None:
    db = load_component_database()
    assert "methane" in db


def test_get_component_record_is_case_insensitive() -> None:
    record = get_component_record("methane")
    assert record["name"] == "Methane"


def test_component_units_are_si() -> None:
    record = get_component_record("Methane")
    props = record["properties"]
    assert props["Tc_K"] == pytest.approx(190.6)
    assert props["Pc_Pa"] == pytest.approx(4.6e6, rel=1e-4)
    assert props["MW_kg_per_mol"] == pytest.approx(0.016042, rel=1e-6)


def test_missing_component_raises() -> None:
    with pytest.raises(PropertyNotFoundError):
        get_component_record("not-a-component")


def test_list_component_names_includes_methane() -> None:
    names = list_component_names()
    assert "Methane" in names
