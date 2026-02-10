import shutil
from pathlib import Path

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
    # New schema: properties are at top level as Parameter objects (dicts)
    assert record["Tc"]["value"] == pytest.approx(190.6)
    assert record["Pc"]["value"] == pytest.approx(4.6e6, rel=1e-4)
    assert record["MW"]["value"] == pytest.approx(0.016042, rel=1e-6)


def test_missing_component_raises() -> None:
    with pytest.raises(PropertyNotFoundError):
        get_component_record("not-a-component")


def test_list_component_names_includes_methane() -> None:
    names = list_component_names()
    assert "Methane" in names


def test_runtime_loader_does_not_depend_on_legacy_database_copy(tmp_path: Path) -> None:
    legacy_copy = Path("database/components.json")
    if not legacy_copy.exists():
        pytest.skip("Legacy authoring copy is absent.")

    backup = tmp_path / "components_legacy_backup.json"
    shutil.move(str(legacy_copy), backup)
    try:
        load_component_database.cache_clear()
        record = get_component_record("Methane")
        assert record["name"] == "Methane"
    finally:
        shutil.move(str(backup), legacy_copy)
        load_component_database.cache_clear()
