"""Tests for the new Pydantic schema."""

import json
from pathlib import Path

from chemthermo.schemas import ComponentData, Database, Parameter


def test_parameter_model():
    """Test Parameter validation."""
    p = Parameter(value=100.0, units="K", source_key="test_ref")
    assert p.value == 100.0
    assert p.units == "K"
    assert p.source_key == "test_ref"


def test_component_model():
    """Test ComponentData validation."""
    c = ComponentData(
        name="Test",
        formula="C",
        MW=Parameter(value=12.01, units="g/mol"),
        Tc=Parameter(value=300.0, units="K"),
        Pc=Parameter(value=10.0, units="bar"),
        omega=Parameter(value=0.1, units="-"),
    )
    assert c.name == "Test"
    assert c.MW.value == 12.01


def test_database_model():
    """Test Database validation."""
    db = Database(components=[])
    assert db.schema_version == 1
    assert db.components == []


def test_load_template_db():
    """Test loading the template components.json."""
    db_path = Path("database/components.json")
    assert db_path.exists()

    with open(db_path) as f:
        data = json.load(f)

    # process components list
    db = Database(**data)
    assert db.schema_version == 1
