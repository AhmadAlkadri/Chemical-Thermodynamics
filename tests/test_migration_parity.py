"""Verify database loading and minimal parity check."""

import pytest
from chemthermo.core.component import Component
from chemthermo.data import list_component_names
from chemthermo import cite

def test_list_components():
    names = list_component_names()
    assert len(names) > 0
    assert "Methane" in names
    assert "Water" in names

def test_component_loading():
    comp = Component.from_database("Methane")
    assert comp.name == "Methane"
    assert comp.formula == "CH4"
    # value is roughly 16.04 g/mol -> 0.01604 kg/mol
    assert abs(comp.mw_kg_per_mol - 0.016042) < 1e-5
    
    # Check citation
    cit = cite("Methane", "Tc")
    # Should resolve to something, not "Unknown citation" if we did it right
    # But wait, migrate_legacy_db uses "legacy_db" key.
    # And references.bib has "legacy_db".
    assert "Legacy Chemical Thermodynamics Database" in cit

def test_legacy_parity_spot_check():
    """Spot check values against known legacy values."""
    # Water: Tc=647.3 K, Pc=220.48 bar
    h2o = Component.from_database("Water")
    assert abs(h2o.tc_k - 647.3) < 1e-1
    assert abs(h2o.pc_pa - 220.48e5) < 1e3  # 1e5 Pa = 1 bar tolerance? No, 220.48 bar is 22MPa.
    # 220.48 * 1e5 = 2.2048e7. Tolerance 1e3 is 1kPa.
