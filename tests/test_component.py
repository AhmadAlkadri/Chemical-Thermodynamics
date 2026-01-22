import pytest

import chemthermo as ct


def test_component_from_database_si_units() -> None:
    comp = ct.Component.from_database("Methane")
    assert comp.mw_kg_per_mol == pytest.approx(0.016042)
    assert comp.tc_k == pytest.approx(190.6)
    assert comp.pc_pa == pytest.approx(4.6e6, rel=1e-4)


def test_component_missing_property_raises() -> None:
    comp = ct.Component(name="dummy", formula="X", properties={})
    with pytest.raises(ct.PropertyNotFoundError):
        _ = comp.pc_pa
