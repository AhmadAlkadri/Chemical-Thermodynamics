import pytest

import chemthermo as ct


def test_component_from_database_si_units() -> None:
    comp = ct.Component.from_database("Methane")
    assert comp.mw_kg_per_mol == pytest.approx(0.016042)
    assert comp.tc_k == pytest.approx(190.6)
    assert comp.pc_pa == pytest.approx(4.6e6, rel=1e-4)


def test_component_loading_works_without_bibliography(monkeypatch) -> None:
    import chemthermo.core.component as component_module

    class BrokenBibliography:
        def __init__(self, path=None):
            raise FileNotFoundError("bibliography not available")

    monkeypatch.setattr(component_module, "BibliographicDatabase", BrokenBibliography)

    comp = ct.Component.from_database("Methane")
    assert comp.name == "Methane"
    assert comp.tc_k == pytest.approx(190.6)


def test_cite_returns_reference_when_bibliography_present() -> None:
    citation = ct.cite("Methane", "Tc")
    assert "Koretsky" in citation
