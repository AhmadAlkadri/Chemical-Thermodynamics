import pytest

import chemthermo as ct


def test_flash_tp_peng_robinson_example_smoke() -> None:
    component_names = ("Methane", "Ethane", "Propane")
    components = tuple(ct.Component.from_database(name) for name in component_names)
    z = (0.50, 0.30, 0.20)

    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=z, basis="mole", normalize=False),
    )
    eos = ct.PengRobinsonEOS()

    result = ct.flash_tp(
        mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=eos,
    )

    assert result.vapor_fraction is not None
    assert 0.0 <= result.vapor_fraction <= 1.0
    assert set(result.phase_names()) == {"liquid", "vapor"}
    assert sum(result.phases["liquid"].composition.fractions) == pytest.approx(1.0)
    assert sum(result.phases["vapor"].composition.fractions) == pytest.approx(1.0)
    assert result.diagnostics.get("converged") is True
