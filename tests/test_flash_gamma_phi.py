import pytest

import chemthermo as ct


def test_flash_tp_gamma_phi_two_phase() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
    )

    result = ct.flash_tp(
        mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=ct.PengRobinsonEOS(),
        activity_model=ct.NRTL(),
        flash_mode="gamma-phi",
    )

    assert result.vapor_fraction is not None
    assert 0.0 < result.vapor_fraction < 1.0
    assert set(result.phase_names()) == {"liquid", "vapor"}
    assert sum(result.phases["liquid"].composition.fractions) == pytest.approx(1.0)
    assert sum(result.phases["vapor"].composition.fractions) == pytest.approx(1.0)
    assert result.diagnostics.get("converged") is True
    assert result.diagnostics.get("phase_count") == 2
    assert result.diagnostics.get("phase_state") == "two_phase"


def test_flash_tp_gamma_phi_single_phase() -> None:
    component = ct.Component.from_database("Methane")
    mixture = ct.Mixture(
        components=(component,),
        composition=ct.Composition(fractions=(1.0,), basis="mole", normalize=False),
    )

    result = ct.flash_tp(
        mixture,
        temperature_K=350.0,
        pressure_Pa=101325.0,
        eos=ct.PengRobinsonEOS(),
        activity_model=ct.NRTL(),
        flash_mode="gamma-phi",
    )

    assert result.phase_names() == ["vapor"]
    assert result.diagnostics.get("phase_count") == 1
    assert result.diagnostics.get("phase_state") == "vapor"
