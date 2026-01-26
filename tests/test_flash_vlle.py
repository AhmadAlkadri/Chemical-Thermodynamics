import pytest

import chemthermo as ct


def test_flash_tp_vlle_scaffold() -> None:
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
        activity_model=ct.NRTL(parameters=ct.ActivityParameters.load("NRTL")),
        flash_mode="vlle",
    )

    assert result.diagnostics.get("phase_regime") == "VLLE"
    assert set(result.phase_names()) == {"vapor", "liquid1", "liquid2"}
    assert result.vapor_fraction is not None
    assert sum(result.phase_fractions.values()) == pytest.approx(1.0)
    for phase in result.phases.values():
        assert sum(phase.composition.fractions) == pytest.approx(1.0)
    assert result.diagnostics.get("converged") is True
    lle_delta = result.diagnostics.get("lle_max_delta_x")
    assert isinstance(lle_delta, float)
    assert lle_delta < 1e-6
