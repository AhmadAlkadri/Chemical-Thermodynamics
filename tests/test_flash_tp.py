import pytest

import chemthermo as ct


def test_flash_tp_single_phase_vapor_detection() -> None:
    mix = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    eos = ct.PengRobinsonEOS()

    result = ct.flash_tp(
        mix,
        temperature_K=450.0,
        pressure_Pa=1.0e5,
        eos=eos,
    )

    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == pytest.approx(1.0)
    assert result.diagnostics.get("iterations") == 0
    assert result.diagnostics.get("converged") is True


def test_flash_tp_two_phase_split() -> None:
    mix = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    eos = ct.PengRobinsonEOS()

    result = ct.flash_tp(
        mix,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=eos,
    )

    assert set(result.phase_names()) == {"liquid", "vapor"}
    assert 0.0 < (result.vapor_fraction or 0.0) < 1.0
    assert result.diagnostics.get("converged") is True
    assert result.vapor_fraction == pytest.approx(0.67451818, rel=1e-6)


def test_flash_tp_convergence_failure() -> None:
    mix = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    eos = ct.PengRobinsonEOS()
    settings = ct.FlashSettings(max_iter=1, tol=1e-12)

    with pytest.raises(ct.ConvergenceError):
        ct.flash_tp(
            mix,
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=eos,
            settings=settings,
        )
