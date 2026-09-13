import math

import numpy as np
import pytest

import chemthermo as ct


def test_pr_eos_pure_methane_phi_regression() -> None:
    methane = ct.Component.from_database("Methane")
    mixture = ct.Mixture.from_components([methane], [1.0])
    eos = ct.PengRobinsonEOS()

    phi = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=300.0,
        pressure_Pa=101325.0,
        composition=[1.0],
        phase="vapor",
    )

    assert phi == pytest.approx([0.9977890849], rel=1e-10)


def test_pr_eos_fugacity_coefficients_shape() -> None:
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.25, 0.75], normalize=True)
    eos = ct.PengRobinsonEOS()

    phi = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=280.0,
        pressure_Pa=5.0e6,
        composition=[0.25, 0.75],
        phase="vapor",
    )

    assert len(phi) == 2
    assert all(math.isfinite(value) and value > 0.0 for value in phi)


def test_pr_eos_default_kij_matches_explicit_zero() -> None:
    """kij=0.0 given explicitly must be bit-identical to the implicit default."""
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)

    phi_default = ct.PengRobinsonEOS().fugacity_coefficients(
        mixture=mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        composition=[0.5, 0.5],
        phase="vapor",
    )
    phi_explicit = ct.PengRobinsonEOS(kij=0.0).fugacity_coefficients(
        mixture=mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        composition=[0.5, 0.5],
        phase="vapor",
    )

    assert phi_default == phi_explicit


@pytest.mark.parametrize(
    "kij",
    [0.0, 0.1, {("Methane", "n-Decane"): 0.05}],
    ids=["kij=0.0", "kij=0.1 (scalar)", "kij mapping (pair not present)"],
)
def test_pr_eos_pure_component_invariant_to_kij(kij: object) -> None:
    """A single-component mixture has no i != j cross term, so phi and Z must
    be identical regardless of kij (scalar, nonzero scalar, or a mapping)."""
    methane = ct.Component.from_database("Methane")
    mixture = ct.Mixture.from_components([methane], [1.0])
    baseline = ct.PengRobinsonEOS(kij=0.0)
    candidate = ct.PengRobinsonEOS(kij=kij)  # type: ignore[arg-type]

    for phase in ("vapor", "liquid"):
        phi_base = baseline.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase=phase,
        )
        phi_candidate = candidate.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase=phase,
        )
        assert phi_candidate == pytest.approx(phi_base, rel=1e-14)

        z_base = baseline.compressibility_factor(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase=phase,
        )
        z_candidate = candidate.compressibility_factor(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase=phase,
        )
        assert z_candidate == pytest.approx(z_base, rel=1e-14)


def test_pr_eos_kij_never_corrupts_the_diagonal() -> None:
    """This is the pre-fix bug, made concrete: at the pure-component limit
    within a binary mixture (y = [1, 0]), a_mix must equal a_ii of the present
    component for every kij value, because the (1, 1) entry of the kij matrix
    is always 0 -- never (1 - kij) as the pre-fix code computed it."""
    mixture = ct.Mixture.from_database(["Methane", "n-Decane"], [1.0, 0.0], normalize=False)
    y = np.array([1.0, 0.0])

    a_i_reference, _ = ct.PengRobinsonEOS._component_parameters(mixture, 300.0)

    for kij in (0.0, 0.3, {("Methane", "n-Decane"): 0.0411}):
        eos = ct.PengRobinsonEOS(kij=kij)  # type: ignore[arg-type]
        a_i, _b_i, _aij, a_mix, _b_mix = eos._mixture_parameters(mixture, 300.0, y)
        assert np.allclose(a_i, a_i_reference)
        assert a_mix == pytest.approx(float(a_i_reference[0]), rel=1e-14)


def test_pr_eos_kij_scalar_matches_equivalent_mapping() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Decane"], [0.4, 0.6], normalize=True)
    scalar_eos = ct.PengRobinsonEOS(kij=0.05)
    mapping_eos = ct.PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.05})

    for phase in ("vapor", "liquid"):
        phi_scalar = scalar_eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=350.0,
            pressure_Pa=2.0e6,
            composition=[0.4, 0.6],
            phase=phase,
        )
        phi_mapping = mapping_eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=350.0,
            pressure_Pa=2.0e6,
            composition=[0.4, 0.6],
            phase=phase,
        )
        assert phi_scalar == pytest.approx(phi_mapping, rel=1e-14)


def test_pr_eos_kij_name_normalization_and_symmetry() -> None:
    reference = ct.PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.0411})
    reordered_and_recased = ct.PengRobinsonEOS(kij={("n-decane", "METHANE"): 0.0411})

    assert reference == reordered_and_recased

    mixture = ct.Mixture.from_database(["Methane", "n-Decane"], [0.5, 0.5], normalize=True)
    phi_reference = reference.fugacity_coefficients(
        mixture=mixture,
        temperature_K=350.0,
        pressure_Pa=2.0e6,
        composition=[0.5, 0.5],
        phase="liquid",
    )
    phi_reordered = reordered_and_recased.fugacity_coefficients(
        mixture=mixture,
        temperature_K=350.0,
        pressure_Pa=2.0e6,
        composition=[0.5, 0.5],
        phase="liquid",
    )
    assert phi_reference == phi_reordered


def test_pr_eos_kij_agreeing_duplicate_order_is_not_an_error() -> None:
    eos = ct.PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.0411, ("n-Decane", "Methane"): 0.0411})
    assert eos.kij == ((("methane", "n-decane"), 0.0411),)


def test_pr_eos_kij_conflicting_pair_order_raises_model_error() -> None:
    with pytest.raises(ct.ModelError):
        ct.PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.0411, ("n-Decane", "Methane"): 0.05})


def test_pr_eos_kij_identical_component_pair_raises_model_error() -> None:
    with pytest.raises(ct.ModelError):
        ct.PengRobinsonEOS(kij={("Methane", "Methane"): 0.01})


def test_pr_eos_kij_unknown_pair_is_silently_unused() -> None:
    """A mapping naming components absent from a given mixture is not an
    error: the EOS does not know the mixture at construction time, so an
    unrelated pair is simply never looked up (documented behavior)."""
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    baseline = ct.PengRobinsonEOS(kij=0.0)
    with_unrelated_pair = ct.PengRobinsonEOS(kij={("Propane", "n-Decane"): 0.09})

    phi_baseline = baseline.fugacity_coefficients(
        mixture=mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        composition=[0.5, 0.5],
        phase="vapor",
    )
    phi_with_unrelated_pair = with_unrelated_pair.fugacity_coefficients(
        mixture=mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        composition=[0.5, 0.5],
        phase="vapor",
    )
    assert phi_baseline == phi_with_unrelated_pair


def test_pr_eos_kij_permutation_invariance() -> None:
    """Reordering mixture components under the mapping form must permute phi
    exactly; this is the property a positional matrix aligned to mixture
    order would break (see ADR-0006)."""
    names_forward = ["Methane", "Ethane", "n-Decane"]
    z_forward = [0.5, 0.2, 0.3]
    names_reversed = list(reversed(names_forward))
    z_reversed = list(reversed(z_forward))

    kij_map = {
        ("Methane", "n-Decane"): 0.0411,
        ("Ethane", "n-Decane"): 0.0170,
        ("Methane", "Ethane"): -0.0026,
    }
    eos = ct.PengRobinsonEOS(kij=kij_map)

    mixture_forward = ct.Mixture.from_database(names_forward, z_forward, normalize=True)
    mixture_reversed = ct.Mixture.from_database(names_reversed, z_reversed, normalize=True)

    phi_forward = eos.fugacity_coefficients(
        mixture=mixture_forward,
        temperature_K=350.0,
        pressure_Pa=3.0e6,
        composition=z_forward,
        phase="liquid",
    )
    phi_reversed = eos.fugacity_coefficients(
        mixture=mixture_reversed,
        temperature_K=350.0,
        pressure_Pa=3.0e6,
        composition=z_reversed,
        phase="liquid",
    )

    assert list(phi_forward) == pytest.approx(list(reversed(phi_reversed)), rel=1e-12)


def test_pr_eos_flash_two_phase_split_regression_unchanged() -> None:
    """kij=0.0 (the default) must keep giving the exact pre-existing
    two-phase split; guards against accidental behavior drift from this
    slice's refactor of the mixing-rule code path."""
    mix = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    result = ct.flash_tp(
        mix,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    assert result.vapor_fraction == pytest.approx(0.67451818, rel=1e-6)
