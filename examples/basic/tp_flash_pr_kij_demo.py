"""TP flash demo: PengRobinsonEOS with a per-pair kij mapping.

Shows the mapping form of `kij` (`{(name_a, name_b): value}`) introduced by
the `pr-kij-matrix` slice, and contrasts it with the same feed run at kij = 0
so the effect of the binary interaction parameter is visible in the printed
output. The kij value used (0.0411 for Methane/n-Decane) is an illustrative,
literature-order-of-magnitude value for a light-heavy alkane pair; it is not
sourced from a specific publication and should not be read as a validated
physical parameter -- see `.agents/brain/validation-cases.md` (Case K-1).
"""

from __future__ import annotations

import chemthermo as ct


def _format_composition(
    label: str, components: tuple[ct.Component, ...], fractions: tuple[float, ...]
) -> str:
    lines = [label]
    for component, fraction in zip(components, fractions):
        lines.append(f"  {component.name:<12} {fraction: .6f}")
    return "\n".join(lines)


def _run(
    label: str,
    eos: ct.PengRobinsonEOS,
    mixture: ct.Mixture,
    components: tuple[ct.Component, ...],
    temperature_K: float,
    pressure_Pa: float,
) -> None:
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=eos,
    )

    print(f"--- {label} (kij={eos.kij!r}) ---")
    print(f"Phase names: {', '.join(result.phase_names())}")
    if result.vapor_fraction is not None:
        print(f"Vapor fraction beta: {result.vapor_fraction:.6f}")

    liquid = result.phases.get("liquid")
    vapor = result.phases.get("vapor")
    if liquid is not None:
        print(_format_composition("x (liquid):", components, liquid.composition.fractions))
    if vapor is not None:
        print(_format_composition("y (vapor):", components, vapor.composition.fractions))


def main() -> None:
    temperature_K = 350.0
    pressure_Pa = 3.0e6

    component_names = ("Methane", "n-Decane")
    components = tuple(ct.Component.from_database(name) for name in component_names)
    z = (0.50, 0.50)

    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=z, basis="mole", normalize=False),
    )

    print("TP flash (Peng-Robinson EOS, per-pair kij mapping)")
    print(f"T [K]: {temperature_K:.2f}")
    print(f"P [Pa]: {pressure_Pa:.3e}")
    print(_format_composition("Feed z (mole):", components, mixture.fractions))
    print()

    _run(
        "Default (kij = 0.0)",
        ct.PengRobinsonEOS(),
        mixture,
        components,
        temperature_K,
        pressure_Pa,
    )
    print()
    _run(
        "Per-pair kij mapping",
        ct.PengRobinsonEOS(kij={("Methane", "n-Decane"): 0.0411}),
        mixture,
        components,
        temperature_K,
        pressure_Pa,
    )


if __name__ == "__main__":
    main()
