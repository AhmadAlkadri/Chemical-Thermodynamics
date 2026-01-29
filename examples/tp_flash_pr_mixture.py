"""TP flash demo using the Peng-Robinson EOS (mixture)."""

from __future__ import annotations

import chemthermo as ct


def _format_composition(
    label: str, components: tuple[ct.Component, ...], fractions: tuple[float, ...]
) -> str:
    lines = [label]
    for component, fraction in zip(components, fractions):
        lines.append(f"  {component.name:<12} {fraction: .6f}")
    return "\n".join(lines)


def _k_values(liquid: ct.PhaseResult, vapor: ct.PhaseResult) -> list[float]:
    ks = []
    for x_i, y_i in zip(liquid.composition.fractions, vapor.composition.fractions):
        ks.append(y_i / x_i if x_i > 0.0 else float("inf"))
    return ks


def main() -> None:
    temperature_K = 240.0
    pressure_Pa = 3.0e6

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
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=eos,
    )

    print("TP flash (Peng-Robinson EOS, mixture)")
    print(f"T [K]: {temperature_K:.2f}")
    print(f"P [Pa]: {pressure_Pa:.3e}")
    print(_format_composition("Feed z (mole):", components, mixture.fractions))

    if result.vapor_fraction is not None:
        print(f"Vapor fraction beta: {result.vapor_fraction:.6f}")

    liquid = result.phases.get("liquid")
    vapor = result.phases.get("vapor")

    if liquid is not None:
        print(_format_composition("x (liquid):", components, liquid.composition.fractions))
    if vapor is not None:
        print(_format_composition("y (vapor):", components, vapor.composition.fractions))

    if liquid is not None and vapor is not None:
        ks = _k_values(liquid, vapor)
        print("K-values (y/x):")
        for component, value in zip(components, ks):
            print(f"  {component.name:<12} {value: .6f}")

    if result.diagnostics:
        print("Diagnostics:")
        for key, value in result.diagnostics.items():
            print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
