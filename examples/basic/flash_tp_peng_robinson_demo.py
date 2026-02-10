"""TP flash demo using the Peng-Robinson EOS (SI units)."""

from __future__ import annotations

from typing import Iterable, Mapping

import chemthermo as ct


def _format_composition(
    label: str, components: Iterable[ct.Component], fractions: Iterable[float]
) -> str:
    lines = [f"{label}"]
    for component, fraction in zip(components, fractions):
        lines.append(f"  {component.name:<12} {fraction: .6f}")
    return "\n".join(lines)


def _format_diagnostics(diag: Mapping[str, float | int | str | bool]) -> str:
    keys = [
        "iterations",
        "converged",
        "rr_status",
        "max_delta_k",
        "k_min",
        "k_max",
        "rr_f0",
        "rr_f1",
    ]
    lines = ["Diagnostics:"]
    for key in keys:
        if key not in diag:
            continue
        value = diag[key]
        if isinstance(value, float):
            if key in {"max_delta_k", "rr_f0", "rr_f1"}:
                formatted = f"{value:.3e}"
            else:
                formatted = f"{value:.6g}"
        else:
            formatted = str(value)
        lines.append(f"  {key}: {formatted}")
    if len(lines) == 1:
        lines.append("  (none)")
    return "\n".join(lines)


def main() -> None:
    # All inputs are SI: temperature in K, pressure in Pa.
    temperature_K = 240.0
    pressure_Pa = 3.0e6

    component_names = ("Methane", "Ethane", "Propane")
    components = tuple(ct.Component.from_database(name) for name in component_names)
    z = (0.50, 0.30, 0.20)  # Mole fractions, sum to 1.0.

    composition = ct.Composition(fractions=z, basis="mole", normalize=False)
    mixture = ct.Mixture(components=components, composition=composition)
    eos = ct.PengRobinsonEOS()

    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=eos,
    )

    print("TP flash (Peng-Robinson EOS)")
    print(f"T [K]: {temperature_K:.2f}")
    print(f"P [Pa]: {pressure_Pa:.3e}")
    print(_format_composition("Feed z (mole):", components, mixture.fractions))
    print(f"Vapor fraction beta: {result.vapor_fraction:.6f}")

    liquid = result.phases.get("liquid")
    vapor = result.phases.get("vapor")

    if liquid is not None:
        print(_format_composition("x (liquid):", components, liquid.composition.fractions))
    else:
        print("x (liquid): (not present)")

    if vapor is not None:
        print(_format_composition("y (vapor):", components, vapor.composition.fractions))
    else:
        print("y (vapor): (not present)")

    print(_format_diagnostics(result.diagnostics))


if __name__ == "__main__":
    main()
