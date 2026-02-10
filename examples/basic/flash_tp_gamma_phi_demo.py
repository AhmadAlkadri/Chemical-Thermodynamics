"""TP flash demo using gamma-phi (NRTL liquid + Peng-Robinson vapor)."""

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
        "convergence_reason",
        "phase_state",
        "phase_count",
        "flash_mode",
        "max_delta_k",
        "k_min",
        "k_max",
    ]
    lines = ["Diagnostics:"]
    for key in keys:
        if key not in diag:
            continue
        value = diag[key]
        if isinstance(value, float):
            if key == "max_delta_k":
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

    component_names = ("Methane", "Ethane")
    components = tuple(ct.Component.from_database(name) for name in component_names)
    z = (0.50, 0.50)  # Mole fractions, sum to 1.0.

    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=z, basis="mole", normalize=False),
    )

    eos = ct.PengRobinsonEOS()
    activity_params = ct.ActivityParameters.load("NRTL")
    activity_model = ct.NRTL(parameters=activity_params)

    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=eos,
        activity_model=activity_model,
        flash_mode="gamma-phi",
    )

    print("TP flash (gamma-phi, NRTL + Peng-Robinson EOS)")
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
