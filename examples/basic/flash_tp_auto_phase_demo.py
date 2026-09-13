"""Automatic 1-vs-2 phase detection in `flash_tp` (ADR-0008).

Three phi-phi states, all Peng-Robinson with `kij = 0` and the packaged
databank constants:

1. A stable feed: `flash_tp` returns one phase because Michelsen's
   tangent-plane test found no negative tangent-plane distance, not because a
   Wilson K-value bound said so.
2. An unstable feed: the tangent-plane minimizer seeds the split, and the
   converged two-phase result is verified (material balance, equal fugacities,
   Gibbs energy below the single-phase feed).
3. The state where the legacy heuristic and the tangent-plane path *disagree*.
   Methane(0.6) / n-Pentane(0.4) at 175 K and 1.778 MPa: the Wilson K-values
   straddle 1 so the K-bound test does not fire, Rachford-Rice finds no root
   for them, and the legacy path reports a single liquid. The feed is in fact
   unstable. `thermo`'s `FlashVL` with the same critical constants gives
   VF = 0.0791276 (validation Case F-2); this script needs no optional
   dependency and prints the Gibbs-energy evidence instead.

Run:

    python examples/basic/flash_tp_auto_phase_demo.py
"""

from __future__ import annotations

from typing import Mapping, Sequence

import chemthermo as ct

LEGACY = ct.FlashSettings(phase_detection="wilson-heuristic")

_DIAGNOSTIC_KEYS = (
    "phase_detection",
    "stability_status",
    "tpd_min",
    "feed_branch",
    "stability_trials",
    "k_seed",
    "incipient_phase",
    "iterations",
    "converged",
    "termination_reason",
    "max_delta_k",
    "k_min",
    "k_max",
    "mass_balance_residual",
    "fugacity_residual",
    "delta_g_split_rt",
    "rr_status",
    "rr_f0",
    "rr_f1",
)
_EXPONENTIAL_KEYS = frozenset(
    {"max_delta_k", "rr_f0", "rr_f1", "mass_balance_residual", "fugacity_residual"}
)


def _print_diagnostics(diagnostics: Mapping[str, float | int | str | bool]) -> None:
    for key in _DIAGNOSTIC_KEYS:
        if key not in diagnostics:
            continue
        value = diagnostics[key]
        if isinstance(value, float):
            formatted = f"{value:.3e}" if key in _EXPONENTIAL_KEYS else f"{value:.6g}"
        else:
            formatted = str(value)
        print(f"    {key}: {formatted}")


def _print_phases(names: Sequence[str], result: ct.FlashResult) -> None:
    beta = result.vapor_fraction
    print(f"    phases: {', '.join(result.phase_names())}")
    if beta is not None:
        print(f"    vapor fraction beta: {beta:.6f}")
    for phase_name in result.phase_names():
        fractions = result.phases[phase_name].composition.fractions
        rendered = ", ".join(f"{name}={fraction:.6f}" for name, fraction in zip(names, fractions))
        print(f"    {phase_name}: {rendered}")


def _run(
    title: str,
    names: Sequence[str],
    z: Sequence[float],
    temperature_K: float,
    pressure_Pa: float,
    *,
    also_legacy: bool = False,
) -> None:
    print(title)
    print(
        f"  z = ({', '.join(f'{name} {value:.2f}' for name, value in zip(names, z))}), "
        f"T = {temperature_K:.1f} K, P = {pressure_Pa:.4e} Pa"
    )
    mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)
    eos = ct.PengRobinsonEOS()

    print("  tangent-plane phase detection (default):")
    result = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=eos)
    _print_phases(names, result)
    _print_diagnostics(result.diagnostics)

    if also_legacy:
        print("  legacy Wilson K-bound heuristic (phase_detection='wilson-heuristic'):")
        try:
            legacy = ct.flash_tp(
                mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=eos,
                settings=LEGACY,
            )
        except ct.ConvergenceError as exc:
            print(f"    ConvergenceError: {exc}")
        else:
            _print_phases(names, legacy)
            _print_diagnostics(legacy.diagnostics)
    print()


def main() -> None:
    print("Automatic phase detection in flash_tp (Peng-Robinson, kij = 0)")
    print()

    _run(
        "1) Stable feed -> one phase, because the tangent-plane test says so",
        ("Methane", "Ethane"),
        (0.5, 0.5),
        450.0,
        1.0e5,
    )

    _run(
        "2) Unstable feed -> split seeded from the tangent-plane minimizer",
        ("Methane", "Ethane", "Propane"),
        (0.5, 0.3, 0.2),
        240.0,
        3.0e6,
    )

    _run(
        "3) Where the two paths DISAGREE (validation Case F-2)",
        ("Methane", "n-Pentane"),
        (0.6, 0.4),
        175.0,
        1.778e6,
        also_legacy=True,
    )

    print("Reading of case 3:")
    print("  The legacy path returns a single liquid with termination_reason")
    print("  'rr_no_root': the Wilson K-values straddle 1, so the K-bound test")
    print("  does not fire, but Rachford-Rice cannot bracket a root for them.")
    print("  The tangent-plane test finds tpd_min = -4.266e-02 < 0, which is a")
    print("  proof that the feed is not one phase, and the split it seeds has")
    print("  delta_g_split_rt < 0, i.e. a lower Gibbs energy than the feed.")
    print("  thermo's FlashVL with the same Tc/Pc/omega and kij = 0 reports")
    print("  VF = 0.0791276 at this state (validation Case F-2).")
    print()
    print("Limits of this release: at most two phases, and the converged phases")
    print("are not themselves re-tested for stability (next slice). 'stable'")
    print("means no negative tangent-plane distance was found from the")
    print("deterministic trial set, not a global proof.")


if __name__ == "__main__":
    main()
