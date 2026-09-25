"""Command-line interface for chemthermo."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Sequence
from typing import Any

from .core import Mixture
from .eos import PCSAFTEOS
from .exceptions import ConvergenceError, ThermoError
from .flash import FlashSettings, flash_tp
from .models import NRTL, PengRobinsonEOS
from .models.base import EquationOfState
from .stability import StabilitySettings, stability_tp

CLI_SCHEMA_VERSION = 1

#: `--eos` choice -> the `solver.eos` label written to the payload (ADR-0033).
EOS_LABELS = {"peng-robinson": "peng_robinson", "pc-saft": "pc_saft"}

#: What "stable" means here (ADR-0005 honesty note, ADR-0033 decision 3).
STABILITY_SCOPE = "bounded-trial-set"
STABLE_MEANING = (
    "no negative tangent-plane distance found from the deterministic trial set (not a global proof)"
)


class _UsageError(Exception):
    """A flag combination argparse cannot express; mapped to exit code 2."""


def _make_eos(choice: str) -> EquationOfState:
    """The packaged-parameter model behind an `--eos` choice."""
    if choice == "pc-saft":
        return PCSAFTEOS()
    return PengRobinsonEOS()


def _json_safe(value: Any) -> Any:
    """Non-finite floats as strings, so stdout stays valid JSON (ADR-0033)."""
    if isinstance(value, float) and not math.isfinite(value):
        return "nan" if math.isnan(value) else ("inf" if value > 0 else "-inf")
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _parse_csv_strings(raw: str, *, label: str) -> list[str]:
    values = [item.strip() for item in raw.split(",")]
    filtered = [item for item in values if item]
    if not filtered:
        raise ValueError(f"{label} must include at least one comma-separated value.")
    return filtered


def _parse_csv_floats(raw: str, *, label: str) -> list[float]:
    tokens = _parse_csv_strings(raw, label=label)
    values: list[float] = []
    for token in tokens:
        try:
            values.append(float(token))
        except ValueError as exc:
            raise ValueError(f"{label} contains a non-numeric value: {token!r}.") from exc
    return values


def _algorithm_label(*, flash_mode: str, settings: FlashSettings) -> str:
    """Name the path the solver actually took.

    The phase split is unchanged (Rachford-Rice plus fixed-point K updates);
    only the initialization and the one-versus-two-phase decision differ, so the
    first token names the phase-detection mode. Gamma-phi always uses the
    Wilson heuristic (ADR-0008).
    """
    if flash_mode == "phi-phi" and settings.phase_detection == "tangent-plane":
        return "tangent-plane+rachford-rice+fixed-point"
    return "wilson+rachford-rice+fixed-point"


def _tp_flash_payload(
    *,
    component_names: list[str],
    feed_fractions: list[float],
    temperature_K: float,
    pressure_Pa: float,
    normalize: bool,
    flash_mode: str,
    settings: FlashSettings,
    result: Any,
    eos_label: str = "peng_robinson",
    max_phases_given: bool = False,
) -> dict[str, Any]:
    phase_names = result.phase_names()

    phases: dict[str, dict[str, list[float]]] = {}
    for phase_name in phase_names:
        phase = result.phases[phase_name]
        phases[phase_name] = {
            "fractions": list(phase.composition.fractions),
        }

    solver_settings: dict[str, Any] = {
        "max_iter": settings.max_iter,
        "tol": settings.tol,
        "damping": settings.damping,
    }
    if max_phases_given:
        # Additive, and only when asked for, so existing invocations keep
        # their exact bytes (ADR-0033 decision 1).
        solver_settings["max_phases"] = settings.max_phases

    return {
        "cli_schema_version": CLI_SCHEMA_VERSION,
        "command": "tp-flash",
        "inputs": {
            "components": component_names,
            "z_mole": feed_fractions,
            "temperature_K": temperature_K,
            "pressure_Pa": pressure_Pa,
            "normalize": normalize,
        },
        "solver": {
            "eos": eos_label,
            "method": flash_mode,
            "algorithm": _algorithm_label(flash_mode=flash_mode, settings=settings),
            "settings": solver_settings,
        },
        "result": {
            "component_order": component_names,
            "phase_names": phase_names,
            "vapor_fraction": result.vapor_fraction,
            "phase_fractions": dict(result.phase_fractions),
            "phases": phases,
        },
        "diagnostics": dict(result.diagnostics),
    }


def _format_tp_flash_text(payload: dict[str, Any]) -> str:
    lines = ["TP flash (chemthermo CLI)"]
    inputs = payload["inputs"]
    result = payload["result"]

    lines.append(f"T [K]: {inputs['temperature_K']:.6f}")
    lines.append(f"P [Pa]: {inputs['pressure_Pa']:.6e}")
    lines.append("Feed z (mole):")
    for name, z_i in zip(inputs["components"], inputs["z_mole"]):
        lines.append(f"  {name:<12} {z_i: .6f}")

    vapor_fraction = result["vapor_fraction"]
    if vapor_fraction is not None:
        lines.append(f"Vapor fraction beta: {vapor_fraction:.6f}")

    for phase_name in result["phase_names"]:
        lines.append(f"{phase_name} composition:")
        fractions = result["phases"][phase_name]["fractions"]
        for name, fraction in zip(result["component_order"], fractions):
            lines.append(f"  {name:<12} {fraction: .6f}")

    diagnostics = payload["diagnostics"]
    lines.append("Diagnostics:")
    if diagnostics:
        for key in sorted(diagnostics):
            value = diagnostics[key]
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.6g}")
            else:
                lines.append(f"  {key}: {value}")
    else:
        lines.append("  (none)")

    return "\n".join(lines)


def _run_tp_flash(args: argparse.Namespace) -> int:
    component_names = _parse_csv_strings(args.components, label="--components")
    feed_fractions = _parse_csv_floats(args.z, label="--z")

    if len(component_names) != len(feed_fractions):
        raise ValueError(
            "--components and --z must contain the same number of comma-separated values."
        )

    flash_mode = args.flash_mode
    if args.eos == "pc-saft" and flash_mode == "gamma-phi":
        raise _UsageError("--eos pc-saft supports --flash-mode phi-phi only.")

    mixture = Mixture.from_database(component_names, feed_fractions, normalize=args.normalize)
    settings_kwargs: dict[str, Any] = {
        "max_iter": args.max_iter,
        "tol": args.tol,
        "damping": args.damping,
    }
    if args.max_phases is not None:
        settings_kwargs["max_phases"] = args.max_phases
    settings = FlashSettings(**settings_kwargs)
    activity_model = NRTL() if flash_mode == "gamma-phi" else None

    result = flash_tp(
        mixture,
        temperature_K=args.temperature_k,
        pressure_Pa=args.pressure_pa,
        eos=_make_eos(args.eos),
        activity_model=activity_model,
        flash_mode=flash_mode,
        settings=settings,
    )

    payload = _tp_flash_payload(
        component_names=component_names,
        feed_fractions=feed_fractions,
        temperature_K=args.temperature_k,
        pressure_Pa=args.pressure_pa,
        normalize=args.normalize,
        flash_mode=flash_mode,
        settings=settings,
        result=result,
        eos_label=EOS_LABELS[args.eos],
        max_phases_given=args.max_phases is not None,
    )

    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True))
    else:
        print(_format_tp_flash_text(payload))

    return 0


def _stability_payload(
    *,
    component_names: list[str],
    feed_fractions: list[float],
    temperature_K: float,
    pressure_Pa: float,
    normalize: bool,
    eos_label: str,
    settings: StabilitySettings,
    result: Any,
) -> dict[str, Any]:
    diagnostics = dict(result.diagnostics)
    return _json_safe(
        {
            "cli_schema_version": CLI_SCHEMA_VERSION,
            "command": "stability-tp",
            "inputs": {
                "components": component_names,
                "z_mole": feed_fractions,
                "temperature_K": temperature_K,
                "pressure_Pa": pressure_Pa,
                "normalize": normalize,
            },
            "solver": {
                "eos": eos_label,
                "method": "tangent-plane",
                "settings": {
                    "max_iter": settings.max_iter,
                    "tol": settings.tol,
                    "tpd_tol": settings.tpd_tol,
                    "trivial_tol": settings.trivial_tol,
                },
            },
            "result": {
                "status": result.status,
                "stable": result.stable,
                "stability_scope": STABILITY_SCOPE,
                "tpd_min": result.tpd_min,
                "component_order": component_names,
                "trial_composition": (
                    None if result.trial_composition is None else list(result.trial_composition)
                ),
                "k_values": None if result.k_values is None else list(result.k_values),
                "feed_branch": result.feed_branch,
                "phase_branch": result.phase_branch,
                "trial_count": len(result.trials),
                "converged_trial_count": sum(1 for trial in result.trials if trial.converged),
                "non_trivial_trial_count": sum(
                    1 for trial in result.trials if trial.converged and not trial.trivial
                ),
            },
            "diagnostics": diagnostics,
        }
    )


def _format_stability_text(payload: dict[str, Any]) -> str:
    inputs = payload["inputs"]
    result = payload["result"]
    lines = ["TP stability (chemthermo CLI)"]
    lines.append(f"EOS: {payload['solver']['eos']}")
    lines.append(f"T [K]: {inputs['temperature_K']:.6f}")
    lines.append(f"P [Pa]: {inputs['pressure_Pa']:.6e}")
    lines.append("Feed z (mole):")
    for name, z_i in zip(inputs["components"], inputs["z_mole"]):
        lines.append(f"  {name:<12} {z_i: .6f}")
    lines.append(f"Status: {result['status']}")
    if result["status"] == "stable":
        lines.append(f"  meaning: {STABLE_MEANING}")
    lines.append(f"tpd_min [RT]: {result['tpd_min']:.6g}")
    lines.append(f"Feed branch: {result['feed_branch']}")
    if result["trial_composition"] is not None:
        lines.append(f"Incipient phase ({result['phase_branch']}) composition:")
        for name, w_i in zip(result["component_order"], result["trial_composition"]):
            lines.append(f"  {name:<12} {w_i: .6f}")
    lines.append(
        "Trials: "
        f"{result['trial_count']} run, {result['converged_trial_count']} converged, "
        f"{result['non_trivial_trial_count']} non-trivial"
    )
    return "\n".join(lines)


def _run_stability_tp(args: argparse.Namespace) -> int:
    component_names = _parse_csv_strings(args.components, label="--components")
    feed_fractions = _parse_csv_floats(args.z, label="--z")
    if len(component_names) != len(feed_fractions):
        raise ValueError(
            "--components and --z must contain the same number of comma-separated values."
        )

    mixture = Mixture.from_database(component_names, feed_fractions, normalize=args.normalize)
    settings = StabilitySettings()
    result = stability_tp(
        mixture,
        temperature_K=args.temperature_k,
        pressure_Pa=args.pressure_pa,
        eos=_make_eos(args.eos),
        settings=settings,
    )
    payload = _stability_payload(
        component_names=component_names,
        feed_fractions=feed_fractions,
        temperature_K=args.temperature_k,
        pressure_Pa=args.pressure_pa,
        normalize=args.normalize,
        eos_label=EOS_LABELS[args.eos],
        settings=settings,
        result=result,
    )

    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True))
    else:
        print(_format_stability_text(payload))

    # ADR-0033 decision 4: no verdict is not an answer, but show the evidence.
    return 3 if result.status == "inconclusive" else 0


def _positive_int(raw: str) -> int:
    try:
        value = int(raw)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected an integer, got {raw!r}") from exc
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {value}")
    return value


def _add_state_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--components",
        required=True,
        help="Comma-separated component names, e.g. Methane,Ethane,Propane",
    )
    parser.add_argument(
        "--z",
        required=True,
        help="Comma-separated feed mole fractions in component order.",
    )
    parser.add_argument("--temperature-k", type=float, required=True, help="Temperature in K.")
    parser.add_argument("--pressure-pa", type=float, required=True, help="Pressure in Pa.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="chemthermo", description="chemthermo command-line tools")
    subparsers = parser.add_subparsers(dest="command", required=True)

    tp_flash = subparsers.add_parser(
        "tp-flash",
        help="Run TP flash: Peng-Robinson (phi-phi or gamma-phi) or PC-SAFT (phi-phi).",
    )
    tp_flash.add_argument(
        "--components",
        required=True,
        help="Comma-separated component names, e.g. Methane,Ethane,Propane",
    )
    tp_flash.add_argument(
        "--z",
        required=True,
        help="Comma-separated feed mole fractions in component order.",
    )
    tp_flash.add_argument("--temperature-k", type=float, required=True, help="Temperature in K.")
    tp_flash.add_argument("--pressure-pa", type=float, required=True, help="Pressure in Pa.")
    tp_flash.add_argument(
        "--flash-mode",
        choices=("phi-phi", "gamma-phi"),
        default="phi-phi",
        help="Flash method selection.",
    )
    tp_flash.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize feed composition if fractions do not sum exactly to 1.",
    )
    tp_flash.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format.",
    )
    tp_flash.add_argument(
        "--max-iter",
        type=int,
        default=100,
        help="Maximum flash iterations.",
    )
    tp_flash.add_argument(
        "--tol",
        type=float,
        default=1e-8,
        help="Flash convergence tolerance.",
    )
    tp_flash.add_argument(
        "--damping",
        type=float,
        default=None,
        help="Optional damping factor for K-value updates.",
    )
    tp_flash.add_argument(
        "--eos",
        choices=tuple(EOS_LABELS),
        default="peng-robinson",
        help=(
            "Equation of state (packaged parameters; pc-saft uses kij = 0 and "
            "supports --flash-mode phi-phi only)."
        ),
    )
    tp_flash.add_argument(
        "--max-phases",
        type=_positive_int,
        default=None,
        help="Largest number of phases the flash may return (library default: 3).",
    )
    tp_flash.set_defaults(handler=_run_tp_flash)

    stability = subparsers.add_parser(
        "stability-tp",
        help="Tangent-plane stability test of a feed (Peng-Robinson or PC-SAFT).",
        description=(
            "Michelsen tangent-plane stability test at fixed T, P and feed. "
            f"'stable' means {STABLE_MEANING}. Exit code 3 when no trial "
            "converged (status 'inconclusive'); the JSON is still printed."
        ),
    )
    _add_state_arguments(stability)
    stability.add_argument(
        "--eos",
        choices=tuple(EOS_LABELS),
        default="peng-robinson",
        help="Equation of state (packaged parameters; pc-saft uses kij = 0).",
    )
    stability.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize feed composition if fractions do not sum exactly to 1.",
    )
    stability.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Output format.",
    )
    stability.set_defaults(handler=_run_stability_tp)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        if isinstance(exc.code, int):
            return exc.code
        return 2

    try:
        return int(args.handler(args))
    except _UsageError as exc:
        parser.print_usage(sys.stderr)
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ConvergenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 3
    except (ThermoError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
