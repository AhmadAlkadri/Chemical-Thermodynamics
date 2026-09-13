"""Command-line interface for chemthermo."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from typing import Any

from .core import Mixture
from .exceptions import ConvergenceError, ThermoError
from .flash import FlashSettings, flash_tp
from .models import NRTL, PengRobinsonEOS

CLI_SCHEMA_VERSION = 1


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
) -> dict[str, Any]:
    phase_names = result.phase_names()

    phases: dict[str, dict[str, list[float]]] = {}
    for phase_name in phase_names:
        phase = result.phases[phase_name]
        phases[phase_name] = {
            "fractions": list(phase.composition.fractions),
        }

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
            "eos": "peng_robinson",
            "method": flash_mode,
            "algorithm": _algorithm_label(flash_mode=flash_mode, settings=settings),
            "settings": {
                "max_iter": settings.max_iter,
                "tol": settings.tol,
                "damping": settings.damping,
            },
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

    mixture = Mixture.from_database(component_names, feed_fractions, normalize=args.normalize)
    settings = FlashSettings(max_iter=args.max_iter, tol=args.tol, damping=args.damping)
    flash_mode = args.flash_mode
    activity_model = NRTL() if flash_mode == "gamma-phi" else None

    result = flash_tp(
        mixture,
        temperature_K=args.temperature_k,
        pressure_Pa=args.pressure_pa,
        eos=PengRobinsonEOS(),
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
    )

    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True))
    else:
        print(_format_tp_flash_text(payload))

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="chemthermo", description="chemthermo command-line tools")
    subparsers = parser.add_subparsers(dest="command", required=True)

    tp_flash = subparsers.add_parser(
        "tp-flash", help="Run TP flash using Peng-Robinson EOS (phi-phi or gamma-phi)."
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
    tp_flash.set_defaults(handler=_run_tp_flash)

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
    except ConvergenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 3
    except (ThermoError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
