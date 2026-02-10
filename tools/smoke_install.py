#!/usr/bin/env python3
"""Create an isolated venv, install chemthermo, and run install smoke checks."""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import textwrap
import venv
from pathlib import Path


def _venv_python(venv_dir: Path) -> Path:
    if sys.platform.startswith("win"):
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Install a package into a temporary virtualenv and run chemthermo smoke checks."
    )
    parser.add_argument(
        "--package",
        default=".",
        help="Package specifier to install (default: current directory).",
    )
    args = parser.parse_args()

    package = Path(args.package)
    package_spec = str(package.resolve()) if package.exists() else args.package

    with tempfile.TemporaryDirectory(prefix="chemthermo-install-smoke-") as tmp_dir:
        venv_dir = Path(tmp_dir) / ".venv"
        venv.EnvBuilder(with_pip=True).create(venv_dir)
        python = _venv_python(venv_dir)

        _run([str(python), "-m", "pip", "install", "--upgrade", "pip"])
        _run([str(python), "-m", "pip", "install", package_spec])

        smoke_code = textwrap.dedent(
            """
            from chemthermo.data import get_component_record
            from chemthermo.parameters import NRTLParameters

            record = get_component_record("Methane")
            if record.get("name") != "Methane":
                raise SystemExit("Unexpected component payload for Methane.")

            tau, alpha = NRTLParameters.load().for_components(["Methane", "Ethane"])
            if tau.shape != (2, 2) or alpha.shape != (2, 2):
                raise SystemExit("Unexpected NRTL matrix shape in smoke check.")

            print("chemthermo install smoke passed")
            """
        )
        _run([str(python), "-c", smoke_code])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
