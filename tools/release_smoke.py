#!/usr/bin/env python3
"""Smoke-test an *installed* chemthermo: packaged data, stability, 1/2/3-phase flash.

Run it with the interpreter of the environment under test, from a directory
outside the source tree, so the import cannot resolve to ``src/``::

    cd /tmp && /path/to/venv/bin/python /path/to/tools/release_smoke.py

It checks phase counts, stability verdicts and mass-balance residuals of four
fixed states (ADR-0031 release gate). It is a smoke test, not validation: the
numbers are checked against the published solver's own invariants, not against
a reference. Pass ``--expect-version X`` to also pin ``chemthermo.__version__``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ATMOSPHERE_PA = 101325.0


def main() -> int:
    parser = argparse.ArgumentParser(description="Smoke-test an installed chemthermo.")
    parser.add_argument("--expect-version", default=None)
    args = parser.parse_args()

    import chemthermo as ct
    from chemthermo.eos import PCSAFTEOS
    from chemthermo.parameters import NRTLParameters, PCSAFTParameters

    location = Path(ct.__file__).resolve()
    print(f"chemthermo {ct.__version__} from {location}")
    if "src" in location.parts and (location.parents[2] / "pyproject.toml").exists():
        print("FAIL: imported from a source tree, not an installed distribution")
        return 1
    if args.expect_version is not None and ct.__version__ != args.expect_version:
        print(f"FAIL: version {ct.__version__} != expected {args.expect_version}")
        return 1

    # Packaged runtime data: databank, bibliography, NRTL and PC-SAFT sets.
    ct.Component.from_database("Methane")
    NRTLParameters.load().for_components(["Methane", "Ethane"])
    PCSAFTParameters.load().for_components(["Methane", "n-Hexane"])

    cases = [
        # label, components, feed, T [K], P [Pa], eos, stability, phases
        ("PR single phase", ["Methane", "Ethane"], [0.5, 0.5], 300.0, 1.0e5,
         ct.PengRobinsonEOS(), "stable", ["vapor"]),
        ("PR VLE", ["Methane", "n-Hexane"], [0.5, 0.5], 300.0, 2.0e6,
         ct.PengRobinsonEOS(), "unstable", ["liquid", "vapor"]),
        ("PR three liquids", ["Water", "Ethanol", "n-Hexane"], [0.2, 0.4, 0.4], 280.0,
         ATMOSPHERE_PA, ct.PengRobinsonEOS(), "unstable", ["liquid1", "liquid2", "liquid3"]),
        ("PC-SAFT VLE", ["Methane", "n-Hexane"], [0.5, 0.5], 300.0, 2.0e6,
         PCSAFTEOS(), "unstable", ["liquid", "vapor"]),
    ]  # fmt: skip
    failures = 0
    for label, names, feed, temperature, pressure, eos, status, phases in cases:
        mixture = ct.Mixture.from_database(names, feed, normalize=True)
        stability = ct.stability_tp(
            mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos
        )
        result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
        balance = result.diagnostics.get("mass_balance_residual")
        ok = (
            stability.status == status
            and sorted(result.phases) == phases
            and abs(sum(result.phase_fractions.values()) - 1.0) < 1e-10
            and (balance is None or float(balance) < 1e-10)
        )
        failures += not ok
        print(
            f"{'ok  ' if ok else 'FAIL'} {label:<17} stability={stability.status:<9} "
            f"phases={sorted(result.phases)} mass_balance={balance}"
        )
    if failures:
        print(f"FAIL: {failures} case(s)")
        return 1
    print("chemthermo release smoke passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
