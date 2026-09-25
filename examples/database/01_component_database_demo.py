"""Inspect the packaged runtime component database."""

from __future__ import annotations

from pathlib import Path

import chemthermo as ct


def main() -> None:
    runtime_db = Path("src/chemthermo/data/components.json")

    names = ct.list_component_names()
    methane = ct.Component.from_database("Methane")

    print("Runtime database demo")
    print(f"Canonical runtime DB path: {runtime_db}")
    print(f"Total components: {len(names)}")
    print(f"First 5 components: {', '.join(names[:5])}")
    print(f"Methane critical point: Tc={methane.tc_k:.3f} K, Pc={methane.pc_pa:.3e} Pa")


if __name__ == "__main__":
    main()
