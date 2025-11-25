"""Utilities to build the component property HDF5 database.

This script parses the raw text files shipped in ``database/`` and writes a
consolidated ``database.h5`` that contains:

- ``organics``: critical and Antoine data for organic compounds
- ``inorganics``: critical and Antoine data for inorganic compounds
- ``all``: combined organic + inorganic data with Antoine constants grouped
- ``/Cp/gases``: ideal-gas heat capacity coefficients

Run directly to rebuild the HDF5 store:

    python sort_database.py --data-dir . --output database.h5
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


DATA_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = DATA_DIR / "database.h5"


def _parse_fixed_tail(path: Path) -> pd.DataFrame:
    """Parse a whitespace-delimited file with two text columns then numeric columns.

    The first line supplies the column headers. The first two columns are treated
    as strings (``Formula`` and ``Name``), and all remaining columns are converted
    to floats. Any additional words between the ``Formula`` token and the numeric
    tail are folded into the ``Name`` field, allowing multi-word component names.
    """

    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty.")

    header = lines[0].split()
    numeric_tail = len(header) - 2
    if numeric_tail < 1:
        raise ValueError(f"{path} header must include at least one numeric column.")

    rows: list[list[object]] = []
    for idx, line in enumerate(lines[1:], start=2):
        tokens = line.split()
        if len(tokens) < numeric_tail + 2:
            raise ValueError(
                f"Line {idx} in {path} has too few fields (expected at least "
                f"{numeric_tail + 2}, got {len(tokens)})."
            )

        formula = tokens[0]
        name_tokens = tokens[1 : len(tokens) - numeric_tail]
        name = " ".join(name_tokens)
        numeric_tokens = tokens[-numeric_tail:]
        try:
            numeric = [float(value) for value in numeric_tokens]
        except ValueError as exc:
            raise ValueError(
                f"Failed to parse numeric values on line {idx} in {path}."
            ) from exc

        rows.append([formula, name, *numeric])

    return pd.DataFrame(rows, columns=header)


def _build_antoine_column(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    """Group Antoine coefficients into a single list-valued column."""

    cols = list(columns)
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"Antoine columns missing from dataframe: {missing}")
    return df[cols].apply(lambda row: [float(row[c]) for c in cols], axis=1)


def build_database(
    organics_path: Path,
    inorganics_path: Path,
    cp_gas_path: Path,
    output_path: Path = DEFAULT_OUTPUT,
) -> Path:
    """Read source text files and write the consolidated HDF5 database."""

    organics = _parse_fixed_tail(organics_path)
    inorganics = _parse_fixed_tail(inorganics_path)
    cp_gases = _parse_fixed_tail(cp_gas_path)

    comps = pd.concat([organics, inorganics], ignore_index=True)

    keep_cols = ["Formula", "Name", "MW[g/mol]", "Tc[K]", "Pc[bar]", "omega"]
    ant_cols = ["A", "B", "C", "Tmin", "Tmax"]

    ant_params = _build_antoine_column(comps, ant_cols).rename("Antoine Constants")
    all_components = pd.concat([comps[keep_cols], ant_params], axis=1)

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.HDFStore(output_path) as store:
        store["organics"] = organics
        store["inorganics"] = inorganics
        store["all"] = all_components
        store["/Cp/gases"] = cp_gases

    return output_path


def parse_args() -> argparse.Namespace:
    """Define and parse command-line arguments."""

    parser = argparse.ArgumentParser(description="Rebuild database.h5 from raw text inputs.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DATA_DIR,
        help="Directory containing organics.txt, inorganics.txt, and Cp_gas.txt.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"Output HDF5 path (default: {DEFAULT_OUTPUT}).",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for CLI usage."""

    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    output_path = (args.output or DEFAULT_OUTPUT).expanduser()

    organics_path = data_dir / "organics.txt"
    inorganics_path = data_dir / "inorganics.txt"
    cp_gas_path = data_dir / "Cp_gas.txt"

    for required in (organics_path, inorganics_path, cp_gas_path):
        if not required.exists():
            raise FileNotFoundError(f"Required input not found: {required}")

    built_path = build_database(
        organics_path=organics_path,
        inorganics_path=inorganics_path,
        cp_gas_path=cp_gas_path,
        output_path=output_path,
    )
    print(f"Wrote {built_path}")


if __name__ == "__main__":
    main()
