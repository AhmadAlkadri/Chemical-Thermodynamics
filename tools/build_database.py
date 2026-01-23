"""Build the runtime JSON component database from raw text sources."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

PA_PER_BAR = 1e5
G_PER_KG = 1000.0
SCHEMA_VERSION = 1


def normalize_name(name: str) -> str:
    return " ".join(name.casefold().split())


def parse_fixed_tail(path: Path) -> list[dict[str, object]]:
    """Parse a whitespace-delimited file with two text columns then numeric columns.

    The first line provides column headers. The first two columns are treated as
    strings (Formula, Name), and the trailing columns are parsed as floats.
    """

    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty.")

    header = lines[0].split()
    numeric_tail = len(header) - 2
    if numeric_tail < 1:
        raise ValueError(f"{path} header must include at least one numeric column.")

    rows: list[dict[str, object]] = []
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
            raise ValueError(f"Failed to parse numeric values on line {idx} in {path}.") from exc

        rows.append(dict(zip(header, [formula, name, *numeric])))

    return rows


def build_records(
    rows: Iterable[dict[str, object]],
    *,
    name_filter: set[str] | None = None,
) -> list[dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for row in rows:
        name = str(row["Name"])
        canonical = normalize_name(name)
        if name_filter is not None and canonical not in name_filter:
            continue
        if canonical in records:
            raise ValueError(f"Duplicate component name after normalization: {name!r}")

        mw_kg_per_mol = float(row["MW[g/mol]"]) / G_PER_KG
        tc_k = float(row["Tc[K]"])
        pc_pa = float(row["Pc[bar]"]) * PA_PER_BAR
        omega = float(row["omega"])

        record = {
            "name": name,
            "formula": str(row["Formula"]),
            "properties": {
                "MW_kg_per_mol": mw_kg_per_mol,
                "Tc_K": tc_k,
                "Pc_Pa": pc_pa,
                "omega": omega,
            },
            "antoine": {
                "A": float(row["A"]),
                "B": float(row["B"]),
                "C": float(row["C"]),
                "Tmin_K": float(row["Tmin"]),
                "Tmax_K": float(row["Tmax"]),
            },
        }
        records[canonical] = record

    return [records[key] for key in sorted(records)]


def build_database(
    data_dir: Path,
    output_path: Path,
    *,
    names: Iterable[str] | None = None,
) -> Path:
    organics = parse_fixed_tail(data_dir / "organics.txt")
    inorganics = parse_fixed_tail(data_dir / "inorganics.txt")

    name_filter = {normalize_name(name) for name in names} if names else None
    records = build_records([*organics, *inorganics], name_filter=name_filter)

    payload = {
        "schema_version": SCHEMA_VERSION,
        "components": records,
    }

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=True)
        handle.write("\n")

    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build components.json from raw text sources.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("database"),
        help="Directory containing organics.txt and inorganics.txt.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("src/chemthermo/data/components.json"),
        help="Output JSON path.",
    )
    parser.add_argument(
        "--names",
        nargs="*",
        default=None,
        help="Optional component names to include (case-insensitive).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_database(args.data_dir, args.output, names=args.names)


if __name__ == "__main__":
    main()
