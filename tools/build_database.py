"""Build the canonical packaged component database from raw text sources."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Iterable

PA_PER_BAR = 1e5
G_PER_KG = 1000.0
SCHEMA_VERSION = 1

CANONICAL_RUNTIME_DB_PATH = Path("src/chemthermo/data/components.json")
DEFAULT_MIRROR_PATH = Path("database/components.mirror.json")


def normalize_name(name: str) -> str:
    return " ".join(name.casefold().split())


def parse_fixed_tail(path: Path) -> list[dict[str, object]]:
    """Parse a whitespace-delimited file with two text columns then numeric columns."""

    lines = path.read_text(encoding="utf-8").splitlines()
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


def _parameter(value: float, units: str, source_key: str) -> dict[str, Any]:
    return {
        "value": value,
        "units": units,
        "source_key": source_key,
        "uncertainty": None,
        "method": None,
    }


def build_records(
    rows: Iterable[dict[str, object]],
    *,
    source_key: str,
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

        record: dict[str, object] = {
            "name": name,
            "formula": str(row["Formula"]),
            "CAS": None,
            "MW": _parameter(mw_kg_per_mol, "kg/mol", source_key),
            "Tc": _parameter(tc_k, "K", source_key),
            "Pc": _parameter(pc_pa, "Pa", source_key),
            "omega": _parameter(omega, "-", source_key),
            "antoine": {
                "A": float(row["A"]),
                "B": float(row["B"]),
                "C": float(row["C"]),
                "Tmin_K": float(row["Tmin"]),
                "Tmax_K": float(row["Tmax"]),
                "units": "bar",
                "source_key": source_key,
            },
        }
        records[canonical] = record

    return [records[key] for key in sorted(records)]


def build_payload(
    data_dir: Path,
    *,
    source_key: str,
    names: Iterable[str] | None = None,
) -> dict[str, object]:
    organics = parse_fixed_tail(data_dir / "organics.txt")
    inorganics = parse_fixed_tail(data_dir / "inorganics.txt")

    name_filter = {normalize_name(name) for name in names} if names else None
    records = build_records([*organics, *inorganics], source_key=source_key, name_filter=name_filter)

    return {
        "schema_version": SCHEMA_VERSION,
        "components": records,
    }


def write_payload(payload: dict[str, object], output_path: Path) -> Path:
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=True)
        handle.write("\n")
    return output_path


def check_canonical_sync(payload: dict[str, object], canonical_path: Path) -> bool:
    canonical_path = canonical_path.expanduser().resolve()
    if not canonical_path.exists():
        raise FileNotFoundError(f"Canonical runtime database file not found: {canonical_path}")

    current_payload = json.loads(canonical_path.read_text(encoding="utf-8"))
    return payload == current_payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the canonical packaged runtime database at "
            "src/chemthermo/data/components.json from database/organics.txt and "
            "database/inorganics.txt."
        )
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("database"),
        help="Directory containing organics.txt and inorganics.txt.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=CANONICAL_RUNTIME_DB_PATH,
        help="Output JSON path when writing payload (default: canonical runtime path).",
    )
    parser.add_argument(
        "--source-key",
        default="koretsky2012engineering",
        help="BibTeX source key used for generated parameter provenance.",
    )
    parser.add_argument(
        "--names",
        nargs="*",
        default=None,
        help="Optional component names to include (case-insensitive).",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Validate generated payload against canonical runtime path "
            "src/chemthermo/data/components.json without writing files."
        ),
    )
    parser.add_argument(
        "--write-mirror",
        action="store_true",
        help=(
            "Also write a non-runtime mirror file to database/components.mirror.json "
            "(or --mirror-output)."
        ),
    )
    parser.add_argument(
        "--mirror-output",
        type=Path,
        default=DEFAULT_MIRROR_PATH,
        help="Mirror output path used only with --write-mirror.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_payload(args.data_dir, source_key=args.source_key, names=args.names)

    if args.check:
        if check_canonical_sync(payload, CANONICAL_RUNTIME_DB_PATH):
            print(f"Canonical runtime database is up to date: {CANONICAL_RUNTIME_DB_PATH}")
            return 0

        print(
            "Canonical runtime database is out of sync. "
            "Regenerate and write src/chemthermo/data/components.json.",
            file=sys.stderr,
        )
        return 1

    written = write_payload(payload, args.output)
    print(f"Wrote canonical payload to {written}")

    if args.write_mirror:
        mirror_written = write_payload(payload, args.mirror_output)
        print(f"Wrote non-runtime mirror payload to {mirror_written}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
