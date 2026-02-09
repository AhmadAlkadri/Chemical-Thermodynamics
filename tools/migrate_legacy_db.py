"""Migrate legacy text database to new JSON schema with provenance."""

import json
import sys
from pathlib import Path

# Add src to path to import schemas
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from chemthermo.schemas import Database, ComponentData, Parameter, AntoineCoefficients

# Constants
PA_PER_BAR = 1e5
G_PER_KG = 1000.0
LEGACY_SOURCE_KEY = "legacy_db"
CP_GAS_SOURCE_KEY = "legacy_cp_db"

def parse_fixed_tail(path: Path) -> list[dict]:
    """Parse a whitespace-delimited file with two text columns then numeric columns."""
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty.")

    header = lines[0].split()
    numeric_tail = len(header) - 2
    if numeric_tail < 1:
        raise ValueError(f"{path} header must include at least one numeric column.")

    rows = []
    for idx, line in enumerate(lines[1:], start=2):
        tokens = line.split()
        if len(tokens) < numeric_tail + 2:
            continue

        formula = tokens[0]
        # Allow multi-word names
        name_tokens = tokens[1 : len(tokens) - numeric_tail]
        name = " ".join(name_tokens)
        numeric_tokens = tokens[-numeric_tail:]

        try:
            numeric = [float(value) for value in numeric_tokens]
        except ValueError:
            print(f"Skipping line {idx} in {path}: numeric parse error")
            continue

        rows.append(dict(zip(header, [formula, name, *numeric])))

    return rows

def migrate_database(data_dir: Path, output_file: Path):
    """Convert legacy files to JSON schema."""
    organics_path = data_dir / "organics.txt"
    inorganics_path = data_dir / "inorganics.txt"
    
    if not organics_path.exists() or not inorganics_path.exists():
        print("Legacy files not found.")
        return

    organics = parse_fixed_tail(organics_path)
    inorganics = parse_fixed_tail(inorganics_path)
    
    components = []
    
    for row in organics + inorganics:
        try:
            name = str(row["Name"])
            formula = str(row["Formula"])
            
            # Convert units
            mw = float(row["MW[g/mol]"])
            tc = float(row["Tc[K]"])
            pc_bar = float(row["Pc[bar]"])
            omega = float(row["omega"])
            
            # Antoine
            # Note: Legacy file header names for Antoine might vary slightly or be consistent
            # In database/organics.txt: A B C Tmin Tmax
            # In inorganics.txt: A B C Tmin Tmax
            
            antoine = None
            if "A" in row and "B" in row:
                antoine = AntoineCoefficients(
                    A=float(row["A"]),
                    B=float(row["B"]),
                    C=float(row["C"]),
                    Tmin_K=float(row["Tmin"]),
                    Tmax_K=float(row["Tmax"]),
                    units="bar",
                    source_key=LEGACY_SOURCE_KEY
                )

            comp = ComponentData(
                name=name,
                formula=formula,
                MW=Parameter(value=mw / G_PER_KG, units="kg/mol", source_key=LEGACY_SOURCE_KEY),
                Tc=Parameter(value=tc, units="K", source_key=LEGACY_SOURCE_KEY),
                Pc=Parameter(value=pc_bar * PA_PER_BAR, units="Pa", source_key=LEGACY_SOURCE_KEY),
                omega=Parameter(value=omega, units="-", source_key=LEGACY_SOURCE_KEY),
                antoine=antoine
            )
            components.append(comp)
        except Exception as e:
            print(f"Failed to migrate component {row.get('Name')}: {e}")
            continue

    # Sort checks? Standardize names?
    components.sort(key=lambda c: c.name.casefold())

    db = Database(schema_version=1, components=components)
    
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(db.model_dump_json(indent=2))
    
    print(f"Migrated {len(components)} components to {output_file}")

if __name__ == "__main__":
    db_dir = Path(__file__).resolve().parents[1] / "database"
    out_file = Path(__file__).resolve().parents[1] / "database" / "components.json"
    migrate_database(db_dir, out_file)
