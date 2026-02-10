"""Interactive CLI to add a new component to the canonical runtime database."""

import json
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from chemthermo.data import normalize_name
from chemthermo.schemas import AntoineCoefficients, ComponentData, Database, Parameter

CANONICAL_DB_PATH = Path(__file__).resolve().parents[1] / "src" / "chemthermo" / "data" / "components.json"

def input_or_default(prompt: str, default: str | None = None) -> str:
    """Get input with optional default."""
    if default:
        user_input = input(f"{prompt} [{default}]: ").strip()
        return user_input if user_input else default
    return input(f"{prompt}: ").strip()

def add_component_cli():
    print("=== Add New Component to ChemThermo Database ===")

    db_path = CANONICAL_DB_PATH
    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        return

    with db_path.open("r", encoding="utf-8") as f:
        db_data = json.load(f)

    # Load into Pydantic model
    db = Database(**db_data)

    # Check for existing
    name = input_or_default("Component Name")
    canonical = normalize_name(name)
    existing_names = {normalize_name(c.name) for c in db.components}

    if canonical in existing_names:
        print(f"Warning: Component '{name}' already exists.")
        if input_or_default("Continue anyway? (y/n)", "n").lower() != "y":
            return

    formula = input_or_default("Chemical Formula")

    print("\n--- Properties ---")
    print("Enter values in SI units (kg/mol, K, Pa). Source keys refer to keys in references.bib.")

    source_key = input_or_default("Default Source Key", "manual_entry")

    mw = float(input_or_default("Molecular Weight (kg/mol)"))
    tc = float(input_or_default("Critical Temperature (K)"))
    pc = float(input_or_default("Critical Pressure (Pa)"))
    omega = float(input_or_default("Acentric Factor"))

    comp = ComponentData(
        name=name,
        formula=formula,
        MW=Parameter(value=mw, units="kg/mol", source_key=source_key),
        Tc=Parameter(value=tc, units="K", source_key=source_key),
        Pc=Parameter(value=pc, units="Pa", source_key=source_key),
        omega=Parameter(value=omega, units="-", source_key=source_key),
    )

    # Add Antoine?
    if input_or_default("\nAdd Antoine Coefficients? (y/n)", "n").lower() == "y":
        print("Model: log10(P) = A - B / (T + C)")
        units = input_or_default("Pressure Units (bar/Pa/kPa)", "bar")
        A = float(input_or_default("A"))
        B = float(input_or_default("B"))
        C = float(input_or_default("C"))
        tmin = float(input_or_default("Tmin (K)"))
        tmax = float(input_or_default("Tmax (K)"))
        
        comp.antoine = AntoineCoefficients(
            A=A, B=B, C=C, Tmin_K=tmin, Tmax_K=tmax, units=units, source_key=source_key
        )

    db.components.append(comp)

    # Sort
    db.components.sort(key=lambda c: c.name.casefold())

    # Save
    with db_path.open("w", encoding="utf-8") as f:
        f.write(db.model_dump_json(indent=2))

    print(f"\nSuccessfully added '{name}' to {db_path}")
    print("Remember to add the citation to database/references.bib if needed!")

if __name__ == "__main__":
    add_component_cli()
