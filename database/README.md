# Component database

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

The canonical runtime component database is:

- `src/chemthermo/data/components.json`

Runtime loading uses packaged resources (`chemthermo.data/components.json`) and
does not depend on `database/components.json`.

## Usage

### Using the Database

```python
from chemthermo import Component, cite

# Load a component
methane = Component.from_database("Methane")
print(f"Tc: {methane.tc_k} K")

# Get a citation
print(cite("Methane", "Tc"))
```

### Adding New Components

To add a new component to the database interactively:

```bash
python tools/add_component.py
```

This script will prompt you for properties and automatically update `database/components.json`.
This script now updates the canonical runtime file:

- `src/chemthermo/data/components.json`

### Rebuilding the Package Data

If you modify `src/chemthermo/data/components.json`, reinstall the package if
needed to refresh site-packages in non-editable installs.

`database/components.json` is a deprecated legacy artifact and is not the
runtime source of truth.

The `database.h5` and `sort_database.py` files are **deprecated** and should not be used.

### Schema

The database uses a JSON schema defined in `src/chemthermo/schemas.py`. Each component has:
- `name`, `formula`, `CAS`
- `MW`, `Tc`, `Pc`, `omega` (as Parameter objects with value, units, source_key)
- `antoine` (optional)

## Legacy Data

The original data from `organics.txt` and `inorganics.txt` has been migrated to `components.json`.
