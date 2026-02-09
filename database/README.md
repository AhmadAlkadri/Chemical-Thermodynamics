# Component database

The `flash` module looks for the component database at `chemthermo/data/components.json` by default.

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

### Rebuilding the Package Data

If you modify `database/components.json` manually or via the tool, you must assume the package uses this source directly (in development mode) or reinstall the package to see changes if installed in site-packages.

The `database.h5` and `sort_database.py` files are **deprecated** and should not be used.

### Schema

The database uses a JSON schema defined in `src/chemthermo/schemas.py`. Each component has:
- `name`, `formula`, `CAS`
- `MW`, `Tc`, `Pc`, `omega` (as Parameter objects with value, units, source_key)
- `antoine` (optional)

## Legacy Data

The original data from `organics.txt` and `inorganics.txt` has been migrated to `components.json`.

