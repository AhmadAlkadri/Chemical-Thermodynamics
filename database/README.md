# Component database

## Canonical runtime database path

The canonical packaged runtime component database is:

- `src/chemthermo/data/components.json`

Runtime loading uses packaged resources (`chemthermo.data/components.json`) and does not depend on any file under `database/`.

## Raw source files

Canonical raw tables for regeneration are:

- `database/organics.txt`
- `database/inorganics.txt`

## Rebuild and check workflow

Regenerate canonical payload:

```bash
python tools/build_database.py
```

Verify deterministic sync against canonical runtime path:

```bash
python tools/build_database.py --check
```

Optionally write a non-runtime mirror (generated artifact):

```bash
python tools/build_database.py --write-mirror
```

Optional mirror path:

- `database/components.mirror.json`

`database/components.mirror.json` is generated-only, git-ignored, and never loaded at runtime.

## Using the runtime database

```python
from chemthermo import Component, cite

methane = Component.from_database("Methane")
print(methane.tc_k)
print(cite("Methane", "Tc"))
```

## Adding components interactively

```bash
python tools/add_component.py
```

This tool updates:

- `src/chemthermo/data/components.json`

If you changed packaged data, reinstall before non-editable smoke tests so site-packages reflects the update.

## Schema

The database schema is defined by `src/chemthermo/schemas.py` and loaded by `src/chemthermo/data/__init__.py`.

Each component record includes:

- `name`, `formula`, `CAS`
- `MW`, `Tc`, `Pc`, `omega` (parameter objects)
- `antoine` (optional)

## Legacy / deprecated artifacts

`database/components.json` and `database.h5` are deprecated legacy artifacts and are not runtime sources of truth.
