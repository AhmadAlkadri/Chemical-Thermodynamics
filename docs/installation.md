# Installation and development setup

Installing chemthermo, the optional validation dependencies, and the component databank. For the quick version see the [README](../README.md).

# Local install

macOS/Linux (development install):
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest -q
```

Windows PowerShell (development install):
```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest -q
```

Non-dev install (wheel/non-editable):
```bash
python -m pip install .
```

Tiny import and data check:
```bash
python - <<'PY'
from chemthermo import Component

methane = Component.from_database("Methane")
print(methane.name, methane.tc_k, methane.pc_pa)
PY
```

# Common issues

- Editable install fails with an old pip: run `python -m pip install --upgrade pip` inside the active `.venv`, then retry.
- Build backend errors mentioning setuptools/wheel: run `python -m pip install --upgrade setuptools wheel`, then retry `pip install -e ".[dev]"`.
- `ModuleNotFoundError: No module named 'bibtexparser.bparser'`: your environment resolved `bibtexparser` 2.x; chemthermo's citation loader uses the 1.x API. Reinstall with a `bibtexparser<2` constraint (the pinned `pyproject.toml` requirement already enforces this for new installs).

# Optional validation dependencies

Install the reference libraries used by validation tests - `thermo` (cubic EOS
and activity models) and `teqp` (PC-SAFT, NIST):

```bash
pip install -e ".[validation]"
```

Deterministic single-case validation scripts:

```bash
python examples/validation/00_reference_case.py
python examples/validation/06_stability_vs_thermo.py
python examples/validation/13_pcsaft_vs_teqp.py
python examples/validation/14_pcsaft_flash_vs_teqp.py
```

Every validation test and script skips cleanly when its optional dependency is
missing. `examples/validation/15_flash_split_robustness.py` (phi-phi split
robustness, validation Case F-4) needs no optional dependency at all.

# Database source of truth

Canonical packaged runtime DB path:

- `src/chemthermo/data/components.json`

Regeneration/check command:

```bash
python tools/build_database.py --check
```

Optional non-runtime mirror path (generated, git-ignored):

- `database/components.mirror.json`
