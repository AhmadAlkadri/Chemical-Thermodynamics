# Dev Contract (CI Parity)

This file is the single source of truth for contributor and agent commands.

## Environment bootstrap

From the repo root:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Optional external-validation extras:

```bash
python -m pip install -e ".[validation]"
```

## CI gates

Run the same checks as `.github/workflows/ci.yml`:

```bash
ruff format --check src tests
ruff check src tests
pyright
pytest
```

## Golden path

Run the canonical demo from the repo root:

```bash
python examples/basic/flash_tp_peng_robinson_demo.py
```

## Cheap check subset

Run these quick checks before the full suite when iterating:

```bash
python -m pytest tests/test_import.py -q
python -m pytest tests/test_validation.py -q
python -m pytest tests/test_flash_tp.py -q
python -m pytest tests/test_examples.py -q
```
