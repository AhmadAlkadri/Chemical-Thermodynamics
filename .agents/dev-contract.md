# Dev Contract (CI Parity)

This file is the single source of truth for contributor and agent commands.
Hard-gate policy lives in `AGENTS.md`.

## Environment bootstrap

From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
.venv/bin/pip install -e ".[dev]"
```

Optional external-validation extras:

```bash
.venv/bin/pip install -e ".[validation]"
```

## CI gates

Run the same checks as `.github/workflows/ci.yml`:

```bash
ruff format --check src tests
ruff check src tests
pyright
pytest -q
```

## Installability smoke checks

Run both editable and non-editable install checks:

```bash
python -c "import chemthermo; print(getattr(chemthermo, '__version__', 'no __version__'))"
python tools/smoke_install.py --package .
```

## Golden path

Run the canonical demo from the repo root:

```bash
python examples/basic/flash_tp_peng_robinson_demo.py
```

## Slice evidence

For every thin vertical slice report, include:

- Slice declaration in this format: "After this change, user can X by running Y."
- Golden path command(s) and outcome(s).
- Focused test command(s) for touched behavior and outcome(s).
- CI-equivalent gate command(s) and outcome(s).
- Installability smoke command(s) and outcome(s).

## Proof of Compliance

Policy source: `AGENTS.md`.

Run:

```bash
git status --porcelain
git log --oneline -n 8
git log -n 20 --format=%B | rg "^Slice:"
```

Expected:

- `git status --porcelain` prints nothing.
- `git log --oneline -n 8` shows recent commits.
- `git log -n 20 --format=%B | rg "^Slice:"` shows slice trailers when slices were used.

## Cheap Checks

Run these quick checks before the full suite when iterating:

```bash
python -m pytest tests/test_import.py -q
python -m pytest tests/test_validation.py -q
python -m pytest tests/test_flash_tp.py -q
python -m pytest tests/test_examples.py -q
```
