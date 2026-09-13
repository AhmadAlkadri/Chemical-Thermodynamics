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

## The `slow` marker

`pyproject.toml` registers `slow` and sets `addopts = "-m 'not slow'"`, so a
plain `pytest -q` (what CI runs) **deselects** slow-marked tests rather than
running them - currently just the full 188-state PC-SAFT validation grid,
`tests/validation/test_flash_split_robustness_pcsaft.py::test_the_whole_grid_answers_and_every_answer_is_verified`
(ADR-0017). A representative 16-state subset of the same grid
(`tests/validation/test_flash_split_robustness_pcsaft_subset.py`) runs in the
default `pytest -q`, so CI still exercises the PC-SAFT phi-phi path on every
run; the exhaustive grid is opt-in:

```bash
pytest -q -m slow
```

Command-line `-m` overrides `addopts`' `-m` (standard pytest behavior: the
last `-m` value wins), so this runs exactly the slow-marked tests and nothing
else. CI does not run it automatically - the full grid does not fit the
suite's runtime budget alongside everything else, which is why the trim
above exists; see validation Case F-5 in `.agents/brain/validation-cases.md`
for the measured before/after. Mark a new test `@pytest.mark.slow` only for a
full/exhaustive grid that already has a cheaper representative subset or
golden-path example covering the same code by default.

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
