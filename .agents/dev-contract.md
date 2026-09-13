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

That installs `thermo` (Peng-Robinson / NRTL cross-checks), `teqp`
(non-associating PC-SAFT, NIST) and, since ADR-0018, `feos` (PC-SAFT
**association**, which teqp's `PCSAFT` kind does not implement). `feos` brings
its own units package along: the PyPI name is `si-units` and the import name is
`si_units`, so it is not listed separately in `pyproject.toml`. Every test and
script that needs one of these skips cleanly when it is absent
(`pytest.importorskip` in `tests/validation/`, a printed message and exit 0 in
`examples/validation/`).

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
running them. What is marked:

- the full 188-state PC-SAFT validation grid (ADR-0017), with a representative
  16-state subset (`tests/validation/test_flash_split_robustness_pcsaft_subset.py`)
  in the default run;
- the exhaustive three-phase checks of ADR-0020 - the 41-point temperature
  scans across the water / n-hexane three-phase temperature, the PC-SAFT
  ternary tie triangle, the verdict-boundary bisection;
- a handful of *repetitions*: a further feed on a tie line the default run
  already checks, a further pressure or temperature on the same map.

```bash
pytest -q -m slow
```

Command-line `-m` overrides `addopts`' `-m` (standard pytest behavior: the
last `-m` value wins), so this runs exactly the slow-marked tests and nothing
else. CI does not run it automatically - the marked tests do not fit the
suite's runtime budget alongside everything else; see validation Cases F-5 and
P-9 in `.agents/brain/validation-cases.md` for the measured before/after.

**When to mark a new test `slow`** (widened by ADR-0020): only when the default
run still covers the same capability. That means one of

- a full/exhaustive grid whose representative subset runs by default, or
- a *repetition* of something the default run does - another feed on the same
  tie line, another temperature on the same map, the other side of a
  two-sided check.

Never for the only test of a capability. Put the reason in a comment next to
the marker, naming what covers it by default.

Heavy **examples** are trimmed the other way, not marked: an example that costs
more than a few seconds takes a `--full` flag and defaults to a cheaper subset,
as `examples/validation/15_flash_split_robustness.py` has since ADR-0017 and
six more examples do since ADR-0020. `tests/test_examples.py` runs every
example with its default, so the smoke test stays cheap and the full run is one
flag away.

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

Note that CI invokes the `pytest` **console script**, not `python -m pytest`,
and the two differ in one way that has bitten a test here before: `python -m
pytest` puts the working directory on `sys.path` and the console script does
not. A test that needs to read something out of a sibling test module must
load it by path (`importlib.util.spec_from_file_location`), not with
`from tests.validation... import ...`.
