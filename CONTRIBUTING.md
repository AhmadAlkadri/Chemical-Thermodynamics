# Contributing

Thanks for looking. Bug reports with a reproducing script (components,
composition, temperature, pressure, model) are the most useful thing you can
send.

## Development setup

```bash
python3.11 -m venv .venv
.venv/bin/pip install -e ".[dev]"          # ruff, pyright, pytest, build
.venv/bin/pip install -e ".[validation]"   # optional: teqp, FeOs, thermo cross-checks
```

## Checks (the same ones CI runs)

```bash
.venv/bin/ruff format --check src tests && .venv/bin/ruff check src tests
.venv/bin/pyright
.venv/bin/pytest -q                  # default suite; `-m slow` runs the exhaustive grids
.venv/bin/python tools/smoke_install.py --package .
```

Tests under `tests/validation/` compare against teqp, FeOs and thermo and skip
cleanly when those are not installed.

## How changes are made here

- Hard gates for every change (clean tree, `Slice: <slug>` commit trailers):
  [`AGENTS.md`](AGENTS.md).
- Commands and conventions: [`.agents/dev-contract.md`](.agents/dev-contract.md).
- Design decisions are recorded as ADRs in [`.agents/brain/adr/`](.agents/brain/adr/);
  every validated number has an entry in the ledger
  [`.agents/brain/validation-cases.md`](.agents/brain/validation-cases.md).
  A change to public API or architecture adds or updates an ADR.
- Tolerances are never loosened and tests never deleted to make a change pass.
