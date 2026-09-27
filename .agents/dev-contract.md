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
  already checks, a further pressure or temperature on the same map. ADR-0028
  added 31 of these in one pass to hold the runtime budget - the table naming
  each one and what still covers it by default is in ledger Case P-17;
- the full 2505-state robustness sweep of ADR-0027 (grown from 2110 states by
  its `robustness-map-coverage` amendment), whose 224-state `--quick` subset
  runs the same code over all ten families by default.

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

## Benchmarks (ADR-0023)

Run the fixed workload and write a record:

```bash
python -m chemthermo.bench --out benchmarks/mine.json
# from a source checkout without an install:
python tools/bench.py --out benchmarks/mine.json
```

`--list` prints the case ids, `--case ID` runs one (repeatable), `--repeats N`
sets the timed repeats after the warm-up (default 5). The whole workload takes
about half a minute.

Compare two records:

```bash
python -m chemthermo.bench --compare benchmarks/baseline_<sha>.json benchmarks/after_<sha>.json
```

It prints the per-case median before, after and ratio, and **exits 1 if any
case's `result_hash` differs**.

**Any performance change is reported with a baseline record, an after record,
identical result hashes and a measured ratio**, the two records measured back
to back on one machine - the workload's three activity-model cases are the
drift control. See `benchmarks/README.md` and ADR-0023. A wall time is never
portable between machines; a result hash is.

## The robustness map (ADR-0027)

Sweep every model family over the fixed coverage grids and write a classified
record:

```bash
python -m chemthermo.bench robustness --out benchmarks/robustness_<sha>.json \
                                      --summary-out benchmarks/robustness_<sha>.md
```

2505 states, about 37 minutes. `--family NAME` runs one family (the sweep
partitions and resumes), `--quick` runs the 224-state subset the test suite
uses (~14 s), `--list` prints the grid. Each state ends in one bucket: a phase
verdict, `converged-invariant-violated`, or one of eight refusal classes.

This measures **coverage, not correctness**: nothing in it is compared against
a published number or another implementation. Use it to rank what to fix, not
to claim something is right. See `benchmarks/README.md`, ADR-0027 and ledger
Cases R-MAP-1 and R-MAP-2.

The committed record is `benchmarks/robustness_761fd57.json` (Linux,
2026-09-25), bucket-for-bucket identical to the macOS `f852726` one.
At `f852726` it refuses **5** of 2505 states, all of them the deprecated
`gamma-phi-legacy` path's `rr-no-bracket` behaviour, which is by design
(ADR-0016 decision 8). The nine `multiphase-solver-failure` states the map
found at `74820b8` in `eos-three-phase` were repaired by ADR-0029; the
original 2110-state grid still refuses nothing, unchanged field for field from
`9adf390`. The whole 2505-state record is **identical field for field to the
`64831bd` one**, wall time excluded - ADR-0030 is a performance change, and
re-running the sweep and diffing it is the widest of the three bit-identity
gates a performance slice here is held to. Read every count here as "nothing
in *this* grid refuses" (or "this is what *this* grid refuses"), not as a
coverage claim - and note that regenerating the record at a new commit means a
new file, so the superseded JSON is pruned and its `.md` summary kept.

## Releases (ADR-0031, ADR-0040)

Policy is ADR-0031 (versions, tags, gates) and ADR-0040 (who publishes and
how); this is the procedure. `X.Y.Z` below is the version being released, and
`SHA` is the full commit the gates ran on. Never `git push --tags`, never `-f`.
The only non-fast-forward updates, the 2026-09-26 attribution rewrite, are
recorded in ADR-0039. A PyPI version can never be re-uploaded, so a bad
release is fixed by the next patch version.

1. **Version commit.** Set `version` in `pyproject.toml` and `__version__` in
   `src/chemthermo/__init__.py` to `X.Y.Z`, add its `CHANGELOG.md` section,
   commit (`Slice: release`), push the branch, and get CI green on it.
2. **Gates from a clean clone** fetched from GitHub, never the working
   checkout, on macOS arm64 under CPython 3.11. That is the ADR-0032 capture
   runtime, where the captured-value guards compare with `==`:

```bash
git clone --branch <branch> https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git rel && cd rel
git rev-parse HEAD            # must equal SHA
python3.11 -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/ruff format --check src tests && .venv/bin/ruff check src tests
.venv/bin/pyright && .venv/bin/pytest -q
```

   Also run `pytest -q` under CPython 3.12 and 3.13, where those guards use
   ADR-0032's bounded comparison.
3. **`main` and the tag.** Fast-forward `main` to `SHA`, then tag it. The
   release workflow refuses a commit that is not on `main`:

```bash
git fetch origin && git merge-base --is-ancestor origin/main SHA && git push origin SHA:refs/heads/main
git tag -a vX.Y.Z SHA -m "chemthermo X.Y.Z"
git push origin refs/tags/vX.Y.Z
git ls-remote origin 'refs/tags/vX.Y.Z^{}'   # must print SHA
```

   A pushed tag never moves. Merging a PR from `dev/sprint` deletes the branch
   (`delete_branch_on_merge`), so a direct fast-forward push is preferred.
   After any merge, check that `dev/sprint` still exists.
4. **Publish by one path, never both** (ADR-0040):
   - **Automated (preferred; needs the trusted publishers of ADR-0040 item 2
     registered).** Optionally rehearse first: `gh workflow run release.yml
     -f tag=vX.Y.Z`, which runs preflight, CI, build and smoke and uploads
     nothing. Then publish the Release with no files attached (the workflow
     attaches the frozen ones):
     `gh release create vX.Y.Z --verify-tag --title "chemthermo X.Y.Z" --notes-file notes.md`.
     Watch it with `gh run watch`. The environments `testpypi` and `pypi`
     can be given required reviewers in the repository settings for a manual
     approval before each upload.
   - **Manual (`twine`, the owner's `~/.pypirc`).** From a clean clone of the
     tag, with `build==1.6.1 twine==7.0.0` in a tooling venv:

```bash
git clone --branch vX.Y.Z https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git rel-tag && cd rel-tag
OUT="$(mktemp -d)"                                                  # outside the clone
python tools/release_preflight.py --tag vX.Y.Z --sha SHA --require-on main
python -m build --outdir dist/ && twine check --strict dist/*
python tools/release_artifacts.py check-dist --version X.Y.Z dist
python tools/release_artifacts.py sums dist > $OUT/SHA256SUMS      # the frozen files
#   install dist/*.whl and dist/*.tar.gz in fresh venvs outside the tree (3.11-3.13):
#   pip check; python tools/release_smoke.py --expect-version X.Y.Z
gh release create vX.Y.Z --draft --verify-tag --title "chemthermo X.Y.Z" --notes-file notes.md dist/* $OUT/SHA256SUMS
twine upload --repository testpypi dist/*
python tools/release_artifacts.py verify-index --index testpypi --version X.Y.Z --sums $OUT/SHA256SUMS --download-dir $OUT/tp
#   install $OUT/tp/*.whl in a fresh venv (dependencies from PyPI only); pip check; release_smoke
twine upload dist/*
python tools/release_artifacts.py verify-index --index pypi --version X.Y.Z --sums $OUT/SHA256SUMS --download-dir $OUT/pypi
#   fresh venv: pip install chemthermo==X.Y.Z; pip check; release_smoke
gh release edit vX.Y.Z --draft=false     # the workflow now audits instead of uploading
```

     The owner runs this, or an agent does so under an explicit owner
     authorization for that release (ADR-0040 item 4; 0.4.1 has one).
     `twine` reads `~/.pypirc` itself. Never print, echo, pass on a command
     line, copy or commit a token, and never edit the owner's `~/.pypirc`.
     Without an authorization, an agent prepares
     `.agents/handoffs/release-packet-vX.Y.Z.md` instead.
5. **Record** on `dev/sprint` (docs-only, `Slice: handoff`), clearly after
   the tagged commit: the validation, the tag object and peel, the Release,
   TestPyPI and PyPI URLs, and the SHA-256 of both files. The tag stays where
   it is.

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
