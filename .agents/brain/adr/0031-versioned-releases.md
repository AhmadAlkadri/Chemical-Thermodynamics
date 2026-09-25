# ADR-0031: Versioned releases pin a tested commit; development stays in slices

Status: accepted
Date: 2026-09-25

## Context
Until this ADR the repository had one tag, `v0.1.0`, on the old `main`
(`a7a8ca7`), whose `pyproject.toml` and `chemthermo.__version__` both say
`0.0.0`, and no GitHub release. The 121-commit equilibrium campaign on
`dev/sprint` was published for the first time on 2026-09-25 (`5041dd7`). Users
and agents need fixed, reproducible milestones rather than a moving branch.
CI (`.github/workflows/ci.yml`) runs tests on every push and publishes nothing.

## Decision
1. **Slices, then milestones.** Work continues in thin vertical slices with
   `Slice:` trailers (AGENTS.md). A *release* is cut only at a validated
   capability milestone, never per slice and never retroactively.
2. **Version scheme.** PEP 440 versions, tag `v<version>` (existing `v0.1.0`
   convention). Pre-1.0: a new capability milestone bumps the minor
   (`0.2.0`, `0.3.0`, ...); a fix-only release bumps the patch. A breaking
   change to a public contract (ADR-0001: `chemthermo.__all__`, the CLI JSON
   schema and exit codes of ADR-0003/0004, packaged-data schema versions) is
   allowed pre-1.0 only with a minor bump and a migration note. `1.0` is a
   separate, explicit decision about contract stability, not a commit count.
3. **Prereleases.** A milestone is first published as `X.Y.0bN` (GitHub
   *prerelease*) when its validation is from a single environment or its
   public contract is new; it is promoted to `X.Y.0` (same content, new tag on
   a commit that only changes the version string) once a second environment
   (GitHub Actions, or a cloud session) reproduces the gates.
4. **One version, three places.** `pyproject.toml` `version`,
   `chemthermo.__version__` (pinned together by
   `test_version_matches_installed_metadata` in `tests/test_import.py`) and the tag
   must agree on the tagged commit, and `CHANGELOG.md` has that version's
   section. Between releases the version stays at the last released value.
5. **A tag names one tested commit, forever.** Annotated tag on the exact
   commit the gates ran on; a published tag is never moved, reused or deleted.
   A bad release is superseded by the next version, not rewritten.
6. **Release gates** (all on the tagged commit, from a clean clone fetched from
   GitHub, not the working checkout): `ruff format --check`, `ruff check`,
   `pyright`, `pytest -q`; `python -m build` producing sdist and wheel; the
   wheel installed into a fresh venv *outside* the source tree with the
   packaged data (`components.json`, `references.bib`, `nrtl.json`,
   `pcsaft.json`) loading and a stability/flash smoke passing. The 37-minute
   robustness sweep and `pytest -m slow` are required only when the release
   contains solver/model changes not already covered by a committed record
   (`benchmarks/robustness_<sha>.*`) at or after the last such change.
7. **Release notes** state capabilities, the full commit SHA, the validation
   actually run (environment, commands, results) versus inherited evidence,
   known limitations and compatibility implications, install pins by tag and
   by SHA, and SHA-256 checksums of attached artifacts.
8. **No registry publication.** Releases are GitHub source releases with the
   built sdist/wheel attached. Publishing to PyPI or any registry, or adding
   automation that does, needs a separate explicit decision.
9. **CI's `aglint` step is best-effort.** `agentslint` lives in a private
   repository that the public workflow cannot clone without a credential, which
   made every CI run fail at install before any test ran. The step now runs
   when the tool installs and otherwise emits a warning; every scientific,
   lint, type, install and hygiene gate stays hard. Run `aglint check --repo .`
   locally (where the tool is available) before a release.

## Alternatives considered
- Tag every accepted slice (rejected: dozens of versions nobody can tell apart).
- Call the published checkpoint `1.0` (rejected: CLI does not expose the
  equilibrium work yet, validation is one machine, and ADR-0005's "stable"
  is bounded by a trial set).
- Move or re-cut `v0.1.0` to match its `0.0.0` metadata (rejected: never move
  a published tag; the mismatch is recorded in `CHANGELOG.md` instead).
- Store a CI secret to reach `agentslint` (rejected: credentials and access
  controls are outside this change's authority).

## Consequences
- Every release is reproducible from a tag or an immutable SHA.
- Development continues on `dev/sprint` (or a session branch based on it);
  `main` is not moved by a release.
- The aglint gate is weaker in public CI than before; it is recorded here, not
  silent.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
