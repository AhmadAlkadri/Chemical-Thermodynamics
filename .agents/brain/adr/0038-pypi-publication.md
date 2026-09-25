# ADR-0038: Publish to PyPI, by the owner, by hand, starting with 0.4.0

Status: accepted
Date: 2026-09-25

## Context
ADR-0031 item 8 kept releases to GitHub and required "a separate explicit
decision" for any registry. On 2026-09-25 the owner set the goal: an honest,
tidy public repository merged to `main` and a pinned version published on
PyPI, which the owner publishes personally. The PyPI names `chemthermo` and
`chemical-thermodynamics` were both unregistered that day (HTTP 404).

## Decision
1. The distribution name is **`chemthermo`** (the import name).
2. The first PyPI release is **0.4.0**, `Development Status :: 4 - Beta`, a
   final (non-prerelease) version, built from the commit on `main` that the
   `v0.4.0` tag names.
3. **The owner publishes, by hand**: build from a clean clone of the tag,
   `twine check --strict`, upload to TestPyPI and install from it, then upload
   to PyPI. No CI job publishes and no token is stored in the repository or its
   settings; adding either needs a new decision.
4. The ADR-0031 gates still apply to the tagged commit; CI additionally runs
   Python 3.11 / 3.12 / 3.13 and `twine check --strict` on every push.
5. The sdist ships the library only (no tests, whose fixtures are not
   redistributable); the wheel ships the package and its four data files.

## Alternatives considered
- Trusted publishing from GitHub Actions (rejected for now: automation that
  publishes needs its own decision; the owner chose manual).
- Publish `0.4.0b1` first (rejected by the owner in favour of 0.4.0 Beta).
- A different distribution name (rejected: `chemthermo` is free and matches
  the import name).

## Consequences
- `pip install chemthermo` becomes the documented install; git pins remain.
- Once uploaded, a version can never be re-uploaded: a bad release is fixed by
  0.4.1, mirroring ADR-0031's rule for tags.

## Supersedes (optional)
Amends ADR-0031 item 8.

## Superseded by (optional)
None.
