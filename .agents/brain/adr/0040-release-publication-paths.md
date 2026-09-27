# ADR-0040: Two publication paths, Trusted Publishing from a GitHub Release or an authorized manual upload, one set of gates

Status: accepted
Date: 2026-09-26

## Context
ADR-0038 item 3 made PyPI publication manual and owner-only: "No CI job
publishes and no token is stored", and the dev contract added "agents never
upload". chemthermo 0.4.0 was published that way on 2026-09-26. For the 0.4.1
campaign the owner then authorized an agent to publish to TestPyPI and PyPI
with the credentials the owner keeps in the local `~/.pypirc`, and asked for
automation that needs no long-lived PyPI token in GitHub. PyPI uploads are
immutable, so both paths need the same gates.

## Decision
1. **Automated path (preferred): `.github/workflows/release.yml`.** Publishing
   a GitHub Release (not a prerelease) for tag `vX.Y.Z` runs, in order:
   `tools/release_preflight.py`, the full CI matrix on the tagged commit,
   one build, `twine check --strict` and `tools/release_artifacts.py
   check-dist`, a smoke test of those files on Linux and macOS under
   3.11-3.13, TestPyPI, a check that TestPyPI serves the same bytes, PyPI,
   the same check against PyPI, and finally the same files attached to the
   Release. Pushing a tag publishes nothing: the GitHub Release is the explicit
   production signal, and the CI run inside the workflow removes any race with
   the push-triggered CI. `workflow_dispatch` rehearses without uploading.
2. **Trusted Publishing only.** The two upload jobs are the only jobs with
   `id-token: write`, and they run in the GitHub environments `testpypi` and
   `pypi`. No PyPI token is stored in GitHub. Third-party actions are pinned
   to full commit SHAs because a tag can be moved and this workflow publishes.
   The one-time registration on each index is an owner-side web action. On
   PyPI and TestPyPI, go to Manage `chemthermo`, then Publishing, then Add a
   GitHub publisher, and enter owner `AhmadAlkadri`, repository
   `Chemical-Thermodynamics`, workflow `release.yml`, and environment `pypi`
   (TestPyPI: `testpypi`). Until that is done, the upload jobs fail at the
   token exchange and nothing is published.
3. **Manual path (supported): `twine` with the owner's `~/.pypirc`.** It uses
   the same gates, tools and order: a clean clone of the tag, one build,
   `check-dist`, `twine check --strict`, TestPyPI, `verify-index`, PyPI,
   `verify-index`. The GitHub Release is created as a **draft** and published
   only after the PyPI upload, so the workflow finds the version already on
   PyPI. It then runs its `audit` job, which checks that PyPI and the Release
   carry the same bytes and that they install, and uploads nothing. Use one
   path per version, never both.
4. **Who may publish.** The owner may publish. An agent may prepare a release
   at any time, and it may run the TestPyPI and PyPI uploads **only under an
   explicit owner authorization for that release campaign**. The 0.4.1
   campaign has one (2026-09-26). An authorization to publish is not a
   licence to skip or reorder a gate. If a gate fails, the agent stops before
   the upload and hands over a packet.
5. **Credentials.** `twine` reads `~/.pypirc` itself. A token is never
   printed, echoed, passed on a command line, written to a file or to the
   repository, or copied elsewhere, and the owner's `~/.pypirc` is never
   edited, rotated or deleted by an agent.
6. **Build once.** The files that pass the smoke test are the files uploaded
   to TestPyPI and PyPI and attached to the Release, identified by SHA-256
   (`SHA256SUMS`). Nothing is rebuilt between those steps.

## Alternatives considered
- Publish on tag push (rejected: an irreversible upload racing CI, with no
  deliberate release step).
- A PyPI API token stored as a GitHub secret (rejected: a long-lived
  credential in a second place; OIDC gives short-lived, workflow-scoped ones).
- Automation only, retiring the manual path (rejected: the owner keeps local
  publication as a supported route, and the automated path cannot run until
  the trusted publishers are registered).
- Keep "agents never upload" (rejected by the owner for this campaign; item 4
  keeps agent publication bounded by an explicit authorization).

## Consequences
- A release needs a GitHub Release click (or `gh release create`), not a
  sequence of hand commands, once the publishers are registered.
- The release commit must be on `main`, which preflight enforces. This
  replaces ADR-0031's "main is not moved by a release" for final releases,
  which 0.4.0 had already done by fast-forward.
- Prereleases (`X.Y.ZbN`) stay GitHub-only (ADR-0031 item 3). The workflow
  skips them.

## Supersedes (optional)
Amends ADR-0038 item 3 (manual, owner-only, no CI publication) and ADR-0031
item 8 as ADR-0038 left it.

## Superseded by (optional)
None.
