# Release record: chemthermo 0.4.1 (published 2026-09-26, Pacific)

A post-release hardening patch, published by an agent with the owner's `~/.pypirc`
under the owner's explicit authorization for this campaign (ADR-0040 item 4),
by the manual path of `.agents/dev-contract.md` "Releases" step 4.

| item | value |
| --- | --- |
| release commit | `3ecbbac28257ad2771b13f23a80eac90c8ceafe5` (= `main` = `dev/sprint` at release) |
| tag | `v0.4.1`, annotated, tag object `48620e8b86f2fe0e1ae605966ce2a2b1235ffd07`, peels to `3ecbbac` |
| GitHub release | https://github.com/AhmadAlkadri/Chemical-Thermodynamics/releases/tag/v0.4.1 (created as a draft, published after PyPI) |
| TestPyPI | https://test.pypi.org/project/chemthermo/0.4.1/ |
| PyPI | https://pypi.org/project/chemthermo/0.4.1/ |
| library change since `v0.4.0` | `src/chemthermo/__init__.py` `__version__` only |
| `v0.4.0` | untouched: tag object `e730a2c`, peels to `afc43db` |

The files were built once from a clean clone of `v0.4.1` (build 1.6.1, twine 7.0.0,
`SOURCE_DATE_EPOCH` = commit time). The same bytes are on TestPyPI, PyPI and
the GitHub Release, each checked by `tools/release_artifacts.py verify-index`,
which compares the listed digests and re-hashes the downloads:

```
1178fd1744b6d97591764ed6e375b0a9861cca2993504b368e7a01fdf3f8d508  chemthermo-0.4.1-py3-none-any.whl
b5811e3f7a3ad76c6c978eeab4d473ee7304736f1134a60f311f1adcdde3b7b8  chemthermo-0.4.1.tar.gz
```

The release workflow's rehearsal build on GitHub (run 36287687647) produced a
wheel byte-identical to this one (`1178fd17...`). Its sdist differed
(`d8930b36...`), because sdists are not reproducible here. Only the files
above were published.

## Validation

| where | what | result |
| --- | --- | --- |
| macOS arm64, clean clone of `12f734a` (the code tree of `3ecbbac`; `3ecbbac` adds README badges and a CHANGELOG line) | ruff format/check, pyright | clean, 0 errors |
| same, CPython 3.11.6, numpy 2.4.6 (capture runtime, exact guards) | `pytest -q`; the 7 capture-guard files | 835 passed, 51 skipped; 235 passed |
| same, CPython 3.12.14 and 3.13.7, numpy 2.5.3 (bounded guards) | same | 835 passed, 51 skipped; 235 passed (each) |
| GitHub Actions, 3.11 / 3.12 / 3.13 | CI on `3ecbbac` | green, run 36286803364 |
| GitHub Actions, `release.yml` rehearsals (`workflow_dispatch`, nothing uploaded) | `v0.4.0`: preflight, CI, audit against PyPI; `v0.4.1`: preflight, CI, build, smoke of wheel and sdist on Ubuntu and macOS x 3.11-3.13 | green: runs 36287274404, 36287687647 |
| macOS arm64, frozen files | `twine check --strict`, `check-dist` (wheel: `chemthermo/` and dist-info only; sdist: no `tests/`, `.agents/`, `benchmarks/`, `notebooks/`, `database/`, `tools/`); wheel and sdist in fresh venvs, 3.11 / 3.12 / 3.13: `pip check`, `release_smoke.py --expect-version 0.4.1` (solver cases, citations, CLI) | all pass |
| TestPyPI | `verify-index`; its wheel installed with dependencies from PyPI only, `pip check`, smoke (3.13) | pass |
| PyPI | `verify-index`; fresh `pip install chemthermo==0.4.1` on 3.11 / 3.12 / 3.13: `pip check`, `__version__ == "0.4.1"`, smoke | pass |
| GitHub Actions, `release.yml` on the published Release | preflight, CI, `audit` (0.4.1 already on PyPI, so no upload) | green, run 36288359204 (PyPI and the Release carry the frozen bytes; the PyPI install passes the smoke) |

Not rerun: the 2505-state robustness map, `pytest -m slow` and the validation
extras. No solver or model code changed after the `761fd57` map (ADR-0031
item 6).

## Notes

- **Trailer slip.** Commit `12f734a` carries `Slice: release-0.4.1`. The dot
  breaks the AGENTS.md slug regex, and that CI run failed its cadence step
  (every test step passed). History was not rewritten. The release commit
  `3ecbbac` carries `Slice: release`. The invalid line stays visible in `git log`.
- **Trusted Publishing is not registered yet.** This is the one owner-side step
  left. On https://pypi.org/manage/project/chemthermo/settings/publishing/ and
  on https://test.pypi.org/manage/project/chemthermo/settings/publishing/, add
  a GitHub publisher with owner `AhmadAlkadri`, repository
  `Chemical-Thermodynamics`, workflow `release.yml`, and environment `pypi`
  (TestPyPI: `testpypi`). Optionally, in GitHub Settings -> Environments,
  give `pypi` and `testpypi` a required reviewer and restrict them to tags
  `v*`. Until then, a Release published for a version that is not yet on PyPI
  runs every gate and then fails at the token exchange, and nothing is
  uploaded.
- **Credentials.** `twine` read `~/.pypirc` itself. No token was printed,
  echoed, passed on a command line, written to a file or committed, and
  `~/.pypirc` was not modified (only its `[section]` names were listed).
- **Local checkout note.** The owner's working checkout has an ignored
  leftover `src/chemthermo/vlle/__pycache__/` (from 2026-09-13), which makes
  `chemthermo.vlle` importable as a namespace package there. It fails
  `test_the_removed_vlle_package_is_gone` in that checkout only. Clean
  clones and CI are unaffected. It was left in place, and deleting that
  directory fixes it.
