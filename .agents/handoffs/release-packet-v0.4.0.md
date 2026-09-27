> **Published 2026-09-26.** Tag `v0.4.0` = tag object
> `e730a2ceb4056dc8c7812ce90986043227ad419e`, peeling to the release commit
> `afc43dbf9862d5548a0bcb4d1c26167e7feaa66a` (= `main` at release; PR #1 merged
> by fast-forward). GitHub release
> https://github.com/AhmadAlkadri/Chemical-Thermodynamics/releases/tag/v0.4.0,
> TestPyPI https://test.pypi.org/project/chemthermo/0.4.0/, PyPI
> https://pypi.org/project/chemthermo/0.4.0/, all carrying the same two files,
> built once from a clean clone of the tag:
>
> ```
> 09e26d38cf9802d0c51e7814241d73bf7e0388c085e0efc62ca5b65fc3940fd2  chemthermo-0.4.0-py3-none-any.whl
> d301b30c0a3dacf2e022dd1e144235fccfdea69d18cea0cca00f65322e203b0f  chemthermo-0.4.0.tar.gz
> ```
>
> Fresh validation on `afc43db` (the GitHub release notes have the table): CI
> green on 3.11 / 3.12 / 3.13 (runs 36280528513, 36280529393); macOS arm64
> CPython 3.11.6: `pytest -q` 824 passed / 0 failed, exact-mode guards 199
> passed, validation extras 1021 passed; the wheel from the tag, from TestPyPI
> and from PyPI passed `pip check`, `release_smoke.py`, CLI and citation checks.
> **Open finding, released with the owner's agreement:** on macOS arm64 with
> CPython 3.12 / 3.13, `test_flash_tp_is_bit_identical_to_the_v3_capture` fails
> exact comparison for `gamma-gamma-tessier2000-near-plait` (last-bit moves that
> follow the CPython version, not numpy; passes under ADR-0032's off-platform
> bound). It predates the release and needs a test-only fix in a later slice.
> This commit is a record after the tag; it also removes a stale copy of the
> 2026-09-25 steps 5-8 and licensing note that the 2026-09-26 revision left in
> this file by mistake (the copy at the tagged commit still has it; `.agents/`
> is not in the sdist or wheel).

# Release packet: chemthermo 0.4.0 - merge to main, tag, GitHub release, PyPI

Prepared 2026-09-25 by the cloud session; updated 2026-09-26 by the owner-side
session after the attribution rewrite (ADR-0039,
`.agents/reports/history-rewrite-2026-09-26.md`). Policy: ADR-0031 (tags,
releases), ADR-0038 (PyPI by hand), ADR-0039 (the rewrite).

| item | value |
| --- | --- |
| release commit | the tip of `dev/sprint` that carries this file's 2026-09-26 revision; `v0.4.0` names it, and the tag, the PR and the GitHub release pin its full SHA (a file cannot name the commit that adds it) |
| its library content | identical to `f2165cd` (was `4aa4753` before the rewrite): the commits after it touch `.agents/` only, which neither the sdist (`MANIFEST.in` prunes it) nor the wheel ship |
| version in the code | `0.4.0` (`pyproject.toml`, `chemthermo.__version__`) |
| PR | https://github.com/AhmadAlkadri/Chemical-Thermodynamics/pull/1 (`dev/sprint` -> `main`) |
| `main` before the merge (as prepared) | `a7a8ca7` (`v0.1.0`, untouched by the rewrite), an ancestor of `dev/sprint`: a fast-forward, no merge commit. After the release `main` = `afc43db`, then the record commits |
| PyPI / TestPyPI name | `chemthermo` (no project on either index on 2026-09-26) |

## Steps

Stop at the first failure and report it. Never move `v0.4.0` once pushed and
never upload a version twice (PyPI refuses re-uploads; a fix is 0.4.1). The
only non-fast-forward ref updates this repository has had are the ones ADR-0039
records; a fast-forward of `main` needs no force. Below, `SHA` is the full id
of the release commit.

1. **CI green on PR #1** (Python 3.11 / 3.12 / 3.13) on the PR head `SHA`.
2. **Gates on a clean clone of `SHA`**, macOS arm64 (the ADR-0032 capture
   platform, where the exact-mode guards run as `==`; ADR-0035/0036 were proven
   dormant bit for bit on Linux only):
   ```bash
   git clone https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git ct-040 && cd ct-040
   git checkout SHA
   python3.11 -m venv .venv && .venv/bin/pip install -e ".[dev]"
   .venv/bin/ruff format --check src tests && .venv/bin/ruff check src tests && .venv/bin/pyright
   .venv/bin/pytest -q
   .venv/bin/pytest -q tests/test_flash_refactor_bit_identity.py tests/test_pcsaft_association.py \
     tests/test_stability_eos_surfaces.py tests/test_pcsaft_polymer.py tests/test_stability_candidates.py
   ```
   Plus `pytest -q` under Python 3.12 and 3.13, and `.[validation]` with
   `tests/validation/`. If an exact-mode guard fails on macOS, stop.
3. **Fast-forward `main` to `SHA`** (the PR then shows as merged):
   ```bash
   git fetch origin
   git merge-base --is-ancestor origin/main SHA && test "$(git rev-parse origin/dev/sprint)" = SHA && \
   git push origin SHA:refs/heads/main
   ```
   Do not squash or rebase-merge (both rewrite the commits).
4. **Tag:**
   ```bash
   git tag -a v0.4.0 SHA -m "chemthermo 0.4.0"
   git push origin refs/tags/v0.4.0
   git ls-remote origin 'refs/tags/v0.4.0^{}'   # must print SHA
   ```
5. **Build once from the tag and freeze the files**: `git clone --branch v0.4.0 ...`,
   `python -m build --outdir dist/`, `twine check --strict dist/*`, inspect the
   contents, install the wheel in a venv outside the tree and run
   `tools/release_smoke.py --expect-version 0.4.0`, `chemthermo --help`, and the
   citation check below; `shasum -a 256 dist/*`. The same two files go to
   TestPyPI, the GitHub release and PyPI; never rebuild between them.
   ```python
   from chemthermo import Component, cite
   water = Component.from_database("Water")
   print(water.get_citation("Tc")); print(cite("Water", "antoine"))
   ```
6. **TestPyPI** (TestPyPI account and API token; `twine` reads `~/.pypirc` or
   `TWINE_USERNAME=__token__` / `TWINE_PASSWORD`, never a command-line token):
   ```bash
   twine upload --repository testpypi dist/*
   python3.11 -m venv /tmp/tp && PIP_CONFIG_FILE=/dev/null /tmp/tp/bin/pip install \
     "numpy>=1.24" "pydantic>=2.0" "bibtexparser>=1.4.0,<2"          # dependencies from PyPI only
   PIP_CONFIG_FILE=/dev/null /tmp/tp/bin/pip download --no-deps --index-url https://test.pypi.org/simple/ \
     -d /tmp/tp-dl chemthermo==0.4.0                                   # the package from TestPyPI only
   shasum -a 256 /tmp/tp-dl/*        # must equal the frozen wheel's digest
   /tmp/tp/bin/pip install --no-deps /tmp/tp-dl/chemthermo-0.4.0-py3-none-any.whl && /tmp/tp/bin/pip check
   cd /tmp && /tmp/tp/bin/python <rel-040>/tools/release_smoke.py --expect-version 0.4.0 && /tmp/tp/bin/chemthermo --help
   ```
   Check the rendered project page on test.pypi.org (README, links, classifiers).
7. **GitHub release** (not a prerelease) with the frozen files:
   `gh release create v0.4.0 --verify-tag --title "chemthermo 0.4.0" --notes-file NOTES.md dist/*`,
   `NOTES.md` = the section below with the fresh validation filled in.
8. **PyPI:** `twine upload dist/*`; then a fresh venv, `pip install chemthermo==0.4.0`,
   compare the digests PyPI reports with the frozen ones, and rerun step 5's smoke checks.
9. **Record it** on `dev/sprint` (docs-only, `Slice: handoff`), clearly after
   the tagged commit: validation, tag peel, release URL, PyPI URL, SHA-256.
   `v0.4.0` stays where it is.

**Packaged-data citation: decided by the owner.** The packaged component table
(82 components: critical constants, acentric factors, Antoine coefficients)
cites Koretsky, *Engineering and Chemical Thermodynamics* (Wiley)
(`src/chemthermo/data/references.bib`). The 2026-09-25 packet left
redistribution on PyPI as a call for the owner. On 2026-09-26 the owner
decided to release 0.4.0 with the existing references and data unchanged.
That is the owner's release decision; it is not a legal opinion and no new
source verification was done for it. The PC-SAFT parameters are published
values cross-checked against FeOs / Clapeyron.jl; the two NRTL pairs are
synthetic; the DECHEMA-derived Tessier fixtures are test-only and are **not**
in the sdist or wheel.

---

## chemthermo 0.4.0

**Release commit:** the full SHA `v0.4.0` names (filled in on the GitHub
release). Its library content is that of `f2165cd` (`4aa4753` before the
2026-09-26 attribution rewrite, ADR-0039, which changed no file). First release
on PyPI: `pip install chemthermo`. Beta: the science is checked against
independent implementations; the public API may still change before 1.0.

### What it is
Phase equilibrium for chemical engineering in Python: tangent-plane stability,
multiphase TP flash (1-3 phases, phase count decided by the stability test),
Peng-Robinson, PC-SAFT (2B association, polymers, residual H/S/U/G) and NRTL,
SI units, a `chemthermo` command line. See the README for the capability
table and the limits.

### Since 0.3.0b1
- **Breaking:** `chemthermo.vlle` removed (deprecated since ADR-0013).
- PC-SAFT residual properties now include association (water, alcohols).
- `stability_tp` reports a `tpd` tie from the better-converged trial.
- `flash_tp` no longer accepts the trivial split as converged: 29 of 135
  neighbouring pressures of a polymer ladder state refused before, none now.
- README rewritten as an honest summary; reference material in `docs/`.
- Packaging for PyPI; Python 3.11-3.13.

### Validation for this release
Rows marked *inherited* ran on the original commit named, whose tree the
rewritten commit in parentheses has unchanged; they did not run on the new id.
The fresh rows on the release commit itself are in the header of this file and
the GitHub release notes.

| where | what | result |
| --- | --- | --- |
| *inherited:* clean clone of `4aa4753` (`f2165cd`), Linux x86_64, CPython 3.11.15, numpy 2.4.6 | ruff, pyright, `pytest -q` | clean, 0 errors, **824 passed**, 51 skipped (extras absent) |
| *inherited:* same | `python -m build`, `twine check --strict` | both artifacts PASSED |
| *inherited:* same, wheel in fresh venvs outside the tree, Python 3.11 / 3.12 / 3.13 | `release_smoke.py --expect-version 0.4.0`, CLI `stability-tp --eos pc-saft`, `chemthermo.vlle` absent | pass on all three |
| *inherited:* Linux, Python 3.12.3 and 3.13.12 (numpy 2.5.3) | `pytest -q` at `bc92c5b` (`f9c5f3d`; library code as released) | 824 passed each |
| *inherited:* Linux, validation extras (teqp 0.23.2, FeOs 0.10.1, thermo 0.6.1) | `pytest -q` incl. `tests/validation/` at `4aa4753` (`f2165cd`) | **1021 passed**, 0 failed, 99 deselected |
| *inherited:* Linux, full robustness map at `761fd57` (`41cd3fd`, after the last solver change) | 2505 states | 2500 converge, 5 by-design gamma-phi refusals, 0 invariant violations; identical state by state to the macOS record at `f852726` (`f030711`) |
| *inherited:* GitHub Actions `ubuntu-latest`, Python 3.11 / 3.12 / 3.13 | CI on `4aa4753` (`f2165cd`) | green (run 36186404223) |
| macOS arm64 | step 2 above, on the release commit `afc43db` | CPython 3.11.6: 824 passed / 0 failed, exact-mode guards 199 passed, validation extras 1021 passed. Under CPython 3.12 / 3.13 the exact capture guard failed on last bits. That was a test-contract finding, corrected in 0.4.1 (ADR-0032 amendment) |

**Not claimed:** agreement with experiment. The checks compare code with
independent implementations (teqp, FeOs, `thermo`) and published worked
problems.

### Known limitations
"Stable" is bounded by a deterministic trial set; no packaged `kij`; the two
packaged NRTL pairs are synthetic; residual properties only (no ideal-gas
`Cp`); PC-SAFT has no polar terms and only the 2B scheme is validated;
polymers are monodisperse; a three-phase PC-SAFT flash takes ~35 s;
`flash_mode="gamma-phi"` is deprecated.

### Install
```bash
pip install chemthermo==0.4.0
pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@v0.4.0"
```

### Artifacts (the 2026-09-25 cloud build of `4aa4753`, for reference; the release carries the digests of the files built from the tag and uploaded)
```
5b14fe01bab2ef811b1d6e9a70dacd03e796b7ad98dfc2c041f926003addee59  chemthermo-0.4.0-py3-none-any.whl
fd23392df1bcdaedd64315caf12272c9f654ba3302d69de045e3918ad9c210a3  chemthermo-0.4.0.tar.gz
```
