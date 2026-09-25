# Release packet: chemthermo 0.4.0 - merge to main, tag, GitHub release, PyPI

Prepared 2026-09-25 by the cloud session. Everything up to "tag" was done and
checked from the cloud; the cloud session can push only `dev/sprint`, so the
steps below are the owner's (or an owner-side agent with `gh` and PyPI
credentials). Policy: ADR-0031 (tags, releases), ADR-0038 (PyPI by hand).

| item | value |
| --- | --- |
| release commit (tested, to be tagged) | `4aa475385791b29f7da76e68332ebfdd39870f1e` |
| version in the code at that commit | `0.4.0` (`pyproject.toml`, `chemthermo.__version__`) |
| PR | https://github.com/AhmadAlkadri/Chemical-Thermodynamics/pull/1 (`dev/sprint` -> `main`) |
| `main` today | `a7a8ca7` (`v0.1.0`), an ancestor of `dev/sprint`: a fast-forward, no merge commit |
| PyPI name | `chemthermo` (unregistered on 2026-09-25) |

## Steps for the owner

Stop at the first failure and report it; never force-push, move a tag, or
upload a version twice (PyPI refuses re-uploads; a fix is 0.4.1).

1. **CI green on PR #1** (all three Python versions), on the PR's head.
2. **macOS arm64 exact check on the release commit** (ADR-0032 capture
   platform; ADR-0035/0036 were proven dormant bit for bit on Linux only):
   ```bash
   git clone https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git ct-040 && cd ct-040
   git checkout 4aa475385791b29f7da76e68332ebfdd39870f1e
   python3.11 -m venv .venv && .venv/bin/pip install -e ".[dev]"
   .venv/bin/pytest -q      # expect 824 passed, 51 skipped, 0 failed
   .venv/bin/pytest -q tests/test_flash_refactor_bit_identity.py tests/test_pcsaft_association.py \
     tests/test_stability_eos_surfaces.py tests/test_pcsaft_polymer.py tests/test_stability_candidates.py
   ```
   If an exact-mode guard fails on macOS, do not continue: record the output
   in the handoff and hand it back to a cloud session.
3. **Fast-forward `main`** (the PR then shows as merged):
   ```bash
   git fetch origin
   git merge-base --is-ancestor origin/main origin/dev/sprint && \
   git merge-base --is-ancestor 4aa475385791b29f7da76e68332ebfdd39870f1e origin/dev/sprint && \
   git push origin origin/dev/sprint:refs/heads/main          # a fast-forward; never --force
   ```
   (Equivalently `gh pr merge 1 --merge` would add a merge commit; the
   fast-forward keeps `main` on the tested history. Do **not** squash or
   rebase-merge: both rewrite the 147 commits.)
4. **Tag the tested commit:**
   ```bash
   git tag -a v0.4.0 4aa475385791b29f7da76e68332ebfdd39870f1e -m "chemthermo 0.4.0"
   git push origin refs/tags/v0.4.0
   git ls-remote origin 'refs/tags/v0.4.0^{}'   # must print 4aa475385791b29f7da76e68332ebfdd39870f1e
   ```
5. **Build from the tag, check, GitHub release** (not a prerelease):
   ```bash
   git clone --branch v0.4.0 https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git rel-040 && cd rel-040
   python3.11 -m venv .venv && .venv/bin/pip install build twine
   .venv/bin/python -m build --outdir dist/ && .venv/bin/twine check --strict dist/*
   shasum -a 256 dist/*
   gh release create v0.4.0 --verify-tag --title "chemthermo 0.4.0" --notes-file NOTES.md dist/*
   ```
   `NOTES.md` = this file from "## chemthermo 0.4.0" down, with the
   checksums replaced by the ones just printed and the macOS row filled in.
6. **TestPyPI** (needs a TestPyPI account and API token; `twine` reads
   `TWINE_USERNAME=__token__` / `TWINE_PASSWORD=<token>` or `~/.pypirc`):
   ```bash
   .venv/bin/twine upload --repository testpypi dist/*
   python3.11 -m venv /tmp/tp && /tmp/tp/bin/pip install \
     --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ chemthermo==0.4.0
   cd /tmp && /tmp/tp/bin/python <path-to>/rel-040/tools/release_smoke.py --expect-version 0.4.0
   /tmp/tp/bin/chemthermo --help
   ```
   Check the rendered project page on test.pypi.org (README, links, classifiers).
7. **PyPI:**
   ```bash
   .venv/bin/twine upload dist/*
   python3.11 -m venv /tmp/pp && /tmp/pp/bin/pip install chemthermo==0.4.0
   cd /tmp && /tmp/pp/bin/python <path-to>/rel-040/tools/release_smoke.py --expect-version 0.4.0
   ```
8. **Record it** on `dev/sprint` (docs-only commit, `Slice: handoff`): handoff
   section 7e with the macOS result, tag peel, release URL, PyPI URL and the
   uploaded files' SHA-256; then fast-forward `main` to that commit the same
   way as step 3 if you want `main` to carry the record.

**One licensing note for the owner to judge before uploading.** The packaged
component table (82 components: critical constants, acentric factors, Antoine
coefficients) cites Koretsky, *Engineering and Chemical Thermodynamics*
(Wiley), as its source (`src/chemthermo/data/references.bib`). These are
physical constants rather than creative text, and the table has shipped in the
public GitHub repository since v0.1.0, but redistributing a textbook appendix
on PyPI is a call for the owner. The PC-SAFT parameters are published values
cross-checked against FeOs / Clapeyron.jl; the two NRTL pairs are synthetic;
the DECHEMA-derived Tessier fixtures are test-only and are **not** in the
sdist or wheel.

---

## chemthermo 0.4.0

**Tested commit:** `4aa475385791b29f7da76e68332ebfdd39870f1e`. First release
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
| where | what | result |
| --- | --- | --- |
| clean clone of `4aa4753`, Linux x86_64, CPython 3.11.15, numpy 2.4.6 | ruff, pyright, `pytest -q` | clean, 0 errors, **824 passed**, 51 skipped (extras absent) |
| same | `python -m build`, `twine check --strict` | both artifacts PASSED |
| same, wheel in fresh venvs outside the tree, Python 3.11 / 3.12 / 3.13 | `release_smoke.py --expect-version 0.4.0`, CLI `stability-tp --eos pc-saft`, `chemthermo.vlle` absent | pass on all three |
| Linux, Python 3.12.3 and 3.13.12 (numpy 2.5.3) | `pytest -q` at `bc92c5b` (library code as released) | 824 passed each |
| Linux, validation extras (teqp 0.23.2, FeOs 0.10.1, thermo 0.6.1) | `pytest -q` incl. `tests/validation/` at `4aa4753` | **1021 passed**, 0 failed, 99 deselected |
| Linux, full robustness map at `761fd57` (after the last solver change) | 2505 states | 2500 converge, 5 by-design gamma-phi refusals, 0 invariant violations; identical state by state to the macOS record at `f852726` |
| GitHub Actions `ubuntu-latest`, Python 3.11 / 3.12 / 3.13 | CI on the release commit | green (run 36186404223) |
| macOS arm64 | step 2 above | *(owner fills in)* |

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

### Artifacts (cloud build from the clean clone; replace with the uploaded files' digests)
```
5b14fe01bab2ef811b1d6e9a70dacd03e796b7ad98dfc2c041f926003addee59  chemthermo-0.4.0-py3-none-any.whl
fd23392df1bcdaedd64315caf12272c9f654ba3302d69de045e3918ad9c210a3  chemthermo-0.4.0.tar.gz
```
