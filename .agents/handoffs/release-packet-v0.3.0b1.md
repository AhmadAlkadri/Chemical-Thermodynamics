# Release packet: chemthermo 0.3.0b1 (prerelease)

Prepared 2026-09-25 by the first cloud session. **Neither the tag nor the
GitHub release exists yet.** The session could push the branch but not a tag:
`git push origin refs/tags/v0.3.0b1` was refused three times with "remote end
hung up" and no proxy-side failure logged, which fits a session git policy
that allows only the designated branch. It also has no release-creation or
asset-upload tool. The release gates below all passed on the exact commit, so
the owner only needs to tag it and publish:

```bash
git fetch origin
git tag -a v0.3.0b1 64130a2d5373e88cc65c28fdd047807f04a5daf4 -m "chemthermo 0.3.0b1"
git push origin refs/tags/v0.3.0b1
git ls-remote origin 'refs/tags/v0.3.0b1^{}'   # must print 64130a2d5373e88cc65c28fdd047807f04a5daf4
git clone --branch v0.3.0b1 https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git rel && cd rel
python3.11 -m venv .venv && .venv/bin/pip install build && .venv/bin/python -m build --outdir dist/
shasum -a 256 dist/*
gh release create v0.3.0b1 --verify-tag --prerelease --title "chemthermo 0.3.0b1" \
  --notes-file ../Chemical-Thermodynamics/.agents/handoffs/release-packet-v0.3.0b1.md dist/*
```

(Trim this preamble from the notes, or pass a copy without it.) Rebuilt
artifacts are not byte-identical across machines, so replace the checksums
below with those of the files actually attached.

---

## chemthermo 0.3.0b1

**Tested commit:** `64130a2d5373e88cc65c28fdd047807f04a5daf4` (branch `dev/sprint`).
Not on PyPI.

### What is new
- **CLI `chemthermo stability-tp`**: tangent-plane stability with Peng-Robinson
  or PC-SAFT; `stable` is reported as `stability_scope = "bounded-trial-set"`
  (no negative TPD found from a deterministic trial set, not a global proof);
  `inconclusive` exits 3 with the JSON still printed. ADR-0033.
- **CLI `tp-flash --eos {peng-robinson,pc-saft}` and `--max-phases N`**: PC-SAFT
  flashes (packaged parameters, `kij = 0`, phi-phi) and an explicit phase
  budget; three-phase answers use the existing name-keyed layout.
- **PC-SAFT residual properties** (non-associating):
  `PCSAFTEOS.residual_helmholtz_temperature_derivative`,
  `PCSAFTEOS.residual_properties` (`h_res`, `u_res`, `s_res_tv/tp`,
  `g_res_tv/tp`). ADR-0034.
- **Test suite green on Linux** for the first time. Guards pinned to floats
  captured on macOS arm64 are still exact there and, elsewhere, exact on every
  discrete field and bounded on floats. ADR-0032.

### Compatibility
- Additive only. `cli_schema_version` stays 1, and 11 pre-existing CLI
  invocations print byte-identical output. No solver changed and no existing
  number moved: `refactor_bit_identity_v3.json` was not regenerated.

### Validation actually run for this release
| where | what | result |
| --- | --- | --- |
| clean clone from GitHub at `64130a2`, Linux 6.18 x86_64, CPython 3.11.15, numpy 2.4.6 | ruff format/check, pyright, `pytest -q` | clean, clean, 0 errors, **814 passed**, 51 skipped (extras absent), 74 deselected |
| same | `python -m build`; wheel in a fresh venv outside the tree | `__version__ == 0.3.0b1`; `tools/release_smoke.py` pass; wheel `stability-tp --eos pc-saft`, three-liquid `tp-flash`, `residual_properties` run |
| GitHub Actions ubuntu-latest (Python 3.11.16) | CI on `64130a2`, run 36127239537 | **green**: install, ruff, pyright, pytest, non-editable smoke, trailers, clean tree |
| cloud host working tree at `64130a2`, validation extras (teqp 0.23.2, FeOs 0.10.1, thermo 0.6.1) | `pytest -q` incl. `tests/validation/` | **975 passed**, 0 failed, 98 deselected, 11:36 (first Linux run of the extras) |
| cloud host | `chemthermo.bench robustness --quick` at `5864ab0` (no solver change since) | 224 / 219 / 5 by-design, identical totals to macOS |

**Inherited, not rerun:** the macOS arm64 gates of 0.2.0b1; `pytest -m slow`;
the full 2505-state robustness sweep (no solver change since its record at
`f852726`, so ADR-0031 item 6 does not require it); benchmark timings.

### Known limitations
- New code was exercised on Linux only; nobody has run the ADR-0032 guards on
  macOS arm64 yet. They should be exact there by construction.
- The association term has no temperature derivative, so both new PC-SAFT
  methods raise `ModelError` for an associating mixture. No `Cp^res`, and no
  total caloric properties (no ideal-gas heat capacities are packaged).
- `modified-raoult` / `gamma-gamma` are not reachable from the CLI.
- One polymer ladder state refuses at 1 of 17 pressures within +-8 ULP of its
  grid point (ledger Case P-17 "cross-platform"). `stability_tp` breaks an
  exact `tpd` tie by trial order, so the reported trial is only good to the
  stationarity tolerance (plan items A3/A4).
- Everything listed under 0.2.0b1 in `CHANGELOG.md` still applies.

### Install pins
```bash
# once the tag is pushed:
pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@v0.3.0b1"
# works now:
pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@64130a2d5373e88cc65c28fdd047807f04a5daf4"
```

### Artifacts (cloud build from the clean clone)
```
b1f32972b13ca9bf35f51b0bfb0e561d31a4cf6a35a7e296636485b38436baa6  chemthermo-0.3.0b1-py3-none-any.whl
9a2b366fcb8314903d206eb0c63fa66a8e18afd924c97bfdae96736e17cfd7f5  chemthermo-0.3.0b1.tar.gz
```
