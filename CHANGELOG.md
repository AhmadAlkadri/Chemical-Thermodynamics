# Changelog

Versions follow PEP 440 and are tagged `v<version>` on one tested commit
(policy: `.agents/brain/adr/0031-versioned-releases.md`). Releases up to
0.3.0b1 were GitHub-only; from 0.4.0 the owner also publishes to PyPI by hand
(ADR-0038).

## Unreleased

### Added
- PC-SAFT **association** temperature derivative: `residual_helmholtz_temperature_derivative`
  and `residual_properties` now work for water, alcohols and their mixtures
  (they raised `ModelError` in 0.3.0b1). Checked against FeOs residual entropy
  and enthalpy (2.1e-15 with matched universal constants). ADR-0034 amendment,
  ledger Case P-20.

### Fixed
- `flash_tp` (phi-phi): a split whose two phases are identical (the trivial
  solution) no longer counts as converged, so the stability-seed ladder runs.
  Around the Mw 53000 polyethylene ladder state, 29 of 135 neighbouring
  pressures (+-64 ULP) refused; now none does. The full 2505-state
  robustness map is unchanged field for field. ADR-0036.

### Changed
- `stability_tp`: on a `tpd` tie (within 1e-12) the reported minimizing trial
  is one converged at least 1000x better, when one exists; flagged by
  `diagnostics["minimizing_trial_tie_break"] = "residual"`. Changes only such
  ties (none on any captured state). ADR-0035.

## 0.3.0b1 (2026-09-25) - prerelease

The CLI reaches the equilibrium work, PC-SAFT gains residual caloric
properties, and the test suite is green on Linux for the first time. A
prerelease because the CLI contract (ADR-0033) and the two PC-SAFT methods
(ADR-0034) are new. Validation now spans two platforms (macOS arm64 inherited
at 0.2.0b1; Linux x86_64 on GitHub Actions and a cloud session for this
release), but this release's new code was exercised on Linux only.
There is no 0.2.0b2; its planned content (the cross-platform guards) is here.

### Compatibility
- Library: additive only (two new `PCSAFTEOS` methods). No existing number
  moved: the new derivative is a separate function, and no solver was touched.
- CLI: `cli_schema_version` stays 1; previously valid invocations print
  byte-identical output (checked on 11 invocations).

### Added
- **CLI `chemthermo stability-tp`**: tangent-plane stability of a feed with
  Peng-Robinson or PC-SAFT (`--eos`), JSON or text; `stable` is reported with
  `stability_scope = "bounded-trial-set"`; an `inconclusive` verdict exits 3
  with the payload printed. ADR-0033.
- **CLI `tp-flash --eos {peng-robinson,pc-saft}` and `--max-phases N`**:
  PC-SAFT flashes (packaged parameters, `kij = 0`, phi-phi) and an explicit
  phase budget; three-phase answers use the existing name-keyed layout.
  `cli_schema_version` stays 1 and every previously valid invocation prints
  byte-identical output. `examples/cli/stability_and_multiphase.sh`.

- **PC-SAFT temperature derivative and residual properties** (non-associating):
  `PCSAFTEOS.residual_helmholtz_temperature_derivative` and
  `PCSAFTEOS.residual_properties` (`h_res`, `u_res`, `s_res_tv/tp`,
  `g_res_tv/tp`, reduced). Checked against teqp `get_Ar10` (<= 5.2e-16 relative
  over 14 states), Gibbs-Helmholtz through the fugacity route, and
  `dH_vap = T dS_vap` at saturation. Associating mixtures raise `ModelError`.
  ADR-0034; `examples/basic/pcsaft_residual_properties_demo.py`.

### Changed (tests and policy only; no library code)
- Guards that compare against floats captured on macOS arm64 are exact there
  and, on any other platform, exact on every discrete field and bounded on
  floats (1e-12 flash fixture, 5e-14 PC-SAFT literals); last-bit ties are
  counted as ties. Makes the Linux CI and cloud runs green without weakening
  the capture-platform check. ADR-0032; ledger Cases P-11 and P-17.

### Known limitations found
- One polymer ladder state (Mw 53000, 15 wt%, 8.1 MPa) refuses at 1 of 17
  pressures within +-8 ULP of the grid point, and its solver route depends on
  the last bit; every converged answer is the same tie line.
- `stability_tp` breaks an exact `tpd` tie by trial order, so the reported
  `trial_composition` can be the less-converged of two tied trials (good to
  the stationarity tolerance, 1e-10 residual).

## 0.2.0b1 (2026-09-25) - prerelease

First release of the phase-equilibrium campaign (121 commits and 38
distinct `Slice:` trailers since `v0.1.0`; ADR-0003 to ADR-0031). A prerelease because its
validation comes from one development machine plus the release gates of
ADR-0031, and because the CLI does not yet expose the equilibrium work.

### Added
- **Phase stability** `chemthermo.stability_tp` (Michelsen tangent-plane
  distance) with `StabilityResult` / `StabilitySettings` / `StabilityTrial`,
  incipient-phase compositions and per-trial diagnostics, for an equation of
  state, an activity model alone (liquid-liquid), or an activity liquid against
  an ideal-gas vapour (`vapor="ideal"`). ADR-0005, 0007, 0010, 0012, 0021, 0025.
- **Phase-count discovery in `flash_tp`**: one, two or three phases found by
  stability analysis, a verified split, post-split stability of every phase,
  and phase addition/removal (`FlashSettings.max_phases`, default 3), on the
  EOS (`phi-phi`) and `modified-raoult` paths; liquid-liquid `gamma-gamma`
  mode (two phases). ADR-0008, 0009, 0011, 0016, 0019, 0020, 0024, 0026, 0029.
- **PC-SAFT** `chemthermo.PCSAFTEOS`: Gross-Sadowski (2001) hard chain +
  dispersion with the 2002 association term, analytic composition/density
  derivatives, density roots, log-space fugacities for long chains; packaged
  parameters for 11 non-associating compounds; user records incl. association
  (`PCSAFTRecord`, `PCSAFTAssociationRecord`) and polymers via
  `segments_per_g` + `Component.custom(...)`. ADR-0014, 0015, 0018, 0022.
- Peng-Robinson `kij` matrices (ADR-0006); phase labels by compressibility
  (ADR-0017).
- CLI `chemthermo tp-flash` (`cli_schema_version` 1; Peng-Robinson phi-phi and
  deprecated gamma-phi; exit codes 0/1/2/3). ADR-0003, 0004.
- Internal (not public API) benchmark and robustness harness
  `python -m chemthermo.bench` (ADR-0023, 0027) and a call-local EOS solve memo
  (ADR-0030).

### Changed / compatibility
- `NRTL.activity_coefficients` now implements the standard Renon-Prausnitz
  equation (column sums); values for **asymmetric** parameters differ from
  `v0.1.0`, which violated Gibbs-Duhem. Bug fix, not a model change.
- `flash_tp`: `eos` is optional; `flash_mode=None` infers the mode; phi-phi
  phase detection defaults to `"tangent-plane"` (`"wilson-heuristic"` keeps the
  old behaviour). Single-phase tangent-plane results no longer carry
  `k_min`/`k_max`/`max_delta_k`/`rr_*` diagnostics.
- `flash_mode="gamma-phi"` and the `chemthermo.vlle` plugin boundary are
  **deprecated** (docs only; still work). ADR-0013.
- Package version moves from `0.0.0` to `0.2.0b1`. (The existing `v0.1.0`
  tag points at a commit whose metadata says `0.0.0`; it is left as published.)

### Known limitations
- `stability_tp` "stable" means no negative tangent-plane distance was found
  from a bounded deterministic trial set - not a proof of global stability.
- PC-SAFT parameters and the association equation rest on secondary sources
  (primary papers not read); polymer parameters are a labelled test fixture
  from one secondary source and are **not packaged**; polymers are monodisperse.
- No temperature derivatives for PC-SAFT: no residual enthalpy/entropy or
  caloric properties.
- The robustness map (2505 states, macOS arm64) refuses 5 states, all in the
  deprecated gamma-phi path by design; it measures coverage, not correctness.
- The CLI does not expose stability, multiphase flash or PC-SAFT yet.
- Default test suite takes ~4-5 min; `pytest -m slow` and the full sweep
  (~37 min) are not run by CI.

## 0.1.0

Tag `v0.1.0` on `a7a8ca7` (Feb 2026). No GitHub release was created and its
package metadata reads `0.0.0`.
