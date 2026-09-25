# Continuation plan (autonomous BUILD)

Owner's standing authorization, 2026-09-25. Read
`.agents/handoffs/cloud-continuation.md` first; its section 6 (boundaries) and
ADR-0031 (releases) bind everything here. Do not pause for another product
interview because several defensible next tasks exist: sequencing and release
decisions are delegated. Reordering is allowed when repository evidence
warrants it - record the reason in `steering-brief.md` and continue.

Every slice: state "After this change, user can X by running Y"; acceptance
criteria written before code; one worker at most; review diff and evidence
before accepting; commit with `Slice: <slug>`; update brain/steering/ledger/ADR;
push. Validation cost is risk-based (item D): docs-only slices run the cheap
checks; model/solver slices run the full default suite plus the relevant
ledger/validation tests; the 37-min sweep only at release boundaries after a
solver change or when a change could move a map bucket.

## A. Cloud baseline (first, every new environment)

(Handoff SHA: the one named in your prompt; the release it contains is
`v0.2.0b1` = `bd37068`. Start from the handoff SHA, not the tag.)

1. `git fetch origin && git checkout <handoff SHA>`; verify
   `git merge-base --is-ancestor 5041dd7071c6fd7128456cf643d59cac07f8583b HEAD`.
   Never start from `main` (it is the Feb 2026 baseline, 121+ commits behind).
2. Inspect the permitted push target (`git push --dry-run`); if only a session
   branch is allowed, create it **from the handoff SHA** and later open a PR
   targeting `dev/sprint`, never `main`.
3. Bootstrap per `.agents/dev-contract.md`; run ruff, pyright, `pytest -q`,
   `tools/smoke_install.py`, the golden path, `tools/release_smoke.py` against
   a non-editable install, and `python -m chemthermo.bench robustness --quick`.
4. Record results in section 7 of the handoff as "cloud baseline" with OS,
   Python, numpy versions. A cloud-only failure is recorded as such (with
   evidence), not assumed to invalidate the science.
5. Expect the 3 Linux failures of handoff section 7 on a Linux cloud host;
   anything else failing is new information to record. Then do A2 before B.

## A2. Cross-platform bit-identity guards (first implementation slice)

CI on ubuntu-latest is red on 3 tests (handoff section 7): two bit-identity
fixtures differ in the last ULP and one diagnostic count flips on a near-tied
state. Verdicts agree. Goal: CI green on Linux **without weakening** what the
guards protect.

- Diagnose first: confirm each difference is floating-point platform noise
  (libm/`exp`/`log`, BLAS/`linalg.solve`, summation order), not
  nondeterminism on one platform (run each test twice on Linux, compare).
- Keep exact bit-identity where the fixture's platform matches (record
  platform/machine/numpy in the fixture) and on other platforms compare with
  an explicit, justified bound (e.g. a few ULP / 1e-13 relative on floats)
  plus **exact** equality of every discrete field (phase names and count,
  statuses, verdicts, converged stage, iteration counts only if stable). Write
  this as an ADR amending the bit-identity contract of ADR-0017/0023/0030.
- The surface-count test: pin the verdicts exactly and make the surface
  statement robust to exact ties (for example, count a state for a surface
  only when its minimum beats the other surface's by a stated margin), with
  the per-state reason recorded in ledger Case P-11.
- Acceptance: GitHub Actions green on `ubuntu-latest`; the local/macOS run
  still exact; ledger and ADR record the rule; no test deleted or marked
  `slow`/`xfail` to get green.
- **Release checkpoint:** once CI is green on Linux, `0.2.0b2` (or `0.2.0` if
  nothing else changed and the owner-level policy in ADR-0031 item 3 is met).

Also queued (docs hygiene, cheap, any time): make `aglint check --repo .`
pass - fix the renamed paths it reports and decide whether `file::test`
node ids are written differently or the linter configured; never delete
evidence references to silence it.

## B. CLI exposure of the delivered equilibrium work (default next priority)

Thin adapters over `stability_tp` / `flash_tp`; no new equilibrium code.
Contract today (ADR-0003/0004): `chemthermo tp-flash`, Peng-Robinson phi-phi
and deprecated gamma-phi, `cli_schema_version = 1`, exit codes 0/1/2/3,
golden fixtures `tests/fixtures/cli/tp_flash_v1.json`,
`tp_flash_gamma_phi_v1.json`, tests `tests/test_cli_tp_flash.py`.

- **B0 contract ADR (ADR-0032).** Decide explicitly: subcommand layout
  (recommended default: new `chemthermo stability-tp`; extend `tp-flash` with
  `--eos {peng-robinson,pc-saft}` and `--max-phases`, and a way to reach
  `modified-raoult`/`gamma-gamma` only where packaged parameters make it
  meaningful); schema evolution (existing invocations must produce
  byte-identical v1 output - additive keys only behind new flags, or a v2
  schema if an N-phase `phases` layout cannot be expressed additively);
  deterministic output (sorted keys, stable float formatting, phase names as
  the library returns them); diagnostics subset; error -> exit code mapping
  (inconclusive stability -> which code?); how "stable" is worded
  (`"stable (bounded trial set)"` or a `stability_scope` field).
- **B1 `stability-tp`**: PR and PC-SAFT; JSON with status, tpd_min, incipient
  composition, trials summary. Acceptance: golden tests for a stable feed
  (Methane/Ethane 300 K, 1e5 Pa, PR) and an unstable feed (Methane/n-Hexane
  300 K, 2e6 Pa); exit-code tests; `--help` text states the bounded meaning.
- **B2 multiphase `tp-flash`**: `--eos pc-saft`, `--max-phases`, N-phase
  output. Acceptance: v1 fixtures unchanged; new golden fixtures for
  single-phase, VLE, and the three-liquid PR water/ethanol/n-hexane 280 K
  1 atm case (and PC-SAFT VLLE only if runtime allows; ~30 s per flash);
  mass-balance and phase-fraction invariants asserted from the CLI JSON.
- **B3 docs + examples**: README CLI section, `examples/cli/` shell examples,
  brain section 2 update.
- **Release checkpoint:** `0.3.0b1` (or `0.3.0` if the contract is settled
  and gates pass in two environments) once B1-B3 are accepted.

## C. PC-SAFT temperature derivatives and the properties they support

- **C1** analytic `d(A^res/RT)/dT` at fixed `(rho, x)` for hard chain +
  dispersion. Checks: teqp `get_Ar10` (to ~1e-12), central finite
  differences, consistency with existing derivatives.
- **C2** association contribution (Michelsen-Hendriks form; the first `T`
  derivative is an explicit partial). Checks: FeOs, finite differences.
- **C3** residual properties only: `H^res`, `S^res`, `U^res` (and `G^res`)
  with explicit reference state (ideal gas at same `T`, `rho` or `T`, `P` -
  state which, test both conversions). **Not** total enthalpy/Cp: those need
  ideal-gas heat capacities the databank does not carry - say so in docs and
  errors; a `Cp^res` needs a second `T` derivative (separate slice).
- **C4** `ln phi` temperature derivatives only if a consumer (e.g. the
  Venkatarathnam-Oellrich criterion, ADR-0017) is in the same slice.
- Each slice: ledger case with independent route, golden example, ADR for the
  public method names. **Release checkpoint:** a minor version once C1-C3
  form a coherent "residual caloric properties" milestone.

## D. Performance and validation cost (when measured)

- Default suite exceeds the 4-minute target (~4:20-4:45 measured locally;
  `tests/test_examples.py` ~86 s). Profile (`pytest --durations=30`) before
  changing anything; move tests to `slow` only per the dev-contract rule
  (never the only test of a capability) and record the change in the ledger.
- Any solver speed-up follows ADR-0023/0030: baseline + after records, equal
  result hashes, bit-identity tests, measured ratio.
- A tiered gate (cheap checks per slice, full suite per model slice, sweep per
  release) may be written into `dev-contract.md` - justified and recorded, not
  silent.

## E. Deferred science (do not start without provenance + executable acceptance)

Polydispersity (ADR-0022 roadmap item), water/alcohol `kij` for the R-MAP-1
(iv) finding, exact vs rounded PR constants, full gamma-phi with Poynting,
`gamma-gamma` third liquid. Each needs a cited source, a reference to check
against, and an acceptance test before a slice is opened.

## Session end checklist

`git status --porcelain` empty; `git log --oneline -n 20`;
`git log -n 50 --format=%B | grep '^Slice:'`; branch pushed; handoff section 7
and this plan updated with the new SHA, what was accepted, and what is next.
