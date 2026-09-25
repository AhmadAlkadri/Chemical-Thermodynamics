# Cloud continuation handoff

Written 2026-09-25 by the local orchestrator when the equilibrium campaign was
first published; section 7b and the checkpoint table updated the same day by
the first cloud session. It is self-contained: a fresh session needs this repository
and nothing from the original machine, home directory or session memory.
What to do next is in `.agents/handoffs/continuation-plan.md`.

## 1. Checkpoint

| item | value |
| --- | --- |
| repository | `https://github.com/AhmadAlkadri/Chemical-Thermodynamics` |
| development branch | `dev/sprint` (default branch `main` is the old baseline; never start from it) |
| campaign checkpoint (first push) | `5041dd7071c6fd7128456cf643d59cac07f8583b` |
| old baseline | `main` = tag `v0.1.0` = `a7a8ca742cadf41ad228ffa8c23d95d3252f5ce9` (an ancestor of `dev/sprint`) |
| handoff commits | on top of `5041dd7`: CI fix, ADR-0031, `tools/release_smoke.py`, these files, `CHANGELOG.md`, version `0.2.0b1` |
| release | `v0.2.0b1` (GitHub prerelease, source + wheel; not on PyPI) - see section 7 |
| first cloud session (2026-09-25) | slices `cloud-baseline`, `cross-platform-guards` (ADR-0032), `cli-contract` + `cli-stability-multiphase` (ADR-0033), `pcsaft-temperature-derivative` (ADR-0034), `release`; tag `v0.3.0b1` on `64130a2d5373e88cc65c28fdd047807f04a5daf4` - section 7b |

Verify before editing: `git merge-base --is-ancestor 5041dd7071c6fd7128456cf643d59cac07f8583b HEAD`
must succeed, and `HEAD` must be the SHA your prompt names (or a descendant
you created).

## 2. Reading order

1. `AGENTS.md` - hard gates (clean tree at handoff, `Slice: <slug>` trailers, evidence commands).
2. `.agents/brain/brain.md` - invariants, public API (section 2), architecture, ADR index (section 8), roadmap (section 9), risks (section 10).
3. `.agents/dev-contract.md` - every command: bootstrap, CI gates, `slow` marker, benchmarks, robustness map, releases.
4. `.agents/brain/steering-brief.md` - "Public API status", "Risks / unknowns", "Next 3 recommended actions".
5. ADRs in `.agents/brain/adr/`: 0001 (public API truth source), 0002 (thin vertical slices), **0003/0004 (CLI contract)**, 0005 (stability API), 0014/0018/0022 (PC-SAFT, association, polymers), 0023/0027/0030 (measurement contracts), **0031 (releases)**.
6. `.agents/brain/validation-cases.md` - the ledger (5800 lines; search by case id: S, K, N, F, L, R, V, P, B, R-MAP).
7. `benchmarks/README.md` and `benchmarks/robustness_f852726.md` - latest map summary.
8. `.agents/skills/chemthermo-change-loop/SKILL.md` - the per-slice workflow.
9. `CHANGELOG.md` - what `0.2.0b1` contains and its limitations.

## 3. What is accepted (inherited evidence, 2026-09-13/14, macOS arm64, one machine)

38 distinct `Slice:` slugs in `v0.1.0..5041dd7` (`git log --format=%B v0.1.0..5041dd7 | grep '^Slice:' | sort -u`).
Each capability's evidence is its ADR plus ledger cases:

- Tangent-plane stability with incipient phases, EOS and activity models - ADR-0005/0007/0010/0012/0021/0025; ledger S-*, N-*, K-*.
- `flash_tp` one/two/three-phase discovery (stability -> verified split -> post-split stability -> phase addition/removal) - ADR-0008/0009/0011/0016/0019/0020/0024/0026/0029; ledger F-*, L-*, V-*, P-9, P-10, P-18.
- PC-SAFT (Gross-Sadowski 2001 + 2002 association), analytic derivatives, density roots, log-space fugacities - ADR-0014/0015/0018; checked against teqp (non-associating, to ~1e-14) and FeOs (association) - ledger P-0..P-11, examples `examples/validation/13,14,16,17,18`.
- Labelled polymer fixture (`tests/fixtures/pcsaft/martini2009_polymers.json`), PE/n-pentane VLE and LLE - ADR-0022/0024/0026/0028; ledger P-12..P-17.
- Robustness map: 2505 states, 5 refusals all in the deprecated gamma-phi path by design, 0 invariant violations, record `benchmarks/robustness_f852726.{json,md}` - ADR-0027/0029; ledger R-MAP-1/2. Coverage, not correctness.
- Performance: benchmark records `benchmarks/*.json`, EOS call-local memo - ADR-0023/0030; ledger B-*.

Fresh 2026-09-25 results (clean clone, GitHub Actions) are in section 7; they
are the only results in this file not inherited.

## 4. Limitations, deferred work, failed approaches

- "Stable" = no negative TPD from a bounded deterministic trial set (never a global proof).
- PC-SAFT primary papers (DOI 10.1021/ie0003887, 10.1021/ie010954d) were paywalled (HTTP 403); parameters rest on two agreeing secondary sources (FeOs `gross2001.json`, Clapeyron.jl) and the association `Delta` uses `sigma^3` - decided by exact numerical agreement with FeOs, not by reading the equation.
- Polymer parameters: one secondary source (Martini et al. 2009, CONICET open repository), fixture only, never packaged; monodisperse only.
- Tessier (2000) NRTL fixtures are DECHEMA-derived: test fixtures only, never packaged runtime data (brain.md section 4).
- No PC-SAFT temperature derivatives: no residual H/S, no caloric properties; also blocks the Venkatarathnam-Oellrich `Pi` criterion ADR-0017 rejected.
- Unadjudicated: PC-SAFT 2B water/ethanol at `kij = 0` gives a verified LLE for a miscible pair (ledger R-MAP-1 (iv)) - a parameterization issue needing a cited `kij`.
- Failed / declined approaches worth remembering:
  - Per-iterate min-Gibbs root re-selection inside the split collapses phases; ADR-0019 pins each phase to its stability branch.
  - Caching root solves on the (frozen) model instance was rejected (ADR-0023); the accepted memo is call-local via `ContextVar` (ADR-0030).
  - The structured Hessian speed-up is **not** bit-identical (0/4 builds, 2.5e-04 rel.) and was declined (ADR-0030); shipping it is a re-audit slice.
  - Clamping `ln W` to [-700, 700] caused stability misses; removed in favour of log-space sums (ADR-0025).
- Cross-checking facts not obvious from code:
  - chemthermo's Peng-Robinson uses rounded 0.45724/0.07780; `thermo` uses exact roots, which floors agreement at ~1e-4 in `ln phi` (ledger S-4).
  - FeOs hard-codes 14-digit PC-SAFT universal constants vs the paper's 10; floors FeOs comparisons at ~1e-9 (`Z`) to ~1e-6 (`ln phi`).
  - `thermo` 0.6.0 cannot split two liquids over one excess model; FeOs 0.10.1 cannot take `kij` and its flash fails on PE/n-pentane.
- Cost: default suite ~4:20-4:45 on the dev machine (target was 4 min; `tests/test_examples.py` alone ~86 s); full robustness sweep ~37 min; `pytest -m slow` is not in CI.
- Settled, do not reopen: license stays MIT; this public library is the reference-capability implementation; stability-first architecture where PC-SAFT is a consumer, not the architecture; a model-family abstraction must be earned by at least two families.

## 5. Setup and commands (from `.agents/dev-contract.md`)

Required: Python >= 3.11 (CI uses 3.11), `pip`, `git`. Runtime deps: numpy, pydantic, bibtexparser<2.

```bash
python3.11 -m venv .venv
.venv/bin/pip install -e ".[dev]"                    # ruff, pyright, pytest, build
.venv/bin/ruff format --check src tests && .venv/bin/ruff check src tests
.venv/bin/pyright
.venv/bin/pytest -q                                  # default suite (slow deselected)
.venv/bin/python tools/smoke_install.py --package .  # non-editable install smoke
python examples/basic/flash_tp_peng_robinson_demo.py # golden path
```

Optional:

- `.venv/bin/pip install -e ".[validation]"` - `thermo`, `teqp`, `feos` (public PyPI, MIT/Apache). Without them the validation tests skip cleanly.
- `.venv/bin/pytest -q -m slow` - exhaustive grids (repetitions of default-covered capabilities).
- `python -m chemthermo.bench --out benchmarks/<name>.json` (~30 s); `--compare A B` exits 1 on a result-hash change.
- `python -m chemthermo.bench robustness --quick` (~14 s) or full (`--out benchmarks/robustness_<sha>.json --summary-out ...md`, ~37 min; `--family NAME` partitions it).
- `aglint check --repo .` needs `agentslint` from a **private** repository (`AhmadAlkadri/agentslint@v0.3.0`); CI skips it with a warning when it cannot install (ADR-0031). Not required for a cloud session.
- Installed-distribution smoke: `cd /tmp && <venv>/bin/python <repo>/tools/release_smoke.py --expect-version <v>`.

Local-only resources: none are required. The paywalled primary papers were never available; benchmark wall times are machine-specific (compare result hashes, not times, across machines).

## 6. Operating boundaries (delegated by the owner, 2026-09-25)

Authorized: focused commits; pushing the development/session branch; creating
and pushing deliberately chosen **new** version tags; GitHub releases or
prereleases at validated checkpoints (ADR-0031); PRs for review/CI.

Not authorized: force-push, history rewrite/squash, deleting branches or tags,
moving any published tag; merging into `main` or changing the default branch;
changing visibility, protections, credentials or access; publishing to PyPI or
any registry or adding such automation; committing private material,
credentials or raw session transcripts. Do not bypass sandbox or permission
restrictions or seek broader credentials; if tags/releases are not permitted,
prepare a release packet (version, notes, tested SHA, commands) and continue.

Operating model: one orchestrator, at most one worker at a time; thin
vertical slices with explicit acceptance criteria; review the implementation
and evidence before accepting; commit with `Slice:` trailers; keep
`brain.md`, `steering-brief.md`, the ledger and ADRs current; push accepted
checkpoints; claims proportional to validation; never broaden expected
failures or weaken tolerances to go green; end every session committed,
pushed, and with this file (or a successor) updated.

## 7. Release and fresh verification record (2026-09-25)

**Release `v0.2.0b1`** (GitHub prerelease, not PyPI):
<https://github.com/AhmadAlkadri/Chemical-Thermodynamics/releases/tag/v0.2.0b1>.
Annotated tag peels to `bd370681056baf487f0419a6dd05213b2d83932c` (verified by
`git ls-remote` and the GitHub commits API). Assets and SHA-256:
`chemthermo-0.2.0b1-py3-none-any.whl` `6096deaaabf6590ea0ec0603ca7990b632ae82f6234bdd17e28c9964122d6f55`,
`chemthermo-0.2.0b1.tar.gz` `8f2d6c833740c046961c89d605fd96c3c569fd5220dcfe3ba4e1cf14b205fc1d`.
Pins: `pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@v0.2.0b1"`
or `...@bd370681056baf487f0419a6dd05213b2d83932c`.

**Fresh, clean clone from GitHub at `bd37068`**, macOS 27.0 arm64 (the
development machine, so not independent cross-platform validation), CPython
3.11.6, numpy 2.4.6, validation extras not installed:

| check | result |
| --- | --- |
| ruff format / ruff check / pyright | clean / clean / 0 errors |
| `pytest -q` | 769 passed, 51 skipped, 74 deselected, 2:57 |
| `tools/smoke_install.py --package .`, golden path | pass |
| `chemthermo.bench robustness --quick` | 224 states, 219 verdicts, 5 by-design `rr-no-bracket` refusals |
| `python -m build`; wheel and sdist in fresh venvs outside the tree | 4 data files packaged; `tools/release_smoke.py --expect-version 0.2.0b1` pass |
| `pip install ...@v0.2.0b1` from GitHub + smoke | pass |
| `aglint check --repo .` (agentslint v0.3.0) | **66 P0**, pre-existing: mostly `file::test_name` references and renamed paths in ADRs/ledger/steering |

**GitHub Actions, ubuntu-latest, Python 3.11, run 36121300741 on `bd37068`**:
install, ruff, pyright pass; **pytest 3 failed / 766 passed** (8:06):

- `tests/test_flash_refactor_bit_identity.py::test_flash_tp_is_bit_identical_to_the_v3_capture` (last-ULP, first at `gamma-gamma-binary|z1=0.05`)
- `tests/test_pcsaft_association.py::test_non_associating_values_are_bit_identical[hexane_300_7700-7700.0]` (`-5.783742760059242` vs `-5.783742760059239`)
- `tests/test_stability_eos_surfaces.py::test_the_peng_robinson_grid_still_reaches_a_verdict_everywhere` (verdicts identical 47/97; `minimizing_trial_surface` counts liquid 46 / vapor 31 vs pinned 45 / 32)

No verdict differs. The bit-identity pins were captured on macOS arm64 and
are same-platform refactoring guards. Nothing was relaxed; see plan slice A2.
The CI run before the handoff (`5041dd7`, run 36120860309) never reached the
tests (private `agentslint` install), so this is the branch's first CI test run.

Not rerun in this session (inherited): validation extras (teqp/FeOs/thermo),
`pytest -m slow`, the full 2505-state map, benchmark timings.

### 7a. Cloud baseline (slice A, first cloud session, 2026-09-25)

Claude Code cloud container at handoff SHA
`5864ab0c9fe43a9a792a87b12b590aa1c644438d`: Linux 6.18.44 x86_64, 4 CPUs,
CPython 3.11.15, numpy 2.4.6, validation extras not installed.

| check | result |
| --- | --- |
| ruff format / ruff check / pyright | clean / clean / 0 errors |
| `pytest -q` | **5 failed**, 764 passed, 51 skipped, 74 deselected, 6:43 |
| `tools/smoke_install.py --package .`, golden path | pass |
| `python -m build`; wheel in a fresh venv, `release_smoke.py --expect-version 0.2.0b1` from outside the tree | pass (4 smoke flashes) |
| `chemthermo.bench robustness --quick` | 224 states, 219 verdicts, 5 by-design `rr-no-bracket` - identical totals to macOS |

The 5 failures: the 3 of CI run 36121300741, plus
`test_pcsaft_polymer.py::test_the_band_the_diverged_k_loop_used_to_end[53000.0-0.15-8100000.0]`
(route label `linear-iterate` vs `stability-w`, same tie line to 1.3e-12) and
`test_stability_candidates.py::test_results_are_invariant_under_component_reordering[z0-361.0]`
(8.4e-11 vs 1e-12: two tied trials reported in the two orders). The CI runner
and this host also disagree with each other on the PC-SAFT literals (`a_res`
exact here, 3 ULP off there). Diagnosis and the rule adopted: ADR-0032, ledger
Cases P-11 and P-17 "cross-platform". Slice A2 (`cross-platform-guards`)
followed; its result is below.

### 7b. First cloud session: what was accepted (2026-09-25)

| slice | commit(s) | evidence |
| --- | --- | --- |
| `cloud-baseline` (A) | `67ce78e` | section 7a |
| `cross-platform-guards` (A2, ADR-0032) | `f53bf0c`, `dd821be`, `2279951` | ledger P-11/P-17 "cross-platform"; `tests/test_capture_identity.py` |
| `cli-contract`, `cli-stability-multiphase` (B0-B3, ADR-0033) | `8620fa3`, `a14c5f0` | 11 pre-change CLI invocations byte-identical; `tests/test_cli_stability_multiphase.py` |
| `pcsaft-temperature-derivative` (C1 + C3, ADR-0034) | `7c0c4ef` | ledger P-19 (teqp `Ar10` <= 5.24e-16; Gibbs-Helmholtz; `dH_vap = T dS_vap` to 1.05e-12) |
| `release` 0.3.0b1 | `64130a2d5373e88cc65c28fdd047807f04a5daf4` | below |

**CI (GitHub Actions ubuntu-latest):** red on every earlier commit of this
branch; **first green run at `a14c5f0`** (run 36125595338), green at
`7c0c4ef` (36127206093) and at the release commit `64130a2` (36127239537).
Run 36124342571 at `f53bf0c` failed one float (Tessier near-plait phase
fraction 2.34e-12 vs a 1e-12 bound) and led to the derived `1e-12 / delta`
fraction bound (ADR-0032).

**Release 0.3.0b1: gates passed, tag not pushed.** Clean clone from GitHub
at `64130a2`, Linux x86_64, CPython 3.11.15, numpy 2.4.6: ruff/pyright clean;
`pytest -q` 814 passed, 51 skipped, 74 deselected; `python -m build`; wheel in
a fresh venv outside the tree: `0.3.0b1`, `release_smoke.py` pass, new CLI
commands and `residual_properties` run. With the validation extras (teqp
0.23.2, FeOs 0.10.1, thermo 0.6.1) on the same commit: 975 passed, 0 failed.
`git push origin refs/tags/v0.3.0b1` was refused by this session's git route
(three attempts, "remote end hung up", no proxy failure logged; branch pushes
work), and there is no release-creation tool, so per the boundaries in
section 6 the owner-side commands, notes and checksums are in
`.agents/handoffs/release-packet-v0.3.0b1.md`. Until then, pin by SHA:
`pip install "chemthermo @ git+https://github.com/AhmadAlkadri/Chemical-Thermodynamics.git@64130a2d5373e88cc65c28fdd047807f04a5daf4"`.

**Not run this session:** anything on macOS (the ADR-0032 guards are exact
there by construction, but nobody has run them there yet: the owner's next
local `pytest -q` is that check); `pytest -m slow`; the full 2505-state map
(no solver change since `f852726`); benchmark timings; `aglint` (private).

**Next session starts at** the handoff commit that follows `64130a2` on
`dev/sprint`. Plan order (continuation-plan.md): owner publishes v0.3.0b1 ->
C2 (association temperature derivative; completes the residual-caloric
milestone, then a minor release) -> A3/A4 (stability tie-break by residual;
polymer ladder neighbourhood) -> D -> E.


### 7c. Local publication (owner machine, 2026-09-25)

The owner-side Claude Code session published v0.3.0b1 from a local machine,
following `.agents/handoffs/release-packet-v0.3.0b1.md`.

**Preflight:** `origin/dev/sprint` was at `d4fc0e1`; `64130a2` is an ancestor
of it; no `v0.3.0b1*` tag existed on the remote; `gh` was authenticated with
`repo` scope and admin permission.

**macOS arm64 check at `64130a2` (the capture platform of ADR-0032): pass.**
Fresh clone from GitHub, checked out at
`64130a2d5373e88cc65c28fdd047807f04a5daf4`, `python3.11 -m venv` +
`pip install -e ".[dev]"`. Environment: macOS 27.0 (Darwin), `uname -m` =
arm64, CPython 3.11.6 (conda-forge build), numpy 2.4.6. In this environment
`tests/_capture_identity.on_capture_platform()` returns `True`, so the guards
ran in their exact mode.

| check | result |
| --- | --- |
| ruff format --check / ruff check (src, tests) | 140 files already formatted / all checks passed |
| pyright | 0 errors, 0 warnings |
| `pytest -q` | **814 passed**, 51 skipped, 74 deselected, 0 failed (2:57) |
| `tests/test_flash_refactor_bit_identity.py` | 1 passed |
| `tests/test_stability_eos_surfaces.py::test_the_peng_robinson_grid_still_reaches_a_verdict_everywhere` | 1 passed |
| `tests/test_pcsaft_association.py` | 117 passed |
| `tests/test_pcsaft_polymer.py` | 34 passed, 18 deselected |

The counts match the Linux clean-clone run in 7b exactly, and the
capture-platform exact checks pass.

**Build and wheel smoke (same clone):** `python -m build --outdir dist/`;
the wheel was installed into a fresh venv outside the tree. From outside the
tree, `tools/release_smoke.py --expect-version 0.3.0b1` passed (PR single
phase, PR VLE, PR three liquids, PC-SAFT VLE), and
`chemthermo stability-tp --components Methane,n-Hexane --z 0.5,0.5
--temperature-k 300 --pressure-pa 2e6 --eos pc-saft --format json` exited 0
with `status = unstable` (minimizing trial `wilson-vapor`, tpd -1.293).

**Tag:** annotated `v0.3.0b1`, pushed as `refs/tags/v0.3.0b1`;
`git ls-remote origin 'refs/tags/v0.3.0b1^{}'` peels to
`64130a2d5373e88cc65c28fdd047807f04a5daf4`.

**Release:** GitHub prerelease "chemthermo 0.3.0b1",
https://github.com/AhmadAlkadri/Chemical-Thermodynamics/releases/tag/v0.3.0b1
(`isPrerelease: true`). The notes are the packet from "## chemthermo 0.3.0b1"
on, with the macOS row added to the validation table, the artifact checksums
replaced, and the "nobody has run the ADR-0032 guards on macOS arm64"
limitation removed because it no longer holds. Attached assets, whose SHA-256
matches the GitHub-reported digests:

```
a931ac7166dd7637eeb6484f345e7930378db3c623c1d067438edf56ba0cb758  chemthermo-0.3.0b1-py3-none-any.whl  (292924 bytes)
ce59c250ee049787326e5c125ad24cdf2e9a4947746d6ec6148cd96c647f5ca0  chemthermo-0.3.0b1.tar.gz  (490832 bytes)
```

These differ from the cloud-build checksums in the packet because rebuilt
artifacts are not byte-identical across machines. The attached files are the
reference.

**Failed:** nothing. **Skipped:** the validation extras,
`pytest -m slow`, the robustness map, and benchmark timings were not rerun
on macOS for this release. PyPI publication is out of scope.

### 7d. Cloud session 2 (2026-09-25): C2, A3, A4, and release readiness for 0.4.0

Continued from `a916370` (after 7c). Slices, in order:
`pcsaft-association-temperature-derivative` (C2, ADR-0034 amendment, ledger
P-20), `stability-tie-break` (A3, ADR-0035, S-9), `flash-trivial-split` (A4,
ADR-0036, P-17 "ADR-0036"), `docs` (plan item D re-measured: 2:57 on macOS is
inside budget), then - **after the owner changed the goal to release
readiness** - `release-tidy` (ADR-0037: `chemthermo.vlle` removed; stray
`AGENT/`, `AGENT.md`, `ROADMAP.md` removed; `CONTRIBUTING.md`), `public-docs`
(README 1585 -> ~190 lines, detail in `docs/`), `packaging` (PyPI metadata,
`py.typed`, library-only sdist, CI matrix 3.11-3.13 + `twine check --strict`),
`release-policy` (ADR-0038: owner publishes to PyPI by hand),
`robustness-record` (`benchmarks/robustness_761fd57.*`, first full map on
Linux, identical state by state to macOS `f852726`), and `release` (0.4.0 at
`4aa475385791b29f7da76e68332ebfdd39870f1e`).

Owner decisions (interview, 2026-09-25): 0.4.0 final, Beta; short README +
`docs/`; keep `.agents/` and `benchmarks/`, tidy; remove `chemthermo.vlle`
only (gamma-phi stays for the CLI v1 contract). The open-ended sprint is
over; the target is `main` + PyPI.

**Release state:** gates passed on a clean clone of `4aa4753` (Linux; wheel on
3.11/3.12/3.13); PR https://github.com/AhmadAlkadri/Chemical-Thermodynamics/pull/1
open (`dev/sprint` -> `main`, fast-forward possible). Tag, GitHub release,
macOS check, merge and PyPI upload are the owner's:
`.agents/handoffs/release-packet-v0.4.0.md` has the exact steps and notes.

**Next session starts** only after the owner has worked the packet; its
first action is to record the outcome (7e) if the owner-side agent did not.
There is no open-ended development queue: plan items E (deferred science) and
C4 (`d ln phi / dT`, needs a consumer) wait for an owner decision after 0.4.0.

