# Cloud continuation handoff

Written 2026-09-25 by the local orchestrator when the equilibrium campaign was
first published. It is self-contained: a fresh session needs this repository
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

## 7. Release and fresh verification record

Recorded after publication; see the section appended below.
