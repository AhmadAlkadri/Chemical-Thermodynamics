# ADR-0013: `chemthermo.vlle` plugin boundary deprecation

Status: accepted
Date: 2026-09-13

## Context
`chemthermo.vlle` (`get_vlle_engine`, `VLLEEngine`, `VLLEResult`/`VLLEPhase`/
`VLLEInputs`, `VLLEError` and friends) is a public, documented subpackage
(ADR-0001) added before any equilibrium engine existed in this repository. Its
contract is: install the optional `chemthermo_vlle` package, which must expose
`get_engine()` returning something implementing `VLLEEngine.solve(...)`, and
`flash_tp(..., flash_mode="vlle")` raises `ModelError` telling the caller to do
exactly that.

Every premise of that contract is now stale:
- ADR-0008 through ADR-0012 built a tangent-plane stability test, phase
  candidates for a modified-Raoult (activity-liquid / ideal-vapor) pair, and a
  multiphase Rachford-Rice with phase addition and removal. The result
  (ADR-0011) is that `flash_tp(..., flash_mode="modified-raoult")` with
  `FlashSettings(max_phases=3)` (the default) *discovers* vapor-liquid-liquid
  equilibrium - verified on a published ternary tie-triangle (validation Cases
  V-1..V-5) - with no external engine.
- The plugin boundary's own `VLLEEngine` protocol duplicates that surface with
  a *different* shape (one `solve()` call taking a bare `model: str | None`,
  returning a flat `VLLEResult` with no stability diagnostics, no candidate
  labels, no `phase_set_history`) that nothing in this repository implements
  or tests end to end.
- The one package that names itself as the plugin, the private sibling
  `chemthermo_vlle` (`~/Documents/GitHub/chemthermo_vlle` from this
  orchestrator's seat; not part of this repository or its CI), is not that
  engine. It is a 164-line scaffold: an `EOSProtocol` adapter, a
  `Polymer` model stub, JSON test fixtures, and a `solve_vlle()` that
  unconditionally returns `VLLEOutput(status="not_implemented", phases=None)`.
  It contains no phase-equilibrium physics. Its `pyproject.toml` pins
  `chemthermo @ git+...@eb9cfda7...`, a commit from January 2026, predating
  every stability/flash ADR above (0005 through 0012 are all dated
  2026-09-13). Its own `WORKPACKAGES.md` lists "implement VLLE solver" as
  future, not-yet-started work, gated on upstream PC-SAFT scaffolding that
  also does not exist yet in this repository (`chemthermo.eos.pcsaft` is
  itself an unvalidated placeholder - see brain.md roadmap, `pcsaft-residual-
  helmholtz`).

So today, a user who calls `chemthermo.vlle.get_vlle_engine()` and has not
installed `chemthermo_vlle` gets `VLLEPluginNotInstalledError`; a user who
*has* installed it (at its pinned, five-plus-ADR-stale commit) gets an engine
whose `solve()` always reports `"not_implemented"`. Neither path reaches
working VLLE. Meanwhile `flash_tp(..., flash_mode="vlle")` sends every caller
toward that dead end instead of the capability that has existed, verified,
since ADR-0011. That is actively misleading, and it is a public API question
(ADR-0001), so it needs its own decision rather than a silent code change -
which is exactly what ADR-0011 "What remains" deferred to this ADR.

## Decision

### 1. Deprecate `chemthermo.vlle`, do not delete it now
`chemthermo.vlle` and its exported names (`VLLEEngine`, `VLLEError`,
`VLLEInputs`, `VLLEPhase`, `VLLEResult`, `VLLEPluginError`,
`VLLEPluginNotInstalledError`, `get_vlle_engine`) remain importable for one
deprecation cycle, with their behavior otherwise unchanged (`get_vlle_engine`
still raises `VLLEPluginNotInstalledError` when `chemthermo_vlle` is absent,
`VLLEPluginError` when it is malformed, and returns whatever engine it
supplies otherwise). What changes:
- Importing `chemthermo.vlle` emits `warnings.warn(..., DeprecationWarning,
  stacklevel=2)` at module import time, naming the in-tree replacement.
- Calling `get_vlle_engine()` emits the same category of warning, independent
  of the import-time one (a caller may already hold a reference to the
  function from an earlier import), so the warning is visible at the point of
  use even in a long-running process.
- `flash_tp(..., flash_mode="vlle")` still raises `ModelError` (that string is
  not a member of `FLASH_MODES` and was never meant to select a mode), but the
  message no longer tells the caller to install `chemthermo_vlle`. It says:
  three-phase equilibrium is discovered automatically by `flash_mode=
  "modified-raoult"` with `FlashSettings(max_phases=3)` (the default), points
  to ADR-0011, and states that `chemthermo.vlle` itself is deprecated
  (ADR-0013).

Removal is **not** decided here. It follows in a later ADR, once no internal
user remains (none does today - nothing in `src/chemthermo` imports
`chemthermo.vlle`) and after the one-deprecation-cycle window has been given
to any external user of the boundary.

### 2. The boundary does not become "the way to reach a different engine"
ADR-0011 posed the open question directly: does `chemthermo.vlle` become the
supported route to a *different* solver (a user's own multiphase code,
something PC-SAFT-backed), or is it deprecated? Decision: deprecated, not
repurposed. Reasons:
- A second public multiphase-equilibrium contract, next to the one ADR-0008
  through ADR-0012 just finished validating, is exactly the kind of duplicated
  boundary ADR-0002's thin-vertical-slice rule exists to prevent - it is
  infrastructure with no exercised user (the same reasoning ADR-0011 "What
  remains" already applied to phi-phi/gamma-gamma phase addition).
- `VLLEEngine.solve()`'s shape (a single opaque call, a flat result, a free-
  form `options: Mapping[str, Scalar]`) cannot carry what the in-tree path now
  reports - candidate labels, `phase_set_history`, post-split stability,
  `equilibrium_residual` - without becoming a second, parallel definition of
  the same information that can silently drift from the first.
- Nothing in this repository, or in the one package that names itself as an
  implementer, currently satisfies `VLLEEngine`. Repurposing a protocol that
  has zero implementations is not less speculative than deprecating it.

### 3. Disposition of the private sibling `chemthermo_vlle`
Recorded here **as a recommendation to the maintainer**, not as an action
taken in this repository or this change: `chemthermo_vlle` should be archived
or deprecated rather than developed further, for four reasons -
(a) it contains no phase-equilibrium physics: `solve_vlle()` is a fixed
`status="not_implemented"` stub; (b) it is pinned to a commit (`eb9cfda`,
January 2026) that predates the entire equilibrium engine this repository now
has - every stability and flash ADR from 0005 onward postdates that pin;
(c) its `VLLEEngine`-shaped boundary duplicates what is now public and
verified in-tree (decision 2 above); and (d) the capability it was scaffolded
toward - polymer/solvent VLLE on a private high-performance PC-SAFT backend -
does not need a separate "solve everything" protocol to exist. It needs the
public contracts this repository already publishes: the `EquationOfState`
protocol and phase-candidate mechanism (ADR-0007, ADR-0010), the tangent-plane
stability test (ADR-0005, ADR-0007), and the multiphase flash with addition
and removal (ADR-0011, ADR-0012). A private engine consuming *those* - adding
its own PC-SAFT parameter packs and a fast residual-Helmholtz implementation
behind the existing `EquationOfState`/`_TangentPlaneEvaluator` seams - gets
the validated stability/flash machinery for free and stays a thin, reviewable
adapter instead of a second solver.

What such a private engine could still legitimately differentiate on, so the
boundary stays clear: optimized numerical kernels (vectorized/compiled
residual-Helmholtz and its derivatives for a specific EOS family); non-public
parameter packs (proprietary or licensed pure-component and binary-interaction
data, e.g. for polymers); initialization and continuation strategies for
hard multi-phase regions (arc-length continuation along a binodal, warm starts
across a state grid); and batch/grid infrastructure (solving many `(T, P, z)`
states efficiently, e.g. for a phase diagram or a process simulation inner
loop). None of that requires re-declaring `VLLEEngine`, `VLLEResult` or a
parallel phase-equilibrium contract - it is downstream of the public
`EquationOfState` and `flash_tp`/`stability_tp` contracts, which is where
ADR-0011 already said future engines (PC-SAFT-backed or otherwise) should
enter: "it enters as further `_PhaseCandidate`s ... with no solver change."

## Alternatives considered
- **Delete `chemthermo.vlle` now.** Rejected: it is public API (ADR-0001); an
  unannounced removal breaks any external caller silently. A deprecation
  cycle is the same choice already made for `flash_mode="gamma-phi"`
  (ADR-0010) and the legacy `wilson-heuristic` phase detection (ADR-0008).
- **Repurpose `VLLEEngine` as the supported multi-engine seam** (decision 2's
  rejected option, expanded). Rejected on both the shape mismatch and the
  zero-implementations point above; also rejected because it would invite a
  second "how do I get VLLE" answer to sit next to the one this repository
  just finished validating, contradicting the entire point of ADR-0011.
  Considered narrower forms - e.g. `VLLEEngine` as *only* a PC-SAFT hook -
  and rejected those too: PC-SAFT's public entry point is already decided
  to be an `EquationOfState`/`_PhaseCandidate` (ADR-0010's mechanism,
  ADR-0011 "What remains"), not a new protocol.
- **Silence the warning by default (e.g. only on `-W`).** Rejected:
  `DeprecationWarning` is the standard category for exactly this, other
  deprecations in this repository (`gamma-phi`, `wilson-heuristic`) are
  recorded in docs/docstrings without a runtime warning, and adding the
  runtime signal here is closer to Python convention for something that was
  never actually functional (unlike `gamma-phi`, which does work, just
  inconsistently).
- **Leave the `flash_mode="vlle"` message untouched.** Rejected: it now tells
  every caller to install a package that cannot answer their question even
  when installed (decision above), which is worse than deprecating quietly -
  it actively wastes the caller's time chasing a dead end that ADR-0011
  already resolved a different way.

## Consequences
- Positive: a user who calls `chemthermo.vlle.get_vlle_engine()` or
  `flash_tp(..., flash_mode="vlle")` is now pointed at the capability that
  actually exists and is validated (`flash_mode="modified-raoult"`,
  `FlashSettings(max_phases=3)`, ADR-0011), instead of an out-of-tree package
  that cannot deliver it.
- Positive: no compatibility break. Every exported name, `VLLEPluginError`
  hierarchy member, and the loader's existing `ModuleNotFoundError`/
  `VLLEPluginError` behavior are unchanged; only warnings are added and one
  error string changed.
- Neutral: the private sibling `chemthermo_vlle` is not modified by this ADR
  (it is outside this repository) - decision 3 is a recommendation, not an
  action. If the maintainer instead extends it, its pin will need to move past
  `eb9cfda` to pick up ADR-0008 through ADR-0012 before it can consume any of
  the mechanism decision 3 recommends it build on.
- Tradeoff: `chemthermo.vlle` stays in the package for at least one more
  release, carrying dead-end plumbing (`loader.py`'s `import_module
  ("chemthermo_vlle")` call) that nothing in this repository's CI exercises
  against a real engine. Accepted as the cost of a compatibility-respecting
  deprecation; removal is the next ADR once the cycle has run.

## Supersedes (optional)
None. Discharges the "vlle" bullet of ADR-0011 "What remains" (the promised
follow-up, `vlle-plugin-boundary-disposition`) and amends ADR-0001's listing
of `chemthermo.vlle` as a "documented public subpackage" to "documented
public, deprecated".

## Superseded by (optional)
None.
