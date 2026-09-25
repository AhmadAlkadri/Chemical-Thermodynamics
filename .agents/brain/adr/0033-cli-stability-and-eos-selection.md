# ADR-0033: CLI exposes stability and PC-SAFT; `cli_schema_version` stays 1

Status: accepted
Date: 2026-09-25

## Context
The CLI contract (ADR-0003/0004) is `chemthermo tp-flash` with Peng-Robinson,
`--flash-mode {phi-phi,gamma-phi}`, `cli_schema_version = 1` and exit codes
0/1/2/3. The library has since gained `stability_tp` (ADR-0005), multiphase
`flash_tp` (ADR-0011/0020, `FlashSettings.max_phases`) and `PCSAFTEOS`
(ADR-0014/0018), none reachable from the CLI. Measured before this change:
today's `tp-flash` already returns three liquids for Peng-Robinson water /
ethanol / n-hexane `(0.2, 0.4, 0.4)` at 280 K and 1 atm, because its payload
keys `phases` and `phase_fractions` by phase name and `vapor_fraction` is
`null` when there is no vapour. An N-phase answer therefore needs no new
layout.

## Decision
1. **Schema.** `cli_schema_version` stays **1**. New commands and new flags are
   additive: an invocation that was valid before produces **byte-identical
   stdout** (checked against captures taken before the change, same machine).
   Keys a new flag adds appear only when that flag is given.
2. **`chemthermo stability-tp`** (new). Flags as `tp-flash`
   (`--components`, `--z`, `--temperature-k`, `--pressure-pa`, `--normalize`,
   `--format {text,json}`) plus `--eos {peng-robinson,pc-saft}` (default
   `peng-robinson`). JSON: `cli_schema_version`, `command`, `inputs`,
   `solver` (`eos`, `method = "tangent-plane"`, `settings`: `tol`,
   `tpd_tol`, `trivial_tol`, `max_iter`), `result` (`status`, `stable`,
   `stability_scope = "bounded-trial-set"`, `tpd_min`, `component_order`,
   `trial_composition`, `k_values`, `feed_branch`, `phase_branch`,
   `trial_count`, `converged_trial_count`, `non_trivial_trial_count`) and
   `diagnostics` (the library's mapping, unfiltered, as `tp-flash` does).
3. **"Stable" is worded as bounded.** Every stability payload carries
   `stability_scope = "bounded-trial-set"`, the text output says "no negative
   tangent-plane distance found from the deterministic trial set (not a
   global proof)", and `--help` says the same (ADR-0005 honesty note).
4. **Exit codes, unchanged meaning.** `stability-tp`: `0` for `stable` or
   `unstable`; `3` for `inconclusive` (no trial converged - the analysis did
   not reach a verdict), **with the payload still printed** so the
   diagnostics are visible; `1` for validation/model errors (unknown
   component, missing PC-SAFT parameters); `2` for usage.
5. **`tp-flash --eos {peng-robinson,pc-saft}`** (default `peng-robinson`;
   `solver.eos` becomes `"pc_saft"` only when chosen). PC-SAFT uses the
   packaged parameter table, `kij = 0`, phi-phi only: `--eos pc-saft` with
   `--flash-mode gamma-phi` is a usage error (exit 2).
6. **`tp-flash --max-phases N`** (N >= 1; default: the library's, 3). When
   given, `solver.settings.max_phases` is added to the payload.
7. **Not in this contract:** `modified-raoult` and `gamma-gamma` modes (the
   packaged NRTL table is too small for them to be meaningful from a command
   line), user parameter files, `kij` flags, and PC-SAFT association
   parameters beyond the packaged table. Each is a later additive decision.
8. **Determinism.** JSON is `json.dumps(..., indent=2, sort_keys=True)`; floats
   are Python `repr`; phase names are exactly the library's. A non-finite
   float in a stability payload (e.g. `tm_at_stationary_point = -inf`, ADR-0025)
   is written as the string `"inf"`, `"-inf"` or `"nan"` so stdout stays
   valid JSON.

## Alternatives considered
- A `v2` schema with a phase list (rejected: the v1 name-keyed layout already
  carries N phases; a bump would break consumers for nothing).
- Folding stability into `tp-flash --stability-only` (rejected: a different
  result type under the same `command` value would make `result` ambiguous).
- Exit 0 with `status = "inconclusive"` (rejected: scripts checking `$?`
  would read a non-answer as an answer).
- Always writing `max_phases` into `solver.settings` (rejected: changes the
  bytes of every existing invocation).

## Consequences
- The CLI can reach every equilibrium capability that needs only packaged
  data. Golden fixtures: `tests/fixtures/cli/stability_tp_v1_*.json`,
  `tests/fixtures/cli/tp_flash_v1_*.json` (new states); the ADR-0003/0004
  fixtures are unchanged.
- PC-SAFT from the CLI is limited to the packaged components; the command
  line gives no way to pass a polymer or a `kij`.

## Supersedes (optional)
Extends ADR-0003 and ADR-0004; supersedes neither.

## Superseded by (optional)
None.
