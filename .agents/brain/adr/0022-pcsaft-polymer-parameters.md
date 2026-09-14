# ADR-0022: A polymer as a PC-SAFT component

Status: accepted
Date: 2026-09-13

## Context
Every PC-SAFT component this package could describe was a small molecule. Three
separate things stopped a polymer from being one, and none of them was the
equation of state:

1. **The parameter record.** Polymer PC-SAFT parameters are published per unit
   *mass* - `m/M` in mol/g - because the chain length depends on the molar mass
   of the particular sample, so one parameter set covers every molar mass of
   that polymer. `PCSAFTRecord` had only `m`, so a caller had to multiply by
   hand and the record could not say what it was.
2. **The component.** `Component` could only be built by databank lookup, and
   `ComponentData` requires `Tc`, `Pc` and `omega`. A polymer has no critical
   point. Supplying one to satisfy the schema would put a fabricated number
   where the Wilson K-value estimate would read it.
3. **The exponential.** `EquationOfState.fugacity_coefficients` returns
   `phi = exp(ln phi)`, and a chain of `m = 1393.9` segments (polyethylene,
   `Mw = 53000`) in n-pentane at 453 K and 10 MPa has `ln phi = -1690.6`.
   `exp` of that is an exact `0.0`, so the model's finite answer was destroyed
   by the last step of the accessor and every consumer correctly refused a
   non-positive fugacity coefficient. `stability_tp` raised "No usable
   fugacity-coefficient branch".

`brain.md` listed `pcsaft-polymer-parameters` as the recommended next slice
precisely because it is the last remaining *capability* gap that changes a
public data contract, and a contract change is cheaper before the performance
work of Stage I than after.

## Decision

1. **`PCSAFTRecord` gains `segments_per_g` (mol/g), an alternative spelling of
   `m`.** Given together with `MW_g_mol` (g/mol) it derives
   `m = segments_per_g * MW_g_mol` once, in `__post_init__`. **Exactly one** of
   `m` and `segments_per_g` may be given; both, or neither, raises, because the
   two say the same thing in different units and a record carrying disagreeing
   versions of one parameter is a data defect rather than a choice.
   `MW_g_mol` is promoted from annotation to load-bearing input *only* on this
   path, and it stays the single molar mass on the record - the brief's
   `molar_mass_g_mol` is spelled `MW_g_mol` here deliberately, because a second
   molar-mass field on the same record is a field that can disagree with the
   first. `from_records` accepts either spelling from a mapping. Every existing
   record is unchanged, packaged set included.

2. **`Component.custom(name, *, mw_kg_per_mol, formula=, tc_k=, pc_pa=,
   omega=, volatile=, antoine=, source=)`** builds a component with no databank
   lookup. `mw_kg_per_mol` is required (it is what converts a mass fraction to
   a mole fraction and what a segments-per-mass record multiplies); the three
   critical constants are optional, and `Component.tc_k` / `pc_pa` / `omega`
   raise `PropertyNotFoundError` when they were not given, so a model that
   needs one says so instead of reading a placeholder.

   **The databank schema is not relaxed.** `ComponentData` still *requires*
   `Tc` / `Pc` / `omega` and `Database` still validates exactly as before, so
   `schema_version` stays 1. The relaxed model is a sibling,
   `CustomComponentData`, sharing a new `ComponentCore` base and deliberately
   **not** a `ComponentData`, so it cannot enter `Database`.

3. **`volatile=False` changes exactly one number: the Wilson K-value
   *estimate*.** `chemthermo.flash._common.NON_VOLATILE_WILSON_K = 1e-10`
   replaces the correlation for such a component - "essentially absent from the
   vapour-like trial". The Wilson correlation is written in `Tc`, `Pc` and
   `omega`; a polymer has none of the three, and without a substitute the
   deterministic trial set of `_EOSTangentPlane` could not be built at all.
   `z K` then gives the vapour-like start (polymer at ~1e-13) and `z / K` the
   polymer-dominant liquid-like one, which is exactly the pair a
   polymer/solvent feed needs. It is an initial **estimate**: the trial
   iterates on the model from there, and no verdict, composition or fugacity is
   a function of this number. Nothing else in the package reads `volatile`.
   Every packaged databank component is volatile, so `wilson_k` is unchanged
   for every mixture that could be built before this ADR (pinned with `==`).

4. **A log-space route for fugacity coefficients, as a guard and not as a
   second model.** `EquationOfState.log_fugacity_coefficients(...)` is new,
   optional and returns `None` by default ("this model cannot say"), which
   makes the guard inert and the pre-ADR-0022 failure the outcome for any model
   that does not implement it. `PCSAFTEOS` implements it as
   `fugacity_coefficients` without the final `exp`.

   `chemthermo.flash._common.eos_branch_terms` is the single place that decides:
   ask the model for `phi`, and when every entry is finite and positive return
   `(phi, np.log(phi))` - the doubles every caller used before. Only where that
   fails does it ask for the logarithms. `_PhaseRoot.branch_terms` carries the
   pair, and the phi-phi split's `K` update uses `phi_l / phi_v` exactly as
   before when `phi` exists and `exp(ln phi_l - ln phi_v)` when it does not.

   **Bit-identity is by construction**, the same argument ADR-0016 used for the
   extended Rachford-Rice: the guarded branch is reached only where the
   unguarded expression has no answer at all.

5. **No polymer parameters are packaged, and none should be.** The values used
   by the tests and examples live in
   `tests/fixtures/pcsaft/martini2009_polymers.json` with a provenance block
   that states what they are: **as tabulated by Martini, Cismondi, Barbosa &
   Brignole, *Sep. Sci. Technol.* 44 (2009) (author manuscript, CONICET open
   repository), citing Gross & Sadowski, *IECR* 41 (2002) 1084; not verified
   against that primary table**, which is paywalled and was not read. A search
   for a second open source printing the same three numbers found none: FeOs
   and Clapeyron.jl package no polymer records, and every open hit for the
   values leads back to this one manuscript. That is a weaker footing than the
   packaged 2001 records, which two independent secondary sources confirm digit
   for digit, and the difference is the reason these stay a fixture. The same
   rule NRTL's Tessier (2000) fixture has followed since that slice.

6. **A polymer here is monodisperse.** One chain length, one component. The
   samples the cited `k_ij` values were fitted to have polydispersities of 1.14
   to 2.94. Representing a distribution as several pseudo-components is a
   different capability and is named as a later slice candidate, not delivered.

## Alternatives considered
- **Relax `ComponentData` itself** so `Tc` / `Pc` / `omega` are optional.
  Rejected: the databank's guarantee that every packaged component has critical
  constants is load-bearing for Peng-Robinson and for the Wilson estimates, and
  relaxing the shared model would weaken it for all 16 packaged records in
  order to describe one component the databank does not hold. A sibling class
  costs six repeated field declarations and keeps `schema_version` honest.
- **Give the polymer a fabricated critical point** (extrapolated from an
  n-alkane series, say) so nothing downstream changes. Rejected: it is
  precisely the "do not invent thermodynamic constants" rule of the brain's
  guardrails, and the number would be invisible at the point where the Wilson
  estimate read it.
- **Make the whole phi-phi path work in log space.** Rejected, and this is the
  important one: `exp(ln phi_l - ln phi_v)` and `phi_l / phi_v` are different
  doubles, so a wholesale switch would move every pinned PC-SAFT number and
  require regenerating the 155-state bit-identity fixture for no thermodynamic
  gain. The guard pattern - compute as before, fall back only where the result
  does not exist - keeps the fixture untouched.
- **Clamp `ln phi` into the exponential's range** inside `PCSAFTEOS`. Rejected:
  a clamp there is a wrong answer rather than a missing one, and it would be
  silently active in exactly the states a user cares about.
- **A mole-number floor in the second-order stage**, anticipated in the slice
  design for `n_polymer ~ 1e-300`. Not added: measured, the smallest mole
  number the stage sees on these systems is 1.2e-20 (the `Mw = 53000`
  solvent-rich phase at 5 MPa), which is four orders above the existing
  `1e-300` clamp and 280 orders above underflow. A guard with no measured
  failure behind it is the scaffolding ADR-0002 forbids.
- **Packaging the polymer parameters anyway**, since they are cited. Rejected;
  see decision 5. If the primary table is ever read, or a second independent
  source found, promoting them is a small follow-up with its own ADR note.

## Consequences
- Positive: `flash_tp(mixture, ..., eos=PCSAFTEOS(parameters=...))` returns a
  verified liquid-liquid split for polyethylene / n-pentane at 453 K - two
  liquids, mass balance 2.7e-20, fugacity residual 3.4e-13,
  `delta_g_split_rt = -1.69e-04`, post-split stable - at a segment-number ratio
  of 160 : 1, and the same tie line comes out of an independently written
  two-equation Newton to 1.7e-15. Validation Cases P-12 and P-13.
- Positive: the guard makes `Mw = 53000` (`m = 1393.9`) usable at all. Before
  it, `stability_tp` refused the state; after it, all 1025 branch evaluations
  of that flash run in log space and the answer is verified the same way.
  Instrumented as dormant (0 of 274, 0 of 600) on ordinary states.
- Positive: three new public names - `Component.custom`,
  `PCSAFTRecord.segments_per_g`, `EquationOfState.log_fugacity_coefficients` -
  and one new schema class, `CustomComponentData`.
- Unchanged: the 155-state bit-identity fixture (`refactor_bit_identity_v3`)
  passes without regeneration, and the 144-state Peng-Robinson replay test now
  additionally asserts the guard stayed dormant at every split iterate.
  No packaged parameter, no solver equation and no pinned number moved.
- Tradeoff: **`k_ij` could not be given to FeOs's PC-SAFT** from Python in
  feos 0.10.1 (`EquationOfState.pcsaft` raises "missing field `k_ij`" for every
  serialization tried), so every FeOs comparison in this slice runs at
  `k_ij = 0` on both sides. The fitted `k_ij = -0.006` answer is checked
  against the in-script Newton instead.
- Tradeoff: **FeOs's own `tp_flash` is not a usable reference inside this
  two-phase region.** It raises at 5 and 8 MPa, returns a degenerate pair at
  3 MPa, and converges only at 10 MPa - which is 8 % from the cloud point, so
  both solvers are near a plait point there (chemthermo's equal-fugacity
  residual 1.06e-07, FeOs's 9.07e-06, both measured in chemthermo's model).
  Recorded as a fact about the reference, not as disagreement.
- Tradeoff: **a genuine gap in the phi-phi split is exposed and pinned, not
  patched.** Below about 3 MPa, where n-pentane still has a vapour density root
  at 453 K, the equilibrium is vapour-liquid rather than liquid-liquid;
  `stability_tp` still reports the feed unstable and the split stops after one
  successive-substitution step with an equal-fugacity residual of order 1e+02,
  so `flash_tp` raises `ConvergenceError`. The ternary does the same at 3 MPa
  (`beta = -6.3e+10`). Both are pinned by test.
- Tradeoff: the polymer parameters rest on a single secondary source and are
  never packaged, so a user of this capability must supply their own - which is
  the honest state of the art here, not a deficiency of the code.
- Cost: `pytest -q` goes from 233.2 s for 750 tests to **265.7 s (4:25) for
  796** - about 19 s of new tests and 8 s of the two new example scripts under
  `tests/test_examples.py` - and `pytest -q -m slow` gains about 44 s for six
  more. That is a little over the ~4 min budget; six tests are already marked
  `slow` per the dev-contract rule (each a finer bisection or a repetition of
  something the default run brackets) and the rest is the capability itself, so
  the next lever is Stage I rather than more markers.

## Next slice
`perf-profile-baseline` (Stage I), unchanged from `brain.md`'s recommendation:
this slice was the last public *data* contract change, so a performance
baseline taken now will not have to be retaken. `pcsaft-polydispersity`
(several pseudo-components for one polymer sample) and
`pcsaft-polymer-vle-ethylene` (the vapour-liquid regime this slice's
low-pressure limitation names) are the later candidates.

## Supersedes (optional)
None. Discharges the `pcsaft-polymer-parameters` (Stage H) roadmap item in
`brain.md`.

## Superseded by (optional)
None.
