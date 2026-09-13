# ADR-0015: PC-SAFT density roots, and PC-SAFT in stability_tp / flash_tp

Status: accepted
Date: 2026-09-13

## Context
ADR-0014 shipped a real, teqp-validated non-associating PC-SAFT, but only at a
state fixed by `(T, rho, x)`. It deliberately left `PCSAFTEOS` outside the
`EquationOfState` hierarchy, so the model computed published thermodynamics and
could not be used for the one thing the package exists to do:

```python
ct.flash_tp(mixture, temperature_K=..., pressure_Pa=..., eos=PCSAFTEOS())  # TypeError
```

Two things were missing. First, `EquationOfState` fixes the state with a
*pressure*, and PC-SAFT is Helmholtz-explicit: `(T, P, x)` is several states
until a density root is chosen, and unlike a cubic the roots are not the zeros
of a polynomial. Second, someone has to decide which root a phase label means.

ADR-0014's "Next slice" note sketched a different shape: PC-SAFT entering
`stability_tp` as additional `_PhaseCandidate`s, one per density root, which is
the case ADR-0010's candidate abstraction was designed for. This ADR does not
take that route, and says why below.

## Decision

1. **Solve the density roots in an internal module,
   `chemthermo.eos._pcsaft_density`, working in the packing fraction.**
   At fixed `(T, x)` every quantity the model needs is a function of
   `eta = zeta_3` alone, because all four `zeta` moments are proportional to
   the density. That reduction does three jobs at once:
   - the whole scan grid can be evaluated in one vectorized call, so a root
     solve costs about as much as a handful of point evaluations and is
     affordable inside a flash loop;
   - `Z = 1 + eta a'(eta)` and
     `(dP/drho)/RT = 1 + 2 eta a'(eta) + eta^2 a''(eta)` are *ordinary*
     derivatives, so `dP/drho` is analytic;
   - the scan bounds are physical (`0 < eta < 0.7405`, hard-sphere close
     packing) rather than arbitrary densities.

   The algorithm: bracket sign changes of `P_model(eta) - P` on a
   deterministic grid (120 geometric points from `1e-14` to `1e-3`, then a
   uniform step of `5e-4` to `0.7405`), refine each bracket with safeguarded
   Newton, discard every root with `dP/drho <= 0`, return the rest sorted by
   density, and report how many brackets were found so the count of rejected
   spinodal roots is visible. No root at all raises `ModelError` naming the
   state. No scipy; deterministic by construction.

2. **The seam is the existing `EquationOfState.fugacity_coefficients`, not a
   new density-root candidate type.** `PCSAFTEOS` now subclasses
   `chemthermo.models.EquationOfState` while keeping `EOSProtocol`, and
   implements

   ```
   fugacity_coefficients(mixture=, temperature_K=, pressure_Pa=,
                         composition=, phase=)
   ```

   with the Peng-Robinson label semantics **verbatim**: `"vapor"` is the
   lowest-density (largest `Z`) admissible root, `"liquid"` the highest, and
   when there is one root both labels return the same values.

   Why this seam and not the candidate one:
   - `_EOSTangentPlane` already holds the two branches as competing
     `_CubicRootCandidate`s and keeps whichever minimizes
     `sum_i w_i ln phi_i(w)` (ADR-0005, ADR-0012). That rule is *exactly*
     right for PC-SAFT's roots too - the same model on two branches, one of
     which may not exist at a given `w` - so a PC-SAFT-specific candidate would
     be the same class with a different `ln phi` source.
   - Routing PC-SAFT through the label interface means `stability_tp`,
     `flash_tp`, `_detect.py`, `_split.py`, `_verify.py` and the CLI are
     **unchanged by this slice**. Nothing that Peng-Robinson relies on was
     touched, so `tests/test_flash_refactor_bit_identity.py` and every pinned
     PR / NRTL / modified-Raoult / VLLE number is untouched by construction,
     not merely by measurement.
   - It keeps the one place that knows about density roots inside the model
     that has them.

   What a density-root candidate evaluator would look like, if a later slice
   wants one: `_EOSTangentPlane` would take a model that can enumerate its
   roots at `(T, P, w)` and build one `_PhaseCandidate` per root instead of two
   fixed labels, each carrying its density. `_select_min_gibbs` would need no
   change at all - it already minimizes over an arbitrary-length candidate
   sequence. The gain would be (a) a third root, if a model ever has one that
   is not the spinodal branch, and (b) the density travelling with the phase
   into `FlashResult`, which is what would let `PhaseResult.properties` carry a
   molar volume. The cost is a second model-facing protocol (`density_roots`)
   that only PC-SAFT would implement today. That is a real slice with a real
   user-visible payoff (phase densities in the result), and it should be taken
   when that payoff is the point - not as a side effect of wiring PC-SAFT in.

3. **`components` becomes optional on `PCSAFTEOS`.** `PCSAFTEOS()` takes its
   component order from the `Mixture` it is handed, which is what makes
   `eos=PCSAFTEOS()` read the same way as `eos=PengRobinsonEOS()`. An instance
   built *with* `components` still works and is now *checked* against the
   mixture (after `normalize_name`); a mismatch raises rather than silently
   reordering the caller's composition. The `(T, rho, x)` methods still require
   `components`, since no `Mixture` reaches them - so the empty-components
   error moved from construction to the point of use.

4. **Two new public helpers**, `PCSAFTEOS.density_roots(...)` returning the
   admissible molar densities ascending and `PCSAFTEOS.molar_volume(...,
   phase=)`. `density_roots` is what makes a phase density reportable today,
   and it is what the validation script compares against teqp's `rhoL`/`rhoV`.

5. **`chemthermo.__all__` gains `PCSAFTEOS`** (`PCSAFTParameters`,
   `PCSAFTRecord` and `PCSAFTParameterError` were already there from
   ADR-0014). Per ADR-0001 that makes it public API.

6. **Nothing about the Wilson trial estimates changes.** `_EOSTangentPlane`
   builds `wilson-vapor` / `wilson-liquid` starts from the databank's Tc, Pc
   and omega. Those are *starting points* for a tangent-plane trial, not a
   model statement; they stay valid for PC-SAFT, and the validation below
   (seven tie lines to `<= 2e-9` in mole fraction) is the evidence that the
   trial set finds the right stationary points for this model too.

## Alternatives considered
- **A `_PCSAFTRootCandidate` inside `_evaluator.py`.** Deferred, see decision 2.
- **Caching density roots on the `PCSAFTEOS` instance.** The tangent-plane
  evaluator asks the vapour and the liquid candidate for `ln phi` at the same
  `(T, P, w)`, so exactly half of the root solves are repeats. Rejected for
  now: a full flash costs 0.25-0.75 s without it, which is not a problem the
  package has, and a mutable cache on a frozen dataclass is a real invariant
  to defend for a constant-factor win. If PC-SAFT ever enters a sweep, this is
  the first thing to do.
- **Newton from a good initial guess instead of a grid scan.** Rejected: an
  initial guess cannot tell you *how many* roots there are, and the number of
  roots is the question - it is what separates "this is a two-phase candidate"
  from "the label is a convention". The scan also makes the failure mode
  explicit and bounded (see Consequences).
- **Finite-difference `dP/drho`.** Rejected: that derivative decides which
  roots are returned, and it is differenced through zero exactly where it
  matters (at the spinodal).
- **Reusing `chemthermo.eos.pcsaft._evaluate` point by point for the scan.**
  Rejected on cost (about 1e3 times slower for a scan) and because it does not
  supply `a''`. The duplication is paid for by
  `tests/test_pcsaft_density.py`, which pins the two derivations against each
  other to 1e-12 relative.

## Consequences
- Positive: `stability_tp(..., eos=PCSAFTEOS())` and
  `flash_tp(..., eos=PCSAFTEOS())` work for non-associating mixtures, with the
  tangent-plane phase detection, the seeded split, the three verification
  residuals and the post-split stability test all unchanged.
- Positive (validated, Cases P-3 to P-5, against teqp): n-hexane's two
  saturation densities at 300 K to `2.2e-16` / `1.3e-12` relative; seven
  methane / n-hexane tie lines at 300 K to `|dx1| <= 2.0e-9` and
  `|dy1| <= 3.7e-12` with phase densities to `1.3e-9` relative; teqp's own
  fugacity coefficients at chemthermo's converged phases give equal fugacities
  to `3.8e-9` relative; bubble pressures found by bisecting `stability_tp`'s
  *verdict* agree with teqp's `mix_VLE_Tx` to `1.4e-8` relative.
- Tradeoff: **two roots closer together than the grid step in `eta` are not
  resolved**, and the state is then reported with one root fewer. That happens
  only where the isotherm is nearly tangent to the target pressure, i.e. at a
  near-critical or near-spinodal state, where the missed phase is physically
  marginal. The step is a named constant (`ETA_UNIFORM_STEP = 5e-4`) so the
  trade is explicit.
- Tradeoff: the *residual* `|P_model - P|/P` at a returned root is limited by
  cancellation, not by the iteration. On a dense liquid branch at a low
  pressure `Z` is a near-total cancellation (`1.17e-3` at n-hexane's 300 K
  saturation state), so a `1e-14` relative error in `a'` is `1e-11` in `P`.
  Measured: `2.2e-14` at 10 MPa, `6.5e-12` at 21.9 kPa. The *density* is
  unaffected and is what the tests pin;
  `DensityRoots.max_relative_residual` reports the number rather than hiding
  it.
- Tradeoff: a phi-phi flash still stops at **two phases** (ADR-0009,
  ADR-0011). PC-SAFT joining the EOS family does not change that; phase
  addition lives on the modified-Raoult path.
- Tradeoff: the phase densities are computable
  (`PCSAFTEOS.density_roots(...)` at the converged composition) but are not
  carried in `FlashResult`. See decision 2 for the slice that would fix it.
- Tradeoff: the hard-chain and dispersion terms now appear twice in the
  package, once differentiated with respect to the density and composition
  (`pcsaft.py`) and once with respect to `eta` (`_pcsaft_density.py`). The
  complexity receipt is in that module's docstring and the cross-check is a
  test.
- Unchanged: `kij` is the ADR-0006 contract; association and polar terms are
  still out of scope; there is still no temperature derivative.

## Next slice
Recommended: **`pcsaft-association`**. Reason in the report and in `brain.md`;
briefly, `flash-phase-addition-eos` would generalize machinery that no packaged
EOS can currently exercise, while association is what the model needs before it
can be pointed at the alcohol/water systems the modified-Raoult path already
covers - and it is the only remaining gap that makes a *named* class of fluids
unusable rather than making an existing answer nicer.

## Supersedes (optional)
Amends ADR-0014 decision 2 and its "Pressure-based methods" alternative: the
deferral is discharged, and `PCSAFTEOS` now implements both state
specifications. ADR-0014 is otherwise unchanged.

## Superseded by (optional)
None.
