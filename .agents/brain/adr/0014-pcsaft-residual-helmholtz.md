# ADR-0014: PC-SAFT residual Helmholtz energy behind the existing EOS protocol

Status: accepted
Date: 2026-09-13

## Context
`chemthermo.eos.pcsaft` shipped a `PCSAFTEOS` dataclass whose
`residual_helmholtz` looked up parameters that did not exist and then raised
`NotImplementedError`. `chemthermo.parameters.PCSAFTParameterRegistry` held an
empty mapping and raised "PC-SAFT parameters are not available in the
open-source build". `tests/test_pcsaft_sanity.py` asserted that the status was
`"missing_parameters"`.

That is exactly the scaffolding ADR-0002 forbids: a public, registered,
documented capability (`get_eos("pcsaft", ...)` succeeds, `list_eos()` names
it, `EOSProtocol` describes it) that cannot do anything. It was also the last
part of the package with no cited parameters, no independent validation route
and no entry in the validation-case ledger, while every other model family
(NRTL, Peng-Robinson, the modified-Raoult stack) had acquired all three.

The brain and the steering brief both list `pcsaft-residual-helmholtz` as the
next slice, and ADR-0013 decision 3 says a future PC-SAFT-backed multiphase
engine should consume the public `EquationOfState` / `stability_tp` /
`flash_tp` contracts rather than a parallel "solve everything" protocol. That
only works if PC-SAFT actually computes something first.

## Decision
1. **Implement the real physics behind the existing `EOSProtocol`**, not a new
   interface. `PCSAFTEOS.residual_helmholtz(temperature_K=, volume_m3=,
   composition=)` keeps its signature and now returns `A^res/(R T)`;
   `volume_m3` is documented as the **molar** volume in m^3/mol. The model is
   Gross & Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244, **hard chain plus
   dispersion only**. Association (Gross & Sadowski, IECR 41 (2002) 5510) and
   the polar terms are explicitly out of scope and are named as such in the
   module docstring, the packaged parameter file and the README.
2. **Add three density-based property methods** rather than pressure-based
   ones: `compressibility_factor`, `pressure_Pa` and
   `ln_fugacity_coefficients`, each taking `(temperature_K=,
   density_mol_m3=, composition=)`. PC-SAFT is a Helmholtz-explicit model:
   its natural independent variables are `(T, rho, x)`, and a `(T, P)` state is
   *several* states until a density root is chosen. This slice therefore does
   **no** root solving; the caller says which root they are on. A state inside
   the spinodal has `Z <= 0` and no fugacity coefficients at all, and
   `ln_fugacity_coefficients` raises `ModelError` there rather than returning
   a `nan`.
3. **All derivatives are analytic and assembled by chain rule over a handful
   of intermediates**, not transcribed equation by equation from the paper's
   appendix. Two structural facts carry it: every `zeta_n` is proportional to
   the density, so `rho d/drho` of anything built from the moments is
   `sum_n zeta_n (d/d zeta_n)`; and the dispersion term depends on composition
   only through `eta`, `mbar`, `m2es3` and `m2e2s3`. That gives `Z` and
   `partial a_res / partial x_k` from the gradient of `a_hs`, the gradient of
   `g_ii`, and four partials of one scalar function. Fewer hand-written
   formulas means fewer places to mistype (one was in fact caught this way
   during development: a missing `I2` factor in `partial a_disp / partial
   m2e2s3`, invisible in `A^res` and `Z`, visible immediately in `ln phi` of a
   mixture).
4. **`C1` is implemented in the typo-corrected form.** Equation (A.11) as
   printed drops the outer exponent `-1` from its right-hand side; NIST TRC's
   PC-SAFT page states the erratum explicitly. The reciprocal is implemented,
   the module docstring says so, and `tests/test_pcsaft.py` pins the shape
   (`C1 * bracket == 1`) rather than only a value.
5. **Parameters get the same treatment NRTL got**: `PCSAFTParameters` loads a
   packaged `src/chemthermo/parameters/data/eos/pcsaft.json` (schema_version 1,
   model `"PC-SAFT"`, a `provenance` block, per-component `m`, `sigma_A`,
   `epsilon_k_K`, optional `MW_g_mol`, `source`) with the eleven
   non-associating compounds of Gross & Sadowski (2001) Table 1 whose names
   also exist in the packaged databank. `PCSAFTParameters.from_records(...)`
   takes user-supplied values and `for_components(names)` returns arrays in
   the requested order. The stub `PCSAFTParameterRegistry` is **replaced**;
   `PCSAFTParameterError` is kept.
6. **`kij` uses the Peng-Robinson contract verbatim** (ADR-0006): a scalar
   applied to every off-diagonal pair and never the diagonal, or a mapping
   from an unordered pair of normalized component names to a value. The
   helpers moved from `chemthermo/models/peng_robinson.py` to the internal
   `chemthermo/models/_kij.py` **unchanged**, so Peng-Robinson stays
   bit-identical; only the error-message prefix is now a parameter.
7. **teqp is the independent validation route**, added to the `validation`
   extra. It implements the same published model but takes every derivative by
   automatic differentiation, so agreement on `Z` and `ln phi` genuinely tests
   the hand-written analytic derivatives.

## Alternatives considered
- **Take an automatic-differentiation dependency** (jax / autograd / a
  dual-number library) and write only `A^res`. Rejected as a *foundation*: a
  core thermodynamic kernel in this repository should be readable as the
  published equations and should not put a large numeric framework on the
  install path of `import chemthermo`. The cost is paid in the derivative
  code, and it is paid back by the finite-difference and teqp tests, which
  would be needed anyway. (Nothing here forecloses adding an autodiff-backed
  *checker* later.)
- **Wrap teqp and call it the PC-SAFT implementation.** Rejected: the
  reference path has to be readable in-tree - the same reason the validation
  scripts re-derive their own references instead of importing `thermo`'s
  answer. It would also make an optional binary wheel a hard runtime
  dependency of a core model, and it would leave nothing to cross-check
  against.
- **Pressure-based methods (`fugacity_coefficients(..., pressure_Pa=,
  phase=)`), matching `EquationOfState`.** Deferred, not rejected: that
  signature requires a density root solver *and* a root-selection rule, both
  of which are the next slice's subject. Shipping the signature now with a
  half-finished root finder underneath would be the same scaffolding this ADR
  exists to remove. `PCSAFTEOS` deliberately does **not** subclass
  `EquationOfState` yet.
- **Packaging a `kij` dataset.** Rejected for the same reason ADR-0006 gave:
  there is no cited PC-SAFT `kij` table in hand, and a registry over zero real
  entries is scaffolding. `kij` stays a per-instance argument.
- **Storing `MW` in the PC-SAFT parameter file.** Rejected: the component
  databank is already the single source of truth for molar mass, PC-SAFT does
  not use it, and a second copy would drift. The optional key is supported for
  user-supplied records.

## Consequences
- Positive: `chemthermo.eos` is no longer "a registry entry with implementation
  details evolving over time". `get_eos("pcsaft", components=[...])` returns a
  model that computes published thermodynamics, validated to 1e-10 against
  teqp at fourteen states (achieved: worst `|dA^res/RT|` 4.4e-15, worst `|dZ|`
  2.6e-14, worst `max |d ln phi|` 2.7e-14).
- Positive: the public API gains `PCSAFTParameters` / `PCSAFTRecord` /
  `get_pcsaft_parameters`, and `chemthermo.eos.PCSAFTEOS` gains three methods.
- Breaking: `chemthermo.parameters.PCSAFTParameterRegistry` no longer exists.
  It was never in `chemthermo.__all__` and its only behavior was to raise, so
  the blast radius is a caller who imported it to catch that raise; such a
  caller should catch `PCSAFTParameterError`, which is unchanged and still
  exported. `PCSAFTEOS.residual_helmholtz` no longer raises
  `NotImplementedError` for the eleven packaged compounds.
- Tradeoff: `PCSAFTEOS` is not yet an `EquationOfState` and therefore cannot be
  handed to `stability_tp` or `flash_tp`. That is stated in the README rather
  than papered over.
- Tradeoff: no temperature derivative (`Ar10` in teqp's naming) is implemented,
  so residual enthalpy/entropy are not available. Not needed by the flash
  path; it would be a separate, small slice.
- Cost: one more optional dependency in the `validation` extra (`teqp>=0.23`,
  a binary wheel). It is skipped cleanly when absent, like `thermo`.

## Next slice
`pcsaft-density-roots-flash`: a density solver at `(T, P)` that returns the
vapour-like and liquid-like candidate roots, so PC-SAFT can enter
`stability_tp` as further `_PhaseCandidate`s (the case ADR-0010's candidate
abstraction was designed for and ADR-0013 decision 3 points at), followed by
flash validation against teqp VLE - the methane / n-hexane 300 K isotherm is
already traced as a reference.

## Supersedes (optional)
None. Discharges the "Experimental/placeholder: PC-SAFT residual Helmholtz
implementation" entry of the steering brief and the `pcsaft-residual-helmholtz`
roadmap item in `brain.md`.

## Superseded by (optional)
None.
