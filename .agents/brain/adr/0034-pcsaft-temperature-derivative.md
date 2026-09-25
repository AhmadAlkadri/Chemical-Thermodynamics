# ADR-0034: PC-SAFT temperature derivative and residual properties (non-associating)

Status: accepted
Date: 2026-09-25

## Context
`PCSAFTEOS` had density and composition derivatives only, so no residual
enthalpy or entropy, no caloric output, and an empty `get_Ar10` column in the
teqp comparison (ledger Case P-1 "Not compared"). Temperature enters the
non-associating model in exactly two places: the segment diameter
`d_i(T)` (Gross & Sadowski 2001, Eq. A.9) and `eps_ij/kT` in the dispersion
mixing rules (A.12-A.13). The association term (2002) adds a third, through
`Delta(T)` and the site fractions.

## Decision
1. **`PCSAFTEOS.residual_helmholtz_temperature_derivative(*, temperature_K,
   volume_m3, composition) -> float`**: `(d(A^res/RT)/dT)_{V,x}` in 1/K,
   analytic, for hard chain + dispersion. Implemented as a separate pure
   function `_temperature_derivative` beside `_evaluate`, which is not
   touched, so every existing number is unchanged by construction.
2. **`PCSAFTEOS.residual_properties(*, temperature_K, density_mol_m3,
   composition) -> dict[str, float]`**: reduced `a_res`, `z`, `u_res`,
   `h_res`, `s_res_tv`, `s_res_tp`, `g_res_tv`, `g_res_tp`. The reference
   state is named in the key: `_tv` = ideal gas at the same `T` and molar
   volume, `_tp` = at the same `T` and pressure; they differ by `ln Z`, and
   `U`, `H` are reference-independent.
3. **Refuse, do not approximate.** An associating mixture raises `ModelError`
   (the association temperature derivative is plan item C2, not done); a
   state with `Z <= 0` raises for `residual_properties` (`ln Z` undefined).
4. **Residual only.** No total `H`, `S`, `Cp`: those need ideal-gas heat
   capacities the databank does not carry. No `Cp^res` (second `T`
   derivative) and no `ln phi` temperature derivatives (plan item C4, only
   with a consumer). Stated in the docstring and README.

## Alternatives considered
- Add `da_dT` to `_PCSAFTState` inside `_evaluate` (rejected for now: would
  put new arithmetic on the path every flash takes and invite a bit-identity
  audit for a quantity no flash uses).
- Finite differences in `T` (rejected: the whole point is an analytic route
  that can be checked against teqp's automatic differentiation).
- Return SI units (rejected: the reduced forms are what every consumer and
  reference tabulates; the conversion is one multiplication and documented).

## Consequences
- Validation (ledger Case P-19): teqp `get_Ar10` at the 14 Case P-1 states,
  asserted 1e-12 relative (achieved <= 5.24e-16 over all 14); fourth-order
  central differences, 1e-8; Gibbs-Helmholtz through the density solver and
  `ln phi` (`H^res/RT = -T d(sum x ln phi)/dT` at fixed `P`), 1e-7; exact
  identities between the two references.
- Unblocks the Venkatarathnam-Oellrich criterion ADR-0017 rejected only once
  C4 (`d ln phi / dT`) exists.

## Amendment (C2, 2026-09-25): association included
The association term's temperature derivative is now added
(`_pcsaft_association.temperature_derivative`), so both methods work for
associating mixtures and the `ModelError` of decision 3 is removed for them
(the `Z <= 0` refusal stays). Michelsen-Hendriks stationarity makes it the
explicit partial at frozen site fractions,
`-(rho/2) sum w_a w_b X_a X_b dDelta_ab/dT`, with `T` in the Boltzmann factor
`exp(eps_ab/kT) - 1` and, through `d_i(T)`, in the contact value (`c_ij`,
`zeta_2`, `eta`); no sensitivity solve. Validation (ledger Case P-20):
FeOs residual entropy (same `T, V` reference) and enthalpy at the 18 Case P-6
states, **2.1e-15** worst with matched universal constants (asserted 1e-12),
1.7e-11 / 3.2e-10 with the published ones; central differences 1e-8 at five
associating states; Gibbs-Helmholtz through `ln phi` for liquid water and
water / ethanol. Still not provided: `Cp^res`, `d ln phi/dT`, total caloric
properties.

## Supersedes (optional)
None. Discharges the "no temperature derivative" gap of ADR-0014 for
non-associating mixtures.

## Superseded by (optional)
None.
