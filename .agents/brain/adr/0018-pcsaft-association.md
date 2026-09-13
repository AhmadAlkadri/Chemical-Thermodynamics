# ADR-0018: PC-SAFT association (Gross & Sadowski 2002), with Michelsen-Hendriks derivatives

Status: accepted
Date: 2026-09-13

## Context

ADR-0014 and ADR-0015 gave the package a real, validated PC-SAFT: hard chain
plus dispersion, analytic derivatives, density roots at `(T, P, x)`, and the
`EquationOfState` interface that `stability_tp` and `flash_tp` speak. What it
did not give it was **water, alcohols, acids or amines**. Those need the
association term of

> J. Gross and G. Sadowski, "Application of the Perturbed-Chain SAFT Equation
> of State to Associating Systems", *Ind. Eng. Chem. Res.* **41** (2002)
> 5510-5515 (DOI 10.1021/ie010954d),

and without it `PCSAFTEOS` had to be documented as "do not use this for water"
- a whole named class of fluids the package could model at low pressure with
NRTL and not at all with an equation of state. `brain.md`'s roadmap named this
"the only remaining gap that makes a named class of fluids unusable rather
than making an existing answer nicer", and ADR-0011's still-open
`flash-phase-addition-eos` item was blocked on it: no non-associating
hydrocarbon binary in this repository exhibits a third phase.

Two things make the association term different from every other term already
in the module:

1. **It has an inner unknown.** The non-bonded site fractions `X` solve a
   system of non-linear mass-action equations at every state. Nothing else in
   PC-SAFT needs an inner solve, and this one sits inside the density-root
   scan, which evaluates the isotherm at ~1600 packing fractions per call.
2. **Its derivatives are not a chain rule over `zeta`, `mbar`, `m2es3` and
   `m2e2s3`.** `X` depends on the state, so a naive derivative needs implicit
   differentiation of the mass-action equations for every variable.

There was also a reference problem. teqp - the independent implementation the
non-associating slices are validated against (Cases P-0 to P-5) - **does not
implement association** in its `PCSAFT` kind. And the 2002 paper is behind the
ACS paywall: `https://pubs.acs.org/doi/10.1021/ie010954d` returns HTTP 403 from
this environment, exactly as the 2001 paper did, so neither its equations nor
its parameter table could be read directly.

## Decision

### 1. The model: general `(na, nb)` sites, A-B bonding only, `sigma_ij^3`

The site fractions `X_alpha` (fraction of sites of kind `alpha` **not**
hydrogen bonded) solve

    1 / X_alpha = 1 + rho sum_beta w_beta X_beta Delta_{alpha beta}

with `w_alpha = x_{i(alpha)} n_alpha`, and

    a_assoc = sum_alpha w_alpha [ ln X_alpha - X_alpha / 2 + 1/2 ]

Each associating component carries `na` type-A sites and `nb` type-B sites;
only A-B pairs bond. `na = nb = 1` is the **2B** scheme, which is what the 2002
paper uses for water and the 1-alkanols and the only scheme this slice
**validates**. The equations are written for general `(na, nb)` so 3B and 4C
can follow without touching them, and `PCSAFTAssociationRecord` accepts those
counts, but nothing here cross-checks them and the README and the module
docstring say so.

The association strength is

    Delta_{alpha beta} = sigma_ij^3 g_ij^hs(d_ij) kappa^{AB}_ij
                         [ exp(eps^{AB}_ij / kT) - 1 ]

with `g_ij^hs` the **same** Boublik-Mansoori contact value (Eq. A.7 of the
2001 paper) the hard-chain term uses, evaluated at the unlike distance through
`c_ij = d_i d_j / (d_i + d_j)`.

**`sigma_ij^3`, not `d_ij^3`, and how that was settled.** Both spellings are in
circulation - SAFT is usually written with `sigma^3`, several presentations of
the perturbed-chain version print `d^3` - and the difference is not cosmetic:
at 300 K water's `d` is 2.9915 A against `sigma = 3.0007` A, so the cubes differ
by 0.93 % and `a_assoc` by 0.16 % of itself at a liquid density - far above any
tolerance in this package. **The paper was not read**, so the convention was settled
*numerically* against FeOs, which is the stronger evidence anyway. At
(300 K, 55000 mol/m^3) with the packaged water parameters, the `sigma^3` form
reproduces FeOs's association contribution to **0.0** (exactly, in double
precision) and the `d^3` form is **8.9e-03** away - eleven orders apart, so the
verdict cannot be ambiguous. Since the parameters and the convention were
regressed together, using the other spelling *with these parameters* would be
wrong whatever the paper prints. `test_the_association_strength_uses_sigma_cubed_and_not_d_cubed`
in `tests/validation/test_pcsaft_association_vs_feos.py` pins both halves of
that statement, recomputing the term from scratch in the test file.

Cross-association uses the Wolbach-Sandler rules (B. D. Wolbach and
S. I. Sandler, *IECR* **37** (1998) 2917):

    eps^{AB}_ij   = (eps^{AB}_ii + eps^{AB}_jj) / 2
    kappa^{AB}_ij = sqrt(kappa_ii kappa_jj)
                    [ sqrt(sigma_i sigma_j) / (0.5 (sigma_i + sigma_j)) ]^3

These were confirmed independently, not only by matching FeOs: Clapeyron.jl's
`PCSAFT_assoc.csv` stores the **water/ethanol cross pair explicitly** (epsilon
2577.05 K, bond volume 0.03356196748232913, source DOI 10.1021/ie010954d), and
the two rules above reproduce both numbers digit for digit from the two pure
records. `test_the_wolbach_sandler_cross_rules_are_the_ones_implemented` pins
that. **Induced association** - a non-associating component solvating with an
associating one - is *not* modelled: it needs a cross `kappa` that is not a
function of the pure-component ones.

### 2. Derivatives: Michelsen & Hendriks's stationary `Q`, not implicit differentiation

M. L. Michelsen and E. M. Hendriks, "Physical properties from association
models", *Fluid Phase Equilib.* **180** (2001) 165-174, give a function

    Q(X) = sum_alpha w_alpha (ln X_alpha - X_alpha + 1)
           - (rho / 2) sum_{alpha beta} w_alpha w_beta X_alpha X_beta Delta_{alpha beta}

whose stationarity `dQ/dX = 0` **is** the mass-action system, and which equals
`a_assoc` when that system is satisfied. Because it is stationary, the
derivative of `a_assoc` with respect to any state variable is the *explicit*
partial of `Q` at frozen `X`. So

- `rho (d a_assoc / d rho)_{T,x}`, and
- `(partial a_assoc / partial x_k)_{T,rho}` (unconstrained, at fixed total
  density - the convention every other term in `pcsaft.py` uses),

are plain algebraic expressions in the converged `X`, with no implicit
differentiation and no extra linear solve. That is what makes `Z` and
`ln phi_i` cost essentially nothing beyond the site-fraction solve itself.

The **second** density derivative is the one place stationarity is not enough,
because it differentiates a quantity that already depends on `X(state)`.
Writing the mass-action residual as `F(X, eta) = 0`,

    a''(eta) = Q_{eta eta} - (w .* F_eta) . ( J^-1 F_eta ),
    J_{alpha beta} = -delta_{alpha beta} / X_alpha^2 - rho w_beta Delta_{alpha beta}

which is the small linear solve Michelsen & Hendriks describe. **It is
analytic; no finite difference enters the shipped code.** That matters because
`a''` is what `(dP/drho)/RT = 1 + 2 eta a' + eta^2 a''` is built from, and
`dP/drho > 0` is the mechanical-stability filter that decides which density
roots the flash ever sees (ADR-0015). `J` is non-singular for `X > 0` even when
a component's mole fraction is zero (the `w_alpha` factors cancel out of the
solve), which is why a vanishing associating component in a mixture is not a
special case.

### 3. The site-fraction solve

`solve_site_fractions` starts from the closed form of the *decoupled* problem,

    X_alpha^(0) = 2 / (1 + sqrt(1 + 4 rho sum_beta w_beta Delta_{alpha beta}))

- exact for a pure 2B fluid - then takes 12 damped successive-substitution
steps and then Newton steps on

    G_alpha(X) = X_alpha (1 + rho (Delta W X)_alpha) - 1

The residual is written in that O(1) form rather than as `1/X - 1 - ...` so
that one absolute tolerance (`1e-14`) means the same thing in a dilute gas and
in dense liquid water, where `X ~ 0.01`. Everything is vectorised over leading
axes, so the density solver's whole packing-fraction grid is **one** call and
the `(npoints, nsites, nsites)` Newton system is one batched `np.linalg.solve`.
Measured: over the full 1599-point grid the worst mass-action residual is
**4.4e-16** for pure water, water/ethanol and water/n-hexane alike. Iteration
counts are fixed constants and the arithmetic is composition-ordered, so the
result is deterministic.

Successive substitution alone is not enough and is not relied on. Measured, to
a mass-action residual of 1e-14: **plain** (undamped) substitution needs 926
steps for pure water at 300 K / 55 kmol/m^3 and **does not converge at all**
within 200,000 steps for a 0.2/0.8 water/ethanol liquid at 300 K /
40 kmol/m^3, where the undamped map oscillates (500 plain steps there leave
`a_assoc` 4.2e-02 wrong). Damping to 0.5 fixes the oscillation - 14 and 52
steps respectively - but 12 of them plus Newton is faster and, more to the
point, gives a *quadratic* tail instead of a linear one: after the 12 damped
steps the pure-water case is already at 1e-14 (0 Newton steps) and the
water/ethanol case needs 2. Over the whole 1599-point density scan grid, 12
damped steps leave a worst residual of 4.4e-16, so Newton is usually a
no-op there and is there for the states where it is not.

### 4. The guard that keeps ADR-0014 bit-identical

`build_setup` returns `None` when **no component carries sites**, and both call
sites (`pcsaft._evaluate` and `_pcsaft_density.PCSAFTIsotherm`) then return
exactly the pre-slice expression. The guard is on *sites present*, not on the
argument being `None`, so passing `association=(None, None)` is the same code
path as passing nothing. Every ADR-0014 / ADR-0015 number is therefore
unchanged by construction, and `tests/test_pcsaft_association.py` asserts a
handful of them with `==`, not `approx`.

### 5. Parameters

`PCSAFTRecord` gains an optional `association: PCSAFTAssociationRecord | None`,
and `PCSAFTParameters.from_records` accepts it as an object or as a mapping
with `kappa_ab`, `epsilon_ab_k_K` and optional `na` / `nb` / `scheme` /
`source`. `scheme` is **documentation**: `na` and `nb` are what the model
reads, and a `scheme` naming counts that contradict them is rejected for the
names the class knows (`2B`, `3B`, `4C`), so a typo cannot silently disagree
with the numbers. `PCSAFTParameters.association_for_components` is the
accessor; `PCSAFTEOS.association_parameters` and `PCSAFTEOS.associates` expose
it on the model.

Five records were added to `src/chemthermo/parameters/data/eos/pcsaft.json`:
Water, Methanol, Ethanol, 1-Propanol and n-Butanol (the paper's "1-butanol";
the packaged databank spells it `n-Butanol`), all 2B. **Provenance, stated the
same way the 2001 set's is**: the primary citation is the 2002 paper, the paper
**was not read** (HTTP 403), and the values were transcribed from two
independent secondary sources that both cite its DOI and agree digit for digit
- FeOs's `parameters/pcsaft/gross2002.json` and Clapeyron.jl's
`PCSAFT_like.csv` / `PCSAFT_assoc.csv`. No cross-association parameters and no
`k_ij` are packaged.

### 6. Two new inspection methods on `PCSAFTEOS`

`residual_helmholtz_terms(temperature_K=, volume_m3=, composition=)` returns
`{"hard-chain", "dispersion", "association", "total"}`, and
`site_fractions(temperature_K=, density_mol_m3=, composition=)` returns the
converged `X` (empty for a non-associating mixture). They exist because the
association term is the thing this slice validates **term by term** against an
external reference, and without an accessor a caller could only ever see the
sum. They are the only additions to the public surface besides
`PCSAFTAssociationRecord`; no signature changed and no solver was touched.

## Validation (Cases P-6 and P-7)

**FeOs** (feos-org/feos, MIT OR Apache-2.0) is the reference: an independent
Rust implementation of the same model that obtains every derivative by
automatic differentiation (`num-dual`), and that *does* implement association.
It reproduces teqp on non-associating n-hexane, so the two references agree
where they overlap.

**One shared input is deliberately not shared.** chemthermo packages the 2001
paper's 42 universal constants **as printed**, to ten figures, which is also
what teqp uses (the two agree to 4e-15). FeOs hard-codes them to **fourteen**
figures; the two tables differ by up to 4.8e-09. That is an input difference,
not an implementation difference, and it floors any comparison of the
*dispersion* term. So every dispersion-dependent quantity is compared twice.
Worst deviation over 18 states (pure water and pure ethanol at four/five each,
water/ethanol at three compositions and two state types, water/n-hexane at
three, one ternary):

| quantity | shipped constants | FeOs's constants |
|---|---:|---:|
| hard-chain | 7.1e-15 | 7.1e-15 |
| dispersion | 7.1e-10 | 8.9e-15 |
| **association** | **3.1e-15** | **3.1e-15** |
| `A^res/RT` | 7.1e-10 | 7.1e-15 |
| `Z` | 3.9e-09 | 3.4e-14 |
| `ln phi_i` | 1.7e-06 | 1.3e-11 |

The association term does not depend on those constants and matches to
round-off either way, which is the comparison this slice is about.

**Derivative discipline** (no external dependency, ten states): `Z` against a
central difference of `A^res` in the density, `ln phi_i` against a central
difference of `n A^res` in mole numbers, the Euler identity
`sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z`, `a'` and `a''` of the
packing-fraction path against central differences, `dP/drho` against a central
difference of the public `pressure_Pa`, the mass-action residual, and
Michelsen-Hendriks stationarity `dQ/dX = 0` with `Q` written out in the test
file from its definition. Achieved: `|dZ|` <= 4.1e-09, `|d ln phi|` <= 7.3e-09,
Euler <= 1.8e-15, mass action <= 4.4e-16, `a'` <= 4.9e-09, `a''` <= 4.8e-09,
`dP/drho` <= 3.7e-10, stationarity <= 1e-12 (all finite-difference-limited).

**Equilibrium** (Case P-7). Pure water saturation at 373.15 K by bisecting
equal fugacity on chemthermo's own two density roots, against FeOs's
`PhaseEquilibrium.pure`: `Psat` 3.7e-10, `rho_L` 6.5e-11, `rho_V` 3.4e-10
relative, against an asserted 1e-6. A water/ethanol vapour-liquid flash at
351 K and 80 kPa (`z = 0.7/0.3`, `vapor_fraction = 0.6047`) with FeOs's
fugacities evaluated at chemthermo's phases and densities: equal fugacities to
6.1e-10. A water/n-hexane liquid-liquid split at 298.15 K and 1 MPa: equal
fugacities to 1.2e-11.

**Model versus experiment, remarked and not asserted.** PC-SAFT puts water's
saturation pressure at 373.15 K at 100.890 kPa, 0.43 % below the 101.325 kPa
that defines the normal boiling point. For water/n-hexane with `k_ij = 0` the
model gives 1.68e-05 hexane in the water-rich phase and 6.30e-03 water in the
hexane-rich phase, against commonly tabulated experimental figures of about
2e-6 and 5e-4 (not verified against a primary source here) - an order of
magnitude out on both, as expected without a fitted binary parameter. The
liquid-liquid case is a check of this code against another implementation, not
of the model against measurement.

## Limitations, recorded rather than worked around

- **A liquid-liquid EOS split is only reachable where the vapour root is
  gone.** At 298.15 K and 1 atm, `stability_tp` correctly calls the
  `z = 0.5/0.5` water/n-hexane feed unstable (`tpd_min = -9.28e-01`) and then
  `flash_tp` **raises** `ConvergenceError`. The reason is structural: the
  phi-phi split evaluates one phase on the model's vapour root and the other
  on its liquid root (`_split._ln_phi_function`), so it cannot put both phases
  on the liquid branch. It converges on a water-rich liquid against a
  hexane-rich *vapour* whose Gibbs energy is **above** the feed's
  (`delta_g_split_rt = +0.259`), and the post-split stability test refuses it.
  That refusal is correct behaviour for a wrong answer, but the message ("a
  third phase is required") describes the wrong disease: at 1 atm and 298 K the
  true answer is two liquids - water's and n-hexane's vapour pressures sum to
  about 23 kPa - and FeOs's own two-phase flash at that state returns exactly
  that. Above about 0.6 MPa the vapour root no longer exists at any
  composition, both phases sit on the single liquid root, and the same
  unchanged machinery returns the real tie line. **No solver was changed to
  make the 1 atm case work**; the slice brief forbade it and the fix belongs to
  `flash-phase-addition-eos`.
- **The labels a liquid-liquid split gets are wrong, and are recorded as
  produced.** At 1 MPa, `phase_identity` (ADR-0017) measures *both* converged
  phases as liquids, so the two-phase orientation rule falls through to its
  documented Wilson-ranking fallback: the phases come back named `"liquid"` and
  `"vapor"`, `diagnostics["phase_label_method"] == "wilson-ranking"`, and
  `vapor_fraction = 0.503164` is really the hexane-rich **liquid**'s fraction.
  ADR-0017's own consequence list anticipated the near-critical version of this
  ("both phases on the same side of the threshold"); this is a second, more
  common way in. Naming them `liquid1` / `liquid2` needs `_detect` to grow that
  vocabulary for phi-phi, which is a flash-module change and out of scope here.
- **2B only is validated.** General `(na, nb)` runs; nothing cross-checks it.
- **No induced association** (see decision 1).
- **No temperature derivative**, in the association term as everywhere else in
  this package's PC-SAFT, so still no caloric properties.
- **`k_ij = 0` is a poor model for water with a hydrocarbon**, and no binary
  table is packaged.
- Performance: a `density_roots` call for pure water costs ~5.0 ms against
  ~2.5 ms for non-associating n-hexane, i.e. the association term roughly
  doubles the isotherm scan. Not optimised further.

## Alternatives considered

- **Implicit differentiation of the mass-action equations for every
  derivative.** Rejected: it needs a linear solve per composition derivative as
  well as per density derivative, it is the textbook source of sign errors in
  association implementations, and Michelsen & Hendriks's stationary form makes
  every *first* derivative explicit for free. The linear solve survives in
  exactly one place - `a''` - where stationarity genuinely does not help.
- **A finite-difference `a''` inside the density-root solver**, which the slice
  brief allowed if documented and verified to 1e-7. Rejected once the
  Michelsen-Hendriks second derivative turned out to be a four-line addition on
  top of the Jacobian the Newton stage already assembles. `a''` decides
  mechanical stability at the spinodal, where `dP/drho` passes through zero and
  a differenced quantity is worst conditioned - the same argument ADR-0015 made
  for the non-associating terms.
- **Successive substitution alone for the site fractions**, with a large
  iteration cap. Rejected on measurement: 500 damped steps leave `a_assoc`
  4.2e-02 wrong in dense liquid water, and the tolerance this slice needs is
  1e-14.
- **teqp as the reference.** Not available: its `PCSAFT` kind has no
  association. FeOs was chosen because it is permissively licensed (MIT OR
  Apache-2.0), uses automatic differentiation rather than hand derivatives,
  ships the 2002 parameters, exposes per-term contributions
  (`residual_molar_helmholtz_energy_contributions`), and agrees with teqp on
  the non-associating part - so the two references corroborate each other.
- **Copying an existing implementation.** Explicitly not done. `pcsaft` (zmeri)
  is GPL-3 and was not consulted; the equations here are derived from the
  primary literature and settled numerically where the literature could not be
  read.
- **Adopting FeOs's fourteen-figure universal constants** to make the
  comparison tight as shipped. Rejected: this package's stated provenance is
  the paper's table **as printed**, teqp uses the same, and silently swapping
  in a third party's extended-precision variant would break bit-identity with
  Cases P-0 to P-5 for a change no user asked for. The difference is reported
  instead.

## Consequences

- Positive: water and the C1-C4 1-alkanols are usable PC-SAFT components. A
  water/ethanol VLE point and a water/n-hexane liquid-liquid tie line come out
  of the **unchanged** `stability_tp` / `flash_tp`.
- Positive: the association term is validated to 3.1e-15 against an
  independent implementation, and every analytic derivative that touches it is
  checked against a finite difference of the quantity it differentiates.
- Positive: `flash-phase-addition-eos` is now unblocked - water/n-hexane at
  1 atm is the real three-phase-capable EOS state ADR-0011 and Case L-4 said
  the repository did not have.
- Tradeoff: a new private module (`_pcsaft_association.py`) and a second call
  site for it, because the package keeps two differentiation paths (`(T, rho,
  x)` and `eta`) for the reasons ADR-0015 records. The association term is
  therefore written twice, once per variable, like the hard-chain and
  dispersion terms already are. `tests/test_pcsaft_association.py` asserts the
  two paths agree to 1e-12 at every validated state, which is what keeps them
  from drifting.
- Tradeoff: ~2x cost in the density-root scan for an associating mixture.
- Known and documented: the two limitations above (no liquid-liquid split at a
  pressure where the vapour root still exists; wrong phase names when a split
  is two liquids).

## Next slice

`flash-phase-addition-eos`, with **water / n-hexane at 298.15 K and 1 atm** as
the target case. It needs two things this ADR deliberately did not do: teach
the phi-phi path to put two phases on the same density branch (so a
liquid-liquid EOS split exists at a pressure where a vapour root also exists),
and give `_detect` the `liquid1` / `liquid2` vocabulary ADR-0017's measurement
already justifies. ADR-0016 supplied the inner `NP = 2` mechanism; this slice
supplies the physical state that exercises it.

## Supersedes (optional)

Amends ADR-0014's decision that "association and polar terms are out of scope"
- association is now in scope and implemented; the polar terms are not.
ADR-0014's other decisions, and all of ADR-0015, ADR-0016 and ADR-0017, are
unchanged.

## Superseded by (optional)

None.
