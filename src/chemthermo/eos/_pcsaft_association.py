"""PC-SAFT association contribution (Gross & Sadowski 2002), ADR-0018.

This module is private (ADR-0001). It is reached through
:class:`chemthermo.eos.PCSAFTEOS`, which adds its contribution to the
hard-chain and dispersion terms of :mod:`chemthermo.eos.pcsaft` and to the
packing-fraction path of :mod:`chemthermo.eos._pcsaft_density`.

Sources
-------
J. Gross and G. Sadowski, "Application of the Perturbed-Chain SAFT Equation of
State to Associating Systems", Ind. Eng. Chem. Res. 41 (2002) 5510-5515
(DOI 10.1021/ie010954d) - the model and the pure-component parameters.

M. L. Michelsen and E. M. Hendriks, "Physical properties from association
models", Fluid Phase Equilib. 180 (2001) 165-174 - the derivative route used
here.

W. G. Chapman, K. E. Gubbins, G. Jackson and M. Radosz, "New reference
equation of state for associating liquids", Ind. Eng. Chem. Res. 29 (1990)
1709-1721 - the SAFT association term the 2002 paper carries over.

J. R. Elliott and co-workers' ``kappa`` cross rule is *not* used; the cross
rules are those of B. D. Wolbach and S. I. Sandler, "Using Molecular Orbital
Calculations To Describe the Phase Behavior of Cross-Associating Mixtures",
Ind. Eng. Chem. Res. 37 (1998) 2917-2928.

The model
---------
A molecule of component ``i`` carries ``na_i`` sites of type A and ``nb_i``
sites of type B; only A-B pairs bond, which is the 2B scheme when
``na = nb = 1`` and the scheme Gross & Sadowski (2002) use for water and the
1-alkanols. Writing ``X_alpha`` for the fraction of sites of kind ``alpha``
that are **not** bonded, the association contribution to the reduced residual
Helmholtz energy per mole is

    a_assoc = sum_alpha w_alpha [ ln X_alpha - X_alpha / 2 + 1/2 ]        (1)

with the site weight ``w_alpha = x_{i(alpha)} n_alpha`` (mole fraction of the
owning component times how many sites of that kind it has), and the
``X_alpha`` solve the mass-action equations

    1 / X_alpha = 1 + rho sum_beta w_beta X_beta Delta_{alpha beta}       (2)

where ``rho`` is the number density of *molecules* in 1/Angstrom^3 and

    Delta_{alpha beta} = sigma_ij^3 g_ij^hs(d_ij) kappa^{AB}_ij
                         [ exp(eps^{AB}_ij / (k T)) - 1 ]                 (3)

for a bonding (A-B) pair of sites on components ``i`` and ``j``, and zero for a
non-bonding (A-A or B-B) pair. ``g_ij^hs`` is the **same** Boublik-Mansoori
contact value the hard-chain term uses (Eq. A.7 of the 2001 paper), evaluated
at the unlike distance through ``c_ij = d_i d_j / (d_i + d_j)``:

    g_ij = 1/U + 3 c_ij zeta_2 / U^2 + 2 c_ij^2 zeta_2^2 / U^3,  U = 1 - eta

The cross parameters are the Wolbach-Sandler rules

    eps^{AB}_ij   = (eps^{AB}_ii + eps^{AB}_jj) / 2                       (4)
    kappa^{AB}_ij = sqrt(kappa_ii kappa_jj)
                    [ sqrt(sigma_i sigma_j) / (0.5 (sigma_i + sigma_j)) ]^3

and ``sigma_ij = 0.5 (sigma_i + sigma_j)`` as everywhere else in PC-SAFT.

``sigma_ij^3``, not ``d_ij^3`` - and how that was settled
---------------------------------------------------------
Both spellings of Eq. (3) are in circulation: the SAFT association term is
usually written with ``sigma_ij^3`` while several presentations of the
perturbed-chain version print ``d_ij^3`` (the temperature-dependent diameter).
The difference is not cosmetic - at 300 K water's ``d`` is 2.9915 A against
``sigma = 3.0007`` A, so the cubes differ by ~0.9 % and the association term
by ~0.16 % of itself at a liquid density - far above any tolerance here.

**The 2002 paper itself was not read**: pubs.acs.org returned HTTP 403 for
DOI 10.1021/ie010954d from the environment that wrote this module, exactly as
for the 2001 paper. The convention was therefore settled *numerically*, which
is the stronger check anyway: with the packaged water parameters and the 2B
scheme, ``sigma^3`` reproduces FeOs's association contribution at
(300 K, 55000 mol/m^3) **exactly** (0.0 in double precision) while ``d^3`` is
off by 8.9e-03 - so the comparison cannot be ambiguous. ``sigma_ij^3`` is what is
implemented, and ``tests/validation/test_pcsaft_association_vs_feos.py``
pins both halves of that statement. Since the parameters and the convention
were regressed together, using the other one with these parameters would be
wrong regardless of what the paper prints.

Derivatives: Michelsen & Hendriks
---------------------------------
Equations (1) and (2) hide the ``X`` dependence on the state, and
differentiating (2) implicitly for every variable would be both tedious and
slow. Michelsen & Hendriks (2001) remove the need. Define, at fixed state,

    Q(X) = sum_alpha w_alpha (ln X_alpha - X_alpha + 1)
           - (rho / 2) sum_{alpha beta} w_alpha w_beta X_alpha X_beta
                                        Delta_{alpha beta}               (5)

Then ``dQ/dX_alpha = 0`` is exactly Eq. (2), and substituting (2) back into (5)
gives ``Q = a_assoc`` of Eq. (1) - so ``Q`` is a *stationary* representation of
the same number. Because it is stationary, the derivative of ``a_assoc`` with
respect to **any** state variable ``v`` is the explicit partial of ``Q`` at
frozen ``X``:

    d a_assoc / d v = (partial Q / partial v)_X                          (6)

That is what makes every first derivative here a plain algebraic expression.
Concretely, with ``S = sum w w X X Delta`` and the two moment-derivative sums
``S_2`` and ``S_3`` formed the same way but with ``dDelta/d zeta_2`` and
``dDelta/d zeta_3`` in place of ``Delta``:

    rho (d a_assoc / d rho)_{T,x} = -(rho/2) (S + zeta_2 S_2 + zeta_3 S_3)

    (partial a_assoc / partial x_k)_{T,rho}
        = sum_{alpha on k} n_alpha (ln X_alpha - X_alpha + 1)
          - (rho/2) [ 2 sum_{alpha on k} n_alpha X_alpha (Delta W X)_alpha
                      + S_2 (d zeta_2 / d x_k) + S_3 (d zeta_3 / d x_k) ]

(the composition derivative is unconstrained, at fixed total density, which is
the convention :mod:`chemthermo.eos.pcsaft` uses for every other term).

The **second** density derivative is the one place Eq. (6) is not enough: it is
the derivative of a quantity that already depends on ``X(state)``. Writing the
mass-action residual as ``F(X, eta) = 0`` and differentiating,

    a''(eta) = Q_{eta eta} - (w .* F_eta) . ( J^-1 F_eta ),
    J_{alpha beta} = -delta_{alpha beta} / X_alpha^2
                     - rho w_beta Delta_{alpha beta}                      (7)

which is the small linear solve Michelsen & Hendriks describe. It is exact and
analytic - no finite difference enters the density-root solver's
``dP/drho``, which is what decides mechanical stability. ``J`` is the Jacobian
of ``F``, is non-singular for ``X > 0`` even when a component's mole fraction
is zero, and is already assembled for the Newton stage of the site-fraction
solve, so Eq. (7) costs one extra back-substitution.

Solving the site fractions
--------------------------
:func:`solve_site_fractions` starts from the closed-form solution of the
*decoupled* problem (which is exact for a pure 2B fluid),

    X_alpha^(0) = 2 / (1 + sqrt(1 + 4 rho sum_beta w_beta Delta_{alpha beta}))

then takes damped successive substitution steps, then Newton steps on
``G_alpha(X) = X_alpha (1 + rho (Delta W X)_alpha) - 1`` - the mass-action
equation written so that its residual is O(1) rather than O(1/X), which is what
lets the same tolerance mean the same thing at a dilute gas and in dense
liquid water. Everything is vectorised over leading axes so the whole
packing-fraction scan grid of :mod:`chemthermo.eos._pcsaft_density` is one
call. The iteration counts are fixed constants and the arithmetic is
composition-ordered, so the result is deterministic.

Limits stated, not hidden
-------------------------
- Only ``na = nb = 1`` (2B) is **validated**. The equations are general in the
  site counts and 3B / 4C will run, but nothing here cross-checks them.
- No induced association (a non-associating component solvating with an
  associating one), because that needs a cross ``kappa`` that is not a function
  of the pure-component ones.
- No temperature derivative, in line with the rest of the PC-SAFT package.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ..exceptions import ModelError
from ..parameters.pcsaft import PCSAFTAssociationRecord

#: Damped successive-substitution steps taken before Newton. Successive
#: substitution alone needs tens of thousands of steps in dense liquid water
#: (measured: ``X ~ 8e-3`` at 300 K and 40 kmol/m^3); these are here to put
#: Newton inside its basin, not to converge on their own.
_SUBSTITUTION_STEPS = 12
#: Damping factor of those steps: ``X <- (1 - w) X + w X_fixed_point``.
_SUBSTITUTION_DAMPING = 0.5
#: Newton steps on the O(1) mass-action residual. Newton is quadratic here, so
#: this is a generous ceiling rather than a working budget.
_NEWTON_STEPS = 60
#: Convergence tolerance on ``max |X (1 + rho Delta W X) - 1|``.
SITE_FRACTION_TOL = 1e-14
#: Floor applied to a Newton iterate so ``X`` cannot leave ``(0, 1]``.
_X_FLOOR = 1e-300


#: What :meth:`AssociationIsotherm.derivatives` returns in the ``a`` slot when
#: the caller asked for the derivatives only (``value=False``). See the
#: matching constant in :mod:`chemthermo.eos._pcsaft_density`.
_UNREQUESTED = np.array([float("nan")])


class AssociationParameterError(ModelError):
    """Raised when a set of association parameters cannot be used."""


@dataclass(frozen=True)
class AssociationTopology:
    """Which association sites exist and which of them may bond.

    Built once per component-name order. ``site_component`` maps each site to
    the index of the component that carries it **in the full component list**
    (so a non-associating component simply owns no site); ``multiplicity`` is
    ``na`` or ``nb``; ``bond`` is 1.0 for an A-B site pair and 0.0 otherwise.
    ``pair_component`` is the ``(nsites, nsites)`` pair of component indices,
    kept so the component-level constants can be spread over the site grid with
    a single fancy-index.
    """

    site_component: np.ndarray
    multiplicity: np.ndarray
    bond: np.ndarray
    associating: np.ndarray

    @property
    def site_count(self) -> int:
        return int(self.site_component.size)


@dataclass(frozen=True)
class AssociationSetup:
    """Everything about the association term that is fixed at ``(T, names)``.

    Attributes:
        topology: The site bookkeeping.
        pair_constant: ``sigma_ij^3 kappa_ij [exp(eps_ij/kT) - 1]`` per site
            pair, in Angstrom^3, with the non-bonding pairs already zeroed.
            Multiplying it by ``g_ij`` gives ``Delta`` of Eq. (3).
        pair_c: ``c_ij = d_i d_j / (d_i + d_j)`` per site pair, in Angstrom.
    """

    topology: AssociationTopology
    pair_constant: np.ndarray
    pair_c: np.ndarray

    def weights(self, x: np.ndarray) -> np.ndarray:
        """Return the site weights ``w_alpha = x_{i(alpha)} n_alpha``."""
        return x[self.topology.site_component] * self.topology.multiplicity


def build_setup(
    *,
    temperature_K: float,
    sigma_A: np.ndarray,
    d: np.ndarray,
    association: Sequence[PCSAFTAssociationRecord | None],
) -> AssociationSetup | None:
    """Return the association setup for one ``(T, component order)``.

    Returns ``None`` when no component in ``association`` carries sites, which
    is the guard that keeps the non-associating path of ADR-0014 bit-identical:
    the caller then never touches this module again.

    Args:
        temperature_K: Temperature in K.
        sigma_A: Segment diameters in Angstrom, one per component.
        d: Temperature-dependent segment diameters (Eq. A.9) in Angstrom, one
            per component.
        association: One :class:`PCSAFTAssociationRecord` (or ``None``) per
            component, in the same order.
    """
    associating = np.array([record is not None for record in association], dtype=bool)
    if not associating.any():
        return None

    owners = np.flatnonzero(associating)
    records = [association[int(index)] for index in owners]

    kappa = np.array([record.kappa_ab for record in records], dtype=float)  # type: ignore[union-attr]
    epsilon = np.array(
        [record.epsilon_ab_k_K for record in records],  # type: ignore[union-attr]
        dtype=float,
    )
    na = np.array([record.na for record in records], dtype=float)  # type: ignore[union-attr]
    nb = np.array([record.nb for record in records], dtype=float)  # type: ignore[union-attr]

    # Two sites per associating component: an A then a B, interleaved so the
    # site order is a deterministic function of the component order.
    site_owner = np.repeat(np.arange(owners.size), 2)
    site_component = owners[site_owner]
    multiplicity = np.empty(2 * owners.size, dtype=float)
    multiplicity[0::2] = na
    multiplicity[1::2] = nb
    is_type_a = np.zeros(2 * owners.size, dtype=bool)
    is_type_a[0::2] = True
    bond = (is_type_a[:, None] != is_type_a[None, :]).astype(float)

    # Component-level (local index) cross parameters: Eq. (4) and sigma_ij.
    sigma_local = sigma_A[owners]
    sigma_ij = 0.5 * (sigma_local[:, None] + sigma_local[None, :])
    kappa_ij = (
        np.sqrt(np.outer(kappa, kappa))
        * (np.sqrt(np.outer(sigma_local, sigma_local)) / sigma_ij) ** 3
    )
    epsilon_ij = 0.5 * (epsilon[:, None] + epsilon[None, :])
    boltzmann = np.expm1(epsilon_ij / float(temperature_K))
    constant_local = sigma_ij**3 * kappa_ij * boltzmann
    if not np.all(np.isfinite(constant_local)):
        raise AssociationParameterError(
            "PC-SAFT association strength is not finite at "
            f"T = {float(temperature_K)!r} K; check epsilon_ab_k_K."
        )

    d_local = d[owners]
    c_local = np.outer(d_local, d_local) / (d_local[:, None] + d_local[None, :])

    spread = np.ix_(site_owner, site_owner)
    return AssociationSetup(
        topology=AssociationTopology(
            site_component=site_component,
            multiplicity=multiplicity,
            bond=bond,
            associating=associating,
        ),
        pair_constant=constant_local[spread] * bond,
        pair_c=c_local[spread],
    )


def solve_site_fractions(
    *,
    rho: np.ndarray | float,
    weights: np.ndarray,
    delta: np.ndarray,
) -> np.ndarray:
    """Return the non-bonded site fractions solving Eq. (2).

    Args:
        rho: Number density of molecules in 1/Angstrom^3. A scalar, or an array
            broadcasting against ``delta``'s leading axes.
        weights: ``(nsites,)`` site weights from
            :meth:`AssociationSetup.weights`.
        delta: ``(..., nsites, nsites)`` association strengths of Eq. (3).

    Returns:
        ``(..., nsites)`` site fractions in ``(0, 1]``.

    Raises:
        AssociationParameterError: If the iteration does not reach
            :data:`SITE_FRACTION_TOL`. The message reports the residual, so a
            failure names a number rather than a fixed point.
    """
    delta = np.asarray(delta, dtype=float)
    # One shape for both callers: a scalar density becomes a ``(1,)`` array,
    # which broadcasts against ``(nsites,)`` exactly as a ``(npoints, 1)``
    # column broadcasts against ``(npoints, nsites)``.
    rho_column = np.asarray(rho, dtype=float)[..., None]
    weighted = delta * weights  # Delta_{ab} w_b, the matrix that acts on X

    # Decoupled start: exact for a pure fluid whose sites are all equivalent.
    diagonal_strength = weighted.sum(axis=-1)
    x_sites = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * rho_column * diagonal_strength))

    for _ in range(_SUBSTITUTION_STEPS):
        target = 1.0 / (1.0 + rho_column * _apply(weighted, x_sites))
        x_sites = (1.0 - _SUBSTITUTION_DAMPING) * x_sites + _SUBSTITUTION_DAMPING * target

    identity = np.eye(weights.size)
    for _ in range(_NEWTON_STEPS):
        bonded = rho_column * _apply(weighted, x_sites)
        residual = x_sites * (1.0 + bonded) - 1.0
        if float(np.max(np.abs(residual))) <= SITE_FRACTION_TOL:
            return x_sites
        jacobian = (
            identity * (1.0 + bonded)[..., None]
            + x_sites[..., :, None] * rho_column[..., None] * weighted
        )
        step = np.linalg.solve(jacobian, residual[..., None])[..., 0]
        # Never leave (0, 1]: halve any step that would, which is what keeps
        # the first Newton step safe on a strongly associating liquid.
        proposed = x_sites - step
        bad = ~(proposed > 0.0)
        while bool(np.any(bad)):
            step = np.where(bad, 0.5 * step, step)
            proposed = x_sites - step
            new_bad = ~(proposed > 0.0)
            if bool(np.all(new_bad == bad)) and bool(np.all(np.abs(step) < _X_FLOOR)):
                break
            bad = new_bad
        x_sites = np.minimum(np.maximum(proposed, _X_FLOOR), 1.0)

    residual = x_sites * (1.0 + rho_column * _apply(weighted, x_sites)) - 1.0
    worst = float(np.max(np.abs(residual)))
    if worst > SITE_FRACTION_TOL:
        raise AssociationParameterError(
            "PC-SAFT association site fractions did not converge: worst mass-action "
            f"residual {worst!r} after {_SUBSTITUTION_STEPS} damped substitutions and "
            f"{_NEWTON_STEPS} Newton steps (tolerance {SITE_FRACTION_TOL!r})."
        )
    return x_sites


def mass_action_residual(
    *,
    rho: np.ndarray | float,
    weights: np.ndarray,
    delta: np.ndarray,
    x_sites: np.ndarray,
) -> np.ndarray:
    """Return ``X_a (1 + rho (Delta W X)_a) - 1``, zero at the solution."""
    rho_column = np.asarray(rho, dtype=float)[..., None]
    weighted = np.asarray(delta, dtype=float) * weights
    return x_sites * (1.0 + rho_column * _apply(weighted, x_sites)) - 1.0


def _apply(weighted: np.ndarray, x_sites: np.ndarray) -> np.ndarray:
    """Return ``(Delta W X)_alpha`` for a batch of matrices and vectors."""
    return np.einsum("...ab,...b->...a", weighted, x_sites)


def helmholtz(weights: np.ndarray, x_sites: np.ndarray) -> np.ndarray:
    """Return ``a_assoc`` of Eq. (1) from converged site fractions."""
    return (weights * (np.log(x_sites) - 0.5 * x_sites + 0.5)).sum(axis=-1)


def q_site_term(weights: np.ndarray, x_sites: np.ndarray) -> np.ndarray:
    """Return ``sum_alpha w_alpha (ln X - X + 1)``, the first sum of Eq. (5)."""
    return (weights * (np.log(x_sites) - x_sites + 1.0)).sum(axis=-1)


def contact_values(
    *, pair_c: np.ndarray, zeta_2: float, eta: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(g_ij, dg/d zeta_2, dg/d zeta_3)`` at unlike contact.

    The Boublik-Mansoori contact value of Eq. (A.7) of the 2001 paper, with
    ``c_ij = d_i d_j / (d_i + d_j)`` instead of the like-pair ``d_i / 2``.
    """
    u = 1.0 - eta
    g = 1.0 / u + 3.0 * pair_c * zeta_2 / u**2 + 2.0 * pair_c**2 * zeta_2**2 / u**3
    dg_dz2 = 3.0 * pair_c / u**2 + 4.0 * pair_c**2 * zeta_2 / u**3
    dg_dz3 = 1.0 / u**2 + 6.0 * pair_c * zeta_2 / u**3 + 6.0 * pair_c**2 * zeta_2**2 / u**4
    return g, dg_dz2, dg_dz3


class AssociationState:
    """The association term and its derivatives at one ``(T, rho, x)``.

    This is the ``(T, rho, x)`` path that :mod:`chemthermo.eos.pcsaft` uses.
    """

    def __init__(
        self,
        setup: AssociationSetup,
        *,
        rho_a3: float,
        x: np.ndarray,
        zeta_2: float,
        eta: float,
        dzeta_dx: np.ndarray,
    ) -> None:
        weights = setup.weights(x)
        g, dg_dz2, dg_dz3 = contact_values(pair_c=setup.pair_c, zeta_2=zeta_2, eta=eta)
        delta = setup.pair_constant * g
        x_sites = solve_site_fractions(rho=rho_a3, weights=weights, delta=delta)

        weighted_sites = weights * x_sites
        pair_sum = float(weighted_sites @ delta @ weighted_sites)
        pair_sum_z2 = float(weighted_sites @ (setup.pair_constant * dg_dz2) @ weighted_sites)
        pair_sum_z3 = float(weighted_sites @ (setup.pair_constant * dg_dz3) @ weighted_sites)

        self.site_fractions = x_sites
        self.a_assoc = float(helmholtz(weights, x_sites))
        # Eq. (6) in the density: every zeta moment is proportional to rho.
        self.z_assoc = -0.5 * rho_a3 * (pair_sum + zeta_2 * pair_sum_z2 + eta * pair_sum_z3)

        # Eq. (6) in the mole fractions, unconstrained and at fixed density.
        component_count = x.size
        owners = setup.topology.site_component
        multiplicity = setup.topology.multiplicity
        explicit = np.bincount(
            owners,
            weights=multiplicity * (np.log(x_sites) - x_sites + 1.0),
            minlength=component_count,
        )
        bonded = delta @ weighted_sites
        cross = np.bincount(
            owners,
            weights=multiplicity * x_sites * bonded,
            minlength=component_count,
        )
        self.da_dx = explicit - 0.5 * rho_a3 * (
            2.0 * cross + pair_sum_z2 * dzeta_dx[2] + pair_sum_z3 * dzeta_dx[3]
        )


class AssociationIsotherm:
    """The association term as a function of ``eta`` at fixed ``(T, x)``.

    This is the packing-fraction path that
    :mod:`chemthermo.eos._pcsaft_density` uses. ``ratio_2`` is ``M_2 / M_3``
    and ``m3`` is ``M_3`` in that module's notation, so ``zeta_2 = ratio_2 eta``
    and ``rho = eta / M_3``.
    """

    def __init__(
        self,
        setup: AssociationSetup,
        *,
        x: np.ndarray,
        ratio_2: float,
        m3: float,
    ) -> None:
        self._setup = setup
        self._weights = setup.weights(x)
        self._m3 = float(m3)
        # g_ij(eta) = 1/U + B eta/U^2 + C eta^2/U^3 - the same shape the like
        # pairs take in _pcsaft_density, with c_ij for the unlike pair.
        self._b = 3.0 * setup.pair_c * ratio_2
        self._c = 2.0 * setup.pair_c**2 * ratio_2**2

    def _pair_terms(
        self, eta: np.ndarray, *, second: bool
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(Delta, dDelta/deta, d2Delta/deta2)`` on a batch of ``eta``.

        The second derivative is :data:`_UNREQUESTED` (``nan``) when ``second``
        is False: nothing reads it there, and on the 1599-point scan grid it
        was a ``(points, sites, sites)`` allocation per call.
        """
        column = eta[..., None, None]
        u = 1.0 - column
        g = 1.0 / u + self._b * column / u**2 + self._c * column**2 / u**3
        g1 = (
            1.0 / u**2
            + self._b * (1.0 + column) / u**3
            + self._c * (2.0 * column + column**2) / u**4
        )
        constant = self._setup.pair_constant
        if not second:
            return constant * g, constant * g1, _UNREQUESTED
        g2 = (
            2.0 / u**3
            + self._b * (4.0 + 2.0 * column) / u**4
            + self._c * (2.0 + 8.0 * column + 2.0 * column**2) / u**5
        )
        return constant * g, constant * g1, constant * g2

    def derivatives(
        self, eta: np.ndarray, *, second: bool, value: bool = True
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(a, a', a'')`` in ``eta``; ``a''`` is zero when not asked.

        ``a'`` is Eq. (6); ``a''`` is Eq. (7), the one derivative that needs the
        site-fraction sensitivity and therefore a linear solve.

        ``value=False`` skips Eq. (1) itself - the caller wants only the
        derivatives - and ``a`` comes back as ``nan``. Nothing else changes:
        Eq. (1) is a function of the converged site fractions alone and feeds
        nothing below it, so ``a'`` and ``a''`` are the same doubles either way
        (ADR-0023).
        """
        eta = np.asarray(eta, dtype=float)
        rho = eta / self._m3
        delta, delta_1, delta_2 = self._pair_terms(eta, second=second)
        weights = self._weights
        x_sites = solve_site_fractions(rho=rho, weights=weights, delta=delta)

        weighted_sites = weights * x_sites
        half = 0.5 / self._m3
        pair_sum = _quadratic(weighted_sites, delta)
        pair_sum_1 = _quadratic(weighted_sites, delta_1)

        a = helmholtz(weights, x_sites) if value else _UNREQUESTED
        a1 = -half * (pair_sum + eta * pair_sum_1)
        if not second:
            return a, a1, np.zeros_like(a1)

        pair_sum_2 = _quadratic(weighted_sites, delta_2)
        q_eta_eta = -half * (2.0 * pair_sum_1 + eta * pair_sum_2)

        # F_alpha(X, eta) = 1/X_a - 1 - rho (Delta W X)_a ; see Eq. (7).
        combined = delta + eta[..., None, None] * delta_1
        f_eta = -(1.0 / self._m3) * _apply(combined * weights, x_sites)
        jacobian = -np.eye(weights.size) / (x_sites**2)[..., None] - (
            rho[..., None, None] * delta * weights
        )
        sensitivity = np.linalg.solve(jacobian, f_eta[..., None])[..., 0]
        a2 = q_eta_eta - (weights * f_eta * sensitivity).sum(axis=-1)
        return a, a1, a2


def _quadratic(weighted_sites: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """Return ``(W X) . M . (W X)`` over a batch of matrices."""
    return np.einsum("...a,...ab,...b->...", weighted_sites, matrix, weighted_sites)


__all__ = [
    "AssociationIsotherm",
    "AssociationParameterError",
    "AssociationSetup",
    "AssociationState",
    "AssociationTopology",
    "SITE_FRACTION_TOL",
    "build_setup",
    "contact_values",
    "helmholtz",
    "mass_action_residual",
    "q_site_term",
    "solve_site_fractions",
]
