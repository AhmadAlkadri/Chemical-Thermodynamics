"""PC-SAFT density roots at a specified ``(T, P, x)`` (ADR-0015).

This module is private (ADR-0001). Its one public-facing result is reached
through :meth:`chemthermo.eos.PCSAFTEOS.density_roots` and, indirectly,
through :meth:`chemthermo.eos.PCSAFTEOS.fugacity_coefficients`.

What has to be solved
---------------------
:mod:`chemthermo.eos.pcsaft` evaluates the model at a state fixed by
``(T, rho, x)``. The ``EquationOfState`` interface the flash and the stability
test speak instead fixes ``(T, P, x)`` and asks for a *phase*. Going from one
to the other means solving

    P_model(T, rho, x) = P

for ``rho`` and keeping only the mechanically stable solutions, those with

    (dP / d rho)_{T,x} > 0 .

A cubic does this by factoring a polynomial. PC-SAFT is not a polynomial in
the density, so the roots have to be found numerically. Unlike a cubic it also
has no a-priori bound on how many roots it has: this module finds all of the
ones its grid resolves and says how many that was, rather than assuming three.

Reduction to one variable: the packing fraction
-----------------------------------------------
At fixed ``T`` and ``x`` every quantity the model needs is a function of the
packing fraction ``eta = zeta_3`` alone, because all four moments are
proportional to the density:

    zeta_n = M_n * rho_A3 ,  M_n = (pi/6) sum_i x_i m_i d_i^n ,  eta = M_3 rho_A3

so ``zeta_n = (M_n / M_3) eta`` and the molar density is ``rho = K eta`` with
``K = 1 / (M_3 N_A 1e-30)``. Writing ``a(eta) = A^res/(R T)`` in that single
variable turns the two derivatives the solver needs into ordinary derivatives:

    Z(eta)       = 1 + eta a'(eta)
    (dP/drho)/RT = 1 + 2 eta a'(eta) + eta^2 a''(eta)

(the second follows from ``P = Z rho R T`` and ``rho d/drho = eta d/deta``).
``a''`` is the one derivative :mod:`chemthermo.eos.pcsaft` does not have -
that module differentiates once with respect to the density and once with
respect to composition - so it is derived here, in the same notation.

Complexity receipt for the re-derivation
----------------------------------------
Why: the scan below needs thousands of pressures per call, and it needs
``dP/drho`` for both the Newton refinement and the mechanical-stability
filter. Evaluating :func:`chemthermo.eos.pcsaft._evaluate` point by point
would be ~1e3 times slower and would still not supply ``a''``.

Bug prevented: a finite-difference ``dP/drho`` would decide mechanical
stability - which root is returned and which is discarded - from a
differenced quantity that is ill-conditioned exactly where it matters, at the
spinodal where ``dP/drho`` passes through zero.

Cost: the hard-chain and dispersion terms appear twice in the package, once
per differentiation variable.

What keeps them honest: ``tests/test_pcsaft_density.py`` asserts that this
module's ``a``, ``Z`` and ``P`` equal
:meth:`~chemthermo.eos.PCSAFTEOS.residual_helmholtz`,
:meth:`~chemthermo.eos.PCSAFTEOS.compressibility_factor` and
:meth:`~chemthermo.eos.PCSAFTEOS.pressure_Pa` to 1e-12 relative over a grid of
states, and that ``dP/drho`` matches a central difference of the *other*
module's pressure. If the two ever drift apart, that test goes red.

The derivatives, term by term
-----------------------------
With ``U = 1 - eta``, ``r_n = M_n / M_3`` and the constants
``A1 = 3 r_1 r_2 / r_0``, ``A2 = r_2^3 / r_0``, ``A3 = r_2^3 / r_0 - 1``, the
Boublik-Mansoori term (Eq. A.6) collapses to

    a_hs  = A1 eta/U + A2 eta/U^2 + A3 ln U
    a_hs' = A1/U^2 + A2 (1+eta)/U^3 - A3/U
    a_hs" = 2 A1/U^3 + A2 (4+2 eta)/U^4 - A3/U^2

and the contact value (Eq. A.7), with ``B_i = 3 c_i r_2`` and
``C_i = 2 c_i^2 r_2^2`` for ``c_i = d_i/2``, to

    g_i  = 1/U + B_i eta/U^2 + C_i eta^2/U^3
    g_i' = 1/U^2 + B_i (1+eta)/U^3 + C_i (2 eta + eta^2)/U^4
    g_i" = 2/U^3 + B_i (4+2 eta)/U^4 + C_i (2 + 8 eta + 2 eta^2)/U^5

so that ``a_hc = mbar a_hs - sum_i x_i (m_i-1) ln g_i`` differentiates by the
chain rule, the second derivative picking up the usual
``g"/g - (g'/g)^2``.

The dispersion term (Eq. A.10) is ``a_disp = (eta/M_3) f(eta)`` with
``f = -2 pi I1 m2es3 - pi mbar C1 I2 m2e2s3``, hence

    a_disp'  = (f + eta f') / M_3
    a_disp" = (2 f' + eta f") / M_3

``I1`` and ``I2`` are degree-six polynomials in ``eta``, so their derivatives
are term-by-term. ``C1 = 1/D`` (the typo-corrected Eq. A.11; see the
:mod:`chemthermo.eos.pcsaft` docstring) gives ``C1' = -D'/D^2`` and
``C1" = -D"/D^2 + 2 D'^2/D^3``; ``D'`` is the bracket of Eq. A.31 and ``D"``
is one more differentiation of it:

    D  = 1 + mbar (8 eta - 2 eta^2)/U^4
           + (1-mbar) (20 eta - 27 eta^2 + 12 eta^3 - 2 eta^4)/gap^2
    D' = mbar (8 + 20 eta - 4 eta^2)/U^5 + (1-mbar) Q/gap^3
    D" = mbar (60 + 72 eta - 12 eta^2)/U^6
           + (1-mbar) (Q' gap - 3 Q gap') / gap^4

with ``gap = (1-eta)(2-eta)``, ``gap' = 2 eta - 3``,
``Q = 2 eta^3 + 12 eta^2 - 48 eta + 40`` and ``Q' = 6 eta^2 + 24 eta - 48``.

How the roots are found
-----------------------
1. **Bracket.** Evaluate ``P_model(eta) - P`` on the deterministic grid
   :data:`ETA_GRID_DESCRIPTION`: 120 geometrically spaced points from
   ``1e-14`` to ``1e-3``, then a uniform step of ``5e-4`` up to ``0.7405``
   (the hard-sphere close-packing limit ``pi/(3 sqrt 2)``, above which the
   Boublik-Mansoori term has no meaning). The low end is geometric because a
   dilute vapour root sits at ``eta ~ 1e-8`` at low pressure while a liquid
   root sits near ``0.4``; the model is monotone there (``P -> rho R T``), so
   spacing cannot cost a root. Above ``1e-3`` it is uniform because that is
   where the non-monotone part lives and a uniform step is what bounds how
   narrow a pair of roots may be before it is missed.

2. **Refine.** Every sign change is refined by a safeguarded Newton iteration
   (Newton step when it stays inside the bracket and makes progress, bisection
   otherwise) on ``g(eta) = P_model(eta) - P`` with the analytic
   ``dP/deta = K dP/drho``. It stops when the bracket is at machine width or
   the residual is below ``1e-14 P``.

   **How small the residual can get is a property of the state, not of the
   iteration.** On a dense liquid branch at a low pressure the pressure is a
   near-total cancellation: at n-hexane's 300 K saturation state
   ``Z = 1 + eta a'(eta)`` is ``1.17e-3``, so a relative error of ``1e-14`` in
   ``a'`` is a relative error of ``1e-11`` in ``P``. The measured worst
   ``|P_model - P|/P`` over the returned roots is therefore ``2e-14`` at
   10 MPa but ``6.5e-12`` at 21.9 kPa on the liquid root and ``1.1e-7`` at
   1 Pa. The *root* is unaffected - it agrees with teqp's own saturation
   density to 2.2e-16 relative - and the density, not the residual, is what
   the tests pin. :attr:`DensityRoots.max_relative_residual` reports the
   number rather than hiding it.

3. **Filter.** A refined root is kept only if ``dP/drho > 0``. The spinodal
   branch - the middle root, where ``dP/drho < 0`` - is discarded, never
   returned.

4. **Report.** The kept roots are returned sorted by increasing density. The
   lowest is the vapour-like candidate, the highest the liquid-like one. One
   root means the two labels name the same state; that is a fact about the
   fluid, not a failure, and is exactly how a cubic behaves outside its
   three-root region. No root at all raises ``ModelError`` naming the state.

Known limit
-----------
A pair of roots closer together than the grid step in ``eta`` is invisible to
step 1 and the state is then reported with one root fewer. That happens where
the isotherm is nearly tangent to the target pressure, i.e. at a near-critical
or near-spinodal state, and there the missed pair is a physically marginal
phase. The grid step is a named constant so the trade is explicit rather than
hidden.
"""

from __future__ import annotations

import math
from typing import NamedTuple, Sequence

import numpy as np

from ..exceptions import ModelError
from ..parameters.pcsaft import PCSAFTAssociationRecord
from ._pcsaft_association import AssociationIsotherm, build_setup

#: Exact SI definitions, mirrored from :mod:`chemthermo.eos.pcsaft` (importing
#: them from there would be circular: that module imports this one).
BOLTZMANN_J_PER_K = 1.380649e-23
AVOGADRO_PER_MOL = 6.02214076e23
R_J_PER_MOL_K = AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K

#: Hard-sphere close packing, ``pi / (3 sqrt 2)``. Above it the
#: Boublik-Mansoori term has no physical meaning, so no root is looked for.
ETA_CLOSE_PACKING = math.pi / (3.0 * math.sqrt(2.0))

#: Upper end of the scan, just below close packing.
ETA_MAX = 0.7405
#: Lower end of the scan. A dilute vapour at 1 Pa and 300 K sits near 1e-11.
ETA_MIN = 1e-14
#: Where the geometric low-density part of the grid hands over to the uniform
#: part. Below it the isotherm is monotone (ideal-gas-like).
ETA_GEOMETRIC_LIMIT = 1e-3
#: Number of geometrically spaced points below :data:`ETA_GEOMETRIC_LIMIT`.
ETA_GEOMETRIC_POINTS = 120
#: Uniform step in ``eta`` above :data:`ETA_GEOMETRIC_LIMIT`. Two roots closer
#: together than this are not resolved; see "Known limit" above.
ETA_UNIFORM_STEP = 5.0e-4

ETA_GRID_DESCRIPTION = (
    f"{ETA_GEOMETRIC_POINTS} geometric points on [{ETA_MIN:g}, {ETA_GEOMETRIC_LIMIT:g}] "
    f"then a uniform step of {ETA_UNIFORM_STEP:g} up to {ETA_MAX:g}"
)

#: Relative residual at which the Newton refinement stops.
_ROOT_RTOL = 1e-14
#: Maximum safeguarded-Newton iterations per bracket. Bisection alone would
#: need ~60 to reach machine width, so this never binds in practice.
_MAX_REFINE_ITERATIONS = 200


def _eta_grid() -> np.ndarray:
    """Return the deterministic scan grid; see :data:`ETA_GRID_DESCRIPTION`."""
    low = np.geomspace(ETA_MIN, ETA_GEOMETRIC_LIMIT, ETA_GEOMETRIC_POINTS)
    count = int(math.ceil((ETA_MAX - ETA_GEOMETRIC_LIMIT) / ETA_UNIFORM_STEP)) + 1
    high = np.linspace(ETA_GEOMETRIC_LIMIT, ETA_MAX, count)
    return np.concatenate((low, high[1:]))


#: Built once: it depends on nothing but the constants above.
_ETA_GRID = _eta_grid()

_POWERS = np.arange(7)


class _EtaDerivatives(NamedTuple):
    """``a(eta)`` and its first two ``eta`` derivatives at fixed ``(T, x)``."""

    a: np.ndarray
    a1: np.ndarray
    a2: np.ndarray


class PCSAFTIsotherm:
    """The PC-SAFT isotherm at one ``(T, x)``, as a function of ``eta``.

    Internal helper (ADR-0001). Built once per ``(T, P, x)`` request and then
    evaluated on the whole scan grid at once, which is what makes the scan
    affordable inside a flash loop.
    """

    def __init__(
        self,
        *,
        temperature_K: float,
        composition: np.ndarray,
        m: np.ndarray,
        sigma_A: np.ndarray,
        epsilon_k_K: np.ndarray,
        kij: np.ndarray,
        association: Sequence[PCSAFTAssociationRecord | None] | None = None,
    ) -> None:
        self.temperature = float(temperature_K)
        x = np.asarray(composition, dtype=float)

        # Eq. A.9, and the four moments of Eq. A.8 divided by the density.
        d = sigma_A * (1.0 - 0.12 * np.exp(-3.0 * epsilon_k_K / self.temperature))
        xm = x * m
        moments = (math.pi / 6.0) * (xm[None, :] * d[None, :] ** np.arange(4)[:, None]).sum(axis=1)
        m3 = float(moments[3])
        if not math.isfinite(m3) or m3 <= 0.0:
            raise ModelError("PC-SAFT segment-size moment M_3 is not positive.")
        self._m3 = m3
        ratios = moments / m3

        #: ``rho_molar = density_per_eta * eta``.
        self.density_per_eta = 1.0 / (m3 * AVOGADRO_PER_MOL * 1e-30)

        self._mbar = float(xm.sum())
        self._a1_const = 3.0 * ratios[1] * ratios[2] / ratios[0]
        self._a2_const = ratios[2] ** 3 / ratios[0]
        self._a3_const = self._a2_const - 1.0

        c_ii = 0.5 * d
        self._b_i = 3.0 * c_ii * ratios[2]
        self._c_i = 2.0 * c_ii**2 * ratios[2] ** 2
        self._chain_weight = x * (m - 1.0)

        # Eq. A.12-A.14.
        sigma_ij3 = (0.5 * (sigma_A[:, None] + sigma_A[None, :])) ** 3
        eps_over_kt = np.sqrt(np.outer(epsilon_k_K, epsilon_k_K)) * (1.0 - kij) / self.temperature
        weights = np.outer(xm, xm)
        self._m2es3 = float(np.sum(weights * eps_over_kt * sigma_ij3))
        self._m2e2s3 = float(np.sum(weights * eps_over_kt**2 * sigma_ij3))

        # Eq. A.18-A.19.
        from .pcsaft import A_UNIVERSAL, B_UNIVERSAL

        ratio_1 = (self._mbar - 1.0) / self._mbar
        ratio_2 = ratio_1 * (self._mbar - 2.0) / self._mbar
        self._a_bar = A_UNIVERSAL[0] + ratio_1 * A_UNIVERSAL[1] + ratio_2 * A_UNIVERSAL[2]
        self._b_bar = B_UNIVERSAL[0] + ratio_1 * B_UNIVERSAL[1] + ratio_2 * B_UNIVERSAL[2]

        # Association (ADR-0018), in the same single variable. ``None`` unless
        # some component carries sites, and then no association code runs at
        # all - which is what keeps the ADR-0014/ADR-0015 numbers bit-identical.
        setup = (
            None
            if association is None
            else build_setup(
                temperature_K=self.temperature, sigma_A=sigma_A, d=d, association=association
            )
        )
        self.association: AssociationIsotherm | None = (
            None
            if setup is None
            else AssociationIsotherm(setup, x=x, ratio_2=float(ratios[2]), m3=m3)
        )

    # -- the model in one variable -----------------------------------------

    def derivatives(self, eta: np.ndarray, *, second: bool) -> _EtaDerivatives:
        """Return ``(a, a', a'')`` at every ``eta``; ``a''`` is zero if not asked."""
        eta = np.asarray(eta, dtype=float)
        u = 1.0 - eta
        column = eta[..., None]
        u_column = u[..., None]

        a_hs = self._a1_const * eta / u + self._a2_const * eta / u**2 + self._a3_const * np.log(u)
        a_hs1 = self._a1_const / u**2 + self._a2_const * (1.0 + eta) / u**3 - self._a3_const / u

        g = 1.0 / u_column + self._b_i * column / u_column**2 + self._c_i * column**2 / u_column**3
        g1 = (
            1.0 / u_column**2
            + self._b_i * (1.0 + column) / u_column**3
            + self._c_i * (2.0 * column + column**2) / u_column**4
        )
        g_ratio = g1 / g

        a_hc = self._mbar * a_hs - (self._chain_weight * np.log(g)).sum(axis=-1)
        a_hc1 = self._mbar * a_hs1 - (self._chain_weight * g_ratio).sum(axis=-1)

        powers = column**_POWERS
        i1 = (self._a_bar * powers).sum(axis=-1)
        i2 = (self._b_bar * powers).sum(axis=-1)
        shifted = np.concatenate((np.zeros_like(column), column ** _POWERS[:-1]), axis=-1)
        i1_1 = (self._a_bar * _POWERS * shifted).sum(axis=-1)
        i2_1 = (self._b_bar * _POWERS * shifted).sum(axis=-1)

        mbar = self._mbar
        gap = u * (2.0 - eta)
        gap_d = 2.0 * eta - 3.0
        q_poly = 2.0 * eta**3 + 12.0 * eta**2 - 48.0 * eta + 40.0

        denominator = (
            1.0
            + mbar * (8.0 * eta - 2.0 * eta**2) / u**4
            + (1.0 - mbar) * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4) / gap**2
        )
        denominator_1 = (
            mbar * (8.0 + 20.0 * eta - 4.0 * eta**2) / u**5 + (1.0 - mbar) * q_poly / gap**3
        )
        c1 = 1.0 / denominator
        c1_1 = -denominator_1 / denominator**2

        f = -2.0 * math.pi * i1 * self._m2es3 - math.pi * mbar * c1 * i2 * self._m2e2s3
        f1 = (
            -2.0 * math.pi * i1_1 * self._m2es3
            - math.pi * mbar * (c1_1 * i2 + c1 * i2_1) * self._m2e2s3
        )

        a = a_hc + eta * f / self._m3
        a1 = a_hc1 + (f + eta * f1) / self._m3

        if not second:
            if self.association is not None:
                a_assoc, a1_assoc, _ = self.association.derivatives(eta, second=False)
                a = a + a_assoc
                a1 = a1 + a1_assoc
            return _EtaDerivatives(a=a, a1=a1, a2=np.zeros_like(a1))

        a_hs2 = (
            2.0 * self._a1_const / u**3
            + self._a2_const * (4.0 + 2.0 * eta) / u**4
            - self._a3_const / u**2
        )
        g2 = (
            2.0 / u_column**3
            + self._b_i * (4.0 + 2.0 * column) / u_column**4
            + self._c_i * (2.0 + 8.0 * column + 2.0 * column**2) / u_column**5
        )
        a_hc2 = self._mbar * a_hs2 - (self._chain_weight * (g2 / g - g_ratio**2)).sum(axis=-1)

        twice_shifted = np.concatenate(
            (np.zeros_like(column), np.zeros_like(column), column ** _POWERS[:-2]), axis=-1
        )
        i1_2 = (self._a_bar * _POWERS * (_POWERS - 1.0) * twice_shifted).sum(axis=-1)
        i2_2 = (self._b_bar * _POWERS * (_POWERS - 1.0) * twice_shifted).sum(axis=-1)

        denominator_2 = (
            mbar * (60.0 + 72.0 * eta - 12.0 * eta**2) / u**6
            + (1.0 - mbar)
            * ((6.0 * eta**2 + 24.0 * eta - 48.0) * gap - 3.0 * q_poly * gap_d)
            / gap**4
        )
        c1_2 = -denominator_2 / denominator**2 + 2.0 * denominator_1**2 / denominator**3

        f2 = (
            -2.0 * math.pi * i1_2 * self._m2es3
            - math.pi * mbar * (c1_2 * i2 + 2.0 * c1_1 * i2_1 + c1 * i2_2) * self._m2e2s3
        )
        a2 = a_hc2 + (2.0 * f1 + eta * f2) / self._m3
        if self.association is not None:
            a_assoc, a1_assoc, a2_assoc = self.association.derivatives(eta, second=True)
            a = a + a_assoc
            a1 = a1 + a1_assoc
            a2 = a2 + a2_assoc
        return _EtaDerivatives(a=a, a1=a1, a2=a2)

    def pressure(self, eta: np.ndarray) -> np.ndarray:
        """Return ``P = Z rho R T`` in Pa at every ``eta``."""
        derivatives = self.derivatives(eta, second=False)
        z_factor = 1.0 + eta * derivatives.a1
        return z_factor * (self.density_per_eta * eta) * R_J_PER_MOL_K * self.temperature

    def pressure_and_slope(self, eta: float) -> tuple[float, float]:
        """Return ``(P, dP/drho)`` at one ``eta``.

        ``dP/drho = R T (1 + 2 eta a' + eta^2 a'')`` - see the module
        docstring; ``dP/deta`` is ``dP/drho`` times
        :attr:`density_per_eta`.
        """
        point = np.array([eta], dtype=float)
        derivatives = self.derivatives(point, second=True)
        a1 = float(derivatives.a1[0])
        a2 = float(derivatives.a2[0])
        z_factor = 1.0 + eta * a1
        pressure = z_factor * (self.density_per_eta * eta) * R_J_PER_MOL_K * self.temperature
        slope = R_J_PER_MOL_K * self.temperature * (1.0 + 2.0 * eta * a1 + eta**2 * a2)
        return pressure, slope


class DensityRoots(NamedTuple):
    """Result of one density-root solve.

    Attributes:
        densities: Admissible molar densities in mol/m^3, sorted ascending.
            The first is the vapour-like root, the last the liquid-like one.
        bracket_count: How many sign changes the scan grid resolved, including
            the mechanically unstable one that was then discarded. The
            difference between this and ``len(densities)`` is the number of
            spinodal-branch roots rejected.
        max_relative_residual: Worst ``|P_model - P| / P`` over the returned
            roots.
    """

    densities: tuple[float, ...]
    bracket_count: int
    max_relative_residual: float


def solve_density_roots(
    isotherm: PCSAFTIsotherm,
    pressure_Pa: float,
    *,
    state_description: str = "",
) -> DensityRoots:
    """Return the mechanically stable density roots of ``isotherm`` at ``pressure_Pa``.

    Args:
        isotherm: The PC-SAFT isotherm at the requested ``(T, x)``.
        pressure_Pa: Target pressure in Pa; must be positive.
        state_description: Text appended to the error raised when no root
            exists, so the caller sees which state failed.

    Returns:
        A :class:`DensityRoots` with at least one density.

    Raises:
        ModelError: If the scan finds no sign change, or finds only
            mechanically unstable roots.
    """
    target = float(pressure_Pa)
    grid = _ETA_GRID
    residual = isotherm.pressure(grid) - target

    finite = np.isfinite(residual)
    brackets: list[tuple[float, float]] = []
    exact: list[float] = []
    for index in range(grid.size - 1):
        if not (finite[index] and finite[index + 1]):
            continue
        left = float(residual[index])
        right = float(residual[index + 1])
        if left == 0.0:
            exact.append(float(grid[index]))
            continue
        if right == 0.0:
            continue
        if left < 0.0 < right or right < 0.0 < left:
            brackets.append((float(grid[index]), float(grid[index + 1])))
    if finite[-1] and float(residual[-1]) == 0.0:
        exact.append(float(grid[-1]))

    bracket_count = len(brackets) + len(exact)
    if bracket_count == 0:
        raise ModelError(
            "PC-SAFT found no density root: the isotherm never crosses "
            f"P = {target!r} Pa on the scan grid ({ETA_GRID_DESCRIPTION})"
            f"{state_description}."
        )

    candidates = [_refine(isotherm, low, high, target) for low, high in brackets]
    candidates.extend(exact)

    densities: list[float] = []
    worst_residual = 0.0
    for eta in sorted(candidates):
        pressure, slope = isotherm.pressure_and_slope(eta)
        if not math.isfinite(slope) or slope <= 0.0:
            continue
        densities.append(isotherm.density_per_eta * eta)
        worst_residual = max(worst_residual, abs(pressure - target) / target)

    if not densities:
        raise ModelError(
            "PC-SAFT found only mechanically unstable density roots (dP/drho <= 0) at "
            f"P = {target!r} Pa{state_description}."
        )

    return DensityRoots(
        densities=tuple(densities),
        bracket_count=bracket_count,
        max_relative_residual=worst_residual,
    )


def _refine(isotherm: PCSAFTIsotherm, low: float, high: float, target: float) -> float:
    """Safeguarded Newton on ``P_model(eta) - P`` inside a sign-changing bracket.

    Newton is taken when the step stays strictly inside the current bracket;
    otherwise the step is a bisection. The bracket therefore never widens and
    the iteration cannot leave the root, which plain Newton can do on the steep
    liquid branch.

    The iterate with the smallest ``|residual|`` is what is returned, not the
    last one. On the steep liquid branch one ``ulp`` of ``eta`` is already
    worth a few times ``1e-12 P`` at a low target pressure, so the last two
    iterates straddle the root with different residuals and there is no reason
    to keep the worse one.
    """
    lo, hi = low, high
    f_lo = float(isotherm.pressure(np.array([lo]))[0]) - target
    eta = 0.5 * (lo + hi)
    best_eta = eta
    best_residual = math.inf
    for _ in range(_MAX_REFINE_ITERATIONS):
        value, slope = isotherm.pressure_and_slope(eta)
        residual = value - target
        if abs(residual) < best_residual:
            best_eta, best_residual = eta, abs(residual)
        if residual == 0.0:
            return eta
        if (f_lo < 0.0) == (residual < 0.0):
            lo, f_lo = eta, residual
        else:
            hi = eta
        if abs(residual) <= _ROOT_RTOL * target:
            return eta
        width = hi - lo
        if width <= 2.0 * np.finfo(float).eps * max(abs(eta), 1e-300):
            break
        derivative = slope * isotherm.density_per_eta
        step = eta - residual / derivative if derivative != 0.0 else math.nan
        if math.isfinite(step) and lo < step < hi:
            eta = step
        else:
            eta = 0.5 * (lo + hi)
    return best_eta


def build_isotherm(
    *,
    temperature_K: float,
    composition: Sequence[float] | np.ndarray,
    m: np.ndarray,
    sigma_A: np.ndarray,
    epsilon_k_K: np.ndarray,
    kij: np.ndarray,
    association: Sequence[PCSAFTAssociationRecord | None] | None = None,
) -> PCSAFTIsotherm:
    """Build a :class:`PCSAFTIsotherm` from already-validated inputs."""
    return PCSAFTIsotherm(
        temperature_K=temperature_K,
        composition=np.asarray(composition, dtype=float),
        m=m,
        sigma_A=sigma_A,
        epsilon_k_K=epsilon_k_K,
        kij=kij,
        association=association,
    )
