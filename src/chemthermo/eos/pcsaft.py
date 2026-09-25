"""PC-SAFT equation of state (Gross & Sadowski 2001, 2002).

Source
------
J. Gross and G. Sadowski, "Perturbed-Chain SAFT: An Equation of State Based on
a Perturbation Theory for Chain Molecules", Ind. Eng. Chem. Res. 40 (2001)
1244-1260 (DOI 10.1021/ie0003887). Equation numbers below are that paper's
appendix numbering, and the hard-chain and dispersion contributions written
out below are its Eq. (A.4) and (A.10).

The **association** contribution of J. Gross and G. Sadowski, "Application of
the Perturbed-Chain SAFT Equation of State to Associating Systems", Ind. Eng.
Chem. Res. 41 (2002) 5510-5515 (DOI 10.1021/ie010954d) is implemented too, in
the sibling module :mod:`chemthermo.eos._pcsaft_association`, which carries
its equations and its derivative route (ADR-0018). It is evaluated only when
some component in the mixture carries association sites; when none does, the
association code does not run at all and this module's numbers are those of
ADR-0014 down to the last bit.

Still out of scope: the **polar** (dipolar / quadrupolar) terms, so a strongly
polar non-associating compound must not be modelled here.

Units and conventions
---------------------
- ``sigma`` is in Angstrom and ``epsilon/k_B`` in K, as published.
- Densities passed in are **molar**, in mol/m^3; molar volumes are in m^3/mol.
  Internally the model works with the number density of *molecules*,
  ``rho_A3 = rho_molar * N_A * 1e-30`` in molecules per Angstrom^3, which is
  what makes the Angstrom-based sigma and d consistent.
- ``k_B = 1.380649e-23`` J/K and ``N_A = 6.02214076e23`` /mol are the exact SI
  definitions; ``R = N_A k_B = 8.31446261815324`` J/(mol K) is used for
  pressure. (``chemthermo.units.R_J_PER_MOL_K`` is the same constant rounded
  to 8.314462618, i.e. 2.2e-11 relative away; the exact product is used here so
  that PC-SAFT pressures agree with reference implementations that do the
  same.)
- Index convention: ``i``, ``j`` run over components, ``n`` over the four
  ``zeta`` moments (0..3) and over the seven powers of ``eta`` (0..6).
- Everything returned is dimensionless and *reduced*: ``A^res/(R T)`` per mole,
  the compressibility factor ``Z``, and natural logarithms of fugacity
  coefficients.

The equations, in the notation used below
-----------------------------------------
Reduced residual Helmholtz energy is the sum of a hard-chain and a dispersion
contribution (Eq. A.3):

    a_res = A^res / (N k T) = A^res / (R T) per mole = a_hc + a_disp [+ a_assoc]

*Temperature-dependent segment diameter* (Eq. A.9):

    d_i = sigma_i * [1 - 0.12 * exp(-3 eps_i / (k T))]

*Moments of the segment-size distribution* (Eq. A.8), with ``rho`` the number
density of molecules:

    zeta_n = (pi / 6) * rho * sum_i x_i m_i d_i^n ,   n = 0, 1, 2, 3

``eta = zeta_3`` is the packing fraction. Write ``U = 1 - eta``.

*Boublik-Mansoori hard-sphere mixture term* (Eq. A.6):

    a_hs = (1 / zeta_0) * [ 3 zeta_1 zeta_2 / U
                            + zeta_2^3 / (zeta_3 U^2)
                            + (zeta_2^3 / zeta_3^2 - zeta_0) ln U ]

*Radial distribution function of the hard-sphere reference at contact*
(Eq. A.7), with ``c_ij = d_i d_j / (d_i + d_j)``:

    g_ij = 1/U + 3 c_ij zeta_2 / U^2 + 2 c_ij^2 zeta_2^2 / U^3

*Hard chain* (Eq. A.4), with ``mbar = sum_i x_i m_i``:

    a_hc = mbar * a_hs - sum_i x_i (m_i - 1) ln g_ii

*Dispersion* (Eq. A.10), with the one-fluid mixing rules of Eq. A.12-A.14
(``sigma_ij = (sigma_i + sigma_j)/2``,
``eps_ij = sqrt(eps_i eps_j) (1 - k_ij)``):

    m2es3  = sum_i sum_j x_i x_j m_i m_j (eps_ij / kT)   sigma_ij^3
    m2e2s3 = sum_i sum_j x_i x_j m_i m_j (eps_ij / kT)^2 sigma_ij^3

    a_disp = -2 pi rho I1 m2es3 - pi rho mbar C1 I2 m2e2s3

*Perturbation integrals* (Eq. A.16-A.19), with the 21 + 21 universal constants
``a_ni``, ``b_ni`` of the paper's Table 1:

    I1 = sum_{n=0..6} abar_n eta^n ,  I2 = sum_{n=0..6} bbar_n eta^n
    abar_n = a_0n + (mbar-1)/mbar a_1n
                  + (mbar-1)/mbar (mbar-2)/mbar a_2n      (same for bbar_n)

*Compressibility of the dispersion term* (Eq. A.11). **The paper's printing of
A.11 has a known typographical error** - the outer exponent ``-1`` is missing
on the right-hand side, so as printed it says ``C1 = (...)^-1 = (...)``. The
correct form, confirmed by NIST TRC's PC-SAFT page
(https://trc.nist.gov/TDE/TDE_Help/eos-PC-SAFT.htm, which states the erratum
explicitly) and by every reference implementation consulted, is the reciprocal:

    C1 = [ 1 + mbar (8 eta - 2 eta^2) / (1 - eta)^4
             + (1 - mbar) (20 eta - 27 eta^2 + 12 eta^3 - 2 eta^4)
               / ((1 - eta)(2 - eta))^2 ]^(-1)

That is what is implemented (see :func:`_c1_terms`), and
``tests/test_pcsaft.py`` pins both the reciprocal shape (``C1 -> 1`` as
``eta -> 0``) and its eta-derivative against finite differences.

Derivatives
-----------
All derivatives are analytic, and all of them are assembled by one chain rule
over a small set of intermediates rather than by transcribing the appendix's
individual derivative equations. Two facts do the work:

1. Every ``zeta_n`` is proportional to the density at fixed composition, so
   ``rho d/drho`` of any function of the ``zeta`` moments is
   ``sum_n zeta_n * (partial / partial zeta_n)``. Hence

       Z = 1 + rho (d a_res / d rho)_{T,x}

   needs only the gradient of ``a_hs`` and of ``g_ii`` with respect to the
   ``zeta`` moments, plus the explicit density prefactor in ``a_disp``.

2. ``d zeta_n / d x_k = (pi/6) rho m_k d_k^n`` at fixed T and fixed total
   density, and ``a_disp`` depends on the composition only through
   ``eta``, ``mbar``, ``m2es3`` and ``m2e2s3``. So the unconstrained
   composition derivative ``partial a_res / partial x_k`` is again one chain
   rule.

The residual chemical potential then follows from the definition
``mu_i^res/kT = (partial (n a_res) / partial n_i)_{T,V}`` (Eq. A.33):

    mu_i^res / kT = a_res + (Z - 1)
                    + (partial a_res / partial x_i)
                    - sum_j x_j (partial a_res / partial x_j)

    ln phi_i = mu_i^res / kT - ln Z                              (Eq. A.32)

Two state specifications, one model (ADR-0015)
----------------------------------------------
The methods above take the state as ``(T, rho, x)`` or ``(T, v, x)``: that is
the natural specification for a residual-Helmholtz model and it is what
``EOSProtocol`` asks for. The flash and the stability test speak the other
one - ``(T, P, x)`` plus a *phase label* - through
``chemthermo.models.EquationOfState``. Since the ``pcsaft-density-roots-flash``
slice this class implements **both**:

    fugacity_coefficients(mixture=, temperature_K=, pressure_Pa=,
                          composition=, phase=)

solves ``P_model(T, rho, x) = P`` for its mechanically stable roots
(:mod:`chemthermo.eos._pcsaft_density`), takes the lowest-density root for
``phase="vapor"`` and the highest for ``phase="liquid"``, and returns
``exp(ln phi)`` there. Those are the Peng-Robinson label semantics verbatim
(largest ``Z`` is the vapour), including the single-root case where both
labels name the same state. :meth:`density_roots` exposes the root set itself.

That is the whole wiring: ``stability_tp`` and ``flash_tp`` are untouched by
this slice. The tangent-plane evaluator already holds the two compressibility
branches as competing phase candidates and keeps the lowest-Gibbs one
(ADR-0005, ADR-0012), which is exactly the right rule for PC-SAFT's roots too.
See ADR-0015 for why the fugacity interface, and not a new density-root
candidate type, is the seam this slice uses.

Still out of scope: the polar terms, and any temperature derivative (so no
caloric properties) - in the association contribution as in the rest.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple, Sequence

import numpy as np

from .._eos_memo import MISS, active_memo, composition_key
from ..core import Mixture
from ..data import normalize_name
from ..exceptions import CompositionError, InputRangeError, ModelError
from ..models._kij import KijInput, KijPairs, canonicalize_kij, kij_matrix
from ..models.base import KAPPA_LIQUID_THRESHOLD, EquationOfState
from ..parameters.pcsaft import (
    PCSAFTAssociationRecord,
    PCSAFTParameterError,
    PCSAFTParameters,
    get_pcsaft_parameters,
)
from ..validation import (
    COMPOSITION_SUM_TOL,
    validate_fractions,
    validate_pressure,
    validate_temperature,
)
from ._pcsaft_association import AssociationState, build_setup
from ._pcsaft_association import temperature_derivative as assoc_temperature_derivative
from ._pcsaft_density import DensityRoots, PCSAFTIsotherm, build_isotherm, solve_density_roots
from .api import EOSProtocol
from .registry import register_eos

#: Phase labels accepted by :meth:`PCSAFTEOS.fugacity_coefficients`, with the
#: Peng-Robinson semantics: "vapor" is the largest ``Z`` (lowest density) root,
#: "liquid" the smallest ``Z`` (highest density) one.
_VAPOR = "vapor"
_LIQUID = "liquid"

#: Exact SI definitions (2019 redefinition).
BOLTZMANN_J_PER_K = 1.380649e-23
AVOGADRO_PER_MOL = 6.02214076e23
#: R = N_A k_B, used for pressure. See the module docstring on why this is the
#: exact product rather than ``chemthermo.units.R_J_PER_MOL_K``.
R_J_PER_MOL_K = AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K

#: Universal model constants of Gross & Sadowski (2001), Table 1: rows are
#: ``a_0n`` / ``a_1n`` / ``a_2n`` and columns are ``n = 0..6``.
#:
#: Provenance: the paper is paywalled and was not read directly; these 42
#: values were taken from two independent sources that agree digit for digit -
#: teqp (NIST, MIT) ``src/data/PCSAFT.cpp``, namespace
#: ``PCSAFTMatrices::GrossSadowski2001``, and the table in the Wikipedia
#: article "PC-SAFT" (section "Dispersion Term"). See validation Case P-0.
A_UNIVERSAL = np.array(
    [
        [
            0.9105631445,
            0.6361281449,
            2.6861347891,
            -26.547362491,
            97.759208784,
            -159.59154087,
            91.297774084,
        ],
        [
            -0.3084016918,
            0.1860531159,
            -2.5030047259,
            21.419793629,
            -65.255885330,
            83.318680481,
            -33.746922930,
        ],
        [
            -0.0906148351,
            0.4527842806,
            0.5962700728,
            -1.7241829131,
            -4.1302112531,
            13.776631870,
            -8.6728470368,
        ],
    ],
    dtype=float,
)

#: Universal model constants ``b_0n`` / ``b_1n`` / ``b_2n``; same provenance.
B_UNIVERSAL = np.array(
    [
        [
            0.7240946941,
            2.2382791861,
            -4.0025849485,
            -21.003576815,
            26.855641363,
            206.55133841,
            -355.60235612,
        ],
        [
            -0.5755498075,
            0.6995095521,
            3.8925673390,
            -17.215471648,
            192.67226447,
            -161.82646165,
            -165.20769346,
        ],
        [
            0.0976883116,
            -0.2557574982,
            -9.1558561530,
            20.642075974,
            -38.804430052,
            93.626774077,
            -29.666905585,
        ],
    ],
    dtype=float,
)

_POWERS_ZETA = np.arange(4)
_POWERS_ETA = np.arange(7)


class _PCSAFTState(NamedTuple):
    """Everything one evaluation of the model produces at a state.

    ``a_res`` is ``A^res/(R T)`` per mole, ``z_minus_one`` is
    ``rho (d a_res / d rho)_{T,x} = Z - 1`` and ``da_dx`` holds the
    unconstrained composition derivatives ``partial a_res / partial x_k`` at
    fixed temperature and total density. ``z_hc``, ``z_disp`` and ``z_assoc``
    are the three contributions to ``z_minus_one``; they are kept apart so each
    term's density derivative can be finite-difference tested on its own.

    ``a_assoc`` and ``z_assoc`` are ``0.0`` and ``site_fractions`` is ``None``
    when no component associates - the association code never runs then
    (ADR-0018), which is what keeps every ADR-0014 number bit-identical.
    """

    a_res: float
    z_minus_one: float
    da_dx: np.ndarray
    eta: float
    a_hc: float
    a_disp: float
    z_hc: float
    z_disp: float
    a_assoc: float = 0.0
    z_assoc: float = 0.0
    site_fractions: np.ndarray | None = None


def _c1_terms(eta: float, mbar: float) -> tuple[float, float, float]:
    """Return ``(C1, dC1/deta, dC1/dmbar)`` for the typo-corrected Eq. A.11.

    ``C1`` is the *reciprocal* of the bracket (see the module docstring on the
    known error in the paper's printing of A.11). ``dC1/deta`` is Eq. A.31.
    """
    one_minus = 1.0 - eta
    gap = one_minus * (2.0 - eta)
    term_chain = (8.0 * eta - 2.0 * eta**2) / one_minus**4
    term_ring = (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4) / gap**2

    c1 = 1.0 / (1.0 + mbar * term_chain + (1.0 - mbar) * term_ring)
    c1_deta = -(c1**2) * (
        mbar * (-4.0 * eta**2 + 20.0 * eta + 8.0) / one_minus**5
        + (1.0 - mbar) * (2.0 * eta**3 + 12.0 * eta**2 - 48.0 * eta + 40.0) / gap**3
    )
    c1_dmbar = -(c1**2) * (term_chain - term_ring)
    return c1, c1_deta, c1_dmbar


def _evaluate(
    *,
    temperature_K: float,
    density_mol_m3: float,
    x: np.ndarray,
    m: np.ndarray,
    sigma_A: np.ndarray,
    epsilon_k_K: np.ndarray,
    kij: np.ndarray,
    association: Sequence[PCSAFTAssociationRecord | None] | None = None,
) -> _PCSAFTState:
    """Evaluate the PC-SAFT model at one state.

    Pure function of arrays; see the module docstring for every equation and
    for the two chain rules that give the derivatives. ``association`` is one
    record (or ``None``) per component; ``None`` for the argument itself, or a
    sequence whose entries are all ``None``, skips the association term
    entirely (ADR-0018).
    """
    rho_a3 = density_mol_m3 * AVOGADRO_PER_MOL * 1e-30  # molecules / Angstrom^3

    # Eq. A.9: temperature-dependent segment diameter, in Angstrom.
    d = sigma_A * (1.0 - 0.12 * np.exp(-3.0 * epsilon_k_K / temperature_K))
    mbar = float(x @ m)
    m_minus_1 = m - 1.0

    # Eq. A.8: the four moments. Each is linear in the density.
    moments = (math.pi / 6.0) * ((x * m)[None, :] * d[None, :] ** _POWERS_ZETA[:, None]).sum(axis=1)
    zeta = moments * rho_a3
    z0, z1, z2, z3 = (float(value) for value in zeta)
    eta = z3
    if not (0.0 < eta < 1.0):
        raise ModelError(
            f"PC-SAFT packing fraction eta = {eta!r} is outside (0, 1); the model has no "
            "meaning at this density."
        )
    u = 1.0 - eta
    ln_u = math.log(u)

    # Eq. A.6 and its gradient with respect to the four moments.
    a_hs = (3.0 * z1 * z2 / u + z2**3 / (z3 * u * u) + (z2**3 / z3**2 - z0) * ln_u) / z0
    dahs_dzeta = np.array(
        [
            -a_hs / z0 - ln_u / z0,
            3.0 * z2 / u / z0,
            (3.0 * z1 / u + 3.0 * z2**2 / (z3 * u * u) + 3.0 * z2**2 / z3**2 * ln_u) / z0,
            (
                3.0 * z1 * z2 / (u * u)
                + z2**3 * (2.0 / (z3 * u**3) - 1.0 / (z3**2 * u * u))
                - 2.0 * z2**3 * ln_u / z3**3
                - (z2**3 / z3**2 - z0) / u
            )
            / z0,
        ],
        dtype=float,
    )

    # Eq. A.7 at contact for the like pairs (c_ii = d_i / 2) and its gradient.
    c_ii = 0.5 * d
    g_ii = 1.0 / u + c_ii * 3.0 * z2 / u**2 + c_ii**2 * 2.0 * z2**2 / u**3
    dg_dz2 = 3.0 * c_ii / u**2 + 4.0 * c_ii**2 * z2 / u**3
    dg_dz3 = 1.0 / u**2 + 6.0 * c_ii * z2 / u**3 + 6.0 * c_ii**2 * z2**2 / u**4

    # Eq. A.4.
    a_hc = mbar * a_hs - float(np.sum(x * m_minus_1 * np.log(g_ii)))

    # Eq. A.12-A.14: one-fluid mixing rules for the dispersion term.
    sigma_ij3 = (0.5 * (sigma_A[:, None] + sigma_A[None, :])) ** 3
    eps_over_kt = np.sqrt(np.outer(epsilon_k_K, epsilon_k_K)) * (1.0 - kij) / temperature_K
    xm = x * m
    weights = np.outer(xm, xm)
    m2es3 = float(np.sum(weights * eps_over_kt * sigma_ij3))
    m2e2s3 = float(np.sum(weights * eps_over_kt**2 * sigma_ij3))

    # Eq. A.18-A.19: the eta-polynomial coefficients and their mbar derivative.
    ratio_1 = (mbar - 1.0) / mbar
    ratio_2 = ratio_1 * (mbar - 2.0) / mbar
    a_bar = A_UNIVERSAL[0] + ratio_1 * A_UNIVERSAL[1] + ratio_2 * A_UNIVERSAL[2]
    b_bar = B_UNIVERSAL[0] + ratio_1 * B_UNIVERSAL[1] + ratio_2 * B_UNIVERSAL[2]
    dratio_1 = 1.0 / mbar**2
    dratio_2 = 3.0 / mbar**2 - 4.0 / mbar**3
    da_bar = dratio_1 * A_UNIVERSAL[1] + dratio_2 * A_UNIVERSAL[2]
    db_bar = dratio_1 * B_UNIVERSAL[1] + dratio_2 * B_UNIVERSAL[2]

    eta_powers = eta**_POWERS_ETA
    eta_powers_shifted = np.concatenate(([0.0], eta ** _POWERS_ETA[:-1]))
    i1 = float(a_bar @ eta_powers)
    i2 = float(b_bar @ eta_powers)
    i1_deta = float((a_bar * _POWERS_ETA) @ eta_powers_shifted)
    i2_deta = float((b_bar * _POWERS_ETA) @ eta_powers_shifted)
    i1_dmbar = float(da_bar @ eta_powers)
    i2_dmbar = float(db_bar @ eta_powers)

    c1, c1_deta, c1_dmbar = _c1_terms(eta, mbar)

    # Eq. A.10, written as a_disp = rho * f(eta, mbar, m2es3, m2e2s3).
    f = -2.0 * math.pi * i1 * m2es3 - math.pi * mbar * c1 * i2 * m2e2s3
    a_disp = rho_a3 * f
    f_deta = (
        -2.0 * math.pi * i1_deta * m2es3 - math.pi * mbar * (c1_deta * i2 + c1 * i2_deta) * m2e2s3
    )
    f_dmbar = (
        -2.0 * math.pi * i1_dmbar * m2es3
        - math.pi * (c1 * i2 + mbar * c1_dmbar * i2 + mbar * c1 * i2_dmbar) * m2e2s3
    )
    f_dm2es3 = -2.0 * math.pi * i1
    f_dm2e2s3 = -math.pi * mbar * c1 * i2

    a_res = a_hc + a_disp

    # Density derivative: every zeta_n is proportional to rho, and a_disp has
    # one explicit factor of rho on top of its eta dependence.
    dahs_drho = float(zeta @ dahs_dzeta)
    dlng_drho = (z2 * dg_dz2 + z3 * dg_dz3) / g_ii
    dahc_drho = mbar * dahs_drho - float(np.sum(x * m_minus_1 * dlng_drho))
    dadisp_drho = a_disp + rho_a3 * f_deta * eta
    z_minus_one = dahc_drho + dadisp_drho

    # Unconstrained composition derivatives at fixed T and total density.
    dzeta_dx = (math.pi / 6.0) * rho_a3 * (m[None, :] * d[None, :] ** _POWERS_ZETA[:, None])
    dahs_dx = dahs_dzeta @ dzeta_dx
    dg_dx = dg_dz2[:, None] * dzeta_dx[2][None, :] + dg_dz3[:, None] * dzeta_dx[3][None, :]
    dahc_dx = (
        m * a_hs
        + mbar * dahs_dx
        - m_minus_1 * np.log(g_ii)
        - ((x * m_minus_1 / g_ii)[:, None] * dg_dx).sum(axis=0)
    )
    dm2es3_dx = 2.0 * m * (xm[None, :] * eps_over_kt * sigma_ij3).sum(axis=1)
    dm2e2s3_dx = 2.0 * m * (xm[None, :] * eps_over_kt**2 * sigma_ij3).sum(axis=1)
    dadisp_dx = rho_a3 * (
        f_deta * dzeta_dx[3] + f_dmbar * m + f_dm2es3 * dm2es3_dx + f_dm2e2s3 * dm2e2s3_dx
    )

    # Association (Gross & Sadowski 2002), ADR-0018. The setup is ``None``
    # unless some component carries sites, and then nothing below runs.
    setup = (
        None
        if association is None
        else build_setup(temperature_K=temperature_K, sigma_A=sigma_A, d=d, association=association)
    )
    if setup is None:
        return _PCSAFTState(
            a_res=float(a_res),
            z_minus_one=float(z_minus_one),
            da_dx=dahc_dx + dadisp_dx,
            eta=eta,
            a_hc=float(a_hc),
            a_disp=float(a_disp),
            z_hc=float(dahc_drho),
            z_disp=float(dadisp_drho),
        )

    assoc = AssociationState(setup, rho_a3=rho_a3, x=x, zeta_2=z2, eta=eta, dzeta_dx=dzeta_dx)
    return _PCSAFTState(
        a_res=float(a_res + assoc.a_assoc),
        z_minus_one=float(z_minus_one + assoc.z_assoc),
        da_dx=dahc_dx + dadisp_dx + assoc.da_dx,
        eta=eta,
        a_hc=float(a_hc),
        a_disp=float(a_disp),
        z_hc=float(dahc_drho),
        z_disp=float(dadisp_drho),
        a_assoc=float(assoc.a_assoc),
        z_assoc=float(assoc.z_assoc),
        site_fractions=assoc.site_fractions,
    )


def _temperature_derivative(
    *,
    temperature_K: float,
    density_mol_m3: float,
    x: np.ndarray,
    m: np.ndarray,
    sigma_A: np.ndarray,
    epsilon_k_K: np.ndarray,
    kij: np.ndarray,
    association: Sequence[PCSAFTAssociationRecord | None] | None = None,
) -> float:
    """``(d a_res / d T)_{rho, x}`` in 1/K (ADR-0034 and its C2 amendment).

    A separate function rather than a new output of :func:`_evaluate`, so that
    no number :func:`_evaluate` produces can move. Temperature enters the
    non-associating model in exactly two places:

    - the segment diameter, Eq. A.9: ``d_i = sigma_i (1 - 0.12 exp(-3 eps_i/T))``,
      so ``d d_i / dT = -0.36 sigma_i (eps_i / T^2) exp(-3 eps_i / T)``, which
      moves every ``zeta_n`` (Eq. A.8) and the contact values ``g_ii`` (A.7);
    - ``eps_ij / (k T)`` in the dispersion mixing rules, A.12-A.13, so
      ``m2es3 ~ 1/T`` and ``m2e2s3 ~ 1/T^2``.

    Everything else is the chain rule through the same gradients
    :func:`_evaluate` forms (``d a_hs / d zeta``, ``d g / d zeta_2,3``,
    ``d f / d eta``). When a component carries association sites, the
    association term's derivative is added by
    :func:`chemthermo.eos._pcsaft_association.temperature_derivative`
    (Michelsen-Hendriks: an explicit partial at frozen site fractions).
    """
    t = temperature_K
    rho_a3 = density_mol_m3 * AVOGADRO_PER_MOL * 1e-30
    boltz = np.exp(-3.0 * epsilon_k_K / t)
    d = sigma_A * (1.0 - 0.12 * boltz)
    dd_dt = -0.36 * sigma_A * epsilon_k_K / (t * t) * boltz
    mbar = float(x @ m)
    m_minus_1 = m - 1.0
    xm = x * m

    powers = _POWERS_ZETA[:, None]
    zeta = (math.pi / 6.0) * rho_a3 * (xm[None, :] * d[None, :] ** powers).sum(axis=1)
    # d(d^n)/dT = n d^(n-1) dd/dT; the n = 0 row is identically zero.
    dpow = powers * np.where(powers > 0, d[None, :] ** np.maximum(powers - 1, 0), 0.0)
    dzeta_dt = (math.pi / 6.0) * rho_a3 * (xm[None, :] * dpow * dd_dt[None, :]).sum(axis=1)
    z0, z1, z2, z3 = (float(value) for value in zeta)
    eta = z3
    if not (0.0 < eta < 1.0):
        raise ModelError(
            f"PC-SAFT packing fraction eta = {eta!r} is outside (0, 1); the model has no "
            "meaning at this density."
        )
    u = 1.0 - eta
    ln_u = math.log(u)

    a_hs = (3.0 * z1 * z2 / u + z2**3 / (z3 * u * u) + (z2**3 / z3**2 - z0) * ln_u) / z0
    dahs_dzeta = np.array(
        [
            -a_hs / z0 - ln_u / z0,
            3.0 * z2 / u / z0,
            (3.0 * z1 / u + 3.0 * z2**2 / (z3 * u * u) + 3.0 * z2**2 / z3**2 * ln_u) / z0,
            (
                3.0 * z1 * z2 / (u * u)
                + z2**3 * (2.0 / (z3 * u**3) - 1.0 / (z3**2 * u * u))
                - 2.0 * z2**3 * ln_u / z3**3
                - (z2**3 / z3**2 - z0) / u
            )
            / z0,
        ],
        dtype=float,
    )
    dahs_dt = float(dahs_dzeta @ dzeta_dt)

    c_ii = 0.5 * d
    dc_dt = 0.5 * dd_dt
    g_ii = 1.0 / u + c_ii * 3.0 * z2 / u**2 + c_ii**2 * 2.0 * z2**2 / u**3
    dg_dc = 3.0 * z2 / u**2 + 4.0 * c_ii * z2**2 / u**3
    dg_dz2 = 3.0 * c_ii / u**2 + 4.0 * c_ii**2 * z2 / u**3
    dg_dz3 = 1.0 / u**2 + 6.0 * c_ii * z2 / u**3 + 6.0 * c_ii**2 * z2**2 / u**4
    dg_dt = dg_dc * dc_dt + dg_dz2 * dzeta_dt[2] + dg_dz3 * dzeta_dt[3]
    dahc_dt = mbar * dahs_dt - float(np.sum(x * m_minus_1 * dg_dt / g_ii))

    sigma_ij3 = (0.5 * (sigma_A[:, None] + sigma_A[None, :])) ** 3
    eps_over_kt = np.sqrt(np.outer(epsilon_k_K, epsilon_k_K)) * (1.0 - kij) / t
    weights = np.outer(xm, xm)
    m2es3 = float(np.sum(weights * eps_over_kt * sigma_ij3))
    m2e2s3 = float(np.sum(weights * eps_over_kt**2 * sigma_ij3))

    ratio_1 = (mbar - 1.0) / mbar
    ratio_2 = ratio_1 * (mbar - 2.0) / mbar
    a_bar = A_UNIVERSAL[0] + ratio_1 * A_UNIVERSAL[1] + ratio_2 * A_UNIVERSAL[2]
    b_bar = B_UNIVERSAL[0] + ratio_1 * B_UNIVERSAL[1] + ratio_2 * B_UNIVERSAL[2]
    eta_powers = eta**_POWERS_ETA
    eta_powers_shifted = np.concatenate(([0.0], eta ** _POWERS_ETA[:-1]))
    i1 = float(a_bar @ eta_powers)
    i2 = float(b_bar @ eta_powers)
    i1_deta = float((a_bar * _POWERS_ETA) @ eta_powers_shifted)
    i2_deta = float((b_bar * _POWERS_ETA) @ eta_powers_shifted)
    c1, c1_deta, _ = _c1_terms(eta, mbar)

    f_deta = (
        -2.0 * math.pi * i1_deta * m2es3 - math.pi * mbar * (c1_deta * i2 + c1 * i2_deta) * m2e2s3
    )
    # f = -2 pi I1 m2es3 - pi mbar C1 I2 m2e2s3, with m2es3 ~ 1/T, m2e2s3 ~ 1/T^2.
    df_dt = (
        f_deta * float(dzeta_dt[3])
        + 2.0 * math.pi * i1 * m2es3 / t
        + 2.0 * math.pi * mbar * c1 * i2 * m2e2s3 / t
    )
    total = float(dahc_dt + rho_a3 * df_dt)

    setup = (
        None
        if association is None
        else build_setup(temperature_K=t, sigma_A=sigma_A, d=d, association=association)
    )
    if setup is None:
        return total
    epsilon_ab = np.array(
        [0.0 if record is None else record.epsilon_ab_k_K for record in association or ()],
        dtype=float,
    )
    return total + assoc_temperature_derivative(
        setup,
        temperature_K=t,
        epsilon_ab_k_K=epsilon_ab,
        d=d,
        dd_dt=dd_dt,
        rho_a3=rho_a3,
        x=x,
        zeta_2=z2,
        eta=eta,
        dzeta2_dt=float(dzeta_dt[2]),
        dzeta3_dt=float(dzeta_dt[3]),
    )


@dataclass(frozen=True)
class PCSAFTEOS(EquationOfState, EOSProtocol):
    """PC-SAFT equation of state (ADR-0014, ADR-0015, ADR-0018).

    Hard chain plus dispersion (Gross & Sadowski 2001) and, for a component
    whose parameter record carries an ``association`` block, the association
    term of Gross & Sadowski (2002). :meth:`associates` reports whether the
    term is active; :meth:`residual_helmholtz_terms` and :meth:`site_fractions`
    expose it.

    Args:
        components: Component names, in the order every composition argument
            uses. Names are matched against the parameter set (and therefore
            the packaged databank) via ``chemthermo.data.normalize_name``.
            **Optional since ADR-0015**: leave it empty (``PCSAFTEOS()``) and
            the ``(T, P, x)`` methods take the order from the ``Mixture`` they
            are handed, which is what lets an instance be passed straight to
            ``stability_tp`` / ``flash_tp``. Give it, and a ``Mixture`` whose
            names differ is rejected rather than silently reordered. The
            ``(T, rho, x)`` methods always need it, since no ``Mixture``
            reaches them.
        parameters: Pure-component parameters. ``None`` (the default) uses the
            packaged Gross & Sadowski (2001) Table 1 set; supply a
            :class:`~chemthermo.parameters.PCSAFTParameters` built with
            ``from_records`` to override or to add a component.
        kij: Binary interaction parameter(s), same contract as
            ``PengRobinsonEOS.kij`` (ADR-0006): a scalar applied to every
            off-diagonal pair (never to the diagonal), or a mapping from an
            unordered pair of component names to a value, e.g.
            ``{("Methane", "n-Decane"): 0.03}``. Missing pairs default to
            ``0.0``. It enters through ``eps_ij = sqrt(eps_i eps_j)(1 - k_ij)``.

    There are two families of methods and they differ in how the state is
    specified, not in the model:

    - :meth:`residual_helmholtz`, :meth:`compressibility_factor`,
      :meth:`pressure_Pa` and :meth:`ln_fugacity_coefficients` take
      ``(T, molar density, x)`` (or molar volume). The caller chooses which
      root they are on. A state inside the mechanically unstable region can
      have ``Z <= 0``, for which fugacity coefficients do not exist;
      :meth:`ln_fugacity_coefficients` raises ``ModelError`` there rather than
      returning a ``nan``.
    - :meth:`fugacity_coefficients`, :meth:`density_roots` and
      :meth:`molar_volume` take ``(T, P, x)`` and solve for the density
      themselves. This is the ``chemthermo.models.EquationOfState`` interface
      the flash and the stability test use.
    """

    components: tuple[str, ...] = ()
    parameters: PCSAFTParameters | None = None
    kij: KijInput | KijPairs = 0.0
    name: str = "PC-SAFT"

    def __post_init__(self) -> None:
        object.__setattr__(self, "components", tuple(self.components))
        object.__setattr__(self, "kij", canonicalize_kij(self.kij, model="PCSAFTEOS"))

    def num_components(self) -> int:
        """Return the number of components the EOS instance was configured for.

        Zero for an instance built without ``components``, which takes its
        component order from the ``Mixture`` it is handed instead.
        """
        return len(self.components)

    def residual_helmholtz(
        self,
        *,
        temperature_K: float,
        volume_m3: float,
        composition: Sequence[float],
    ) -> float:
        """Return ``A^res / (R T)`` at ``(T, molar volume, x)``.

        Args:
            temperature_K: Temperature in K.
            volume_m3: **Molar** volume in m^3/mol (the reciprocal of the molar
                density). Must be positive.
            composition: Mole fractions, summing to 1 within
                ``COMPOSITION_SUM_TOL``.

        Returns:
            The dimensionless reduced residual Helmholtz energy.
        """
        volume = float(volume_m3)
        if not math.isfinite(volume) or volume <= 0.0:
            raise InputRangeError(f"Molar volume must be positive and finite (got {volume_m3!r}).")
        return self._state(temperature_K, 1.0 / volume, composition).a_res

    def residual_helmholtz_temperature_derivative(
        self,
        *,
        temperature_K: float,
        volume_m3: float,
        composition: Sequence[float],
    ) -> float:
        """Return ``d(A^res / (R T)) / dT`` at fixed molar volume and composition.

        Analytic, in 1/K (ADR-0034), including the association term for a
        mixture with association sites (Michelsen-Hendriks, no extra solve).
        teqp's ``get_Ar10`` is ``-T`` times this number.
        """
        volume = float(volume_m3)
        if not math.isfinite(volume) or volume <= 0.0:
            raise InputRangeError(f"Molar volume must be positive and finite (got {volume_m3!r}).")
        temperature = validate_temperature(temperature_K)
        m, sigma_A, epsilon_k_K = self.component_parameters()
        x = self._composition(composition, None)
        return _temperature_derivative(
            temperature_K=temperature,
            density_mol_m3=1.0 / volume,
            x=x,
            m=m,
            sigma_A=sigma_A,
            epsilon_k_K=epsilon_k_K,
            kij=self.kij_matrix(),
            association=self.association_parameters(),
        )

    def residual_properties(
        self,
        *,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> dict[str, float]:
        """Reduced residual properties at ``(T, molar density, x)`` (ADR-0034).

        All values are dimensionless; multiply by ``R T`` (energies) or ``R``
        (entropy) for SI. A residual property is the real fluid's value minus
        an ideal gas's, and **which** ideal gas matters for ``S`` and ``G``:

        - ``*_tv`` keys: ideal gas at the same ``T`` and molar volume;
        - ``*_tp`` keys: ideal gas at the same ``T`` and pressure (the usual
          tabulated departure function). They differ by ``ln Z``.

        ``U`` and ``H`` are the same for both references (an ideal gas's
        energy depends on ``T`` only).

        Keys: ``"a_res"`` (``A^res/RT``), ``"z"``, ``"u_res"``
        (``-T da/dT``), ``"h_res"`` (``u_res + Z - 1``), ``"s_res_tv"``
        (``-a - T da/dT``), ``"s_res_tp"`` (``s_res_tv + ln Z``),
        ``"g_res_tv"`` (``a + Z - 1``), ``"g_res_tp"`` (``g_res_tv - ln Z``,
        equal to ``sum_i x_i ln phi_i``).

        These are **residual** properties only. Total enthalpy, entropy or
        heat capacity also need ideal-gas heat capacities, which the databank
        does not carry, and ``Cp^res`` needs a second temperature derivative;
        neither is provided.

        Raises:
            ModelError: where ``Z <= 0`` (inside the spinodal), since the
                ``*_tp`` values need ``ln Z``.
        """
        density = _validated_density(density_mol_m3)
        state = self._state(temperature_K, density, composition)
        z = 1.0 + state.z_minus_one
        if z <= 0.0:
            raise ModelError(
                f"Z = {z!r} <= 0 at this density (inside the spinodal); the constant-pressure "
                "residual properties need ln Z and do not exist here."
            )
        t_da_dt = float(temperature_K) * self.residual_helmholtz_temperature_derivative(
            temperature_K=temperature_K, volume_m3=1.0 / density, composition=composition
        )
        a = state.a_res
        ln_z = math.log(z)
        s_tv = -a - t_da_dt
        g_tv = a + state.z_minus_one
        return {
            "a_res": a,
            "z": z,
            "u_res": -t_da_dt,
            "h_res": -t_da_dt + state.z_minus_one,
            "s_res_tv": s_tv,
            "s_res_tp": s_tv + ln_z,
            "g_res_tv": g_tv,
            "g_res_tp": g_tv - ln_z,
        }

    def residual_helmholtz_terms(
        self,
        *,
        temperature_K: float,
        volume_m3: float,
        composition: Sequence[float],
    ) -> dict[str, float]:
        """Return ``A^res/(R T)`` split into the terms that make it up (ADR-0018).

        Keys: ``"hard-chain"`` (Eq. A.4, the Boublik-Mansoori hard-sphere
        reference and the chain correction together, which is how the 2001
        paper groups them), ``"dispersion"`` (Eq. A.10), ``"association"``
        (Gross & Sadowski 2002; **0.0** for a non-associating mixture, where
        the term is not evaluated at all) and ``"total"``, which is exactly
        what :meth:`residual_helmholtz` returns.

        This exists so the association contribution is *visible*: it is the
        term this package validates against an external reference term by
        term, and without an accessor a caller could only see the sum.

        Args:
            temperature_K: Temperature in K.
            volume_m3: **Molar** volume in m^3/mol.
            composition: Mole fractions, summing to 1 within
                ``COMPOSITION_SUM_TOL``.
        """
        volume = float(volume_m3)
        if not math.isfinite(volume) or volume <= 0.0:
            raise InputRangeError(f"Molar volume must be positive and finite (got {volume_m3!r}).")
        state = self._state(temperature_K, 1.0 / volume, composition)
        return {
            "hard-chain": state.a_hc,
            "dispersion": state.a_disp,
            "association": state.a_assoc,
            "total": state.a_res,
        }

    def site_fractions(
        self,
        *,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> list[float]:
        """Return the non-bonded association site fractions ``X`` at this state.

        Two entries per associating component - its type-A site then its
        type-B site - in the instance's component order, and an **empty list**
        for a non-associating mixture. ``X`` is the fraction of sites of that
        kind that are *not* hydrogen bonded: 1 in the ideal-gas limit and a few
        per cent in liquid water. It is the inner unknown of the association
        term (see :mod:`chemthermo.eos._pcsaft_association`), and it is solved
        to a mass-action residual of ``1e-14``.
        """
        state = self._state(temperature_K, density_mol_m3, composition)
        if state.site_fractions is None:
            return []
        return [float(value) for value in state.site_fractions]

    def compressibility_factor(
        self,
        *,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> float:
        """Return ``Z = 1 + rho (d a_res / d rho)_{T,x}`` at the given state."""
        state = self._state(temperature_K, density_mol_m3, composition)
        return 1.0 + state.z_minus_one

    def pressure_Pa(
        self,
        *,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> float:
        """Return ``P = Z rho R T`` in Pa at the given state.

        Negative pressures are returned as computed: inside the spinodal the
        model genuinely has them, and hiding that would hide the root
        structure the next slice has to solve.
        """
        density = _validated_density(density_mol_m3)
        temperature = validate_temperature(temperature_K)
        state = self._state(temperature_K, density_mol_m3, composition)
        return (1.0 + state.z_minus_one) * density * R_J_PER_MOL_K * temperature

    def ln_fugacity_coefficients(
        self,
        *,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> list[float]:
        """Return ``ln phi_i`` (natural logs), one per component, at this state.

        The state is fixed by ``(T, rho, x)``; the pressure implied by it is
        :meth:`pressure_Pa`. Raises ``ModelError`` when ``Z <= 0``, where
        ``ln phi`` does not exist.
        """
        return self._ln_fugacity_coefficients(
            names=self.components,
            temperature_K=temperature_K,
            density_mol_m3=density_mol_m3,
            composition=composition,
        )

    def _ln_fugacity_coefficients(
        self,
        *,
        names: tuple[str, ...],
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> list[float]:
        state = self._state(temperature_K, density_mol_m3, composition, names)
        z_factor = 1.0 + state.z_minus_one
        if z_factor <= 0.0:
            raise ModelError(
                f"PC-SAFT compressibility factor Z = {z_factor!r} is not positive at "
                f"T = {temperature_K!r} K, rho = {density_mol_m3!r} mol/m^3; fugacity "
                "coefficients do not exist inside the mechanically unstable region."
            )
        x = self._composition(composition, names)
        mu_res = state.a_res + state.z_minus_one + state.da_dx - float(x @ state.da_dx)
        return (mu_res - math.log(z_factor)).tolist()

    # -- the (T, P, x) interface: EquationOfState (ADR-0015) ----------------

    def fugacity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> list[float]:
        """Return fugacity coefficients on the named density root.

        This is the ``chemthermo.models.EquationOfState`` interface, with the
        Peng-Robinson label semantics: ``"vapor"`` selects the largest ``Z``
        (lowest density) admissible root and ``"liquid"`` the smallest ``Z``
        (highest density) one. When the solver finds only one root the two
        labels return the same values, exactly as a cubic does outside its
        three-root region - the tangent-plane evaluator's minimum-Gibbs
        selection then has nothing to choose between and keeps the first
        candidate (ADR-0005).

        Args:
            mixture: Mixture providing the component names (and, for the
                stability test's Wilson estimates, their critical properties).
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, summing to 1 within
                ``COMPOSITION_SUM_TOL``.
            phase: ``"vapor"`` or ``"liquid"``.

        Returns:
            Fugacity coefficients (dimensionless), one per component.

        Raises:
            CompositionError: If the composition length does not match.
            ModelError: If the instance's components disagree with the
                mixture's, if no admissible density root exists at this state,
                or if ``ln phi`` does not exist on the selected root.
            ValueError: If ``phase`` is neither ``"vapor"`` nor ``"liquid"``.
        """
        if phase not in (_VAPOR, _LIQUID):
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        names = self._resolve_components(mixture)
        density = self._root_for_phase(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
            phase=phase,
        )
        ln_phi = self._ln_fugacity_coefficients(
            names=names,
            temperature_K=temperature_K,
            density_mol_m3=density,
            composition=composition,
        )
        return np.exp(np.asarray(ln_phi, dtype=float)).tolist()

    def log_fugacity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> list[float]:
        """Return ``ln phi_i`` on the named density root (ADR-0022's numerical guard).

        Exactly :meth:`fugacity_coefficients` without the final ``exp``: same
        root, same ``ln phi``, one fewer lossy step. It exists because a long
        chain makes ``exp(ln phi)`` underflow - ``ln phi_PE = -1690.6`` for a
        53 000 g/mol polyethylene in n-pentane at 453 K and 10 MPa, where
        ``exp`` is ``0.0`` and every consumer refuses a non-positive fugacity
        coefficient. See ``chemthermo.models.EquationOfState`` for the contract
        and for why callers must reach for it only when the ``phi`` they
        already have is unusable.
        """
        if phase not in (_VAPOR, _LIQUID):
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        names = self._resolve_components(mixture)
        density = self._root_for_phase(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
            phase=phase,
        )
        return self._ln_fugacity_coefficients(
            names=names,
            temperature_K=temperature_K,
            density_mol_m3=density,
            composition=composition,
        )

    def ln_fugacity_branches(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
    ) -> dict[str, list[float]]:
        """Return ``ln phi_i`` on both density roots from **one** root solve (ADR-0023).

        The ``EquationOfState`` capability, and for PC-SAFT it is the one that
        pays: :meth:`fugacity_coefficients` spends essentially all of its time
        in :func:`chemthermo.eos._pcsaft_density.solve_density_roots` - a
        1599-point isotherm scan plus a safeguarded Newton per bracket - and a
        caller that needs both branches to compare them used to pay for that
        twice at the same ``(T, P, x)``, for the same roots.

        The returned values are ``ln phi`` *before* the exponential
        :meth:`fugacity_coefficients` finishes with, which is what makes the
        reuse bit-identical rather than merely close: ``exp`` of what comes
        back here is the same double that method returns, on the same root.

        When the state has a single admissible root both labels carry that
        root's values - the documented behaviour of ``phase="vapor"`` and
        ``phase="liquid"`` there - and the evaluation is done once, not twice.

        Raises:
            ModelError: If no admissible density root exists at this state, or
                if ``ln phi`` does not exist on a root that was found. The
                caller then falls back to the per-branch route, which raises
                the same refusal with the same message.
        """
        names = self._resolve_components(mixture)
        roots = self._density_roots(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
        ).densities
        vapor = self._ln_fugacity_coefficients(
            names=names,
            temperature_K=temperature_K,
            density_mol_m3=roots[0],
            composition=composition,
        )
        if len(roots) == 1:
            return {_VAPOR: vapor, _LIQUID: list(vapor)}
        liquid = self._ln_fugacity_coefficients(
            names=names,
            temperature_K=temperature_K,
            density_mol_m3=roots[-1],
            composition=composition,
        )
        return {_VAPOR: vapor, _LIQUID: liquid}

    def density_roots(
        self,
        *,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        mixture: Mixture | None = None,
    ) -> tuple[float, ...]:
        """Return the mechanically stable molar densities at ``(T, P, x)``.

        The roots of ``P_model(T, rho, x) = P`` that satisfy
        ``(dP/drho)_{T,x} > 0``, in mol/m^3, sorted ascending: the first entry
        is the vapour-like root and the last the liquid-like one. A
        single-element result means the fluid has one phase-like state at this
        ``(T, P, x)`` - a dense liquid, a dilute gas or a supercritical fluid -
        and is not an error. The mechanically unstable (spinodal-branch) root
        is never returned.

        Args:
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, summing to 1 within
                ``COMPOSITION_SUM_TOL``.
            mixture: Supplies the component names when the instance was built
                without ``components``; otherwise optional, and checked
                against them.

        Returns:
            Admissible molar densities in mol/m^3, ascending.

        Raises:
            ModelError: If the scan finds no admissible root at this state.
                The message names the state and the scan grid.
        """
        names = self._resolve_components(mixture)
        return self._density_roots(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
        ).densities

    def molar_volume(
        self,
        *,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
        mixture: Mixture | None = None,
    ) -> float:
        """Return the molar volume in m^3/mol on the named density root.

        The reciprocal of the root :meth:`fugacity_coefficients` would use for
        the same ``phase``; see :meth:`density_roots` for the root semantics.
        """
        if phase not in (_VAPOR, _LIQUID):
            raise ValueError("phase must be 'vapor' or 'liquid'.")
        names = self._resolve_components(mixture)
        density = self._root_for_phase(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
            phase=phase,
        )
        return 1.0 / density

    def phase_identity(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> str:
        """Return "liquid" or "vapor" from the compressibility criterion (ADR-0017).

        ``kappa = P / (rho * dP/drho)`` is the dimensionless isothermal
        compressibility times pressure, evaluated at the density root
        ``phase`` selects (see :meth:`density_roots`). ``dP/drho`` is the
        analytic derivative already computed by
        :meth:`chemthermo.eos._pcsaft_density.PCSAFTIsotherm.pressure_and_slope`
        for the Newton refinement and the mechanical-stability filter of the
        density-root solver, so this needs no new derivative and no finite
        difference. See ``chemthermo.models.base.KAPPA_LIQUID_THRESHOLD`` for
        the threshold and ADR-0017 for the measured separation between liquid
        and vapor roots.

        When the isotherm has a single admissible root at this composition,
        both ``phase="liquid"`` and ``phase="vapor"`` name that same root and
        this method returns the same identity either way - which is exactly
        the case the historical vapor-first Gibbs tie-break (ADR-0008
        decision 3, superseded by ADR-0017) could not tell apart.
        """
        if phase not in (_VAPOR, _LIQUID):
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        names = self._resolve_components(mixture)
        density = self._root_for_phase(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
            phase=phase,
        )
        isotherm = self._isotherm(names=names, temperature_K=temperature_K, composition=composition)
        eta = density / isotherm.density_per_eta
        model_pressure, slope = isotherm.pressure_and_slope(eta)
        if not math.isfinite(slope) or slope <= 0.0:
            raise ModelError(
                f"PC-SAFT root at rho = {density!r} mol/m^3 is mechanically unstable "
                "(dP/drho <= 0); phase identity is undefined there."
            )
        kappa = model_pressure / (density * slope)
        return _LIQUID if kappa < KAPPA_LIQUID_THRESHOLD else _VAPOR

    # -- parameter and component plumbing -----------------------------------

    def component_parameters(
        self, components: Sequence[str] | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(m, sigma_A, epsilon_k_K)`` for this instance's components.

        Args:
            components: Optional explicit name order; defaults to the
                instance's own ``components``.
        """
        names = self.components if components is None else tuple(components)
        if not names:
            raise ModelError(
                "PC-SAFT requires at least one component: this instance was built without "
                "'components', so the component order has to come from a Mixture."
            )
        return self._source().for_components(names)

    def association_parameters(
        self, components: Sequence[str] | None = None
    ) -> tuple[PCSAFTAssociationRecord | None, ...]:
        """Return one association record (or ``None``) per component (ADR-0018).

        An all-``None`` result means this instance is non-associating and the
        association term is not evaluated at all, so its numbers are those of
        ADR-0014 down to the last bit.

        Args:
            components: Optional explicit name order; defaults to the
                instance's own ``components``.
        """
        names = self.components if components is None else tuple(components)
        if not names:
            raise ModelError(
                "PC-SAFT requires at least one component: this instance was built without "
                "'components', so the component order has to come from a Mixture."
            )
        return self._source().association_for_components(names)

    def associates(self, components: Sequence[str] | None = None) -> bool:
        """Return whether any of the named components carries association sites."""
        return any(record is not None for record in self.association_parameters(components))

    def _source(self) -> PCSAFTParameters:
        return self.parameters if self.parameters is not None else get_pcsaft_parameters()

    def kij_matrix(self, components: Sequence[str] | None = None) -> np.ndarray:
        """Return the dense ``n x n`` kij matrix in the given (or instance) order."""
        names = self.components if components is None else tuple(components)
        if not names:
            raise ModelError(
                "PC-SAFT requires at least one component: this instance was built without "
                "'components', so the component order has to come from a Mixture."
            )
        return kij_matrix(self.kij, names)

    def _resolve_components(self, mixture: Mixture | None) -> tuple[str, ...]:
        """Return the component-name order to use for a ``(T, P, x)`` request.

        An instance built without ``components`` takes the mixture's order.
        An instance built *with* them must agree with the mixture, compared
        after ``normalize_name``: reordering the caller's composition silently
        would be a wrong answer rather than a convenience.
        """
        if not self.components:
            if mixture is None:
                raise ModelError(
                    "PCSAFTEOS was built without 'components' and no Mixture was supplied, "
                    "so the component order is unknown."
                )
            return tuple(mixture.component_names)

        if mixture is not None:
            expected = [normalize_name(name) for name in self.components]
            actual = [normalize_name(name) for name in mixture.component_names]
            if expected != actual:
                raise ModelError(
                    "PCSAFTEOS components do not match the mixture: the model was built for "
                    f"{tuple(self.components)!r} but the mixture is {tuple(mixture.component_names)!r}. "
                    "Build PCSAFTEOS() without 'components' to take the mixture's order."
                )
        return self.components

    def _isotherm(
        self,
        *,
        names: tuple[str, ...],
        temperature_K: float,
        composition: Sequence[float],
    ) -> PCSAFTIsotherm:
        temperature = validate_temperature(temperature_K)
        x = self._composition(composition, names)
        m, sigma_A, epsilon_k_K = self.component_parameters(names)
        return build_isotherm(
            temperature_K=temperature,
            composition=x,
            m=m,
            sigma_A=sigma_A,
            epsilon_k_K=epsilon_k_K,
            kij=self.kij_matrix(names),
            association=self.association_parameters(names),
        )

    def _density_roots(
        self,
        *,
        names: tuple[str, ...],
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
    ) -> DensityRoots:
        """The admissible density roots at ``(T, P, x)``, through the call memo (ADR-0030).

        Every ``(T, P, x)`` entry point of this class funnels through here -
        :meth:`fugacity_coefficients`, :meth:`log_fugacity_coefficients`,
        :meth:`ln_fugacity_branches`, :meth:`density_roots`,
        :meth:`molar_volume` and :meth:`phase_identity` - and it is where
        essentially all of a PC-SAFT evaluation's time goes: a 1599-point
        isotherm scan plus a safeguarded Newton per bracket. Inside a
        ``flash_tp`` / ``stability_tp`` call a memo is active
        (:mod:`chemthermo._eos_memo`) and a repeat of an argument list solved
        earlier in *that call* returns the identical :class:`DensityRoots`
        object rather than an equal one, so nothing downstream can see a
        different double. Outside such a call there is no memo and this runs
        exactly as it did before.

        The key carries ``id(self)`` because two ``PCSAFTEOS`` instances may
        differ in parameters or ``kij`` while agreeing on every argument here;
        the memo keeps the instance alive for as long as the entry lives, so
        that id cannot be reused underneath it.
        """
        memo = active_memo()
        key: tuple[object, ...] | None = None
        if memo is not None:
            key = (
                "pcsaft-density-roots",
                id(self),
                names,
                temperature_K,
                pressure_Pa,
                composition_key(composition),
            )
            cached = memo.lookup(key)
            if cached is not MISS:
                return cached

        pressure = validate_pressure(pressure_Pa)
        isotherm = self._isotherm(names=names, temperature_K=temperature_K, composition=composition)
        description = (
            f" at T = {float(temperature_K)!r} K, x = {list(map(float, composition))!r} "
            f"for {names!r}"
        )
        roots = solve_density_roots(isotherm, pressure, state_description=description)
        if memo is not None and key is not None:
            memo.store(key, roots, self)
        return roots

    def _root_for_phase(
        self,
        *,
        names: tuple[str, ...],
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> float:
        roots = self._density_roots(
            names=names,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=composition,
        ).densities
        return roots[0] if phase == _VAPOR else roots[-1]

    def _composition(
        self, composition: Sequence[float], names: Sequence[str] | None = None
    ) -> np.ndarray:
        expected = len(self.components) if names is None else len(names)
        if len(composition) != expected:
            raise CompositionError(
                "Composition length must match number of PC-SAFT components "
                f"({len(composition)} != {expected})."
            )
        fractions = validate_fractions(composition, normalize=False, tol=COMPOSITION_SUM_TOL)
        return np.array(fractions, dtype=float)

    def _state(
        self,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
        names: Sequence[str] | None = None,
    ) -> _PCSAFTState:
        temperature = validate_temperature(temperature_K)
        density = _validated_density(density_mol_m3)
        # Resolved first: an instance built without ``components`` has no
        # component order at all, and "length 1 != 0" would be a poor way to
        # say so.
        m, sigma_A, epsilon_k_K = self.component_parameters(names)
        x = self._composition(composition, names)
        return _evaluate(
            temperature_K=temperature,
            density_mol_m3=density,
            x=x,
            m=m,
            sigma_A=sigma_A,
            epsilon_k_K=epsilon_k_K,
            kij=self.kij_matrix(names),
            association=self.association_parameters(names),
        )


def _validated_density(density_mol_m3: float) -> float:
    density = float(density_mol_m3)
    if not math.isfinite(density) or density <= 0.0:
        raise InputRangeError(
            f"Molar density must be positive and finite (got {density_mol_m3!r})."
        )
    return density


register_eos("pcsaft", PCSAFTEOS)

__all__ = [
    "AVOGADRO_PER_MOL",
    "A_UNIVERSAL",
    "BOLTZMANN_J_PER_K",
    "B_UNIVERSAL",
    "PCSAFTAssociationRecord",
    "PCSAFTEOS",
    "PCSAFTParameterError",
    "PCSAFTParameters",
    "R_J_PER_MOL_K",
]
