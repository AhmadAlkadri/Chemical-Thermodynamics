"""PC-SAFT equation of state (Gross & Sadowski 2001), non-associating.

Source
------
J. Gross and G. Sadowski, "Perturbed-Chain SAFT: An Equation of State Based on
a Perturbation Theory for Chain Molecules", Ind. Eng. Chem. Res. 40 (2001)
1244-1260 (DOI 10.1021/ie0003887). Equation numbers below are that paper's
appendix numbering. Only the hard-chain and dispersion contributions are
implemented here; **association** (Gross & Sadowski, Ind. Eng. Chem. Res. 41
(2002) 5510) and the **polar** terms are deliberately out of scope for this
slice, so this module must not be used for associating or strongly polar
fluids.

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

    a_res = A^res / (N k T) = A^res / (R T) per mole = a_hc + a_disp

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

Scope of this slice
-------------------
No density root solving happens here: every property method takes the density
(or molar volume) as an input. Finding the vapour/liquid density roots at a
given ``(T, P)``, and wiring PC-SAFT into ``stability_tp`` / ``flash_tp``, is
the next slice (``pcsaft-density-roots-flash``). See ADR-0014.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple, Sequence

import numpy as np

from ..exceptions import CompositionError, InputRangeError, ModelError
from ..models._kij import KijInput, KijPairs, canonicalize_kij, kij_matrix
from ..parameters.pcsaft import (
    PCSAFTParameterError,
    PCSAFTParameters,
    get_pcsaft_parameters,
)
from ..validation import COMPOSITION_SUM_TOL, validate_fractions, validate_temperature
from .api import EOSProtocol
from .registry import register_eos

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
    fixed temperature and total density. ``z_hc`` and ``z_disp`` are the two
    contributions to ``z_minus_one``; they are kept apart so each term's
    density derivative can be finite-difference tested on its own.
    """

    a_res: float
    z_minus_one: float
    da_dx: np.ndarray
    eta: float
    a_hc: float
    a_disp: float
    z_hc: float
    z_disp: float


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
) -> _PCSAFTState:
    """Evaluate the non-associating PC-SAFT model at one state.

    Pure function of arrays; see the module docstring for every equation and
    for the two chain rules that give the derivatives.
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


@dataclass(frozen=True)
class PCSAFTEOS(EOSProtocol):
    """PC-SAFT equation of state for non-associating fluids (ADR-0014).

    Args:
        components: Component names, in the order every composition argument
            uses. Names are matched against the parameter set (and therefore
            the packaged databank) via ``chemthermo.data.normalize_name``.
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

    Every property method takes the state as ``(T, molar density)`` or
    ``(T, molar volume)``: **this slice does no density root solving**, so the
    caller chooses which root they are on. A state inside the mechanically
    unstable region can have ``Z <= 0``, for which fugacity coefficients do not
    exist; :meth:`ln_fugacity_coefficients` raises ``ModelError`` there rather
    than returning a ``nan``.
    """

    components: tuple[str, ...]
    parameters: PCSAFTParameters | None = None
    kij: KijInput | KijPairs = 0.0
    name: str = "PC-SAFT"

    def __post_init__(self) -> None:
        if not self.components:
            raise ModelError("PC-SAFT requires at least one component.")
        object.__setattr__(self, "components", tuple(self.components))
        object.__setattr__(self, "kij", canonicalize_kij(self.kij, model="PCSAFTEOS"))

    def num_components(self) -> int:
        """Return the number of components the EOS instance was configured for."""
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
        state = self._state(temperature_K, density_mol_m3, composition)
        z_factor = 1.0 + state.z_minus_one
        if z_factor <= 0.0:
            raise ModelError(
                f"PC-SAFT compressibility factor Z = {z_factor!r} is not positive at "
                f"T = {temperature_K!r} K, rho = {density_mol_m3!r} mol/m^3; fugacity "
                "coefficients do not exist inside the mechanically unstable region."
            )
        x = self._composition(composition)
        mu_res = state.a_res + state.z_minus_one + state.da_dx - float(x @ state.da_dx)
        return (mu_res - math.log(z_factor)).tolist()

    def component_parameters(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(m, sigma_A, epsilon_k_K)`` for this instance's components."""
        source = self.parameters if self.parameters is not None else get_pcsaft_parameters()
        return source.for_components(self.components)

    def kij_matrix(self) -> np.ndarray:
        """Return the dense ``n x n`` kij matrix in this instance's order."""
        return kij_matrix(self.kij, self.components)

    def _composition(self, composition: Sequence[float]) -> np.ndarray:
        if len(composition) != len(self.components):
            raise CompositionError(
                "Composition length must match number of PC-SAFT components "
                f"({len(composition)} != {len(self.components)})."
            )
        fractions = validate_fractions(composition, normalize=False, tol=COMPOSITION_SUM_TOL)
        return np.array(fractions, dtype=float)

    def _state(
        self,
        temperature_K: float,
        density_mol_m3: float,
        composition: Sequence[float],
    ) -> _PCSAFTState:
        temperature = validate_temperature(temperature_K)
        density = _validated_density(density_mol_m3)
        x = self._composition(composition)
        m, sigma_A, epsilon_k_K = self.component_parameters()
        return _evaluate(
            temperature_K=temperature,
            density_mol_m3=density,
            x=x,
            m=m,
            sigma_A=sigma_A,
            epsilon_k_K=epsilon_k_K,
            kij=self.kij_matrix(),
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
    "PCSAFTEOS",
    "PCSAFTParameterError",
    "PCSAFTParameters",
    "R_J_PER_MOL_K",
]
