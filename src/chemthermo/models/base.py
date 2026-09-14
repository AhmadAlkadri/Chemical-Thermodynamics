"""Interfaces for thermodynamic models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Mapping, Sequence

from ..core import Mixture

#: Compressibility-based liquid/vapor threshold (ADR-0017). An implementation
#: of :meth:`EquationOfState.phase_identity` reports ``"liquid"`` when the
#: dimensionless isothermal-compressibility ratio
#: ``kappa = (1 / rho) * (d rho / dP)_T * P`` (equivalently
#: ``kappa = -P / (V * (dP/dV)_T)``) is below this value at the evaluated
#: root, and ``"vapor"`` otherwise. ``kappa`` is exactly 1 for an ideal gas
#: (``P = rho R T``) and well below 1 for a liquid; the value below was chosen
#: so that it separates every liquid root from every vapor root, with a wide
#: margin, over both the Peng-Robinson and the PC-SAFT validation grids - see
#: ``.agents/brain/adr/0017-phase-identity-by-compressibility.md`` for the
#: measured kappa tables.
KAPPA_LIQUID_THRESHOLD: float = 0.5


class EquationOfState(ABC):
    """Interface for equation-of-state models.

    Implementations should accept SI units (K, Pa) and return dimensionless
    fugacity coefficients for a specified phase.
    """

    name: str = "generic-eos"

    @abstractmethod
    def fugacity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> Sequence[float]:
        """Return fugacity coefficients for a phase composition.

        Args:
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: Phase label (e.g., "vapor" or "liquid").

        Returns:
            Fugacity coefficients (dimensionless), one per component.
        """

    def log_fugacity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> Sequence[float] | None:
        """Return ``ln phi_i`` on the same root :meth:`fugacity_coefficients` would use.

        **Optional, and a numerical guard rather than a second model** (ADR-0022).
        ``phi_i = exp(ln phi_i)`` is not representable as a double once
        ``|ln phi_i|`` passes about 709, and a chain molecule makes that
        ordinary: a polyethylene of ``Mw = 53 000`` (``m = 1393.9``) dissolved
        in n-pentane at 453 K and 10 MPa has ``ln phi_PE = -1690.6``, so
        ``fugacity_coefficients`` returns an exact ``0.0`` for it and every
        caller correctly refuses a non-positive fugacity coefficient. The model
        is perfectly well defined there; only the exponential is not.

        Callers must use this **only where the ``phi`` they already computed is
        unusable**, and must otherwise keep taking ``log`` of that ``phi``, so
        that no converging number moves. The default implementation returns
        ``None`` ("this model cannot say"), which makes the guard inert and the
        pre-ADR-0022 failure the outcome, exactly as before.

        Args:
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: Phase label, with the same semantics
                :meth:`fugacity_coefficients` uses.

        Returns:
            ``ln phi_i``, one per component, or ``None`` when this
            implementation cannot produce them without exponentiating.
        """
        return None

    def ln_fugacity_branches(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
    ) -> Mapping[str, Sequence[float]] | None:
        """Return ``ln phi_i`` on **every** density root, from one root solve (ADR-0023).

        **Optional, and a performance capability rather than a second model.**
        The tangent-plane evaluator's minimum-Gibbs rule (ADR-0005) and the
        split's lowest-Gibbs fallback (ADR-0019) both have to compare the
        branches against each other, so they ask the model for
        ``phase="vapor"`` and then for ``phase="liquid"`` at the same
        ``(T, P, x)``. Each of those calls solves for the roots from scratch -
        a cubic for Peng-Robinson, a 1599-point isotherm scan plus a
        safeguarded Newton for PC-SAFT - and the second solve reproduces the
        first exactly. An implementation of this method does the root solve
        once and evaluates every root it found.

        **The values are the model's own ``ln phi``**: the numbers
        :meth:`log_fugacity_coefficients` would return on the same root, which
        is to say ``fugacity_coefficients`` *without* the final ``exp``.
        Callers reconstruct what they used before by exponentiating - see
        :func:`chemthermo.flash._common.eos_branch_terms_all`, which is the
        only caller and which keeps the ADR-0022 guard unchanged. That is what
        makes the capability bit-identical by construction rather than by
        measurement: ``exp`` of the returned value **is** the ``phi`` the
        per-branch call returns, on the same root.

        The keys are the ``phase`` labels of :meth:`fugacity_coefficients`. A
        label the model cannot evaluate here is **omitted**, and a caller that
        finds a label missing falls back to the per-branch route for that
        composition, so the error a refusal produces is the per-branch one
        rather than a paraphrase. Where the model has a single admissible root
        both labels are present and carry that one root's values, exactly as
        the two per-branch calls would.

        The default implementation returns ``None`` ("this model has no such
        capability"), and every caller then keeps its per-branch behaviour.

        Args:
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.

        Returns:
            A mapping from phase label to ``ln phi_i``, or ``None`` when this
            implementation does not offer the capability.

        Raises:
            ModelError: If no root exists at this state at all. A caller
                treats that the way it treats any refusal from the model.
        """
        return None

    def phase_identity(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> str | None:
        """Return ``"liquid"`` or ``"vapor"`` from a compressibility criterion (ADR-0017).

        This identifies the *root* ``phase`` selects (the same root
        :meth:`fugacity_coefficients` would use) from the model's own
        pressure-volume behavior, rather than from which of the two labels
        happens to have lower Gibbs energy - the two agree whenever the model
        has two distinct roots, and only the compressibility criterion is
        meaningful when it has one (see
        ``chemthermo.models.base.KAPPA_LIQUID_THRESHOLD``).

        Args:
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: Phase label (``"vapor"`` or ``"liquid"``) selecting which
                root to identify, with the same semantics
                :meth:`fugacity_coefficients` uses.

        Returns:
            ``"liquid"`` or ``"vapor"``, or ``None`` when this implementation
            does not support compressibility-based identification. **The
            default implementation always returns ``None``**, which is the
            documented fallback: a caller that receives ``None`` must keep its
            pre-ADR-0017 behavior (the min-Gibbs tie-break for a single-phase
            result, the Wilson volatility ranking for a phi-phi split)
            unchanged.
        """
        return None


class ActivityModel(ABC):
    """Interface for activity coefficient models.

    Implementations should accept SI units (K) and return dimensionless
    activity coefficients for the liquid phase.
    """

    name: str = "generic-activity"

    @abstractmethod
    def activity_coefficients(
        self, *, mixture: Mixture, temperature_K: float, composition: Sequence[float]
    ) -> Sequence[float]:
        """Return activity coefficients for a phase composition.

        Args:
            mixture: Mixture providing component identities.
            temperature_K: Temperature in K.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.

        Returns:
            Activity coefficients (dimensionless), one per component.
        """
