"""Internal tangent-plane evaluator contract for :mod:`chemthermo.stability`.

This module is private (ADR-0007). Nothing here is exported from
``chemthermo`` or from ``chemthermo.stability``.

Why it exists
-------------
Michelsen's tangent-plane machinery -- the successive-substitution map, the
second-order stage, trivial-solution detection, the trial summary and the
result types -- is identical for an equation of state and for an
activity-coefficient model. Only two things differ:

1. what the "fugacity term" of component ``i`` at a trial composition ``w`` is
   (``ln phi_i(w)`` on the minimum-Gibbs compressibility root for an EOS,
   ``ln gamma_i(w)`` for an activity model), and
2. which deterministic initial estimates make sense (Wilson K-value estimates
   are a vapor-liquid device and are meaningless for an activity-only model).

Those two differences are the whole contract:

    ln_fugacity_terms(w) -> (ndarray, str | None)
    initial_estimates(z, active) -> list[(label, w0)]

The solver in :mod:`chemthermo.stability.tp` sees nothing else, so it does not
know which model family it is serving.
"""

from __future__ import annotations

import math
from typing import Protocol, Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import ModelError
from ..flash._common import as_float_array, wilson_k
from ..models import ActivityModel, EquationOfState

# Trace amount kept on the non-dominant components of a pure-component-dominant
# initial estimate. A hard zero would pin those components at W_i = 0 forever.
_PURE_TRIAL_TRACE = 1e-3


class _TangentPlaneEvaluator(Protocol):
    """What the tangent-plane solver needs from a thermodynamic model.

    Attributes:
        model_family: ``"eos"`` or ``"activity"``; recorded in diagnostics.
        pressure_dependent: True when the returned terms depend on pressure.
            False for an activity-only model, where ``pressure_Pa`` is still
            validated for API uniformity but does not affect the result.
    """

    model_family: str
    pressure_dependent: bool

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        """Return the tangent-plane fugacity terms and an optional branch label.

        Args:
            composition: Normalized mole fractions ``w``.

        Returns:
            ``(terms, branch)`` where ``terms[i]`` is ``ln phi_i(w)`` (EOS, on
            the minimum-Gibbs root) or ``ln gamma_i(w)`` (activity model), and
            ``branch`` names the selected compressibility root for an EOS or is
            None when the model has a single branch.

        Raises:
            ModelError: If the model returns unusable values.
        """
        ...

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[tuple[str, np.ndarray]]:
        """Return deterministic ``(label, w0)`` trial-phase initial estimates."""
        ...


class _EOSTangentPlane:
    """Evaluator backed by an equation of state.

    ``ln phi`` is taken from the minimum-Gibbs compressibility root, selected
    generically over the existing ``EquationOfState`` interface (ADR-0005).
    """

    model_family = "eos"
    pressure_dependent = True

    def __init__(
        self,
        eos: EquationOfState,
        *,
        mixture: Mixture,
        temperature: float,
        pressure: float,
    ) -> None:
        self._eos = eos
        self._mixture = mixture
        self._temperature = temperature
        self._pressure = pressure

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        return _ln_phi_min_gibbs(
            self._eos,
            mixture=self._mixture,
            temperature=self._temperature,
            pressure=self._pressure,
            composition=composition,
        )

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[tuple[str, np.ndarray]]:
        """Two Wilson estimates plus one pure-component-dominant estimate each."""
        estimates: list[tuple[str, np.ndarray]] = []

        k = wilson_k(self._mixture, self._temperature, self._pressure)
        estimates.append(("wilson-vapor", _normalized(k * z, active)))
        estimates.append(("wilson-liquid", _normalized(z / k, active)))
        estimates.extend(_pure_component_estimates(self._mixture, z, active))
        return estimates


class _ActivityTangentPlane:
    """Evaluator backed by an activity-coefficient model (liquid-liquid).

    Both phases are liquids with the same pure-liquid reference state, so the
    reference fugacities cancel from the tangent-plane distance and
    ``ln gamma_i`` takes the place of ``ln phi_i`` exactly (see the module
    docstring of :mod:`chemthermo.stability.tp`). There is a single branch, so
    no root selection is needed and the branch label is None.
    """

    model_family = "activity"
    pressure_dependent = False

    def __init__(
        self,
        activity_model: ActivityModel,
        *,
        mixture: Mixture,
        temperature: float,
    ) -> None:
        self._model = activity_model
        self._mixture = mixture
        self._temperature = temperature

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        try:
            values = as_float_array(
                self._model.activity_coefficients(
                    mixture=self._mixture,
                    temperature_K=self._temperature,
                    composition=composition.tolist(),
                )
            )
        except ModelError:
            raise
        except Exception as exc:  # noqa: BLE001 - surfaced as a ModelError below
            raise ModelError(f"Activity model failed during stability analysis: {exc}") from exc

        if values.shape != composition.shape:
            raise ModelError("Activity model returned an inconsistent number of coefficients.")
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ModelError("Non-finite or non-positive activity coefficients.")
        return np.log(values), None

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[tuple[str, np.ndarray]]:
        """Pure-component-dominant estimates only.

        Wilson K-values are a vapor-liquid construction built from Tc, Pc and
        omega; for a liquid-liquid test driven by an activity model they carry
        no information about the split and are therefore not used. Michelsen
        (1982) recommends pure-component-dominant estimates for liquid-liquid
        stability, and one per component is what this evaluator supplies. A
        single-component feed has no composition degree of freedom, so the feed
        itself is the only admissible trial.
        """
        if int(np.count_nonzero(active)) <= 1:
            names = self._mixture.component_names
            index = int(np.argmax(active))
            return [(f"pure-{names[index]}", _normalized(z, active))]
        return _pure_component_estimates(self._mixture, z, active)


def _ln_phi_min_gibbs(
    eos: EquationOfState,
    *,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    composition: np.ndarray,
) -> tuple[np.ndarray, str]:
    """Return ``ln phi`` on the minimum-Gibbs root branch and the branch name.

    The branch minimizing ``sum_i w_i ln phi_i(w)`` minimizes the molar Gibbs
    energy at fixed ``(T, P, w)`` because that sum is the reduced residual Gibbs
    energy and the ideal-mixing contribution is root independent.
    """
    best_ln_phi: np.ndarray | None = None
    best_branch = ""
    best_g = math.inf
    failures: list[str] = []

    for branch in ("vapor", "liquid"):
        try:
            values = as_float_array(
                eos.fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=temperature,
                    pressure_Pa=pressure,
                    composition=composition.tolist(),
                    phase=branch,
                )
            )
        except Exception as exc:  # noqa: BLE001 - branch may be unavailable
            failures.append(f"{branch}: {exc}")
            continue

        if values.shape != composition.shape:
            failures.append(f"{branch}: inconsistent fugacity coefficient shape")
            continue
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            failures.append(f"{branch}: non-finite or non-positive fugacity coefficients")
            continue

        ln_phi = np.log(values)
        g_res = float(np.sum(composition * ln_phi))
        if not math.isfinite(g_res):
            failures.append(f"{branch}: non-finite reduced residual Gibbs energy")
            continue
        if g_res < best_g:
            best_g = g_res
            best_ln_phi = ln_phi
            best_branch = branch

    if best_ln_phi is None:
        raise ModelError(
            "No usable fugacity-coefficient branch for stability analysis ("
            + "; ".join(failures)
            + ")."
        )
    return best_ln_phi, best_branch


def _pure_component_estimates(
    mixture: Mixture, z: np.ndarray, active: np.ndarray
) -> list[tuple[str, np.ndarray]]:
    """One pure-component-dominant estimate per active component."""
    n_active = int(np.count_nonzero(active))
    if n_active <= 1:
        return []

    trace = _PURE_TRIAL_TRACE / (n_active - 1)
    names: Sequence[str] = mixture.component_names
    estimates: list[tuple[str, np.ndarray]] = []
    for index in range(z.size):
        if not active[index]:
            continue
        w = np.where(active, trace, 0.0)
        w[index] = 1.0 - _PURE_TRIAL_TRACE
        estimates.append((f"pure-{names[index]}", _normalized(w, active)))
    return estimates


def _normalized(values: np.ndarray, active: np.ndarray) -> np.ndarray:
    """Zero out inactive components and renormalize to sum to one."""
    w = np.where(active, np.asarray(values, dtype=float), 0.0)
    w = np.where(w > 0.0, w, 0.0)
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0.0:
        raise ModelError("Stability trial initial estimate has a non-positive total.")
    return w / total
