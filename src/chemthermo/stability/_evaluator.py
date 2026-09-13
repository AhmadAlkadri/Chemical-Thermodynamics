"""Internal tangent-plane evaluator contract for :mod:`chemthermo.stability`.

This module is private (ADR-0007, ADR-0010). Nothing here is exported from
``chemthermo`` or from ``chemthermo.stability``.

Why it exists
-------------
Michelsen's tangent-plane machinery -- the successive-substitution map, the
second-order stage, trivial-solution detection, the trial summary and the
result types -- is identical for an equation of state, for an
activity-coefficient model, and for a set of *competing candidate phases*. Only
two things differ:

1. what the "fugacity term" of component ``i`` at a trial composition ``w`` is,
   and
2. which deterministic initial estimates make sense.

Those two differences are the whole contract:

    ln_fugacity_terms(w) -> (ndarray, str | None)
    ln_terms_on_surface(w, surface) -> (ndarray, str | None, bool)
    initial_estimates(z, active) -> list[_InitialEstimate]

The solver in :mod:`chemthermo.stability.tp` sees nothing else, so it does not
know which model family it is serving.

Phase candidates (ADR-0010)
---------------------------
An evaluator holds one or more **phase candidates**. A candidate is a
thermodynamic description of what the mixture could be at a composition ``w``,
and it contributes the term ``ln( f_i(w) / (x_i P) )`` of the shared reference:

===========================  ==========================  ======================
family                       candidates                  term of candidate
===========================  ==========================  ======================
``"eos"``                    the cubic's vapor and        ``ln phi_i(w)``
                             liquid compressibility
                             roots
``"activity"``               one activity-model liquid    ``ln gamma_i(w)``
``"modified-raoult"``        an activity-model liquid     ``ln gamma_i(w)
                             and an ideal gas             + ln(Psat_i/P)``
                                                          and ``0``
===========================  ==========================  ======================

``ln_fugacity_terms(w)`` returns the terms of the candidate with the **lowest
Gibbs energy** at ``w`` together with that candidate's label. At fixed
``(T, P, w)``

    G/RT = sum_i w_i [ g_i^0/RT + ln(w_i P / P^0) ] + sum_i w_i term_i(w)

and only the last sum depends on the candidate (the ideal-mixing part is
candidate independent), so minimizing ``sum_i w_i term_i(w)`` selects the
lowest-Gibbs candidate. That is one rule, not three: for a cubic EOS it is the
minimum-Gibbs root selection of ADR-0005, for a single activity-model liquid it
is a no-op, and for the modified-Raoult pair it decides whether ``w`` is a
liquid or a vapor.

The one-candidate case reports ``None`` as its label: there was no choice to
make, so there is nothing to report.

Trial surfaces (ADR-0012)
-------------------------
Re-selecting the lowest-Gibbs candidate *inside* a trial iteration is right for
the cubic roots and wrong for the heterogeneous modified-Raoult pair, so an
initial estimate may name the candidate its trial belongs to:

    initial_estimates(z, active) -> list[_InitialEstimate(label, w0, surface)]

``surface`` is None for the EOS and activity-only evaluators, which keeps their
iteration exactly the min-Gibbs one they have always used. When it is a
candidate label the solver calls :meth:`ln_terms_on_surface` at every iteration
instead, so the trial walks one fixed Gibbs surface. ADR-0012 gives the reason:
a missing cubic root is the *same* model failing to exist at that composition,
whereas the liquid and the ideal vapor are two different models whose surfaces
both exist everywhere, so swapping between them mid-iteration makes the
successive-substitution map discontinuous and non-monotone.

The tangent-plane distance *reported* for a trial is always the min-Gibbs one
at the converged composition: the true distance to the tangent plane is the
minimum over candidates, and a trial that iterated on the vapor surface must
not claim a vapor distance if the liquid lies lower there.

A future PC-SAFT model (several density roots) or a user-supplied Gibbs-energy
phase model would enter as further candidates behind the same two methods,
without a solver change.
"""

from __future__ import annotations

import math
from typing import Mapping, NamedTuple, Protocol, Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import ModelError
from ..flash._common import as_float_array, wilson_k
from ..models import ActivityModel, EquationOfState
from ..models._antoine import antoine_saturation_pressures, antoine_temperature_range

# Trace amount kept on the non-dominant components of a pure-component-dominant
# initial estimate. A hard zero would pin those components at W_i = 0 forever.
_PURE_TRIAL_TRACE = 1e-3

#: Candidate labels used by the modified-Raoult pair and by the cubic roots.
_LIQUID = "liquid"
_VAPOR = "vapor"


class _InitialEstimate(NamedTuple):
    """One deterministic trial-phase start.

    Attributes:
        label: Deterministic identifier of the estimate, reported as
            ``StabilityTrial.label``.
        composition: Normalized initial trial composition ``w0``.
        surface: Label of the phase candidate the trial is pinned to, or None
            to iterate on the lowest-Gibbs candidate re-selected at every
            iterate (the pre-ADR-0012 behavior, kept for the EOS and
            activity-only families).
    """

    label: str
    composition: np.ndarray
    surface: str | None = None


class _PhaseCandidate(Protocol):
    """One thermodynamic description the mixture could take at a composition.

    Attributes:
        label: Reported as the branch / candidate label of a stationary point.
        optional: True when the candidate may legitimately be absent at some
            compositions (a cubic root that does not exist there), so the
            selector should record the failure and carry on. False when a
            failure is a real error that must reach the caller: silently
            answering with the *other* candidate would be a wrong phase, not a
            missing one.
    """

    label: str
    optional: bool

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
        """Return ``ln( f_i(w) / (w_i P) )`` for this candidate at ``w``.

        Raises:
            ModelError: If the candidate is not evaluable at ``w`` (an absent
                compressibility root, a model failure, unusable values). The
                selector records the reason and tries the other candidates.
        """
        ...


class _TangentPlaneEvaluator(Protocol):
    """What the tangent-plane solver needs from a thermodynamic model.

    Attributes:
        model_family: ``"eos"``, ``"activity"`` or ``"modified-raoult"``;
            recorded in diagnostics.
        pressure_dependent: True when the returned terms depend on pressure.
            False for an activity-only model, where ``pressure_Pa`` is still
            validated for API uniformity but does not affect the result.
        diagnostics: Extra diagnostics keys contributed by the evaluator
            (empty for the EOS and activity-only families).
    """

    model_family: str
    pressure_dependent: bool
    diagnostics: Mapping[str, float | int | str | bool]

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        """Return the tangent-plane fugacity terms and the candidate label.

        Args:
            composition: Normalized mole fractions ``w``.

        Returns:
            ``(terms, label)`` where ``terms[i]`` is the term of the
            lowest-Gibbs candidate at ``w`` and ``label`` names that candidate.
            ``label`` is None when the evaluator holds a single candidate, so
            no selection was made.

        Raises:
            ModelError: If no candidate is usable at ``w``.
        """
        ...

    def ln_terms_on_surface(
        self, composition: np.ndarray, surface: str
    ) -> tuple[np.ndarray, str | None, bool]:
        """Return the terms of the *named* candidate at ``w`` (ADR-0012).

        Args:
            composition: Normalized mole fractions ``w``.
            surface: Label of the candidate to evaluate.

        Returns:
            ``(terms, label, fell_back)``. ``fell_back`` is True when the named
            candidate was optional and unavailable at ``w``, in which case the
            lowest-Gibbs candidate was used instead and ``label`` names it.

        Raises:
            ModelError: If the named candidate is unknown, or if it is
                mandatory and unusable at ``w``, or if no candidate is usable.
        """
        ...

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[_InitialEstimate]:
        """Return deterministic trial-phase initial estimates."""
        ...

    def identity_label(self, composition: np.ndarray, label: str | None) -> str | None:
        """Replace a min-Gibbs ``label`` with a compressibility identity, if possible.

        ``label`` is a candidate label already selected by
        :func:`_select_min_gibbs` / :func:`_select_surface` (a min-Gibbs tie-break
        when two candidates coincide, e.g. a cubic's single real root). This
        gives the evaluator a chance to replace it with a model-measured
        identity instead (ADR-0017); it must never change ``composition`` or
        any numeric term, only the label reported alongside it.

        The default the EOS family relies on is ``EquationOfState.phase_identity``
        returning ``None``; families with no such measurement (activity-only,
        modified-Raoult) return ``label`` unchanged.
        """
        ...


# ---------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------


class _CubicRootCandidate:
    """One compressibility branch of an equation of state (``ln phi_i``)."""

    optional = True

    def __init__(
        self,
        eos: EquationOfState,
        *,
        label: str,
        mixture: Mixture,
        temperature: float,
        pressure: float,
    ) -> None:
        self.label = label
        self._eos = eos
        self._mixture = mixture
        self._temperature = temperature
        self._pressure = pressure

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
        try:
            values = as_float_array(
                self._eos.fugacity_coefficients(
                    mixture=self._mixture,
                    temperature_K=self._temperature,
                    pressure_Pa=self._pressure,
                    composition=composition.tolist(),
                    phase=self.label,
                )
            )
        except Exception as exc:  # noqa: BLE001 - branch may be unavailable
            raise ModelError(str(exc)) from exc

        if values.shape != composition.shape:
            raise ModelError("inconsistent fugacity coefficient shape")
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ModelError("non-finite or non-positive fugacity coefficients")
        return np.log(values)


class _ActivityLiquidCandidate:
    """An activity-coefficient liquid (``ln gamma_i``, plus an optional offset).

    The offset is the reduced pure-liquid reference fugacity
    ``ln(f_i^0 / P)``. It is None for a liquid-liquid problem, where both
    phases share the reference and it cancels, and
    ``ln(Psat_i(T) / P)`` for the modified-Raoult pair, where the liquid has to
    be compared against a vapor and the reference cannot cancel.
    """

    optional = False

    def __init__(
        self,
        activity_model: ActivityModel,
        *,
        label: str,
        mixture: Mixture,
        temperature: float,
        reference_offset: np.ndarray | None = None,
    ) -> None:
        self.label = label
        self._model = activity_model
        self._mixture = mixture
        self._temperature = temperature
        self._offset = reference_offset

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
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
        if self._offset is None:
            return np.log(values)
        return np.log(values) + self._offset


class _IdealVaporCandidate:
    """An ideal gas: ``phi_i = 1`` and the reference is ``P`` itself, so the term is 0."""

    optional = False

    def __init__(self, *, label: str) -> None:
        self.label = label

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
        return np.zeros_like(composition)


def _select_min_gibbs(
    candidates: Sequence[_PhaseCandidate],
    composition: np.ndarray,
    *,
    failure_message: str,
) -> tuple[np.ndarray, str]:
    """Return the terms and label of the lowest-Gibbs candidate at ``composition``.

    The candidate minimizing ``sum_i w_i term_i(w)`` minimizes the molar Gibbs
    energy at fixed ``(T, P, w)``: that sum is the only candidate-dependent part
    of ``G/RT`` (see the module docstring). Candidates are tried in order and
    ties keep the earlier one, so the selection is deterministic.

    A candidate that raises while ``optional`` is False re-raises: answering
    with the surviving candidate would report the wrong *phase*, not a missing
    branch.
    """
    best_terms: np.ndarray | None = None
    best_label = ""
    best_g = math.inf
    failures: list[str] = []

    for candidate in candidates:
        try:
            terms = candidate.ln_fugacity_terms(composition)
        except ModelError as exc:
            if not candidate.optional:
                raise
            failures.append(f"{candidate.label}: {exc}")
            continue

        g_res = float(np.sum(composition * terms))
        if not math.isfinite(g_res):
            failures.append(f"{candidate.label}: non-finite reduced residual Gibbs energy")
            continue
        if g_res < best_g:
            best_g = g_res
            best_terms = terms
            best_label = candidate.label

    if best_terms is None:
        raise ModelError(failure_message + " (" + "; ".join(failures) + ").")
    return best_terms, best_label


def _select_surface(
    candidates: Sequence[_PhaseCandidate],
    composition: np.ndarray,
    surface: str,
    *,
    failure_message: str,
) -> tuple[np.ndarray, str, bool]:
    """Return the terms of the candidate labelled ``surface`` at ``composition``.

    A trial pinned to one candidate iterates on that candidate's Gibbs surface
    (ADR-0012). Two things can still go wrong and they are treated differently:

    - the label is not one this evaluator holds. That is a programming error in
      the evaluator's own trial set, so it raises rather than guessing;
    - the candidate is ``optional`` and not evaluable at ``composition`` (an
      absent compressibility root). Then there is no surface to walk, and the
      only defined thing left is the lowest-Gibbs candidate; the caller is told
      through the returned flag so the fallback is recorded rather than hidden.

    A *mandatory* candidate that fails re-raises, exactly as in
    :func:`_select_min_gibbs`.
    """
    for candidate in candidates:
        if candidate.label != surface:
            continue
        try:
            terms = candidate.ln_fugacity_terms(composition)
        except ModelError:
            if not candidate.optional:
                raise
            break
        if math.isfinite(float(np.sum(composition * terms))):
            return terms, candidate.label, False
        break
    else:
        raise ModelError(
            f"Unknown phase-candidate surface {surface!r}; this evaluator holds "
            + ", ".join(repr(candidate.label) for candidate in candidates)
            + "."
        )

    terms, label = _select_min_gibbs(candidates, composition, failure_message=failure_message)
    return terms, label, True


# ---------------------------------------------------------------------------
# Evaluators
# ---------------------------------------------------------------------------


class _EOSTangentPlane:
    """Evaluator backed by an equation of state.

    The candidates are the ``"vapor"`` and ``"liquid"`` compressibility branches
    of the cubic; ``ln phi`` is taken from the minimum-Gibbs one, selected
    generically over the existing ``EquationOfState`` interface (ADR-0005).
    """

    model_family = "eos"
    pressure_dependent = True
    diagnostics: Mapping[str, float | int | str | bool] = {}

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
        self._candidates = _cubic_root_candidates(
            eos, mixture=mixture, temperature=temperature, pressure=pressure
        )

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        return _select_min_gibbs(
            self._candidates,
            composition,
            failure_message="No usable fugacity-coefficient branch for stability analysis",
        )

    def ln_terms_on_surface(
        self, composition: np.ndarray, surface: str
    ) -> tuple[np.ndarray, str | None, bool]:
        """Terms of one named compressibility branch.

        Implemented for contract completeness only: this evaluator's trial set
        names no surface (ADR-0012 keeps minimum-Gibbs root selection at every
        iterate for cubics), so the solver never calls it.
        """
        return _select_surface(
            self._candidates,
            composition,
            surface,
            failure_message="No usable fugacity-coefficient branch for stability analysis",
        )

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[_InitialEstimate]:
        """Two Wilson estimates plus one pure-component-dominant estimate each.

        No estimate names a surface: a cubic's trials keep re-selecting the
        minimum-Gibbs root at every iterate (ADR-0005, ADR-0012).
        """
        estimates: list[_InitialEstimate] = []

        k = wilson_k(self._mixture, self._temperature, self._pressure)
        estimates.append(_InitialEstimate("wilson-vapor", _normalized(k * z, active)))
        estimates.append(_InitialEstimate("wilson-liquid", _normalized(z / k, active)))
        estimates.extend(_pure_component_estimates(self._mixture, z, active))
        return estimates

    def identity_label(self, composition: np.ndarray, label: str | None) -> str | None:
        """Ask the model for a compressibility identity of the ``label`` root (ADR-0017).

        ``label`` names which of the two compressibility branches
        :func:`_select_min_gibbs` kept - the branch computed with
        ``eos.fugacity_coefficients(..., phase=label)`` - so asking
        ``eos.phase_identity(..., phase=label)`` measures the identity of that
        *same* root, never a different one. When the model returns ``None``
        (the ``EquationOfState.phase_identity`` default: not implemented) or
        raises evaluating a root that was just evaluated successfully (should
        not happen, defensive only), ``label`` is returned unchanged.
        """
        if label is None:
            return None
        try:
            identity = self._eos.phase_identity(
                mixture=self._mixture,
                temperature_K=self._temperature,
                pressure_Pa=self._pressure,
                composition=composition.tolist(),
                phase=label,
            )
        except ModelError:
            return label
        return identity if identity in (_LIQUID, _VAPOR) else label


class _ActivityTangentPlane:
    """Evaluator backed by an activity-coefficient model (liquid-liquid).

    Both phases are liquids with the same pure-liquid reference state, so the
    reference fugacities cancel from the tangent-plane distance and
    ``ln gamma_i`` takes the place of ``ln phi_i`` exactly (see the module
    docstring of :mod:`chemthermo.stability.tp`). There is a single candidate,
    so no selection is needed and the label is None.
    """

    model_family = "activity"
    pressure_dependent = False
    diagnostics: Mapping[str, float | int | str | bool] = {}

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
        self._candidate = _ActivityLiquidCandidate(
            activity_model, label=_LIQUID, mixture=mixture, temperature=temperature
        )

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        return self._candidate.ln_fugacity_terms(composition), None

    def ln_terms_on_surface(
        self, composition: np.ndarray, surface: str
    ) -> tuple[np.ndarray, str | None, bool]:
        """There is one candidate, so the only admissible surface is that one.

        Implemented for contract completeness only: this evaluator's trial set
        names no surface, so the solver never calls it.
        """
        return _select_surface(
            (self._candidate,),
            composition,
            surface,
            failure_message="No usable activity-model candidate for stability analysis",
        )

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[_InitialEstimate]:
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
            return [_InitialEstimate(f"pure-{names[index]}", _normalized(z, active))]
        return _pure_component_estimates(self._mixture, z, active)

    def identity_label(self, composition: np.ndarray, label: str | None) -> str | None:
        """No compressibility identity for an activity-only liquid: ``label`` (always ``None``)."""
        return label


class _ModifiedRaoultTangentPlane:
    """Evaluator holding a modified-Raoult liquid and an ideal-gas vapor (ADR-0010).

    Both candidates are expressed against the same reference ``ln(f_i / (x_i P))``:

        liquid: ln gamma_i(w) + ln( Psat_i(T) / P )
        vapor:  0

    so one tangent plane covers vapor-liquid *and* liquid-liquid behavior. The
    candidate label of a stationary point is what tells the flash whether the
    incipient phase is a vapor or a second liquid.
    """

    model_family = "modified-raoult"
    pressure_dependent = True

    def __init__(
        self,
        activity_model: ActivityModel,
        *,
        mixture: Mixture,
        temperature: float,
        pressure: float,
    ) -> None:
        self._mixture = mixture
        self._temperature = temperature
        self._pressure = pressure
        self._psat = antoine_saturation_pressures(mixture, temperature)
        self._k_raoult = self._psat / pressure
        lower, upper = antoine_temperature_range(mixture)
        self.diagnostics: Mapping[str, float | int | str | bool] = {
            "antoine_valid_Tmin_K": lower,
            "antoine_valid_Tmax_K": upper,
        }
        self._candidates: tuple[_PhaseCandidate, ...] = modified_raoult_candidates(
            activity_model,
            mixture=mixture,
            temperature=temperature,
            pressure=pressure,
            saturation_pressures=self._psat,
        )

    def ln_fugacity_terms(self, composition: np.ndarray) -> tuple[np.ndarray, str | None]:
        return _select_min_gibbs(
            self._candidates,
            composition,
            failure_message="No usable phase candidate for modified-Raoult stability analysis",
        )

    def ln_terms_on_surface(
        self, composition: np.ndarray, surface: str
    ) -> tuple[np.ndarray, str | None, bool]:
        return _select_surface(
            self._candidates,
            composition,
            surface,
            failure_message="No usable phase candidate for modified-Raoult stability analysis",
        )

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[_InitialEstimate]:
        """The trial set, one fixed candidate surface per trial (ADR-0012).

        For ``n`` active components the set is ``n + 2`` trials:

        ==================  ========  =====================================
        label               surface   initial estimate ``W0``
        ==================  ========  =====================================
        ``raoult-vapor``    vapor     ``z K^Raoult``
        ``raoult-liquid``   liquid    ``z / K^Raoult``
        ``pure-<name>``     liquid    component ``<name>`` dominant
        ==================  ========  =====================================

        ``K_i^Raoult = Psat_i / P`` is the ideal K-value of the same model with
        ``gamma = 1``, so ``W = z K^Raoult`` is a vapor-like estimate and
        ``W = z / K^Raoult`` a liquid-like one (the estimate a vapor feed needs
        in order to find its incipient liquid). The pure-component-dominant
        estimates are what finds a liquid-liquid split, exactly as for the
        activity-only evaluator.

        **The vapor surface needs exactly one trial, and its starting point is
        irrelevant.** The ideal-gas term is identically zero, so the
        successive-substitution map (equation (8) of
        :mod:`chemthermo.stability.tp`) on that surface is the *constant* map
        ``ln W_i <- d_i - 0 = d_i``: one substitution lands on the vapor
        surface's unique stationary point from any start, with an exactly zero
        residual at the next evaluation. Adding the pure-component estimates on
        the vapor surface would therefore add ``n`` trials that all return the
        same point as ``raoult-vapor``. They are deliberately not run, and this
        is why the trial count did not grow when the surfaces were fixed.

        The liquid surface has no such structure - ``ln gamma`` is a genuine
        function of ``w`` - so it keeps the liquid-like and the
        pure-component-dominant starts.
        """
        estimates: list[_InitialEstimate] = []
        if int(np.count_nonzero(active)) > 1:
            estimates.append(
                _InitialEstimate("raoult-vapor", _normalized(self._k_raoult * z, active), _VAPOR)
            )
            estimates.append(
                _InitialEstimate("raoult-liquid", _normalized(z / self._k_raoult, active), _LIQUID)
            )
            estimates.extend(_pure_component_estimates(self._mixture, z, active, surface=_LIQUID))
        else:
            # One active component: there is no composition degree of freedom,
            # so the feed itself is the only admissible trial and it must be
            # evaluated on the candidate the *feed* sits on. Naming a surface
            # here would pin the trial to the wrong one for half the feeds (a
            # pure vapor tested on the liquid surface can never meet the
            # stationarity condition), so this degenerate trial keeps the
            # minimum-Gibbs selection.
            names = self._mixture.component_names
            index = int(np.argmax(active))
            estimates.append(_InitialEstimate(f"pure-{names[index]}", _normalized(z, active)))
        return estimates

    def identity_label(self, composition: np.ndarray, label: str | None) -> str | None:
        """No compressibility identity: an activity liquid and an ideal gas are already
        two different models, not two branches of one equation of state, so there is
        no tie to break and ``label`` is returned unchanged.
        """
        return label


def _cubic_root_candidates(
    eos: EquationOfState,
    *,
    mixture: Mixture,
    temperature: float,
    pressure: float,
) -> tuple[_PhaseCandidate, ...]:
    """The two compressibility branches of a cubic, in the ADR-0005 order."""
    return tuple(
        _CubicRootCandidate(
            eos, label=label, mixture=mixture, temperature=temperature, pressure=pressure
        )
        for label in (_VAPOR, _LIQUID)
    )


def modified_raoult_candidates(
    activity_model: ActivityModel,
    *,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    saturation_pressures: np.ndarray | None = None,
) -> tuple[_PhaseCandidate, _PhaseCandidate]:
    """Return the ``(liquid, vapor)`` candidates of the modified-Raoult model.

    Exposed to :mod:`chemthermo.flash` (internal, not public) so that the split
    can evaluate each converged phase with the candidate the stability test
    assigned to it.
    """
    psat = (
        antoine_saturation_pressures(mixture, temperature)
        if saturation_pressures is None
        else saturation_pressures
    )
    liquid = _ActivityLiquidCandidate(
        activity_model,
        label=_LIQUID,
        mixture=mixture,
        temperature=temperature,
        reference_offset=np.log(psat / pressure),
    )
    vapor = _IdealVaporCandidate(label=_VAPOR)
    return liquid, vapor


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
    energy and the ideal-mixing contribution is root independent. This is the
    two-candidate case of :func:`_select_min_gibbs`.
    """
    return _select_min_gibbs(
        _cubic_root_candidates(eos, mixture=mixture, temperature=temperature, pressure=pressure),
        composition,
        failure_message="No usable fugacity-coefficient branch for stability analysis",
    )


def _pure_component_estimates(
    mixture: Mixture,
    z: np.ndarray,
    active: np.ndarray,
    *,
    surface: str | None = None,
) -> list[_InitialEstimate]:
    """One pure-component-dominant estimate per active component."""
    n_active = int(np.count_nonzero(active))
    if n_active <= 1:
        return []

    trace = _PURE_TRIAL_TRACE / (n_active - 1)
    names: Sequence[str] = mixture.component_names
    estimates: list[_InitialEstimate] = []
    for index in range(z.size):
        if not active[index]:
            continue
        w = np.where(active, trace, 0.0)
        w[index] = 1.0 - _PURE_TRIAL_TRACE
        estimates.append(_InitialEstimate(f"pure-{names[index]}", _normalized(w, active), surface))
    return estimates


def _normalized(values: np.ndarray, active: np.ndarray) -> np.ndarray:
    """Zero out inactive components and renormalize to sum to one."""
    w = np.where(active, np.asarray(values, dtype=float), 0.0)
    w = np.where(w > 0.0, w, 0.0)
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0.0:
        raise ModelError("Stability trial initial estimate has a non-positive total.")
    return w / total
