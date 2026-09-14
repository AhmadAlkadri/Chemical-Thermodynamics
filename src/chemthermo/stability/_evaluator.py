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
from ..flash._common import (
    EosBranches,
    as_float_array,
    eos_branch_terms,
    eos_branch_terms_all,
    wilson_k,
)
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
            iterate (the pre-ADR-0012 behavior, kept for an activity-only
            model and for a single-active-component feed).
    """

    label: str
    composition: np.ndarray
    surface: str | None = None


class _SurfaceTerms(NamedTuple):
    """What a pinned trial gets back from one evaluation at ``w``.

    Attributes:
        terms: The terms the *iteration* uses: those of the pinned candidate,
            or of the lowest-Gibbs one where the pinned candidate was not
            available (``fell_back``).
        label: The candidate ``terms`` came from.
        fell_back: True when the pinned candidate was not available at ``w``.
        min_gibbs: ``(terms, label)`` of the lowest-Gibbs candidate at ``w``
            when the evaluation already produced them, else None. This is a
            *result*, not a cache: an evaluator that has to look at every
            candidate to answer at all (the density roots of an equation of
            state, which have to be compared to see whether the pinned one
            exists) knows the lowest-Gibbs terms as a by-product, and
            :func:`chemthermo.stability.tp._reported_terms` needs exactly
            those to report the tangent-plane distance. Handing them back
            spares a second pass over the model; None simply means "ask".
    """

    terms: np.ndarray
    label: str | None
    fell_back: bool
    min_gibbs: tuple[np.ndarray, str | None] | None = None


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

    def ln_terms_on_surface(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """Return the terms of the *named* candidate at ``w`` (ADR-0012, ADR-0021).

        Args:
            composition: Normalized mole fractions ``w``.
            surface: Label of the candidate to evaluate.

        Returns:
            A :class:`_SurfaceTerms`. ``fell_back`` is True when the named
            candidate was not available at ``w``, in which case the
            lowest-Gibbs candidate was used instead and ``label`` names it.

        Raises:
            ModelError: If the named candidate is unknown, or if it is
                mandatory and unusable at ``w``, or if no candidate is usable.
        """
        ...

    def ln_report_terms(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """Lowest-Gibbs terms at ``w``, for a trial pinned to ``surface``, at its stop.

        Called once per pinned trial, where the solver has to look at the lower
        envelope anyway to report the tangent-plane distance. That is the one
        place a pinned trial can see both candidates for free, so it is also
        where the single-candidate condition is measured (``fell_back``): see
        :func:`_select_density_root_surface` and ADR-0021.

        Args:
            composition: Normalized mole fractions ``w``.
            surface: The surface the trial iterated on.

        Returns:
            A :class:`_SurfaceTerms` whose ``min_gibbs`` is always populated and
            whose ``fell_back`` says whether ``surface`` was unavailable as a
            distinct candidate at ``w``.

        Raises:
            ModelError: If no candidate is usable at ``w``.
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


class _EOSBranchState:
    """The ``(eos, mixture, T, P)`` every density-root candidate of one evaluator shares.

    It exists so that a selector which needs *all* the branches can say so in
    one call (ADR-0023). A candidate on its own can only ask the model for its
    own branch, and two such asks at the same composition make the model solve
    for its density roots twice - a cubic twice for Peng-Robinson, a
    1599-point isotherm scan plus a safeguarded Newton twice for PC-SAFT - for
    roots that are equal to the last bit. Holding the shared state in one
    object lets :func:`_all_branch_terms` recognise a homogeneous candidate
    set and route it through
    :func:`chemthermo.flash._common.eos_branch_terms_all`.

    The *iteration* of a pinned trial still goes one branch at a time
    (ADR-0021 decision 3): it has no use for the other branch, and asking for
    both would undo that saving.
    """

    __slots__ = ("eos", "mixture", "pressure", "temperature")

    def __init__(
        self,
        eos: EquationOfState,
        *,
        mixture: Mixture,
        temperature: float,
        pressure: float,
    ) -> None:
        self.eos = eos
        self.mixture = mixture
        self.temperature = temperature
        self.pressure = pressure

    def branch_terms(self, composition: np.ndarray, labels: Sequence[str]) -> EosBranches:
        """Every branch in ``labels`` at ``composition``, from one root solve where possible."""
        return eos_branch_terms_all(
            self.eos,
            mixture=self.mixture,
            temperature=self.temperature,
            pressure=self.pressure,
            composition=composition,
            phases=labels,
        )


class _CubicRootCandidate:
    """One compressibility branch of an equation of state (``ln phi_i``)."""

    optional = True

    def __init__(self, state: _EOSBranchState, *, label: str) -> None:
        self.label = label
        self.branch_state = state

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
        """``ln phi`` on this branch, asking the model for this branch alone.

        Delegated to :func:`chemthermo.flash._common.eos_branch_terms`, which
        takes ``np.log`` of the model's own ``phi`` whenever that is
        representable - the pre-ADR-0022 double, unchanged - and falls back to
        the model's logarithmic route only where ``exp(ln phi)`` has
        under/overflowed, which is what makes a long-chain polymer usable here
        at all.

        A selector that wants *both* branches does not come through here; see
        :func:`_all_branch_terms`.
        """
        state = self.branch_state
        return eos_branch_terms(
            state.eos,
            mixture=state.mixture,
            temperature=state.temperature,
            pressure=state.pressure,
            composition=composition,
            phase=self.label,
        ).ln_phi


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


def _all_branch_terms(
    candidates: Sequence[_PhaseCandidate], composition: np.ndarray
) -> EosBranches | None:
    """Pre-evaluate a homogeneous density-root candidate set in one root solve.

    Returns None - and the caller then evaluates each candidate the way it
    always did - unless every candidate is a :class:`_CubicRootCandidate`
    standing on the *same* :class:`_EOSBranchState`. That is the only case
    where one model call can answer for all of them (ADR-0023); the
    modified-Raoult pair, for instance, is two different models and has
    nothing to share.
    """
    states = {getattr(candidate, "branch_state", None) for candidate in candidates}
    if len(states) != 1:
        return None
    state = states.pop()
    if state is None:
        return None
    return state.branch_terms(composition, tuple(candidate.label for candidate in candidates))


def _candidate_terms(
    candidate: _PhaseCandidate, composition: np.ndarray, precomputed: EosBranches | None
) -> np.ndarray:
    """One candidate's terms, from the pre-evaluated branches when there are any.

    Raises:
        ModelError: With the message the per-branch route would have raised,
            since :func:`chemthermo.flash._common.eos_branch_terms_all` falls
            back to that route whenever the one-solve capability cannot answer.
    """
    if precomputed is None:
        return candidate.ln_fugacity_terms(composition)
    terms = precomputed.terms.get(candidate.label)
    if terms is not None:
        return terms.ln_phi
    raise ModelError(precomputed.failures[candidate.label])


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
    precomputed = _all_branch_terms(candidates, composition)

    for candidate in candidates:
        try:
            terms = _candidate_terms(candidate, composition, precomputed)
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
) -> _SurfaceTerms:
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
            return _SurfaceTerms(terms, candidate.label, False)
        break
    else:
        raise ModelError(
            f"Unknown phase-candidate surface {surface!r}; this evaluator holds "
            + ", ".join(repr(candidate.label) for candidate in candidates)
            + "."
        )

    terms, label = _select_min_gibbs(candidates, composition, failure_message=failure_message)
    return _SurfaceTerms(terms, label, True, (terms, label))


def _select_density_root_surface(
    candidates: Sequence[_PhaseCandidate],
    composition: np.ndarray,
    surface: str,
    *,
    failure_message: str,
) -> _SurfaceTerms:
    """Both branches at ``composition``: the lower envelope, and the root count.

    This is what an equation of state answers a pinned trial's *reporting*
    call with (:meth:`_EOSTangentPlane.ln_report_terms`), and it exists because
    a missing root is not reported the way a missing candidate is (ADR-0021).
    Both ``PengRobinsonEOS`` and ``PCSAFTEOS`` answer
    ``fugacity_coefficients(..., phase="vapor")`` and ``phase="liquid"`` with
    the *same* numbers where the model has one admissible root -- documented
    behaviour of both, not an accident -- so the pinned trial never sees a
    ``ModelError`` there, and :func:`_select_surface`, which evaluates only the
    named branch, would report ``fell_back = False`` at exactly the
    compositions where there was no surface to choose.

    Here both branches are evaluated, so the single-root case is identified
    exactly -- the two branches come from one density, so their terms are
    bit-identical -- and ``min_gibbs`` comes out of the same pass, which is the
    reporting call's actual purpose. The *iteration* does not pay for this: it
    goes through :func:`_select_surface` and evaluates one branch (ADR-0021
    decision 3), which is why a single-root region met mid-iteration is not
    counted. It cannot change a result: where the model has one root, the
    pinned surface and the lowest-Gibbs surface are the same surface.

    Raises:
        ModelError: If ``surface`` is not a label this evaluator holds, or if
            no branch is usable at ``composition``.
    """
    if not any(candidate.label == surface for candidate in candidates):
        raise ModelError(
            f"Unknown phase-candidate surface {surface!r}; this evaluator holds "
            + ", ".join(repr(candidate.label) for candidate in candidates)
            + "."
        )

    available: list[tuple[str, np.ndarray, float]] = []
    failures: list[str] = []
    precomputed = _all_branch_terms(candidates, composition)
    for candidate in candidates:
        try:
            terms = _candidate_terms(candidate, composition, precomputed)
        except ModelError as exc:
            if not candidate.optional:
                raise
            failures.append(f"{candidate.label}: {exc}")
            continue
        g_res = float(np.sum(composition * terms))
        if not math.isfinite(g_res):
            failures.append(f"{candidate.label}: non-finite reduced residual Gibbs energy")
            continue
        available.append((candidate.label, terms, g_res))

    if not available:
        raise ModelError(failure_message + " (" + "; ".join(failures) + ").")

    # Ties keep the earlier candidate, exactly as `_select_min_gibbs` does.
    best_label, best_terms, _best_g = min(available, key=lambda entry: entry[2])

    for label, terms, _g_res in available:
        if label != surface:
            continue
        degenerate = (
            all(
                np.array_equal(terms, other_terms)
                for other_label, other_terms, _other_g in available
                if other_label != surface
            )
            and len(available) > 1
        )
        if degenerate:
            return _SurfaceTerms(terms, best_label, True, (best_terms, best_label))
        return _SurfaceTerms(terms, label, False, (best_terms, best_label))

    return _SurfaceTerms(best_terms, best_label, True, (best_terms, best_label))


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

    def ln_terms_on_surface(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """Terms of one named density root, for one iterate of a pinned trial (ADR-0021).

        A trial pinned to a surface has no use for the other branch, so only
        the named one is evaluated - half the model calls the minimum-Gibbs
        rule of ADR-0005 made at the same iterate. Where the model has a single
        admissible root it answers *both* labels with that root (documented
        behaviour of ``PengRobinsonEOS`` and ``PCSAFTEOS`` alike), so the trial
        walks the only surface there without a second evaluation to discover
        it; that condition is measured at the trial's stopping point instead,
        by :meth:`ln_report_terms`.
        """
        return _select_surface(
            self._candidates,
            composition,
            surface,
            failure_message="No usable fugacity-coefficient branch for stability analysis",
        )

    def ln_report_terms(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """Both branches at the stopping point: the lower envelope, and the root count."""
        return _select_density_root_surface(
            self._candidates,
            composition,
            surface,
            failure_message="No usable fugacity-coefficient branch for stability analysis",
        )

    def initial_estimates(self, z: np.ndarray, active: np.ndarray) -> list[_InitialEstimate]:
        """Two Wilson estimates plus one pure-component-dominant estimate each.

        Each trial is pinned to the density root its start is an estimate *of*
        (ADR-0021, extending ADR-0012 from the modified-Raoult candidates to
        the roots of an equation of state):

        ==================  ========  =====================================
        label               surface   initial estimate ``W0``
        ==================  ========  =====================================
        ``wilson-vapor``    vapor     ``z K^Wilson``
        ``wilson-liquid``   liquid    ``z / K^Wilson``
        ``pure-<name>``     liquid    component ``<name>`` dominant
        ==================  ========  =====================================

        Re-selecting the minimum-Gibbs root at every iterate - what ADR-0012
        kept here - makes the successive-substitution map discontinuous across
        the composition where the two roots exchange Gibbs energy, and a
        vapour-like trial started above that crossing is dragged onto the
        liquid root and collapses onto the feed's partner liquid. Validation
        Case P-9 (iv) is that failure: from the hexane-rich liquid of a
        water / n-hexane feed above the three-phase temperature, all four
        trials returned ``tpd ~ -3e-09`` while a vapour stationary point with
        ``tpd = -6.5e-03`` sat unvisited.

        The pure-component-dominant estimates are liquid-like by construction
        (a nearly pure component at 1 atm is a dense fluid), and they are the
        estimates Michelsen recommends for finding a second *liquid*; running
        them on the vapour root as well was measured over the Peng-Robinson
        and PC-SAFT grids of Case P-11 and found no stationary point the four
        trials above do not already find, so those ``n`` extra trials are not
        run and the trial count is unchanged.

        A feed with a single active component is the exception, and for
        ADR-0012's reason: there is no composition degree of freedom, both
        Wilson estimates are the feed itself, and pinning them would test a
        root the feed is not on. That degenerate case keeps minimum-Gibbs
        selection.
        """
        estimates: list[_InitialEstimate] = []

        k = wilson_k(self._mixture, self._temperature, self._pressure)
        if int(np.count_nonzero(active)) <= 1:
            estimates.append(_InitialEstimate("wilson-vapor", _normalized(k * z, active)))
            estimates.append(_InitialEstimate("wilson-liquid", _normalized(z / k, active)))
            return estimates

        estimates.append(_InitialEstimate("wilson-vapor", _normalized(k * z, active), _VAPOR))
        estimates.append(_InitialEstimate("wilson-liquid", _normalized(z / k, active), _LIQUID))
        estimates.extend(_pure_component_estimates(self._mixture, z, active, surface=_LIQUID))
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

    def ln_terms_on_surface(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
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

    def ln_report_terms(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """Contract completeness only: one candidate, and no trial names a surface."""
        return self.ln_terms_on_surface(composition, surface)

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

    def ln_terms_on_surface(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        return _select_surface(
            self._candidates,
            composition,
            surface,
            failure_message="No usable phase candidate for modified-Raoult stability analysis",
        )

    def ln_report_terms(self, composition: np.ndarray, surface: str) -> _SurfaceTerms:
        """The lower envelope at the stopping point.

        Both candidates are evaluable at every composition (an activity liquid
        and an ideal gas), so ``fell_back`` here is always False and this is
        exactly the pre-ADR-0021 reporting call, bit for bit.
        """
        terms, label = self.ln_fugacity_terms(composition)
        return _SurfaceTerms(terms, label, False, (terms, label))

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
    """The two compressibility branches of a cubic, in the ADR-0005 order.

    Both stand on one :class:`_EOSBranchState`, which is what lets a selector
    that needs both of them get both from a single root solve (ADR-0023).
    """
    state = _EOSBranchState(eos, mixture=mixture, temperature=temperature, pressure=pressure)
    return tuple(_CubicRootCandidate(state, label=label) for label in (_VAPOR, _LIQUID))


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
