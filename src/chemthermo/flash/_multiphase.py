"""Multiphase split, and the phase addition / removal loop (ADR-0011).

How many phases a feed splits into is an **output**, discovered the same way
the one-versus-two decision already is: by tangent-plane stability. This module
holds the three pieces that turn that into an answer with more than two phases.

1. Multiphase successive substitution
-------------------------------------
A phase set is a list of phase *candidates* (a modified-Raoult activity liquid,
an ideal vapor, ...) each carrying a composition. Choose the first as the
reference ``r``. Equal fugacity of component ``i`` between phase ``j`` and the
reference is

    ln x_i^j + t_i^j(x^j) = ln x_i^r + t_i^r(x^r)                          (1)

where ``t_i^j`` is the phase candidate's tangent-plane fugacity term - the same
``ln( f_i / (x_i P) )`` the stability evaluator uses, so
``ln gamma_i + ln(Psat_i/P)`` for a liquid and ``0`` for an ideal vapor. (1)
rearranges to the successive-substitution update

    ln K_i^j = t_i^r(x^r) - t_i^j(x^j),     K_i^j = x_i^j / x_i^r          (2)

Feeding those ``K`` to the multiphase Rachford-Rice solver of
:mod:`chemthermo.flash._multiphase_rr` gives the phase fractions and, through
``x_i^r = z_i / t_i`` and ``x_i^j = K_i^j z_i / t_i``, the next compositions.
That is exactly the two-phase loop of :mod:`chemthermo.flash._split` with the
scalar Rachford-Rice replaced by its multiphase form; for two phases the two
are the same iteration.

2. A second-order stage
-----------------------
Successive substitution on (1) converges linearly, and slowly: on the ternary
tie-triangle validated in this slice it needs 324 iterations to reach 1e-12,
where the second-order stage below reaches 7e-16 in 3 iterations from a
20-iteration start. The stage generalizes the two-phase Newton minimization of
:mod:`chemthermo.flash._second_order` to any number of phases.

Take one mole of feed and let ``n_i^j`` be the moles of component ``i`` in each
*non-reference* phase ``j``, so the reference phase holds
``n_i^r = z_i - sum_{j != r} n_i^j``. With ``N^j = sum_i n_i^j`` and
``x_i^j = n_i^j / N^j`` the part of the reduced Gibbs energy that depends on
the split is

    g(n) = sum_j sum_i n_i^j [ ln x_i^j + t_i^j(x^j) ]                     (3)

Differentiating, the terms in which the logarithms move cancel phase by phase:
for any phase, ``sum_i n_i d ln(x_i f_i) = (sum_i dn_i - dN) + 0 = 0``, the
first bracket because ``d ln x_i = dn_i / n_i - dN / N`` and the second by
Gibbs-Duhem at fixed ``T, P``. What survives is

    dg / dn_k^j = [ ln x_k^j + t_k^j(x^j) ] - [ ln x_k^r + t_k^r(x^r) ]    (4)

so **the gradient is exactly the equal-fugacity residual of (1)**: a stationary
point of (3) is equilibrium and a *minimum* of (3) is the equilibrium rather
than any other stationary point - in particular not the trivial solution, which
is a stationary ridge at ``g = g(feed)``. The Hessian is built by central
differences of (4), symmetrized and shifted to positive definiteness, and the
step is taken with a backtracking line search inside the box
``0 < n_i^j`` and ``sum_j n_i^j < z_i`` (every phase present, no negative mole
numbers).

3. Phase addition and removal
-----------------------------
    solve the phase set
      -> a phase fraction <= 0, or two phases collapsed?  remove it, re-solve
      -> a converged phase unstable?                      add it, re-solve
      -> all fractions positive and every phase stable?   that is the answer

*Addition* seeds the new phase from the tangent-plane minimizer found by the
post-split stability test on the failing phase, which is Michelsen's incipient
phase: setting ``x^new = w`` and computing ``K`` from (2) reproduces his
``W``-scaled seed ``K_i = W_i / x_i^r`` exactly, because at a stationary point
``ln W_i = ln x_i^r + t_i^r - t_i(w)``.

*Removal* is what the ordinary Rachford-Rice cannot express and the Okuno
formulation can: the feasible region constrains the phase *compositions*, not
the signs of the phase fractions, so a phase that should not be there converges
to a non-positive fraction (the "negative flash") instead of making the solve
fail. That is the resolution of validation Case R-3: in the ~0.135 K window
below the binary three-phase temperature the deepest tangent-plane minimum is a
vapor, so the search starts ``L -> LV``, the pair is unstable towards a second
liquid, ``LV -> LLV``, and the three-phase solve drives the vapor fraction
negative, ``LLV -> LL``. The correct two-liquid answer is reached by *removing*
a phase that addition had to add first.

The loop is bounded by ``FlashSettings.max_phases`` (default 3) and by a round
budget; exceeding either raises :class:`chemthermo.ConvergenceError` rather
than returning an unverified phase set.
"""

from __future__ import annotations

import math
from typing import Callable, Literal, Mapping, Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError, ModelError
from ..models import ActivityModel
from ._assemble import _multi_phase_result
from ._multiphase_rr import _multiphase_rachford_rice, _NoMultiphaseSolution
from ._verify import _is_same_phase, _post_split_report, _PostSplitReport, _reduced_g
from .results import FlashResult
from .settings import FlashSettings

#: Candidate labels. These are also the phase names of a vapor-liquid result.
_LIQUID = "liquid"
_VAPOR = "vapor"

#: Central-difference step for the multiphase second-order Hessian.
_HESSIAN_STEP = 1e-7
#: Smallest accepted backtracking scale in the second-order line search.
_MIN_LINE_SEARCH_SCALE = 1e-14
#: Armijo constant for the line search on the multiphase Gibbs energy.
_ARMIJO_C = 1e-4
#: A converged phase fraction at or below this is treated as "this phase is not
#: there". Zero rather than a positive tolerance: a genuinely tiny phase (Okuno
#: et al.'s Example 3 converges to beta = 2.2e-06) is a real phase, and the
#: post-split stability test - not a magnitude - is what decides existence.
_ABSENT_PHASE_FRACTION = 0.0
#: Extra rounds allowed beyond the additions and removals a successful search
#: needs, so that one add-then-remove detour still terminates cleanly.
_EXTRA_ROUNDS = 4


class _MultiphaseSolution:
    """A converged (or removal-flagged) multiphase split.

    Attributes:
        labels: Phase-candidate label of each phase, reference phase first.
        compositions: Normalized composition of each phase.
        fractions: Phase mole fractions, same order. May contain a non-positive
            entry, which is the removal signal.
        terms: The tangent-plane fugacity terms of each phase at its own
            composition.
        residual: ``max`` over all phase pairs and components of
            ``|ln(x_i^j f_i^j) - ln(x_i^k f_i^k)|``.
        ssi_iterations: Successive-substitution iterations performed.
        second_order_iterations: Second-order iterations performed.
        rr_iterations: Rachford-Rice Newton iterations summed over the solve.
        converged_stage: Which stage met the tolerance, or None.
        removal_index: Index of the phase whose fraction is non-positive, or
            None when every fraction is positive.
    """

    __slots__ = (
        "compositions",
        "converged_stage",
        "fractions",
        "labels",
        "removal_index",
        "residual",
        "rr_iterations",
        "second_order_iterations",
        "ssi_iterations",
        "terms",
    )

    def __init__(
        self,
        *,
        labels: tuple[str, ...],
        compositions: list[np.ndarray],
        fractions: np.ndarray,
        terms: list[np.ndarray],
        residual: float,
        ssi_iterations: int,
        second_order_iterations: int,
        rr_iterations: int,
        converged_stage: str | None,
        removal_index: int | None,
    ) -> None:
        self.labels = labels
        self.compositions = compositions
        self.fractions = fractions
        self.terms = terms
        self.residual = residual
        self.ssi_iterations = ssi_iterations
        self.second_order_iterations = second_order_iterations
        self.rr_iterations = rr_iterations
        self.converged_stage = converged_stage
        self.removal_index = removal_index


def _phase_residual(compositions: Sequence[np.ndarray], terms: Sequence[np.ndarray]) -> float:
    """``max`` over phase pairs and components of the equal-fugacity difference."""
    worst = 0.0
    for first in range(len(compositions)):
        for second in range(first + 1, len(compositions)):
            x, y = compositions[first], compositions[second]
            both = (x > 0.0) & (y > 0.0)
            if not np.any(both):
                continue
            difference = (np.log(x[both]) + terms[first][both]) - (
                np.log(y[both]) + terms[second][both]
            )
            worst = max(worst, float(np.max(np.abs(difference))))
    return worst


class _SsiOutcome:
    """What one run of multiphase successive substitution produced."""

    __slots__ = (
        "compositions",
        "fractions",
        "iterations",
        "removal_index",
        "residual",
        "rr_iterations",
        "terms",
    )

    def __init__(
        self,
        *,
        compositions: list[np.ndarray],
        fractions: np.ndarray,
        terms: list[np.ndarray],
        residual: float,
        iterations: int,
        rr_iterations: int,
        removal_index: int | None,
    ) -> None:
        self.compositions = compositions
        self.fractions = fractions
        self.terms = terms
        self.residual = residual
        self.iterations = iterations
        self.rr_iterations = rr_iterations
        self.removal_index = removal_index


def _multiphase_ssi(
    *,
    z: np.ndarray,
    compositions: list[np.ndarray],
    terms_by_phase: Sequence[Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
    budget: int,
) -> _SsiOutcome:
    """Successive substitution (2) with the multiphase Rachford-Rice inner solve.

    A phase set with no Rachford-Rice solution at all (the feasible region
    recedes, so ``F`` has no minimum - Gibbs' phase rule forbidding, say, three
    phases in a binary away from its three-phase temperature) is *not* an
    error: the recession direction names the phase whose amount runs negative,
    and that index is returned for removal.
    """
    beta: np.ndarray | None = None
    fractions = np.zeros(len(compositions))
    evaluated = [function(x) for function, x in zip(terms_by_phase, compositions)]
    residual = _phase_residual(compositions, evaluated)
    rr_iterations = 0
    iterations = 0

    for iteration in range(1, budget + 1):
        iterations = iteration
        K = np.column_stack(
            [np.exp(evaluated[0] - evaluated[index]) for index in range(1, len(compositions))]
        )
        if np.any(~np.isfinite(K)) or np.any(K <= 0.0):
            raise ConvergenceError("Non-finite or non-positive K-values in the multiphase split.")

        try:
            solution = _multiphase_rachford_rice(z, K, beta0=beta, tol=min(settings.tol, 1e-12))
        except _NoMultiphaseSolution as recession:
            rates = recession.fraction_rates
            departing = int(np.argmin(rates))
            if float(rates[departing]) >= 0.0:  # pragma: no cover - a recession
                raise  # direction always drives some fraction down
            return _SsiOutcome(
                compositions=compositions,
                fractions=fractions,
                terms=evaluated,
                residual=residual,
                iterations=iteration,
                rr_iterations=rr_iterations,
                removal_index=departing,
            )
        beta = solution.beta
        rr_iterations += solution.iterations
        fractions = solution.phase_fractions

        updated = [z / solution.t] + [K[:, index] * z / solution.t for index in range(K.shape[1])]
        compositions = []
        for candidate in updated:
            total = float(np.sum(candidate))
            if not math.isfinite(total) or total <= 0.0:
                raise ConvergenceError("A multiphase split phase collapsed to zero total.")
            compositions.append(candidate / total)

        evaluated = [function(x) for function, x in zip(terms_by_phase, compositions)]
        residual = _phase_residual(compositions, evaluated)
        if residual < settings.tol:
            break
        if float(np.min(fractions)) <= _ABSENT_PHASE_FRACTION:
            # A phase is leaving. Iterating a negative flash further only
            # refines a split the caller is about to discard.
            break

    smallest = int(np.argmin(fractions))
    return _SsiOutcome(
        compositions=compositions,
        fractions=fractions,
        terms=evaluated,
        residual=residual,
        iterations=iterations,
        rr_iterations=rr_iterations,
        removal_index=(smallest if float(fractions[smallest]) <= _ABSENT_PHASE_FRACTION else None),
    )


def _multiphase_second_order(
    *,
    z: np.ndarray,
    compositions: Sequence[np.ndarray],
    fractions: np.ndarray,
    terms_by_phase: Sequence[Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray], float, int] | None:
    """Damped Newton minimization of (3) in the non-reference mole numbers.

    Returns ``None`` when the starting point is not usable (a degenerate phase),
    otherwise ``(compositions, fractions, terms, residual, iterations)``.
    """
    active = z > 0.0
    index = np.flatnonzero(active)
    others = len(compositions) - 1
    width = index.size
    size = others * width
    if size == 0 or others == 0:
        return None

    def unpack(vector: np.ndarray) -> list[np.ndarray]:
        moles = []
        for phase in range(others):
            values = np.zeros_like(z)
            values[index] = vector[phase * width : (phase + 1) * width]
            moles.append(values)
        return moles

    def energy_and_gradient(
        vector: np.ndarray,
    ) -> tuple[float, np.ndarray, list[np.ndarray], list[np.ndarray]]:
        moles = unpack(vector)
        reference = z - sum(moles)
        every = [reference, *moles]
        totals = [float(np.sum(values)) for values in every]
        if any(total <= 0.0 for total in totals):
            raise ValueError("degenerate phase")
        if any(np.any(values[active] <= 0.0) for values in every):
            raise ValueError("negative mole numbers")
        phase_compositions = [values / total for values, total in zip(every, totals)]
        activities = []
        for composition, function in zip(phase_compositions, terms_by_phase):
            value = np.zeros_like(z)
            value[active] = np.log(composition[active]) + function(composition)[active]
            activities.append(value)
        energy = float(
            sum(
                np.sum(values[active] * activity[active])
                for values, activity in zip(every, activities)
            )
        )
        gradient = np.concatenate(
            [(activities[phase + 1] - activities[0])[index] for phase in range(others)]
        )
        return energy, gradient, phase_compositions, activities

    vector = np.concatenate(
        [(float(fractions[phase + 1]) * compositions[phase + 1])[index] for phase in range(others)]
    )
    try:
        energy, gradient, phase_compositions, _ = energy_and_gradient(vector)
    except (ValueError, ModelError):
        return None
    residual = float(np.max(np.abs(gradient)))

    identity = np.eye(size)
    iterations = 0
    for iteration in range(1, settings.second_order_max_iter + 1):
        if residual < settings.second_order_tol:
            break
        iterations = iteration

        hessian = np.zeros((size, size), dtype=float)
        try:
            for column in range(size):
                plus = vector.copy()
                minus = vector.copy()
                plus[column] += _HESSIAN_STEP
                minus[column] -= _HESSIAN_STEP
                _, gradient_plus, _, _ = energy_and_gradient(plus)
                _, gradient_minus, _, _ = energy_and_gradient(minus)
                hessian[:, column] = (gradient_plus - gradient_minus) / (2.0 * _HESSIAN_STEP)
        except (ValueError, ModelError):
            break

        hessian = 0.5 * (hessian + hessian.T)
        direction: np.ndarray
        try:
            smallest = float(np.min(np.linalg.eigvalsh(hessian)))
            shift = 0.0 if smallest > 1e-10 else (1e-10 - smallest)
            direction = np.linalg.solve(hessian + shift * identity, -gradient)
        except np.linalg.LinAlgError:  # pragma: no cover - the shift keeps it solvable
            direction = -gradient
        if not np.all(np.isfinite(direction)) or float(gradient @ direction) >= 0.0:
            direction = -gradient

        slope = float(gradient @ direction)
        scale = 1.0
        accepted = False
        while scale >= _MIN_LINE_SEARCH_SCALE:
            candidate = vector + scale * direction
            try:
                (
                    candidate_energy,
                    candidate_gradient,
                    candidate_compositions,
                    _,
                ) = energy_and_gradient(candidate)
            except (ValueError, ModelError):
                scale *= 0.5
                continue
            candidate_residual = float(np.max(np.abs(candidate_gradient)))
            if (
                candidate_energy < energy + _ARMIJO_C * scale * slope
                or candidate_residual < residual
            ):
                vector = candidate
                energy = candidate_energy
                gradient = candidate_gradient
                residual = candidate_residual
                phase_compositions = candidate_compositions
                accepted = True
                break
            scale *= 0.5

        if not accepted:
            break

    moles = unpack(vector)
    every = [z - sum(moles), *moles]
    refined_fractions = np.array([float(np.sum(values)) for values in every])
    evaluated = [
        function(composition) for function, composition in zip(terms_by_phase, phase_compositions)
    ]
    return (
        list(phase_compositions),
        refined_fractions,
        evaluated,
        _phase_residual(phase_compositions, evaluated),
        iterations,
    )


def _solve_phase_set(
    *,
    z: np.ndarray,
    labels: Sequence[str],
    compositions: Sequence[np.ndarray],
    candidates: Mapping[str, Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
) -> _MultiphaseSolution:
    """Solve one fixed phase set: successive substitution, then the Newton stage."""
    terms_by_phase = [candidates[label] for label in labels]
    budget = (
        min(settings.ssi_iterations, settings.max_iter)
        if settings.second_order
        else settings.max_iter
    )

    outcome = _multiphase_ssi(
        z=z,
        compositions=[np.array(value, dtype=float) for value in compositions],
        terms_by_phase=terms_by_phase,
        settings=settings,
        budget=budget,
    )
    current = outcome.compositions
    fractions = outcome.fractions
    evaluated = outcome.terms
    residual = outcome.residual
    ssi_iterations = outcome.iterations
    rr_iterations = outcome.rr_iterations

    if outcome.removal_index is not None:
        return _MultiphaseSolution(
            labels=tuple(labels),
            compositions=list(current),
            fractions=fractions,
            terms=list(evaluated),
            residual=residual,
            ssi_iterations=ssi_iterations,
            second_order_iterations=0,
            rr_iterations=rr_iterations,
            converged_stage=None,
            removal_index=outcome.removal_index,
        )

    converged_stage = "successive-substitution" if residual < settings.tol else None
    second_order_iterations = 0
    if settings.second_order and residual > settings.second_order_tol:
        refined = _multiphase_second_order(
            z=z,
            compositions=current,
            fractions=fractions,
            terms_by_phase=terms_by_phase,
            settings=settings,
        )
        if refined is not None:
            (
                refined_compositions,
                refined_fractions,
                refined_terms,
                refined_residual,
                second_order_iterations,
            ) = refined
            if refined_residual < residual:
                current = refined_compositions
                fractions = refined_fractions
                evaluated = refined_terms
                residual = refined_residual
                converged_stage = "second-order" if residual < settings.tol else converged_stage

    smallest = int(np.argmin(fractions))
    removal_index = smallest if float(fractions[smallest]) <= _ABSENT_PHASE_FRACTION else None
    if removal_index is None:
        removal_index = _collapsed_phase(current, settings)

    if removal_index is None and residual > settings.tol:
        raise ConvergenceError(
            "The multiphase split did not converge; equal-fugacity residual="
            f"{residual:.3e} after {ssi_iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations."
        )

    return _MultiphaseSolution(
        labels=tuple(labels),
        compositions=list(current),
        fractions=fractions,
        terms=list(evaluated),
        residual=residual,
        ssi_iterations=ssi_iterations,
        second_order_iterations=second_order_iterations,
        rr_iterations=rr_iterations,
        converged_stage=converged_stage,
        removal_index=removal_index,
    )


def _collapsed_phase(compositions: Sequence[np.ndarray], settings: FlashSettings) -> int | None:
    """Index of a phase that has merged onto an earlier one, if any.

    Two phases with the same composition are one phase counted twice: the
    Rachford-Rice system is then singular and the split is the trivial
    solution. The measure is the stability module's trivial-solution metric,
    ``sum_i ln(x_i / x_i')^2``, so no new tolerance is introduced.
    """
    from ..stability import StabilitySettings

    trivial_tol = (settings.stability_settings or StabilitySettings()).trivial_tol
    for later in range(1, len(compositions)):
        for earlier in range(later):
            if _is_same_phase(compositions[later], compositions[earlier], trivial_tol):
                return later
    return None


def _phase_set_label(labels: Sequence[str]) -> str:
    """Compact name of a phase set, e.g. ``"LLV"`` for two liquids and a vapor.

    Liquids are written first so that the same phase set always has the same
    name whatever order the solve happens to hold it in, which is what makes
    ``phase_set_history`` comparable between runs and between feeds.
    """
    return "".join(sorted((label[0].upper() for label in labels), key=lambda c: c != "L"))


def _phase_names(labels: Sequence[str]) -> tuple[str, ...]:
    """Result phase names for a set of candidate labels.

    A single liquid is ``"liquid"``; two or more are ``"liquid1"``,
    ``"liquid2"``, ... in the order the split holds them. The numbering is a
    *role*, not an identity: nothing distinguishes two liquids the way
    volatility distinguishes a vapor from a liquid, so callers must compare the
    phase *set*. A vapor is always ``"vapor"``.
    """
    liquids = sum(1 for label in labels if label != _VAPOR)
    names: list[str] = []
    seen = 0
    for label in labels:
        if label == _VAPOR:
            names.append(_VAPOR)
            continue
        seen += 1
        names.append(_LIQUID if liquids == 1 else f"{_LIQUID}{seen}")
    return tuple(names)


def _phase_regime(labels: Sequence[str]) -> str:
    """``"single-phase"``, ``"VLE"``, ``"LLE"`` or ``"VLLE"``."""
    if len(labels) == 1:
        return "single-phase"
    vapors = sum(1 for label in labels if label == _VAPOR)
    liquids = len(labels) - vapors
    if vapors and liquids >= 2:
        return "VLLE"
    if vapors:
        return "VLE"
    return "LLE"


def _phase_state(count: int) -> str:
    return {1: "single_phase", 2: "two_phase", 3: "three_phase"}.get(count, f"{count}_phase")


def _flash_tp_phase_addition(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    candidates: Mapping[str, Callable[[np.ndarray], np.ndarray]],
    labels: Sequence[str],
    compositions: Sequence[np.ndarray],
    history: list[str],
    ln_f_feed: np.ndarray,
    two_phase_g_rt: float | None,
    activity_model: ActivityModel,
    vapor: Literal["none", "ideal"],
    settings: FlashSettings,
    base: Mapping[str, float | int | str | bool],
    additions: int = 0,
) -> FlashResult:
    """Discover the equilibrium phase count by adding and removing phases.

    Args:
        z: Feed mole fractions (normalized).
        candidates: Tangent-plane term callables keyed by candidate label.
        labels: Candidate labels of the starting phase set (the converged
            two-phase split), reference phase first.
        compositions: Compositions of that starting phase set.
        history: Phase-set labels visited so far, e.g. ``["L", "LV"]``; this
            function appends to it.
        additions: Phases the caller already added before handing over, counted
            into ``diagnostics["phases_added"]``.
        ln_f_feed: Tangent-plane terms of the single-phase feed, for
            ``delta_g_split_rt``.
        two_phase_g_rt: Reduced Gibbs energy of the two-phase candidate that
            started the search, for ``delta_g_vs_two_phase_rt``. None when
            there was none.
        settings: Flash settings; ``max_phases`` bounds the search.
        base: Diagnostics shared with the two-phase paths.

    Raises:
        ConvergenceError: If a phase set is unstable at ``max_phases`` phases,
            if the round budget is exhausted, or if a solve fails.
    """
    current_labels = list(labels)
    current = [np.array(value, dtype=float) for value in compositions]
    rounds = settings.max_phases + _EXTRA_ROUNDS
    removals = 0
    total_ssi = 0
    total_second_order = 0
    total_rr = 0

    for _round in range(rounds):
        solution = _solve_phase_set(
            z=z,
            labels=current_labels,
            compositions=current,
            candidates=candidates,
            settings=settings,
        )
        total_ssi += solution.ssi_iterations
        total_second_order += solution.second_order_iterations
        total_rr += solution.rr_iterations

        if solution.removal_index is not None:
            if len(current_labels) <= 2:
                raise ConvergenceError(
                    "A two-phase set converged to a non-positive phase fraction, which "
                    "would leave no split at all. This is a solver failure, not a phase "
                    "count: the tangent-plane test had already proved the feed unstable."
                )
            index = solution.removal_index
            current_labels.pop(index)
            current = [
                value for position, value in enumerate(solution.compositions) if position != index
            ]
            removals += 1
            history.append(_phase_set_label(current_labels))
            continue

        names = _phase_names(current_labels)
        report = _post_split_report(
            mixture,
            temperature,
            pressure,
            eos=None,
            activity_model=activity_model,
            phases=tuple(zip(names, solution.compositions)),
            settings=settings,
            vapor=vapor,
        )

        if report.status != "stable":
            if report.status == "inconclusive":
                raise ConvergenceError(
                    "A post-split stability test was inconclusive for phase(s) "
                    f"{', '.join(report.inconclusive)}, so flash_tp cannot decide whether "
                    f"the {len(current_labels)}-phase set is the answer."
                )
            if len(current_labels) >= settings.max_phases:
                detail = ", ".join(failure.phase_name for failure in report.instabilities)
                raise ConvergenceError(
                    f"The converged {len(current_labels)}-phase solution is not a stable "
                    f"phase set: the post-split stability test reports 'unstable' for "
                    f"phase(s) {detail} (most negative post-split tpd = "
                    f"{report.tpd_min:.6e}). A further phase is required, and "
                    f"FlashSettings.max_phases = {settings.max_phases} forbids it. Raise "
                    "max_phases, or pass FlashSettings(post_split_stability=False) to "
                    "receive this phase set anyway with the failure in diagnostics."
                )
            failure = report.instabilities[0]
            current_labels.append(failure.branch or _LIQUID)
            current = [*solution.compositions, failure.composition]
            additions += 1
            history.append(_phase_set_label(current_labels))
            continue

        return _assemble(
            mixture,
            temperature,
            pressure,
            z=z,
            solution=solution,
            report=report,
            history=history,
            ln_f_feed=ln_f_feed,
            two_phase_g_rt=two_phase_g_rt,
            base=base,
            additions=additions,
            removals=removals,
            total_ssi=total_ssi,
            total_second_order=total_second_order,
            total_rr=total_rr,
        )

    raise ConvergenceError(
        "The phase addition/removal search did not settle on a stable phase set within "
        f"{rounds} rounds; the sets visited were {' -> '.join(history)}."
    )


def _assemble(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    solution: _MultiphaseSolution,
    report: _PostSplitReport,
    history: list[str],
    ln_f_feed: np.ndarray,
    two_phase_g_rt: float | None,
    base: Mapping[str, float | int | str | bool],
    additions: int,
    removals: int,
    total_ssi: int,
    total_second_order: int,
    total_rr: int,
) -> FlashResult:
    """Verify the converged phase set and build the `FlashResult`."""
    # Report liquids before the vapor, whatever order the solve held them in,
    # so that a three-phase answer always reads liquid1 / liquid2 / vapor. The
    # Rachford-Rice reference phase (index 0 of the solve) is unaffected.
    order = sorted(range(len(solution.labels)), key=lambda index: solution.labels[index] == _VAPOR)
    labels = tuple(solution.labels[index] for index in order)
    names = _phase_names(labels)
    fractions = np.array([solution.fractions[index] for index in order])
    compositions = [solution.compositions[index] for index in order]
    phase_terms = [solution.terms[index] for index in order]

    recombined = sum(
        float(fraction) * composition for fraction, composition in zip(fractions, compositions)
    )
    mass_balance = float(np.max(np.abs(z - recombined)))
    split_g = float(
        sum(
            float(fraction) * _reduced_g(composition, terms)
            for fraction, composition, terms in zip(fractions, compositions, phase_terms)
        )
    )
    feed_g = _reduced_g(z, ln_f_feed)

    diagnostics: dict[str, float | int | str | bool] = {
        **base,
        "iterations": total_ssi + total_second_order,
        "converged": True,
        "termination_reason": "tolerance_met",
        "phase_count": len(names),
        "phase_state": _phase_state(len(names)),
        "phase_regime": _phase_regime(labels),
        "phase_set_history": " -> ".join(history),
        "phases_added": additions,
        "phases_removed": removals,
        "k_seed": "stability",
        "ssi_iterations": total_ssi,
        "second_order_iterations": total_second_order,
        "rachford_rice_iterations": total_rr,
        "converged_stage": solution.converged_stage or "successive-substitution",
        "equilibrium_residual": solution.residual,
        "mass_balance_residual": mass_balance,
        "delta_g_split_rt": split_g - feed_g,
        **report.diagnostics,
    }
    if two_phase_g_rt is not None:
        diagnostics["delta_g_vs_two_phase_rt"] = split_g - two_phase_g_rt

    return _multi_phase_result(
        mixture,
        temperature,
        pressure,
        compositions=compositions,
        fractions=fractions,
        names=names,
        diagnostics=diagnostics,
    )
