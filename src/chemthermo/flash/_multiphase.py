"""Multiphase split, and the phase addition / removal loop (ADR-0011).

How many phases a feed splits into is an **output**, discovered the same way
the one-versus-two decision already is: by tangent-plane stability. This module
holds the three pieces that turn that into an answer with more than two phases.

1. Multiphase successive substitution
-------------------------------------
A phase set is a list of phases, each carrying a composition **and its own
tangent-plane surface** - the thing that turns a composition into fugacities.
What that surface is differs by model family and is the whole of
:class:`_PhaseSetModel` (ADR-0020): a modified-Raoult *phase candidate* (an
activity liquid, an ideal vapor) whose label is also the phase's identity, or,
for an equation of state, one pinned density root
(:class:`chemthermo.flash._split._PhaseRoot`, ADR-0019) whose identity has to
be *measured*. The second case is why the surface travels with the phase
rather than being looked up from its label: a vapour and two liquids is a
phase set in which two phases carry the label ``"liquid"`` and are two
different phases on two different roots.

Choose the first phase as the reference ``r``. Equal fugacity of component
``i`` between phase ``j`` and the reference is

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
``ln W_i = ln x_i^r + t_i^r - t_i(w)``. The new phase's surface is pinned to
the branch that stability test reported for it, which is ADR-0019's rule for
the two phases of a two-phase split, applied to the third.

A minimizer can also *duplicate* a phase already in the set, and then the
solve removes it again before it can move. Adding the same one back would
cycle, so the next stationary point of the same report is tried instead, each
at most once (ADR-0020 decision 3, and the comment on that branch below for
the state that forced it).

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
from typing import Callable, Mapping, Protocol, Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ._assemble import _multi_phase_result
from ._multiphase_log_space import MultiphaseLogSpaceSplit, multiphase_log_space_split
from ._multiphase_rr import _multiphase_rachford_rice, _NoMultiphaseSolution
from ._split import _PhaseRoot
from ._verify import (
    _is_same_phase,
    _PhaseInstability,
    _post_split_report,
    _PostSplitReport,
    _reduced_g,
)
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


class _PhaseSetModel(Protocol):
    """The model side of a phase set: **each phase carries its own surface.**

    The search above is model-agnostic: it only ever needs, per phase, a
    tangent-plane fugacity term ``t_i(x)``. What supplies that term differs by
    family, and the difference is exactly the one ADR-0012 and ADR-0019 already
    drew:

    - modified-Raoult (ADR-0010, ADR-0011): the surface is the phase
      *candidate* the label names - an activity liquid with Antoine reference
      fugacities, or an ideal vapor - and the label **is** the phase identity,
      because those two candidates are two different models rather than two
      roots of one;
    - an equation of state (ADR-0019, generalized here to any number of
      phases): the surface is a :class:`chemthermo.flash._split._PhaseRoot`,
      one density/compressibility branch of one model, pinned for the whole
      solve. A vapour and two liquids then sit on three independent roots,
      which is precisely the bookkeeping a three-phase EOS set needs and the
      reason a single label-to-callable mapping is not enough: two phases may
      carry the *same* label ``"liquid"`` and still be two different phases.

    Which phase is called what is then a separate question from which surface
    it sits on, and only the EOS family has to *measure* it (ADR-0017).

    Complexity receipt: one protocol with four methods replaces a
    ``Mapping[label, callable]``. It buys the three-phase EOS set (two phases
    with one label), the per-phase root pinning ADR-0019 requires, and the
    ADR-0017 measurement of what each converged phase is; without it the loop
    would have to branch on the model family internally, which is what
    ADR-0007 forbids. Cost: two small classes below and one extra
    ``phase_identity`` call per phase per round on the EOS path.
    """

    #: Order two liquid phases by composition (ADR-0019 decision 3) rather than
    #: by the order the search created them (ADR-0011 decision 5).
    composition_ordered_liquids: bool

    def surface(self, label: str) -> Callable[[np.ndarray], np.ndarray]:
        """The tangent-plane term callable of a new phase pinned to ``label``."""
        ...  # pragma: no cover - protocol

    def identity(self, label: str, composition: np.ndarray) -> str | None:
        """What the phase *is* (``"liquid"`` / ``"vapor"``), or None if unmeasurable."""
        ...  # pragma: no cover - protocol

    def report(
        self,
        mixture: Mixture,
        temperature: float,
        pressure: float,
        *,
        phases: Sequence[tuple[str, np.ndarray]],
        settings: FlashSettings,
    ) -> _PostSplitReport:
        """Post-split stability of every converged phase, in this family."""
        ...  # pragma: no cover - protocol

    def label_diagnostics(self, measured: Sequence[str | None]) -> Mapping[str, str]:
        """Extra diagnostics recording how the phase names were decided."""
        ...  # pragma: no cover - protocol


class _ActivityPhaseSet:
    """Modified-Raoult phase set: the label names the candidate (ADR-0010).

    Unchanged behaviour, expressed through the protocol: ``surface`` hands back
    the very same bound method the pre-slice ``candidates[label]`` lookup did,
    so the arithmetic is identical, and ``identity`` is the label itself
    because an activity liquid and an ideal gas are two models, not two roots.
    """

    composition_ordered_liquids = False

    def __init__(
        self,
        activity_model: ActivityModel,
        *,
        candidates: Mapping[str, Callable[[np.ndarray], np.ndarray]],
    ) -> None:
        self._activity_model = activity_model
        self._candidates = candidates

    def surface(self, label: str) -> Callable[[np.ndarray], np.ndarray]:
        return self._candidates[label]

    def identity(self, label: str, composition: np.ndarray) -> str | None:
        return label

    def report(
        self,
        mixture: Mixture,
        temperature: float,
        pressure: float,
        *,
        phases: Sequence[tuple[str, np.ndarray]],
        settings: FlashSettings,
    ) -> _PostSplitReport:
        return _post_split_report(
            mixture,
            temperature,
            pressure,
            eos=None,
            activity_model=self._activity_model,
            phases=phases,
            settings=settings,
            vapor="ideal",
        )

    def label_diagnostics(self, measured: Sequence[str | None]) -> Mapping[str, str]:
        return {}


class _EosPhaseSet:
    """Equation-of-state phase set: the label names a density root (ADR-0019).

    Each phase gets its **own** :class:`chemthermo.flash._split._PhaseRoot`,
    pinned to the branch the tangent-plane test found that phase on and held
    for the whole solve, with the lowest-Gibbs rule as the fallback and the
    post-split stability test as the check - ADR-0019's rule, unchanged, just
    applied to three phases instead of two.

    The names come from ``EquationOfState.phase_identity`` measured on the root
    each phase converged on (ADR-0017), not from the pinned label, because two
    liquids and a vapour cannot be told apart by a label that says ``"liquid"``
    twice.
    """

    composition_ordered_liquids = True

    def __init__(
        self,
        eos: EquationOfState,
        mixture: Mixture,
        temperature: float,
        pressure: float,
    ) -> None:
        self._eos = eos
        self._mixture = mixture
        self._temperature = temperature
        self._pressure = pressure

    def surface(self, label: str) -> Callable[[np.ndarray], np.ndarray]:
        root = _PhaseRoot(self._eos, self._mixture, self._temperature, self._pressure, branch=label)
        return root.ln_fugacity_terms

    def identity(self, label: str, composition: np.ndarray) -> str | None:
        try:
            measured = self._eos.phase_identity(
                mixture=self._mixture,
                temperature_K=self._temperature,
                pressure_Pa=self._pressure,
                composition=np.asarray(composition, dtype=float).tolist(),
                phase=label,
            )
        except ModelError:
            return None
        return measured if measured in (_LIQUID, _VAPOR) else None

    def report(
        self,
        mixture: Mixture,
        temperature: float,
        pressure: float,
        *,
        phases: Sequence[tuple[str, np.ndarray]],
        settings: FlashSettings,
    ) -> _PostSplitReport:
        return _post_split_report(
            mixture,
            temperature,
            pressure,
            eos=self._eos,
            activity_model=None,
            phases=phases,
            settings=settings,
        )

    def label_diagnostics(self, measured: Sequence[str | None]) -> Mapping[str, str]:
        every = all(value in (_LIQUID, _VAPOR) for value in measured)
        return {"phase_label_method": "compressibility" if every else "tie-break"}


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
        removal_order: Every phase this solve names as removable, best first,
            so ``removal_order[0] == removal_index`` whenever there is one.
            The tail is what makes a removal **reversible** (ADR-0029): the
            search can come back and take the next candidate when the first one
            led nowhere. Empty when ``removal_index`` is None.
        log_space_iterations: Iterations spent in the ADR-0029 multiphase
            log-space stage. Zero unless that stage ran, which it does only
            where this solve was about to raise.
    """

    __slots__ = (
        "compositions",
        "converged_stage",
        "fractions",
        "labels",
        "log_space_iterations",
        "removal_index",
        "removal_order",
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
        removal_order: tuple[int, ...] = (),
        log_space_iterations: int = 0,
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
        self.removal_order = removal_order
        self.log_space_iterations = log_space_iterations


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


def _ranked_removals(values: Sequence[float] | np.ndarray) -> tuple[int, ...]:
    """Phases that ``values`` marks as removable, most negative first.

    ``values`` is either the converged phase fractions or, for a receding
    feasible region, the rate of change of each fraction along the recession
    direction. In both cases a non-positive entry says "this phase is not
    there", and ``argmin`` is the entry the search has always acted on; the
    ranking is that same rule read past its first place, so the first element
    is exactly the pre-ADR-0029 choice.
    """
    array = np.asarray(values, dtype=float)
    ranked = [int(index) for index in np.argsort(array, kind="stable")]
    return tuple(index for index in ranked if float(array[index]) <= _ABSENT_PHASE_FRACTION)


class _SsiOutcome:
    """What one run of multiphase successive substitution produced."""

    __slots__ = (
        "compositions",
        "fractions",
        "iterations",
        "removal_index",
        "removal_order",
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
        removal_order: tuple[int, ...] = (),
    ) -> None:
        self.compositions = compositions
        self.fractions = fractions
        self.terms = terms
        self.residual = residual
        self.iterations = iterations
        self.rr_iterations = rr_iterations
        self.removal_index = removal_index
        self.removal_order = removal_order


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

    Both removal signals name a *ranking*, not one phase: several fractions can
    be non-positive at once, and several rates can be negative along one
    recession direction. The first entry is the removal this function has
    always chosen; the rest are recorded because that choice can be wrong and
    the search has to be able to take it back (ADR-0029).
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
                removal_order=_ranked_removals(rates),
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
        removal_order=_ranked_removals(fractions),
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
    terms_by_phase: Sequence[Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
) -> _MultiphaseSolution:
    """Solve one fixed phase set: successive substitution, then the Newton stage.

    ``terms_by_phase`` pairs positionally with ``labels`` and ``compositions``:
    one surface per phase, so two phases may share a label and still be
    evaluated on two different density roots (ADR-0019).
    """
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
            removal_order=outcome.removal_order,
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

    removal_index, removal_order = _removal_signal(current, fractions, settings)

    log_space_iterations = 0
    if removal_index is None and residual > settings.tol:
        # ADR-0029. Reached only on the branch that raises below, so no result
        # that was ever returned can move. The linear stage above cannot be
        # made to work on a phase set holding a component at `x ~ 1e-12` in the
        # phase that carries the mass balance; see `_multiphase_log_space`.
        refined = _multiphase_log_space_stage(
            z=z,
            compositions=current,
            fractions=fractions,
            terms_by_phase=terms_by_phase,
            settings=settings,
        )
        if refined is not None:
            log_space_residual = _phase_residual(refined.compositions, refined.terms)
            if log_space_residual < residual:
                current = refined.compositions
                fractions = refined.fractions
                evaluated = refined.terms
                residual = log_space_residual
                log_space_iterations = refined.iterations
                converged_stage = "second-order-log" if residual < settings.tol else converged_stage
                removal_index, removal_order = _removal_signal(current, fractions, settings)

    if removal_index is None and residual > settings.tol:
        detail = (
            ""
            if log_space_iterations == 0
            else f" and {log_space_iterations} log-space Newton iterations (ADR-0029)"
        )
        raise ConvergenceError(
            "The multiphase split did not converge; equal-fugacity residual="
            f"{residual:.3e} after {ssi_iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations{detail}."
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
        removal_order=removal_order,
        log_space_iterations=log_space_iterations,
    )


def _removal_signal(
    compositions: Sequence[np.ndarray],
    fractions: np.ndarray,
    settings: FlashSettings,
) -> tuple[int | None, tuple[int, ...]]:
    """``(removal_index, removal_order)`` of a converged phase set.

    A non-positive phase fraction comes first (the negative flash), and two
    phases that have merged onto one another come second; the merged case names
    a single phase, because "these two are one" is not a ranking.
    """
    ranked = _ranked_removals(fractions)
    if ranked:
        return ranked[0], ranked
    collapsed = _collapsed_phase(compositions, settings)
    if collapsed is None:
        return None, ()
    return collapsed, (collapsed,)


def _multiphase_log_space_stage(
    *,
    z: np.ndarray,
    compositions: Sequence[np.ndarray],
    fractions: np.ndarray,
    terms_by_phase: Sequence[Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
) -> MultiphaseLogSpaceSplit | None:
    """Walk the ADR-0029 log-space stage, unsafeguarded then safeguarded.

    Two entries for the same reason ADR-0028's ladder has three: the
    unsafeguarded iteration is the cheaper one and is what every state
    measured for ADR-0029 needs, and the ADR-0026 curvature safeguard is what a
    phase set parked next to an indefinite Hessian would need. The best
    residual wins; ``None`` means neither entry could even start.
    """
    best: MultiphaseLogSpaceSplit | None = None
    for curvature_safeguard in (False, True):
        attempt = multiphase_log_space_split(
            z=z,
            compositions=compositions,
            fractions=fractions,
            terms_by_phase=terms_by_phase,
            settings=settings,
            curvature_safeguard=curvature_safeguard,
        )
        if attempt is None:
            return best
        if best is None or attempt.residual < best.residual:
            best = attempt
        if best.residual < settings.tol:
            break
    return best


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


class _RemovalChoice:
    """A phase set as it stood before a removal, and the removals not yet tried.

    Removal is the one step of the ADR-0011 search that was irreversible: the
    solve names the phase whose amount is most negative, the search drops it,
    and there is no way back. ADR-0029 makes it reversible, because the ranking
    the solve produces can be wrong - see the branch in
    :func:`_flash_tp_phase_addition` that pops this.

    Attributes:
        labels: Candidate labels of the pre-removal set.
        surfaces: Its per-phase tangent-plane term callables.
        compositions: The compositions that solve reached, pre-removal.
        remaining: Indices into the pre-removal set that have not been removed
            yet, best first.
    """

    __slots__ = ("compositions", "labels", "remaining", "surfaces")

    def __init__(
        self,
        *,
        labels: list[str],
        surfaces: list[Callable[[np.ndarray], np.ndarray]],
        compositions: list[np.ndarray],
        remaining: list[int],
    ) -> None:
        self.labels = labels
        self.surfaces = surfaces
        self.compositions = compositions
        self.remaining = remaining


def _record_removal(
    undo: list[_RemovalChoice],
    *,
    labels: Sequence[str],
    surfaces: Sequence[Callable[[np.ndarray], np.ndarray]],
    compositions: Sequence[np.ndarray],
    candidates: Sequence[int],
) -> None:
    """Remember a phase set and the removals it offers beyond the one taken.

    ``candidates[0]`` is the removal the caller is about to make, so only the
    tail is recorded. A set that offers no alternative is still pushed, with an
    empty tail, so the stack mirrors the search's own history.
    """
    undo.append(
        _RemovalChoice(
            labels=list(labels),
            surfaces=list(surfaces),
            compositions=[np.array(value, dtype=float) for value in compositions],
            remaining=[int(index) for index in candidates[1:]],
        )
    )


def _take_next_removal(
    undo: list[_RemovalChoice],
) -> tuple[list[str], list[Callable[[np.ndarray], np.ndarray]], list[np.ndarray]] | None:
    """Undo removals until one offers an untried candidate, and take it.

    Returns the restored phase set with that candidate removed, or ``None``
    when every removal on the stack has been exhausted - which is when the
    caller raises, exactly as it always did. Each candidate is taken at most
    once and the stack only shrinks, so the search still terminates.
    """
    while undo:
        choice = undo[-1]
        if not choice.remaining:
            undo.pop()
            continue
        index = choice.remaining.pop(0)
        labels = [value for position, value in enumerate(choice.labels) if position != index]
        surfaces = [value for position, value in enumerate(choice.surfaces) if position != index]
        compositions = [
            value for position, value in enumerate(choice.compositions) if position != index
        ]
        return labels, surfaces, compositions
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


def _ordered_phases(
    model: _PhaseSetModel,
    labels: Sequence[str],
    compositions: Sequence[np.ndarray],
) -> tuple[list[int], tuple[str, ...], tuple[str, ...], list[str | None]]:
    """Measure what each phase is, then order and name the set.

    Returns ``(order, identities, names, measured)`` where ``order`` lists the
    solve-order indices in *report* order (liquids before the vapour) and
    ``identities`` / ``names`` are in **solve order**, so a caller can index
    them alongside ``compositions``.

    Two rules differ by family and both live here rather than in the loop:

    - what a phase *is* comes from ``model.identity``, which is the label
      itself for the modified-Raoult candidates and an ADR-0017 compressibility
      measurement for an equation of state;
    - two liquids are ordered by the first component's mole fraction when the
      model says so (ADR-0019 decision 3, so the labels do not swap between two
      feeds on one tie line) and otherwise by the order the search created
      them (ADR-0011 decision 5, where they are roles).

    The vapour-last sort is stable, so on the modified-Raoult path this
    reproduces the pre-slice ordering and naming exactly.
    """
    measured = [
        model.identity(label, composition) for label, composition in zip(labels, compositions)
    ]
    identities = tuple(
        value if value is not None else label for value, label in zip(measured, labels)
    )
    if model.composition_ordered_liquids:
        order = sorted(
            range(len(identities)),
            key=lambda index: (
                identities[index] == _VAPOR,
                tuple(-value for value in compositions[index].tolist()),
            ),
        )
    else:
        order = sorted(range(len(identities)), key=lambda index: identities[index] == _VAPOR)

    ordered_names = _phase_names([identities[index] for index in order])
    names: list[str] = [""] * len(identities)
    for position, index in enumerate(order):
        names[index] = ordered_names[position]
    return order, identities, tuple(names), measured


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


def _extend(
    labels: list[str],
    surfaces: list[Callable[[np.ndarray], np.ndarray]],
    compositions: list[np.ndarray],
    *,
    model: _PhaseSetModel,
    failure: _PhaseInstability,
) -> int:
    """Append one incipient phase to the set in place; return its index."""
    label = failure.branch or _LIQUID
    labels.append(label)
    surfaces.append(model.surface(label))
    compositions.append(failure.composition)
    return len(labels) - 1


def _flash_tp_phase_addition(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    model: _PhaseSetModel,
    labels: Sequence[str],
    surfaces: Sequence[Callable[[np.ndarray], np.ndarray]],
    compositions: Sequence[np.ndarray],
    history: list[str],
    ln_f_feed: np.ndarray,
    two_phase_g_rt: float | None,
    settings: FlashSettings,
    base: Mapping[str, float | int | str | bool],
    additions: int = 0,
) -> FlashResult:
    """Discover the equilibrium phase count by adding and removing phases.

    Args:
        z: Feed mole fractions (normalized).
        model: The phase-set model (:class:`_PhaseSetModel`): it supplies a new
            phase's surface, the identity of a converged one, and the
            post-split stability report for this model family.
        labels: Candidate label of each phase of the starting set (the
            converged two-phase split), reference phase first. On the EOS path
            a label names a density root, so two phases may carry the same one.
        surfaces: Tangent-plane term callable of each starting phase, pairing
            positionally with ``labels`` and ``compositions``.
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
    current_surfaces = list(surfaces)
    current = [np.array(value, dtype=float) for value in compositions]
    rounds = settings.max_phases + _EXTRA_ROUNDS
    #: Stationary points of the last post-split report that have not been tried
    #: as the incipient phase, deepest first, and the index of the phase the
    #: search added last. See the "duplicate incipient phase" branch below.
    pending: list[_PhaseInstability] = []
    added_index: int | None = len(current_labels) - 1 if additions else None
    #: Removals whose *other* candidates have not been tried, newest last. See
    #: the "the removal was the wrong one" branch below (ADR-0029).
    undo: list[_RemovalChoice] = []
    removals = 0
    total_ssi = 0
    total_second_order = 0
    total_rr = 0
    total_log_space = 0

    for _round in range(rounds):
        solution = _solve_phase_set(
            z=z,
            labels=current_labels,
            compositions=current,
            terms_by_phase=current_surfaces,
            settings=settings,
        )
        total_ssi += solution.ssi_iterations
        total_second_order += solution.second_order_iterations
        total_rr += solution.rr_iterations
        total_log_space += solution.log_space_iterations

        if solution.removal_index is not None:
            if len(current_labels) <= 2:
                # ADR-0029: the removal that produced this pair was a *choice*
                # among the phases the previous solve named, and a two-phase
                # set with a non-positive fraction is that choice turning out
                # to be wrong - the feed is not inside this pair's tie line,
                # while the tangent-plane test has proved it is not one phase
                # either. Undo the removal and take the next candidate. The
                # measured state: PC-SAFT water / n-hexane, z_water = 0.05,
                # just above T3, where LLV recedes, the recession direction
                # names the *hexane-rich* liquid as the phase leaving, and the
                # water-rich pair that is left over then negative-flashes; the
                # answer is the pair the other removal leaves (Case P-18 (i)).
                restored = _take_next_removal(undo)
                if restored is None:
                    raise ConvergenceError(
                        "A two-phase set converged to a non-positive phase fraction, which "
                        "would leave no split at all. This is a solver failure, not a phase "
                        "count: the tangent-plane test had already proved the feed unstable."
                    )
                current_labels, current_surfaces, current = restored
                removals += 1
                added_index = None
                history.append(_phase_set_label(current_labels))
                continue
            index = solution.removal_index
            _record_removal(
                undo,
                labels=current_labels,
                surfaces=current_surfaces,
                compositions=solution.compositions,
                candidates=solution.removal_order,
            )
            # ADR-0020: a phase the search has just added, removed again by the
            # multiphase Rachford-Rice *before it could move*, is a stationary
            # point that duplicates a phase already in the set - not a phase
            # that tried to exist and could not. Adding it again would cycle
            # (measured: PC-SAFT water / n-hexane, z = 0.5/0.5, 1 atm, 330 K,
            # where the deepest minimum found from the vapour is the water-rich
            # liquid already present, and the search runs
            # LV -> LLV -> LV -> ... until the round budget). The *next*
            # stationary point of the same report is tried instead, each one at
            # most once, so the rule terminates with the phase count.
            duplicate = index == added_index and bool(pending)
            current_labels.pop(index)
            current_surfaces.pop(index)
            current = [
                value for position, value in enumerate(solution.compositions) if position != index
            ]
            removals += 1
            added_index = None
            history.append(_phase_set_label(current_labels))
            if duplicate:
                failure = pending.pop(0)
                added_index = _extend(
                    current_labels,
                    current_surfaces,
                    current,
                    model=model,
                    failure=failure,
                )
                additions += 1
                history.append(_phase_set_label(current_labels))
            continue

        order, identities, names, measured = _ordered_phases(
            model, current_labels, solution.compositions
        )
        report = model.report(
            mixture,
            temperature,
            pressure,
            phases=tuple(zip(names, solution.compositions)),
            settings=settings,
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
            pending = list(report.instabilities[1:])
            current = list(solution.compositions)
            added_index = _extend(
                current_labels,
                current_surfaces,
                current,
                model=model,
                failure=report.instabilities[0],
            )
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
            model=model,
            order=order,
            identities=identities,
            measured=measured,
            history=history,
            ln_f_feed=ln_f_feed,
            two_phase_g_rt=two_phase_g_rt,
            base=base,
            additions=additions,
            removals=removals,
            total_ssi=total_ssi,
            total_second_order=total_second_order,
            total_rr=total_rr,
            total_log_space=total_log_space,
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
    model: _PhaseSetModel,
    order: Sequence[int],
    identities: Sequence[str],
    measured: Sequence[str | None],
    history: list[str],
    ln_f_feed: np.ndarray,
    two_phase_g_rt: float | None,
    base: Mapping[str, float | int | str | bool],
    additions: int,
    removals: int,
    total_ssi: int,
    total_second_order: int,
    total_rr: int,
    total_log_space: int,
) -> FlashResult:
    """Verify the converged phase set and build the `FlashResult`.

    ``order`` / ``identities`` / ``measured`` come from :func:`_ordered_phases`
    on the same compositions, so the ``phase_stability_<name>`` keys the
    post-split report just wrote carry the same names this result does.
    """
    # Report liquids before the vapor, whatever order the solve held them in,
    # so that a three-phase answer always reads liquid1 / liquid2 / vapor. The
    # Rachford-Rice reference phase (index 0 of the solve) is unaffected.
    labels = tuple(identities[index] for index in order)
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
        **model.label_diagnostics(measured),
        **report.diagnostics,
    }
    if total_log_space:
        # Conditional, on the ADR-0016 principle: a phase set that converged
        # without the ADR-0029 stage carries the diagnostics mapping it carried
        # before that slice, key for key.
        diagnostics["log_space_iterations"] = total_log_space
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
