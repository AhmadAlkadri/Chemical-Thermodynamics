"""The **multiphase** split in log mole numbers (ADR-0029).

:mod:`chemthermo.flash._log_space` carries the *two-phase* split in
``u = ln n`` (ADR-0024, safeguarded in ADR-0026, seeded in ADR-0028). This
module is that same parametrization for a phase set of any size: it is to
:func:`chemthermo.flash._multiphase._multiphase_second_order` what
:func:`chemthermo.flash._log_space.log_space_split` is to
:func:`chemthermo.flash._second_order._second_order_split`.

Why the linear stage is not enough
----------------------------------
Take one mole of feed and let ``n_i^j`` be the moles of component ``i`` in each
*non-reference* phase; the reference phase holds ``n_i^r = z_i - sum_j n_i^j``.
:func:`chemthermo.flash._multiphase._multiphase_second_order` minimizes the
reduced Gibbs energy in those ``n`` with a **central difference of 1e-7** and
the box ``0 < n_i^j`` and ``sum_j n_i^j < z_i``.

A Peng-Robinson water / ethanol / n-hexane three-liquid set breaks both halves
of that at once (validation Case P-18 (ii)). The water-rich liquid holds
n-hexane at ``x = 1.6e-12`` (``5.3e-14`` at 280 K), so

- a step of ``1e-7`` in any non-reference ``n_i`` is eleven orders of magnitude
  larger than the reference phase's whole inventory of that component: the very
  first Hessian column leaves the box, the stage aborts at its first iteration,
  and successive substitution's residual is reported as a failure;
- and **no** step size repairs it while that component is in the *reference*
  phase, because ``n_i^r`` is there a difference of ``O(1)`` numbers -
  ``0.2 - 0.167 - 0.0326 = 2e-13`` - whose relative accuracy is already spent.
  A logarithmic step on the non-reference phases does not help: the quantity
  that cannot be resolved is not a variable at all.

So this module changes two things, and it needs both:

1. **The reference phase is re-chosen** to the phase whose smallest mole
   fraction (over the components present in the feed) is largest. That is the
   phase in which ``z - sum n`` is a well-conditioned difference, and it moves
   every trace composition into the *variables*, where a logarithm can hold it:
   ``n = 2e-13`` is ``u = -29``, an ordinary double. The reference phase of the
   multiphase Rachford-Rice is a free choice - equation (1) of
   :mod:`chemthermo.flash._multiphase` is symmetric in the phases - so nothing
   thermodynamic is being decided here.
2. **The variables are ``u_i^j = ln n_i^j``**, with the multiplicative
   central-difference step of :mod:`chemthermo.flash._log_space`. The gradient
   of the reduced Gibbs energy in the new variable is ``dg/du_i^j = n_i^j
   r_i^j`` by the chain rule, where ``r`` is the equal-fugacity residual of
   :mod:`chemthermo.flash._multiphase` equation (4), so the fixed points of the
   two parametrizations are the same and only their conditioning differs.

The Newton system, the descent test, the ADR-0026 curvature safeguard and the
line-search acceptance rule are :func:`chemthermo.flash._log_space.log_space_
split`'s, term for term, with the residual and the mole numbers flattened over
``(phase, component)`` instead of over ``component``. The constants are
imported from that module rather than restated, so the two stages cannot drift
apart.

Dormancy
--------
:func:`chemthermo.flash._multiphase._solve_phase_set` calls this **only** where
it was about to raise: the phase set has no phase to remove and the residual is
still above ``FlashSettings.tol`` after successive substitution and the linear
second-order stage. A result that was ever returned therefore cannot move.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np

from ..exceptions import ModelError
from ._log_space import (
    _ARMIJO_C,
    _CURVATURE_EIGENVALUE_FLOOR,
    _HESSIAN_STEP,
    _MIN_LINE_SEARCH_SCALE,
    _log_sum_exp,
)
from .settings import FlashSettings

#: Mole number a component is restarted at, relative to the feed, when the
#: linear stage left it at an exact zero in a phase that must still hold it.
#: The same rule and the same factor as
#: :func:`chemthermo.flash._log_space.seed_from_iterate`.
_RESTART_FLOOR = 1e-8


class MultiphaseLogSpaceSplit:
    """Outcome of the multiphase log-space stage, in the caller's phase order.

    Attributes:
        compositions: Normalized composition of every phase, in the **input**
            order (the internal reference re-ordering is undone before this is
            built), so it pairs positionally with the caller's labels and
            surfaces.
        fractions: Phase mole fractions, same order.
        terms: Each phase's tangent-plane fugacity terms at its composition.
        residual: ``max`` over components and **non-reference phases** of the
            equal-fugacity difference against the reference phase, which is the
            quantity this stage's Newton system drives to zero. It is within a
            factor two of the all-pairs measure of
            :func:`chemthermo.flash._multiphase._phase_residual`; the caller
            re-measures with that function before comparing anything, because
            that is the number a result reports.
        iterations: Newton iterations actually taken.
        reference: Index (in the input order) of the phase used as the
            reference during the solve.
    """

    __slots__ = ("compositions", "fractions", "iterations", "reference", "residual", "terms")

    def __init__(
        self,
        *,
        compositions: list[np.ndarray],
        fractions: np.ndarray,
        terms: list[np.ndarray],
        residual: float,
        iterations: int,
        reference: int,
    ) -> None:
        self.compositions = compositions
        self.fractions = fractions
        self.terms = terms
        self.residual = residual
        self.iterations = iterations
        self.reference = reference


def reference_phase_index(compositions: Sequence[np.ndarray], active: np.ndarray) -> int:
    """Index of the phase that makes ``z - sum n`` best conditioned.

    The reference phase's mole numbers are never variables: they are the
    difference between the feed and everything else. A component that phase
    holds only in trace is therefore a catastrophic cancellation, and the
    remedy is to put a phase that holds *every* component in quantity there.
    "In quantity" is measured as the smallest mole fraction over the components
    present in the feed, and the largest such minimum wins; ties keep the
    earlier phase, so the choice is deterministic.
    """
    best = 0
    best_minimum = -math.inf
    for index, composition in enumerate(compositions):
        values = np.asarray(composition, dtype=float)[active]
        minimum = float(np.min(values)) if values.size else 0.0
        if minimum > best_minimum:
            best, best_minimum = index, minimum
    return best


def _seed(
    *,
    z: np.ndarray,
    compositions: Sequence[np.ndarray],
    fractions: np.ndarray,
    index: np.ndarray,
) -> np.ndarray:
    """``u = ln n`` of the non-reference phases, from a linear iterate.

    ``compositions`` / ``fractions`` are already in solve order (reference
    first). A component the linear iterate left at an exact zero in a phase is
    restarted a factor :data:`_RESTART_FLOOR` below the feed rather than at
    minus infinity, exactly as
    :func:`chemthermo.flash._log_space.seed_from_iterate` does.
    """
    rows = []
    floor = _RESTART_FLOOR * z[index]
    for phase in range(1, len(compositions)):
        moles = float(fractions[phase]) * np.asarray(compositions[phase], dtype=float)[index]
        with np.errstate(divide="ignore"):
            rows.append(np.log(np.where(moles > 0.0, moles, floor)))
    return np.array(rows, dtype=float)


class _Iterate:
    """One admissible point of the multiphase log-space stage."""

    __slots__ = ("compositions", "energy", "moles", "residual", "terms", "totals")

    def __init__(
        self,
        *,
        energy: float,
        residual: np.ndarray,
        compositions: list[np.ndarray],
        terms: list[np.ndarray],
        moles: np.ndarray,
        totals: np.ndarray,
    ) -> None:
        self.energy = energy
        self.residual = residual
        self.compositions = compositions
        self.terms = terms
        self.moles = moles
        self.totals = totals


def _directions(
    *,
    jacobian: np.ndarray | None,
    residual: np.ndarray,
    moles: np.ndarray,
    curvature_safeguard: bool,
) -> list[tuple[np.ndarray, float]]:
    """Descent directions for ``g``, best first.

    The list is :func:`chemthermo.flash._log_space._safeguarded_directions`
    with the flattened ``(phase, component)`` residual in place of the
    two-phase one: the Newton direction on ``r = 0`` while it descends ``g``;
    then, when ``curvature_safeguard`` is set, the modified-Newton direction of
    the Gibbs Hessian ``H = diag(n) J + diag(n r)`` with its eigenvalues
    replaced by their floored magnitudes (Gill & Murray, ADR-0026); then
    ``-r``, the log-space successive-substitution step, which is always a
    descent direction.
    """
    gradient = moles * residual
    candidates: list[tuple[np.ndarray, float]] = []

    def add(direction: np.ndarray) -> None:
        if not np.all(np.isfinite(direction)):
            return
        slope = float(gradient @ direction)
        if slope < 0.0:
            candidates.append((direction, slope))

    if jacobian is not None:
        try:
            add(np.linalg.solve(jacobian, -residual))
        except np.linalg.LinAlgError:  # pragma: no cover - a singular Jacobian
            pass
        if curvature_safeguard:
            with np.errstate(over="ignore", under="ignore"):
                hessian = moles[:, None] * jacobian + np.diag(gradient)
            if np.all(np.isfinite(hessian)):
                symmetric = 0.5 * (hessian + hessian.T)
                eigenvalues, vectors = np.linalg.eigh(symmetric)
                magnitudes = np.abs(eigenvalues)
                floor = _CURVATURE_EIGENVALUE_FLOOR * float(np.max(magnitudes))
                if floor > 0.0:
                    with np.errstate(over="ignore", under="ignore"):
                        add(-(vectors @ ((vectors.T @ gradient) / np.maximum(magnitudes, floor))))

    candidates.append((-residual, float(gradient @ -residual)))
    return candidates


def multiphase_log_space_split(
    *,
    z: np.ndarray,
    compositions: Sequence[np.ndarray],
    fractions: np.ndarray,
    terms_by_phase: Sequence[Callable[[np.ndarray], np.ndarray]],
    settings: FlashSettings,
    curvature_safeguard: bool = False,
) -> MultiphaseLogSpaceSplit | None:
    """Damped Newton on a phase set of any size in ``u = ln n`` (ADR-0029).

    Args:
        z: Feed mole fractions (normalized).
        compositions: The linear stage's compositions, one per phase, in the
            caller's order.
        fractions: The linear stage's phase mole fractions, same order.
        terms_by_phase: Each phase's tangent-plane term callable, same order.
        settings: ``second_order_tol`` and ``second_order_max_iter``.
        curvature_safeguard: Use the ADR-0026 direction list. The caller walks
            ``False`` then ``True``, and reaches either only on a phase set
            that was about to be refused.

    Returns:
        The best point reached, in the caller's phase order, or ``None`` when
        the stage cannot start (fewer than two phases, no component present, or
        a seed that is already outside the box). The caller decides whether
        ``residual`` is good enough.
    """
    z = np.asarray(z, dtype=float)
    active = z > 0.0
    index = np.flatnonzero(active)
    width = index.size
    count = len(compositions)
    others = count - 1
    if width == 0 or others < 1:
        return None

    reference = reference_phase_index(compositions, active)
    order = [reference, *[phase for phase in range(count) if phase != reference]]
    ordered_compositions = [np.asarray(compositions[phase], dtype=float) for phase in order]
    ordered_fractions = np.array([float(fractions[phase]) for phase in order])
    ordered_terms = [terms_by_phase[phase] for phase in order]
    size = others * width
    z_active = z[index]

    def evaluate(u: np.ndarray) -> _Iterate:
        with np.errstate(over="ignore", under="ignore"):
            moles = np.exp(u)
        if not np.all(np.isfinite(moles)):
            raise ValueError("mole numbers left the representable range")
        remainder = z_active - moles.sum(axis=0)
        if np.any(remainder <= 0.0):
            raise ValueError("the reference phase emptied of a component")
        every = np.vstack([remainder[None, :], moles])
        totals = every.sum(axis=1)
        if np.any(totals <= 0.0):
            raise ValueError("a phase is empty")
        phase_compositions: list[np.ndarray] = []
        log_compositions: list[np.ndarray] = []
        reference_log = np.log(remainder) - math.log(float(totals[0]))
        for phase in range(count):
            if phase == 0:
                log_x = reference_log
            else:
                log_x = u[phase - 1] - _log_sum_exp(u[phase - 1])
            log_compositions.append(log_x)
            full = np.zeros_like(z)
            with np.errstate(under="ignore"):
                full[index] = np.exp(log_x)
            phase_compositions.append(full)
        terms = [
            function(composition)
            for function, composition in zip(ordered_terms, phase_compositions)
        ]
        activities = [log_compositions[phase] + terms[phase][index] for phase in range(count)]
        energy = float(sum(every[phase] @ activities[phase] for phase in range(count)))
        if not math.isfinite(energy):
            raise ValueError("non-finite Gibbs energy")
        residual = np.concatenate([activities[phase] - activities[0] for phase in range(1, count)])
        return _Iterate(
            energy=energy,
            residual=residual,
            compositions=phase_compositions,
            terms=terms,
            moles=moles.reshape(size),
            totals=totals,
        )

    u = _seed(z=z, compositions=ordered_compositions, fractions=ordered_fractions, index=index)
    try:
        current = evaluate(u)
    except (ValueError, ModelError):
        return None
    residual = float(np.max(np.abs(current.residual)))
    iterations = 0

    for iteration in range(1, settings.second_order_max_iter + 1):
        if residual < settings.second_order_tol:
            break
        iterations = iteration

        jacobian = np.zeros((size, size), dtype=float)
        usable = True
        for column in range(size):
            phase, component = divmod(column, width)
            plus = u.copy()
            minus = u.copy()
            plus[phase, component] += _HESSIAN_STEP
            minus[phase, component] -= _HESSIAN_STEP
            try:
                forward = evaluate(plus).residual
                backward = evaluate(minus).residual
            except (ValueError, ModelError):
                usable = False
                break
            jacobian[:, column] = (forward - backward) / (2.0 * _HESSIAN_STEP)

        directions = _directions(
            jacobian=jacobian if usable else None,
            residual=current.residual,
            moles=current.moles,
            curvature_safeguard=curvature_safeguard,
        )

        accepted = False
        for direction, slope in directions:
            scale = 1.0
            while scale >= _MIN_LINE_SEARCH_SCALE:
                candidate = u + scale * direction.reshape(others, width)
                try:
                    trial = evaluate(candidate)
                except (ValueError, ModelError):
                    scale *= 0.5
                    continue
                trial_residual = float(np.max(np.abs(trial.residual)))
                if (
                    trial.energy < current.energy + _ARMIJO_C * scale * slope
                    or trial_residual < residual
                ):
                    u = candidate
                    current = trial
                    residual = trial_residual
                    accepted = True
                    break
                scale *= 0.5
            if accepted:
                break

        if not accepted:
            break

    inverse = [order.index(phase) for phase in range(count)]
    return MultiphaseLogSpaceSplit(
        compositions=[current.compositions[position] for position in inverse],
        fractions=np.array([float(current.totals[position]) for position in inverse]),
        terms=[current.terms[position] for position in inverse],
        residual=residual,
        iterations=iterations,
        reference=reference,
    )
