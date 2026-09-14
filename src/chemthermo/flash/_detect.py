"""Tangent-plane phase detection: decide 1-vs-2 phases and seed the split.

Both reference paths - phi-phi (:func:`_flash_tp_tangent_plane`) and
gamma-gamma (:func:`_flash_tp_liquid_liquid`) - decide one phase versus two
from Michelsen's tangent-plane stability criterion (ADR-0008, ADR-0009)
rather than from Wilson K-value bounds:

    stability_tp(feed) -> single phase | seeded split
                        -> post-split stability of every converged phase

``stability_tp`` is run on the feed; a stable feed returns a single-phase
result immediately, an unstable feed seeds the shared split loop
(:mod:`chemthermo.flash._split`) from the stationary point
(:func:`_stability_k_seed`), and the converged split is verified and
post-split-checked (:mod:`chemthermo.flash._verify`) before being assembled
into a `FlashResult` (:mod:`chemthermo.flash._assemble`). An inconclusive
stability result raises rather than guessing.
"""

from __future__ import annotations

import math

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ._assemble import _single_phase_result, _two_phase_result
from ._common import wilson_k
from ._log_space import (
    has_trace_component,
    log_space_seed,
    log_space_split,
    seed_from_iterate,
)
from ._multiphase import (
    _ActivityPhaseSet,
    _EosPhaseSet,
    _flash_tp_phase_addition,
    _phase_set_label,
)
from ._second_order import _second_order_split
from ._split import (
    _ln_gamma_function,
    _PhaseRoot,
    _rachford_rice,
    _solve_k_loop,
    _SplitSolution,
)
from ._verify import (
    _equilibrium_residual,
    _post_split_report,
    _post_split_stability,
    _reduced_g,
    _verify_split,
)
from .results import FlashResult
from .settings import FlashSettings

#: Seed K-value used for components absent from the feed (``z_i == 0``). Those
#: components have ``x_i = y_i = 0`` at every iteration and are rewritten from
#: the model on the first update, so the seed value cannot affect the result.
_INERT_SEED_K = 1.0

#: Phase names of a liquid-liquid result. Roles, not identities; see
#: :func:`chemthermo.flash.tp.flash_tp`.
_LIQUID_I = "liquid1"
_LIQUID_II = "liquid2"

#: Phase-candidate labels of the modified-Raoult pair, which are also the phase
#: names of a vapor-liquid result on that path.
_LIQUID = "liquid"
_VAPOR = "vapor"
_MODIFIED_RAOULT = "modified-raoult"

#: ``diagnostics["phase_label_method"]`` values (ADR-0017, ADR-0019): the phase
#: name(s) came from ``EquationOfState.phase_identity`` ("compressibility"), or
#: from the pre-ADR-0017 convention because the model does not implement it, or
#: because the two phases of a split both measured "vapor" ("wilson-ranking" -
#: the historical name, kept because that is what decided the orientation the
#: label defaults to; see :func:`_name_two_phase_result` and
#: :func:`_phase_label_method`). Two phases that both measure "liquid" are no
#: longer a fallback case since ADR-0019: they are named ``liquid1`` /
#: ``liquid2`` and the method stays "compressibility".
_LABEL_COMPRESSIBILITY = "compressibility"
_LABEL_WILSON_RANKING = "wilson-ranking"
_LABEL_TIE_BREAK = "tie-break"

#: ``diagnostics["k_seed"]`` values whose split is seeded from the
#: tangent-plane stationary point, and whose two phases are therefore pinned to
#: the branches that stationary point named (ADR-0019). ``"stability-log"`` is
#: the ADR-0024 seed: the same stationary point, carried as logarithms because
#: its K-values are not doubles.
_STABILITY_SEEDS = ("stability", "stability-log")


def _phase_label_method(
    eos: EquationOfState,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    composition: np.ndarray,
    phase_name: str,
) -> str:
    """How a single-phase result's ``phase_name`` was decided (ADR-0017).

    ``phase_name`` is already the (possibly compressibility-relabeled) branch
    ``stability_tp`` reported; this asks the model once more, on the same root,
    purely to record *which rule* produced it - "compressibility" when
    ``eos.phase_identity`` is implemented and usable here, "tie-break"
    (the pre-ADR-0017 min-Gibbs convention) otherwise. It cannot change
    ``phase_name`` itself: in the single-real-root case where the label was
    ambiguous, ``phase="liquid"`` and ``phase="vapor"`` name the same root and
    return the same identity either way.
    """
    try:
        identity = eos.phase_identity(
            mixture=mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=composition.tolist(),
            phase=phase_name,
        )
    except ModelError:
        identity = None
    return _LABEL_COMPRESSIBILITY if identity in (_LIQUID, _VAPOR) else _LABEL_TIE_BREAK


class _TwoPhaseNaming:
    """What the two converged phi-phi phases are called, and why.

    Attributes:
        names: Pairs positionally with ``(x, y)``: ``names[0]`` is ``x``'s.
        vapor_fraction: The value reported as ``FlashResult.vapor_fraction`` -
            ``beta``, ``1 - beta``, or None when the phase set holds no vapour.
        regime: ``diagnostics["phase_regime"]``.
        method: ``diagnostics["phase_label_method"]``.
        identities: The measured identity of each phase's converged root.
    """

    __slots__ = ("identities", "method", "names", "regime", "vapor_fraction")

    def __init__(
        self,
        *,
        names: tuple[str, str],
        vapor_fraction: float | None,
        regime: str,
        method: str,
        identities: tuple[str, str],
    ) -> None:
        self.names = names
        self.vapor_fraction = vapor_fraction
        self.regime = regime
        self.method = method
        #: What each phase's converged root *is*, measured by ADR-0017 on that
        #: root, falling back to the candidate label the selector used when the
        #: model cannot measure one. Reported as ``phase_i_branch`` /
        #: ``phase_ii_branch``.
        self.identities = identities


def _phase_identity_on(
    eos: EquationOfState,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    composition: np.ndarray,
    branch: str,
) -> str | None:
    """``eos.phase_identity`` on ``branch``, or None when it cannot be measured."""
    try:
        return eos.phase_identity(
            mixture=mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=composition.tolist(),
            phase=branch,
        )
    except ModelError:
        return None


def _liquid_order(x: np.ndarray, y: np.ndarray) -> tuple[str, str]:
    """Name two liquid phases ``liquid1`` / ``liquid2`` deterministically.

    ADR-0019 decision 3. Unlike the gamma-gamma path - whose ``liquid1`` /
    ``liquid2`` are *roles assigned by the seed* and may swap between two feeds
    on the same tie line - the phi-phi pair is ordered by composition:
    ``liquid1`` is the phase with the **larger mole fraction of the first
    component**, ties broken by the second component and so on, and finally by
    position. Two feeds on one tie line therefore come back with the same
    labels on the same phases, which is what makes a lever-rule check
    meaningful. The order is relative to the mixture's component order, so
    permuting the components permutes which phase is ``liquid1``; the phase
    *set* is unchanged.

    Returns:
        ``(name_of_x, name_of_y)``.
    """
    for value_x, value_y in zip(x.tolist(), y.tolist()):
        if value_x > value_y:
            return (_LIQUID_I, _LIQUID_II)
        if value_y > value_x:
            return (_LIQUID_II, _LIQUID_I)
    return (_LIQUID_I, _LIQUID_II)


def _name_two_phase_result(
    eos: EquationOfState,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    x: np.ndarray,
    y: np.ndarray,
    beta: float,
    branches: tuple[str, str],
) -> _TwoPhaseNaming:
    """Name the two converged phi-phi phases (ADR-0017, extended by ADR-0019).

    Each phase is asked for its own compressibility identity **on the root it
    actually converged on** (``branches``, from
    :class:`chemthermo.flash._split._PhaseRoot`) - which since ADR-0019 need
    not be ``("liquid", "vapor")``. This never touches ``x``, ``y`` or ``beta``
    themselves, only which name and which ``vapor_fraction`` value goes with
    each:

    ==========================  =====================================  ========
    measured identities         names (positional with ``(x, y)``)     regime
    ==========================  =====================================  ========
    one liquid, one vapour      ``"liquid"`` / ``"vapor"``, the         VLE
                                vapour-identified phase carrying
                                ``vapor_fraction``
    both liquid                 ``"liquid1"`` / ``"liquid2"``           LLE
                                (:func:`_liquid_order`),
                                ``vapor_fraction = None``
    both vapour, or either      ``"liquid"`` / ``"vapor"`` in the       VLE
    identity unavailable        historical ``(x, y)`` orientation
    ==========================  =====================================  ========

    The last row is the documented last-resort fallback: a model that does not
    implement ``phase_identity``, or a near-critical split whose phases both
    measure ``"vapor"``. It keeps ADR-0008 decision 3's Wilson-ranking
    orientation and records ``phase_label_method = "wilson-ranking"``. Ranking
    two same-side phases by the *magnitude* of ``kappa`` would need the number
    itself, which ``EquationOfState`` does not expose (only the verdict); that
    is deliberately left to a later slice rather than guessed at here.
    """
    identity_x = _phase_identity_on(eos, mixture, temperature, pressure, x, branches[0])
    identity_y = _phase_identity_on(eos, mixture, temperature, pressure, y, branches[1])
    identities = (identity_x or branches[0], identity_y or branches[1])

    if identity_x in (_LIQUID, _VAPOR) and identity_y in (_LIQUID, _VAPOR):
        if identity_x != identity_y:
            if identity_x == _VAPOR:
                return _TwoPhaseNaming(
                    names=(_VAPOR, _LIQUID),
                    vapor_fraction=1.0 - beta,
                    regime="VLE",
                    method=_LABEL_COMPRESSIBILITY,
                    identities=identities,
                )
            return _TwoPhaseNaming(
                names=(_LIQUID, _VAPOR),
                vapor_fraction=beta,
                regime="VLE",
                method=_LABEL_COMPRESSIBILITY,
                identities=identities,
            )
        if identity_x == _LIQUID:
            return _TwoPhaseNaming(
                names=_liquid_order(x, y),
                vapor_fraction=None,
                regime="LLE",
                method=_LABEL_COMPRESSIBILITY,
                identities=identities,
            )

    return _TwoPhaseNaming(
        names=(_LIQUID, _VAPOR),
        vapor_fraction=beta,
        regime="VLE",
        method=_LABEL_WILSON_RANKING,
        identities=identities,
    )


def _phi_phi_roots(
    eos: EquationOfState,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    feed_branch: str | None,
    incipient_branch: str | None,
    incipient_phase: str,
    seed_label: str,
) -> tuple[_PhaseRoot, _PhaseRoot]:
    """The two per-phase root holders of a phi-phi split, pinned (ADR-0019).

    ``_stability_k_seed`` sets ``K = y / x`` so that ``y`` is the phase the
    Wilson ranking calls vapour-like. Which of ``x`` and ``y`` is therefore the
    *feed-like* phase and which the *incipient* one depends on
    ``incipient_phase``:

    - ``incipient_phase == "vapor"``: ``K = W / z``, so ``x`` is feed-like and
      ``y`` incipient;
    - ``incipient_phase == "liquid"``: ``K = z / W``, so ``x`` is incipient and
      ``y`` feed-like.

    Each phase is pinned to the branch the stability test reported for the
    phase it came from (``feed_branch`` for the feed-like phase,
    ``phase_branch`` for the incipient one) and stays there for the whole
    split; see :class:`chemthermo.flash._split._PhaseRoot` for why the branch
    is held rather than re-selected per iterate. When the Rachford-Rice
    fallback replaced the stability seed with Wilson K-values
    (``seed_label == "wilson"``) those roles do not exist, so the historical
    ``("liquid", "vapor")`` pinning is used instead.
    """
    branches: tuple[str | None, str | None]
    if seed_label not in _STABILITY_SEEDS:
        branches = (_LIQUID, _VAPOR)
    elif incipient_phase == _VAPOR:
        branches = (feed_branch, incipient_branch)
    else:
        branches = (incipient_branch, feed_branch)
    return (
        _PhaseRoot(eos, mixture, temperature, pressure, branch=branches[0]),
        _PhaseRoot(eos, mixture, temperature, pressure, branch=branches[1]),
    )


def _flash_tp_tangent_plane(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Phi-phi TP flash whose 1-vs-2 phase decision is a stability test."""
    # Imported here, not at module scope: chemthermo.stability._evaluator imports
    # chemthermo.flash._common, so a module-level import would make the two
    # packages' import order significant.
    from ..stability import stability_tp
    from ..stability.tp import _ln_phi_min_gibbs

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        eos=eos,
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": mode,
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }
    if stability.feed_branch is not None:
        base["feed_branch"] = stability.feed_branch

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one phase or two. Loosen "
            "FlashSettings.stability_settings, or pass "
            "FlashSettings(phase_detection='wilson-heuristic') to use the legacy "
            "K-bound heuristic instead."
        )

    if stability.status == "stable":
        phase_name = stability.feed_branch or "vapor"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=1.0 if phase_name == "vapor" else 0.0,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
                "phase_label_method": _phase_label_method(
                    eos, mixture, temperature, pressure, z, phase_name
                ),
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )

    trial_w = np.array(trial, dtype=float)
    k_seed, incipient_phase = _stability_k_seed(
        mixture,
        temperature,
        pressure,
        z=z,
        w=trial_w,
        tpd_min=float(stability.tpd_min),
    )

    seed_label = "stability"
    log_seed: np.ndarray | None = None
    vapor_fraction, _f0, _f1 = _rachford_rice(z, k_seed)
    if vapor_fraction is None and has_trace_component(trial_w, z > 0.0):
        # ADR-0024 decision 2. Both conditions hold together only in the
        # geometry this route exists for: a stationary point with a component
        # below `TRACE_MOLE_FRACTION` *and* a K-set that brackets no vapor
        # fraction at all, which is what an essentially pure polymer melt gives
        # (validation Case P-14). Either condition alone would divert states
        # that converge today - a trace component with a workable bracket is
        # ordinary, and a collapsed bracket on an ordinary stationary point is
        # what the Wilson fallback below is for - so the split is only rewritten
        # where both are true and successive substitution has nowhere to start.
        seed_label = "stability-log"
        log_seed = log_space_seed(
            z=z,
            w=trial_w,
            tpd_min=float(stability.tpd_min),
            incipient_vapor=incipient_phase == _VAPOR,
        )
    elif vapor_fraction is None:
        # Documented fallback: the stationary point is a valid starting phase but
        # its K-values need not bracket a Rachford-Rice root in every geometry.
        seed_label = "wilson"
        k_seed = wilson_k(mixture, temperature, pressure)
        vapor_fraction, _f0, _f1 = _rachford_rice(z, k_seed)
        if vapor_fraction is None:
            raise ConvergenceError(
                "Feed is unstable (tpd_min="
                f"{stability.tpd_min:.6e}) but neither the stability-seeded nor the "
                "Wilson K-values bracket a Rachford-Rice root."
            )

    roots_x, roots_y = _phi_phi_roots(
        eos,
        mixture,
        temperature,
        pressure,
        feed_branch=stability.feed_branch,
        incipient_branch=stability.phase_branch,
        incipient_phase=incipient_phase,
        seed_label=seed_label,
    )

    stage_diagnostics: dict[str, float | int | str | bool] = {}
    if log_seed is None:
        # Either the stability seed bracketed a root or the Wilson fallback did;
        # the case where neither does has already raised above.
        assert vapor_fraction is not None
        try:
            split = _solve_k_loop(
                mixture,
                temperature,
                pressure,
                eos=eos,
                activity_model=None,
                mode=mode,
                settings=settings,
                z=z,
                K=k_seed,
                vapor_fraction=vapor_fraction,
                # The full pre-ADR-0016 budget is spent on successive substitution
                # before the second-order stage is allowed to touch anything, so every
                # state that converged before this slice converges identically now.
                # ``FlashSettings.ssi_iterations`` is deliberately *not* consulted here;
                # see ADR-0016.
                max_iter=settings.max_iter,
                allow_unconverged=settings.second_order,
                extended_rachford_rice=True,
                roots_x=roots_x,
                roots_y=roots_y,
            )
        except ModelError:
            if not settings.second_order:
                raise
            # ADR-0024 decision 2, second trigger: successive substitution could
            # not even form a K-value - the two phases' ``ln phi`` differ by more
            # than the exponential's range, which is the Mw = 53 000 chain below
            # the solvent's saturation pressure. There is no iterate to hand
            # over, so the log-space stage starts from the stationary point
            # instead. Only a state that raised before this slice reaches here.
            log_seed = log_space_seed(
                z=z,
                w=trial_w,
                tpd_min=float(stability.tpd_min),
                incipient_vapor=incipient_phase == _VAPOR,
            )
            seed_label = "stability-log"

    if log_seed is not None:
        split, stage_diagnostics = _phi_phi_log_space(
            mixture,
            settings=settings,
            z=z,
            u0=log_seed,
            roots_x=roots_x,
            roots_y=roots_y,
            seed="stability-w",
        )
        x, y, beta = split.x, split.y, split.vapor_fraction
        ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    else:
        x, y, beta = split.x, split.y, split.vapor_fraction
        ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y

        if not split.converged or not 0.0 < beta < 1.0:
            # Successive substitution ran out of budget (ADR-0016), or it
            # converged on the *trivial* solution, whose vapor fraction is
            # outside [0, 1] and which the check below refuses. The second case
            # is ADR-0024 decision 1: before this slice it raised here, so
            # handing it to the second-order stage cannot move a number that
            # was ever returned. The extra keys are added only on this branch: a
            # converged first stage must produce the diagnostics mapping it
            # produced before ADR-0016, bit for bit.
            x, y, beta, ln_f_x, ln_f_y, stage_diagnostics = _phi_phi_second_order(
                mixture,
                settings=settings,
                z=z,
                split=split,
                roots_x=roots_x,
                roots_y=roots_y,
            )

    if not 0.0 < beta < 1.0:
        raise ConvergenceError(
            "The phi-phi split converged to a vapor fraction outside (0, 1) "
            f"(beta={beta:.6e}), i.e. to a single phase, while the tangent-plane "
            f"stability test reports the feed unstable (tpd_min={stability.tpd_min:.6e}). "
            "The extended (negative-flash) Rachford-Rice window admits such a root "
            "during iteration, but it cannot be an answer: it contradicts the "
            "stability verdict. This is a solver failure, not a single-phase state."
        )

    ln_phi_feed, _feed_branch = _ln_phi_min_gibbs(
        eos, mixture=mixture, temperature=temperature, pressure=pressure, composition=z
    )
    checks = _verify_split(
        z=z,
        x=x,
        y=y,
        beta=beta,
        ln_f_x=ln_f_x,
        ln_f_y=ln_f_y,
        ln_f_feed=ln_phi_feed,
        residual_key="fugacity_residual",
    )

    # ADR-0017 / ADR-0019: each phase converged on its *own* density root
    # (`_split._PhaseRoot`), so the pair of branches is read back from the two
    # holders rather than assumed to be ("liquid", "vapor"); only which
    # *name*, which `vapor_fraction` value and which `phase_regime` go with the
    # already-converged pair are decided here. `x`, `y` and `beta` themselves
    # are untouched.
    branches = (roots_x.selected or _LIQUID, roots_y.selected or _VAPOR)
    naming = _name_two_phase_result(eos, mixture, temperature, pressure, x, y, beta, branches)
    names = naming.names

    # Conditional keys, on the same principle as ADR-0016's stage keys: a split
    # whose two phases converged on the historical liquid/vapour pair of roots
    # carries the diagnostics mapping it carried before ADR-0019, bit for bit,
    # and the keys appear exactly when the per-phase selection put the phases
    # somewhere the pre-slice code could not express. Absent means
    # ``("liquid", "vapor")``. What is reported is the ADR-0017 *identity* of
    # each converged root, the same quantity `feed_branch` / `phase_branch`
    # report for the stability test - not the raw candidate label, which on a
    # single-real-root state is a tie-break between two names for one root and
    # would make these keys fire on states where nothing moved.
    branch_diagnostics: dict[str, float | int | str | bool] = {}
    if naming.identities != (_LIQUID, _VAPOR):
        branch_diagnostics = {
            "phase_i_branch": naming.identities[0],
            "phase_ii_branch": naming.identities[1],
        }

    report = _post_split_report(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=None,
        phases=((names[0], x), (names[1], y)),
        settings=settings,
    )

    if report.status != "stable" and settings.post_split_stability:
        # ADR-0020: the phi-phi path joins the ADR-0011 phase addition/removal
        # search. Every state of the bit-identity fixture and of the PR /
        # PC-SAFT grids is post-split *stable*, so this branch is not entered
        # there and those results are unchanged, diagnostics included.
        if settings.max_phases < 3:
            detail = ", ".join(
                [failure.phase_name for failure in report.instabilities] or report.inconclusive
            )
            raise ConvergenceError(
                "The converged two-phase solution is not a stable phase set: the "
                f"post-split stability test reports '{report.status}' for phase(s) "
                f"{detail} (most negative post-split tpd = {report.tpd_min:.6e}). A third "
                "phase is required, and FlashSettings.max_phases = "
                f"{settings.max_phases} forbids it. Raise max_phases (the default 3 "
                "resolves this state), or pass FlashSettings(post_split_stability=False) "
                "to receive the two-phase result anyway, with this failure recorded in "
                "diagnostics."
            )
        if report.status == "inconclusive":
            raise ConvergenceError(
                "A post-split stability test was inconclusive for phase(s) "
                f"{', '.join(report.inconclusive)}, so flash_tp cannot decide whether the "
                "two-phase set is the answer."
            )
        failure = report.instabilities[0]
        third_branch = failure.branch or _LIQUID
        phase_set = _EosPhaseSet(eos, mixture, temperature, pressure)
        history = [
            _phase_set_label((stability.feed_branch or _LIQUID,)),
            _phase_set_label(branches),
            _phase_set_label((*branches, third_branch)),
        ]
        two_phase_g_rt = (1.0 - float(beta)) * _reduced_g(x, ln_f_x) + float(beta) * _reduced_g(
            y, ln_f_y
        )
        return _flash_tp_phase_addition(
            mixture,
            temperature,
            pressure,
            z=z,
            model=phase_set,
            labels=(branches[0], branches[1], third_branch),
            # The two converged phases keep the very holders the split used,
            # so the multiphase solve continues on the roots they are already
            # on (ADR-0019); the third gets a fresh one pinned to the branch
            # the stability test found it on.
            surfaces=(
                roots_x.ln_fugacity_terms,
                roots_y.ln_fugacity_terms,
                phase_set.surface(third_branch),
            ),
            compositions=(x, y, failure.composition),
            history=history,
            ln_f_feed=ln_phi_feed,
            two_phase_g_rt=two_phase_g_rt,
            settings=settings,
            base={**base, "incipient_phase": incipient_phase},
            additions=1,
        )

    post_split = report.diagnostics

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        x,
        y,
        beta,
        names=names,
        vapor_fraction=naming.vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations
            + int(stage_diagnostics.get("second_order_iterations", 0)),
            "converged": True,
            "termination_reason": "tolerance_met",
            "max_delta_k": split.max_delta,
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": naming.regime,
            "k_seed": seed_label,
            "incipient_phase": incipient_phase,
            "phase_label_method": naming.method,
            **branch_diagnostics,
            **stage_diagnostics,
            **checks,
            **post_split,
        },
    )


def _log_space_diagnostics(
    mixture: Mixture,
    *,
    z: np.ndarray,
    ln_x_ii: np.ndarray,
    iterations: int,
    residual: float,
    seed: str,
) -> dict[str, float | int | str | bool]:
    """The ``log_space_*`` keys (ADR-0024 decision 3).

    A phase composition that came out of the log-space stage may carry an exact
    ``0.0`` where the model's ``ln phi`` difference exceeds the exponential's
    range. The number itself is not lost - it is
    ``exp(log_space_ln_x_min)`` for the component named in
    ``log_space_ln_x_min_component`` - so the smallest log mole fraction of
    phase II is reported here, whether or not it underflowed.

    ``log_space_residual`` is reported for the same reason.
    ``fugacity_residual`` is formed by
    :func:`chemthermo.flash._verify._equilibrium_residual` over the components
    present in **both** phases, so a component whose mole fraction underflowed
    to ``0.0`` in one of them drops out of it; the stage's own residual is the
    same quantity taken in log space over every component present in the feed,
    and it is the number that says the polymer's equal-fugacity condition is
    satisfied too.
    """
    active = z > 0.0
    masked = np.where(active, ln_x_ii, math.inf)
    position = int(np.argmin(masked))
    with np.errstate(under="ignore"):
        zeros = int(np.count_nonzero(active & (np.exp(masked) == 0.0)))
    return {
        "log_space_seed": seed,
        "log_space_iterations": iterations,
        "log_space_residual": float(residual),
        "log_space_ln_x_min": float(masked[position]),
        "log_space_ln_x_min_component": mixture.components[position].name,
        "log_space_zero_fractions": zeros,
    }


def _phi_phi_log_space(
    mixture: Mixture,
    *,
    settings: FlashSettings,
    z: np.ndarray,
    u0: np.ndarray,
    roots_x: _PhaseRoot,
    roots_y: _PhaseRoot,
    seed: str,
) -> tuple[_SplitSolution, dict[str, float | int | str | bool]]:
    """Solve a phi-phi split entirely in log mole numbers (ADR-0024).

    The successive-substitution loop is skipped: this route is taken only where
    it has nowhere to start (see :func:`_flash_tp_tangent_plane`), so the split
    is the log-space Newton stage alone, seeded directly in ``u``. The result is
    packed into the same :class:`chemthermo.flash._split._SplitSolution` the
    K-loop returns, so everything downstream - verification, naming, post-split
    stability, assembly - is the code that was already there.

    ``K = x^II / x^I`` is formed from the *logarithms* and may underflow to
    ``0.0`` for a component whose mole fraction is not a double; ``k_min`` is
    then an honest zero and the magnitude is in ``log_space_ln_x_min``.
    """
    refined = log_space_split(
        z=z,
        u0=u0,
        terms_i=roots_x.ln_fugacity_terms,
        terms_ii=roots_y.ln_fugacity_terms,
        settings=settings,
    )
    if refined.residual > settings.tol:
        raise ConvergenceError(
            "flash_tp did not converge the phi-phi split in log mole numbers; "
            f"equal-fugacity residual={refined.residual:.3e} after {refined.iterations} "
            f"log-space Newton iterations from the {seed} seed."
        )

    active = z > 0.0
    with np.errstate(over="ignore", under="ignore"):
        ln_k = np.where(active, refined.ln_x_ii - np.log(np.where(active, refined.x_i, 1.0)), 0.0)
        k_values = np.exp(ln_k)
    split = _SplitSolution(
        x=refined.x_i,
        y=refined.x_ii,
        vapor_fraction=refined.beta,
        K=k_values,
        ln_f_x=refined.ln_f_i,
        ln_f_y=refined.ln_f_ii,
        iterations=0,
        max_delta=refined.max_delta_k,
    )
    diagnostics: dict[str, float | int | str | bool] = {
        "ssi_iterations": 0,
        "second_order_iterations": refined.iterations,
        "converged_stage": "second-order-log",
        "negative_flash_steps": 0,
        **_log_space_diagnostics(
            mixture,
            z=z,
            ln_x_ii=refined.ln_x_ii,
            iterations=refined.iterations,
            residual=refined.residual,
            seed=seed,
        ),
    }
    return split, diagnostics


def _phi_phi_second_order(
    mixture: Mixture,
    *,
    settings: FlashSettings,
    z: np.ndarray,
    split: _SplitSolution,
    roots_x: _PhaseRoot,
    roots_y: _PhaseRoot,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
    np.ndarray,
    np.ndarray,
    dict[str, float | int | str | bool],
]:
    """Finish an unconverged phi-phi split with the ADR-0009 Newton stage.

    The same Gibbs-energy minimization the liquid-liquid and modified-Raoult
    splits use, given the EOS tangent-plane terms: ``ln phi`` on each phase's
    **own** density root, from the same two
    :class:`chemthermo.flash._split._PhaseRoot` holders the
    successive-substitution loop used (ADR-0019), so the stage continues on the
    branches the first stage was on rather than on a fixed liquid/vapour pair.
    The derivation in
    :func:`chemthermo.flash._second_order._second_order_split` is unchanged -
    equation (2) there holds phase by phase, and holds for ``ln phi`` for the
    same reason it holds for ``ln gamma``: the Gibbs-Duhem relation at fixed
    ``T, P``.

    Returns:
        ``(x, y, beta, ln_f_x, ln_f_y, diagnostics)`` for the better of the two
        stages.

    Raises:
        ConvergenceError: If neither stage reached ``settings.tol`` on the
            equal-fugacity residual.
    """
    terms_x = roots_x.ln_fugacity_terms
    terms_y = roots_y.ln_fugacity_terms

    x, y, beta = split.x, split.y, split.vapor_fraction
    ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
    converged_stage = "successive-substitution"
    second_order_iterations = 0

    if settings.second_order and residual > settings.second_order_tol:
        # The stage parameterizes the split by phase-II mole numbers
        # ``n = beta y``, which must satisfy ``0 < n_i < z_i``. A negative-flash
        # iterate does not, so the starting vapor fraction is pulled back into
        # the physical range first; only the *starting point* moves, and the
        # stage is a descent method from wherever it starts.
        seed_beta = beta
        if not 0.0 < seed_beta < 1.0 or np.any(seed_beta * y >= z):
            seed_beta = float(np.min(np.where(y > 0.0, z / np.maximum(y, 1e-300), 1.0))) * 0.5
            seed_beta = min(max(seed_beta, 1e-8), 1.0 - 1e-8)
        refined = _second_order_split(
            z=z,
            x_ii=y,
            beta=seed_beta,
            terms_i=terms_x,
            terms_ii=terms_y,
            settings=settings,
        )
        second_order_iterations = refined.iterations
        if refined.residual < residual:
            x, y, beta = refined.x_i, refined.x_ii, refined.beta
            ln_f_x, ln_f_y = terms_x(x), terms_y(y)
            residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
            converged_stage = "second-order"

    log_space_diagnostics: dict[str, float | int | str | bool] = {}
    if settings.second_order and residual > settings.tol:
        # ADR-0024 decision 1, second entry point: the linear stage failed, so
        # the same minimization is retried in log mole numbers from the best
        # iterate it produced. Only a state that would otherwise raise below
        # reaches this, so no converged number moves.
        try:
            log_refined = log_space_split(
                z=z,
                u0=seed_from_iterate(z=z, x_ii=y, beta=beta),
                terms_i=terms_x,
                terms_ii=terms_y,
                settings=settings,
            )
        except ModelError:
            log_refined = None
        if log_refined is not None and log_refined.residual < residual:
            x, y, beta = log_refined.x_i, log_refined.x_ii, log_refined.beta
            ln_f_x, ln_f_y = log_refined.ln_f_i, log_refined.ln_f_ii
            residual = log_refined.residual
            converged_stage = "second-order-log"
            log_space_diagnostics = _log_space_diagnostics(
                mixture,
                z=z,
                ln_x_ii=log_refined.ln_x_ii,
                iterations=log_refined.iterations,
                residual=log_refined.residual,
                seed="linear-iterate",
            )
            second_order_iterations += log_refined.iterations

    if residual > settings.tol:
        raise ConvergenceError(
            "flash_tp did not converge the phi-phi split; equal-fugacity residual="
            f"{residual:.3e} after {split.iterations} successive-substitution "
            f"(max_delta_k={split.max_delta:.3e}) and {second_order_iterations} "
            "second-order iterations."
        )

    diagnostics: dict[str, float | int | str | bool] = {
        "ssi_iterations": split.iterations,
        "second_order_iterations": second_order_iterations,
        "converged_stage": converged_stage,
        "negative_flash_steps": split.negative_flash_steps,
        **log_space_diagnostics,
    }
    return x, y, float(beta), ln_f_x, ln_f_y, diagnostics


def _flash_tp_liquid_liquid(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    activity_model: ActivityModel,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Liquid-liquid TP flash driven entirely by an activity model (ADR-0009).

    The phase count is an *output*: the feed is tested with
    ``stability_tp(..., activity_model=...)`` and a split is attempted only when
    that test finds a negative tangent-plane distance. The split itself is the
    shared Rachford-Rice / successive-substitution loop with
    ``K_i = gamma_i^I / gamma_i^II``, followed by the second-order stage of
    :func:`chemthermo.flash._second_order._second_order_split` when the
    equal-activity residual is still above ``settings.second_order_tol``.
    """
    from ..stability import stability_tp

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        activity_model=activity_model,
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": "gamma-gamma",
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one liquid phase or two. Loosen "
            "FlashSettings.stability_settings and try again."
        )

    if stability.status == "stable":
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name="liquid",
            vapor_fraction=None,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": "liquid",
                "phase_regime": "single-phase",
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )

    # Same seed as the phi-phi path: Michelsen's unnormalized mole numbers
    # W = w exp(-tpd), so that f_RR(0) = sum_i W_i - 1 > 0 is bracketable.
    # Phase I is feed-like, phase II incipient-like, K_i = x_i^II / x_i^I.
    w = np.array(trial, dtype=float)
    sum_capital_w = math.exp(-float(stability.tpd_min))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (w * sum_capital_w) / z
    usable = np.isfinite(ratio) & (ratio > 0.0)
    k_seed = np.where(usable, ratio, _INERT_SEED_K)

    beta, _f0, _f1 = _rachford_rice(z, k_seed)
    if beta is None:
        raise ConvergenceError(
            f"Feed is unstable (tpd_min={stability.tpd_min:.6e}) but the stability-seeded "
            "K-values do not bracket a Rachford-Rice root, so no liquid-liquid split can "
            "be started."
        )

    ln_gamma = _ln_gamma_function(activity_model, mixture, temperature)

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        mode="gamma-gamma",
        settings=settings,
        z=z,
        K=k_seed,
        vapor_fraction=beta,
        max_iter=(
            min(settings.ssi_iterations, settings.max_iter)
            if settings.second_order
            else settings.max_iter
        ),
        allow_unconverged=settings.second_order,
    )

    x_i, x_ii, beta = split.x, split.y, split.vapor_fraction
    ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    residual = _equilibrium_residual(x_i, x_ii, ln_f_x, ln_f_y)
    second_order_iterations = 0
    converged_stage = "successive-substitution" if split.converged else None

    if settings.second_order and residual > settings.second_order_tol:
        refined = _second_order_split(
            z=z, x_ii=x_ii, beta=beta, terms_i=ln_gamma, settings=settings
        )
        second_order_iterations = refined.iterations
        if refined.residual < residual:
            x_i, x_ii, beta = refined.x_i, refined.x_ii, refined.beta
            ln_f_x, ln_f_y = ln_gamma(x_i), ln_gamma(x_ii)
            residual = _equilibrium_residual(x_i, x_ii, ln_f_x, ln_f_y)
            converged_stage = "second-order"

    if residual > settings.tol and not split.converged:
        raise ConvergenceError(
            "flash_tp did not converge the liquid-liquid split; equal-activity residual="
            f"{residual:.3e} after {split.iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations."
        )

    checks = _verify_split(
        z=z,
        x=x_i,
        y=x_ii,
        beta=beta,
        ln_f_x=ln_f_x,
        ln_f_y=ln_f_y,
        ln_f_feed=ln_gamma(z / float(np.sum(z))),
        residual_key="equilibrium_residual",
    )

    post_split = _post_split_stability(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        phases=((_LIQUID_I, x_i), (_LIQUID_II, x_ii)),
        settings=settings,
    )

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        x_i,
        x_ii,
        beta,
        names=(_LIQUID_I, _LIQUID_II),
        vapor_fraction=None,
        diagnostics={
            **base,
            "iterations": split.iterations + second_order_iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": "LLE",
            "k_seed": "stability",
            "ssi_iterations": split.iterations,
            "second_order_iterations": second_order_iterations,
            "converged_stage": converged_stage or "successive-substitution",
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            **checks,
            **post_split,
        },
    )


def _flash_tp_modified_raoult(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    activity_model: ActivityModel,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Low-pressure gamma-phi TP flash from one tangent plane (ADR-0010).

    The feed is tested against **two** phase candidates at once - an
    activity-coefficient liquid with Antoine pure-liquid reference fugacities
    and an ideal-gas vapor - so the same call returns a vapor-liquid split, a
    liquid-liquid split, or a single phase, and the *candidate label* of the
    stationary point is what decides which:

    ========================  ========================  ====================
    feed candidate            incipient candidate       result
    ========================  ========================  ====================
    liquid                    vapor                     VLE, bubble side
    vapor                     liquid                    VLE, dew side
    liquid                    liquid                    LLE
    ========================  ========================  ====================

    The split loop then evaluates each phase with the candidate it was assigned
    (:func:`chemthermo.flash._split._solve_k_loop`), which makes
    ``K_i = gamma_i Psat_i / P`` for a vapor-liquid pair and
    ``K_i = gamma_i^I / gamma_i^II`` for a liquid-liquid pair without the loop
    knowing the difference. Every converged phase is then re-tested against
    **both** candidates, so a state that needs all three phases raises rather
    than being returned (see :func:`chemthermo.flash._verify._post_split_stability`).
    """
    from ..stability import stability_tp
    from ..stability._evaluator import modified_raoult_candidates

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        activity_model=activity_model,
        vapor="ideal",
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": _MODIFIED_RAOULT,
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }
    for key in ("antoine_valid_Tmin_K", "antoine_valid_Tmax_K"):
        if key in stability.diagnostics:
            base[key] = stability.diagnostics[key]
    if stability.feed_branch is not None:
        base["feed_branch"] = stability.feed_branch

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one phase or two. Loosen "
            "FlashSettings.stability_settings and try again."
        )

    feed_label = stability.feed_branch or _LIQUID
    if stability.status == "stable":
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=feed_label,
            vapor_fraction=1.0 if feed_label == _VAPOR else 0.0,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": feed_label,
                "phase_regime": "single-phase",
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )
    incipient_label = stability.phase_branch or _LIQUID
    if feed_label == _VAPOR and incipient_label == _VAPOR:
        raise ConvergenceError(
            "The tangent-plane minimizer is a second ideal-gas vapor, which cannot "
            "coexist with the first: an ideal gas has no composition range over which "
            "it demixes. This indicates a model or numerical failure, not a phase split."
        )

    liquid, vapor = modified_raoult_candidates(
        activity_model, mixture=mixture, temperature=temperature, pressure=pressure
    )
    candidates = {_LIQUID: liquid, _VAPOR: vapor}
    terms_x = candidates[feed_label].ln_fugacity_terms
    terms_y = candidates[incipient_label].ln_fugacity_terms

    # Same seed as the other tangent-plane paths: Michelsen's unnormalized mole
    # numbers W = w exp(-tpd), with phase x feed-like and phase y incipient, so
    # K_i = y_i / x_i = W_i / z_i and f_RR(0) = sum_i W_i - 1 > 0 is bracketable.
    w = np.array(trial, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (w * math.exp(-float(stability.tpd_min))) / z
    usable = np.isfinite(ratio) & (ratio > 0.0)
    k_seed = np.where(usable, ratio, _INERT_SEED_K)

    beta, _f0, _f1 = _rachford_rice(z, k_seed)
    if beta is None:
        raise ConvergenceError(
            f"Feed is unstable (tpd_min={stability.tpd_min:.6e}) but the stability-seeded "
            "K-values do not bracket a Rachford-Rice root, so no split can be started."
        )

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        mode=_MODIFIED_RAOULT,
        settings=settings,
        z=z,
        K=k_seed,
        vapor_fraction=beta,
        max_iter=(
            min(settings.ssi_iterations, settings.max_iter)
            if settings.second_order
            else settings.max_iter
        ),
        allow_unconverged=settings.second_order,
        terms_x=terms_x,
        terms_y=terms_y,
    )

    x, y, beta = split.x, split.y, split.vapor_fraction
    ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
    second_order_iterations = 0
    converged_stage = "successive-substitution" if split.converged else None

    if settings.second_order and residual > settings.second_order_tol:
        refined = _second_order_split(
            z=z,
            x_ii=y,
            beta=beta,
            terms_i=terms_x,
            terms_ii=terms_y,
            settings=settings,
        )
        second_order_iterations = refined.iterations
        if refined.residual < residual:
            x, y, beta = refined.x_i, refined.x_ii, refined.beta
            ln_f_x, ln_f_y = terms_x(x), terms_y(y)
            residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
            converged_stage = "second-order"

    if residual > settings.tol and not split.converged:
        raise ConvergenceError(
            "flash_tp did not converge the modified-Raoult split; equilibrium residual="
            f"{residual:.3e} after {split.iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations."
        )

    checks = _verify_split(
        z=z,
        x=x,
        y=y,
        beta=beta,
        ln_f_x=ln_f_x,
        ln_f_y=ln_f_y,
        ln_f_feed=terms_x(z / float(np.sum(z))),
        residual_key="equilibrium_residual",
    )

    if feed_label == _LIQUID and incipient_label == _LIQUID:
        names = (_LIQUID_I, _LIQUID_II)
        regime = "LLE"
        vapor_fraction: float | None = None
    else:
        names = (feed_label, incipient_label)
        regime = "VLE"
        vapor_fraction = float(beta) if incipient_label == _VAPOR else 1.0 - float(beta)

    report = _post_split_report(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        phases=((names[0], x), (names[1], y)),
        settings=settings,
        vapor="ideal",
    )

    if report.status != "stable" and settings.post_split_stability:
        if settings.max_phases < 3:
            detail = ", ".join(
                [failure.phase_name for failure in report.instabilities] or report.inconclusive
            )
            raise ConvergenceError(
                "The converged two-phase solution is not a stable phase set: the "
                f"post-split stability test reports '{report.status}' for phase(s) "
                f"{detail} (most negative post-split tpd = {report.tpd_min:.6e}). A third "
                "phase is required, and FlashSettings.max_phases = "
                f"{settings.max_phases} forbids it. Raise max_phases (the default 3 "
                "resolves this state), or pass FlashSettings(post_split_stability=False) "
                "to receive the two-phase result anyway, with this failure recorded in "
                "diagnostics."
            )
        if report.status == "inconclusive":
            raise ConvergenceError(
                "A post-split stability test was inconclusive for phase(s) "
                f"{', '.join(report.inconclusive)}, so flash_tp cannot decide whether the "
                "two-phase set is the answer."
            )
        # The phase set is provably not the answer, and the stability minimizer
        # found on the failing phase is the incipient third phase. Hand over to
        # the phase addition / removal search (ADR-0011).
        failure = report.instabilities[0]
        history = [_phase_set_label((feed_label,)), _phase_set_label(names)]
        history.append(_phase_set_label((*names, failure.branch or _LIQUID)))
        two_phase_g_rt = (1.0 - float(beta)) * _reduced_g(x, ln_f_x) + float(beta) * _reduced_g(
            y, ln_f_y
        )
        surfaces = {
            _LIQUID: candidates[_LIQUID].ln_fugacity_terms,
            _VAPOR: candidates[_VAPOR].ln_fugacity_terms,
        }
        third_label = failure.branch or _LIQUID
        phase_set = _ActivityPhaseSet(activity_model, candidates=surfaces)
        return _flash_tp_phase_addition(
            mixture,
            temperature,
            pressure,
            z=z,
            model=phase_set,
            labels=(feed_label, incipient_label, third_label),
            surfaces=(
                surfaces[feed_label],
                surfaces[incipient_label],
                surfaces[third_label],
            ),
            compositions=(x, y, failure.composition),
            history=history,
            ln_f_feed=terms_x(z / float(np.sum(z))),
            two_phase_g_rt=two_phase_g_rt,
            settings=settings,
            base=base,
            additions=1,
        )

    post_split = report.diagnostics

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        x,
        y,
        beta,
        names=names,
        vapor_fraction=vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations + second_order_iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": regime,
            "k_seed": "stability",
            "incipient_phase": incipient_label,
            "ssi_iterations": split.iterations,
            "second_order_iterations": second_order_iterations,
            "converged_stage": converged_stage or "successive-substitution",
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            **checks,
            **post_split,
        },
    )


def _stability_k_seed(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    w: np.ndarray,
    tpd_min: float,
) -> tuple[np.ndarray, str]:
    """Initial K-values from the tangent-plane minimizer.

    The stability test returns the normalized incipient composition ``w`` and
    the reduced tangent-plane distance at that stationary point. Michelsen's
    unnormalized mole numbers follow from the stationary-point identity
    ``tpd = -ln(sum_i W_i)``, so ``W = w * exp(-tpd)``. Using ``W`` rather than
    ``w`` matters: with ``K_i = w_i / z_i`` the Rachford-Rice function at
    ``beta = 0`` is ``sum_i w_i - 1 = 0`` exactly, a degenerate root that cannot
    be bracketed, whereas ``K_i = W_i / z_i`` gives ``sum_i W_i - 1 > 0`` for an
    unstable feed. Both seeds have the same fixed point; only the bracket
    differs.

    The returned K is ``y / x``, so it is ``W / z`` when the incipient phase is
    the vapor-like one and ``z / W`` when it is the liquid-like one. Which of
    the two converged phases is *named* "vapor" is a labelling convention -
    ``EquationOfState`` exposes no molar volume - and is resolved by volatility
    ordering: see :func:`_incipient_is_vapor_like`.

    Returns:
        ``(K, incipient_phase)`` with ``incipient_phase`` in
        ``{"vapor", "liquid"}``.
    """
    sum_capital_w = math.exp(-tpd_min) if math.isfinite(tpd_min) else 1.0
    capital_w = w * sum_capital_w

    incipient_vapor = _incipient_is_vapor_like(mixture, temperature, pressure, z=z, w=w)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = capital_w / z if incipient_vapor else z / capital_w

    usable = np.isfinite(ratio) & (ratio > 0.0)
    K = np.where(usable, ratio, _INERT_SEED_K)
    return K, "vapor" if incipient_vapor else "liquid"


def _incipient_is_vapor_like(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    w: np.ndarray,
) -> bool:
    """Name the incipient phase "vapor" or "liquid" by volatility ordering.

    The tangent-plane test finds *a* second phase; it does not say which of the
    two converged phases should be called the vapor. The successive-substitution
    loop is symmetric under swapping the two phases (``x <-> y``,
    ``beta <-> 1 - beta``, ``K <-> 1/K``), so an arbitrary orientation would
    return the mirror-labelled solution roughly half the time.

    The label is decided here by the only volatility ordering the package has:
    the Wilson correlation's ranking of the components (built from ``Tc``,
    ``Pc`` and ``omega``). Writing ``hi`` for the component with the largest
    Wilson K among those present in the feed and ``lo`` for the smallest, the
    incipient phase is called vapor-like when it is enriched in ``hi`` relative
    to ``lo`` compared with the feed, i.e. when
    ``ln(w_hi / z_hi) - ln(w_lo / z_lo) >= 0``.

    Only the *ranking* is used, never the magnitudes, and it decides the name
    only - never the one-versus-two-phase verdict, the compositions, or the
    vapor fraction, all of which come from the tangent-plane test and the
    converged split. The minimum-Gibbs branch label reported by the stability
    test is deliberately *not* used for this: whenever the cubic has a single
    real root (common for dense and near-critical states) both branch calls
    return identical fugacity coefficients and the label is only a tie-break.

    There is no analogue for two liquid phases, which is why the liquid-liquid
    path names its phases by role instead (see
    :func:`chemthermo.flash.tp.flash_tp`).
    """
    active = z > 0.0
    if int(np.count_nonzero(active)) < 2:
        return True

    k_wilson = wilson_k(mixture, temperature, pressure)
    masked = np.where(active, k_wilson, np.nan)
    hi = int(np.nanargmax(masked))
    lo = int(np.nanargmin(masked))
    if hi == lo:
        return True

    tiny = np.finfo(float).tiny
    score = math.log(max(float(w[hi]), tiny) / float(z[hi])) - math.log(
        max(float(w[lo]), tiny) / float(z[lo])
    )
    return score >= 0.0
