"""A phase split must be verified, not just converged.

Mass balance, the equal-fugacity/equal-activity residual, the Gibbs-energy
reduction against the single-phase feed (ADR-0008 decision 4), and the
post-split stability re-test of every converged phase (ADR-0009 decision 4,
generalized to any number of phases in ADR-0011).
"""

from __future__ import annotations

import math
from typing import Literal, Sequence

import numpy as np

from ..core import Composition, Mixture
from ..exceptions import ConvergenceError
from ..models import ActivityModel, EquationOfState
from ..validation import COMPOSITION_SUM_TOL
from .settings import FlashSettings


class _PhaseInstability:
    """A converged phase that a post-split stability test found unstable.

    This is what phase *addition* consumes: ``composition`` is the normalized
    tangent-plane minimizer found from that phase, ``tpd`` its reduced
    tangent-plane distance, and ``branch`` the phase-candidate label of the
    incipient phase (``"liquid"`` / ``"vapor"`` for the modified-Raoult pair,
    None for a single-candidate evaluator).
    """

    __slots__ = ("branch", "composition", "phase_name", "tpd")

    def __init__(
        self,
        *,
        phase_name: str,
        composition: np.ndarray,
        tpd: float,
        branch: str | None,
    ) -> None:
        self.phase_name = phase_name
        self.composition = composition
        self.tpd = tpd
        self.branch = branch


class _PostSplitReport:
    """Outcome of re-testing every converged phase for stability.

    Attributes:
        diagnostics: The ``post_split_*`` / ``phase_stability_*`` keys.
        status: ``"stable"``, ``"unstable"`` or ``"inconclusive"``.
        tpd_min: Most negative post-split tangent-plane distance found.
        instabilities: One entry per genuinely unstable phase, deepest first.
        inconclusive: Names of phases whose stability test did not converge.
    """

    __slots__ = ("diagnostics", "inconclusive", "instabilities", "status", "tpd_min")

    def __init__(
        self,
        *,
        diagnostics: dict[str, float | int | str | bool],
        status: str,
        tpd_min: float,
        instabilities: list[_PhaseInstability],
        inconclusive: list[str],
    ) -> None:
        self.diagnostics = diagnostics
        self.status = status
        self.tpd_min = tpd_min
        self.instabilities = instabilities
        self.inconclusive = inconclusive


def _equilibrium_residual(
    x: np.ndarray, y: np.ndarray, ln_f_x: np.ndarray, ln_f_y: np.ndarray
) -> float:
    """``max_i |ln(x_i f_i^x) - ln(y_i f_i^y)|`` over components in both phases."""
    both = (x > 0.0) & (y > 0.0)
    if not np.any(both):
        return float("nan")
    return float(
        np.max(np.abs((np.log(x[both]) + ln_f_x[both]) - (np.log(y[both]) + ln_f_y[both])))
    )


def _reduced_g(fractions: np.ndarray, ln_f: np.ndarray) -> float:
    """``sum_i x_i (ln x_i + t_i)``: the composition-dependent part of ``G/RT``."""
    mask = fractions > 0.0
    return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_f[mask])))


def _verify_split(
    *,
    z: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    beta: float,
    ln_f_x: np.ndarray,
    ln_f_y: np.ndarray,
    ln_f_feed: np.ndarray,
    residual_key: str,
) -> dict[str, float | int | str | bool]:
    """Residuals that a converged two-phase solution must satisfy.

    - ``mass_balance_residual``: ``max_i |z_i - (beta y_i + (1 - beta) x_i)|``.
    - ``residual_key``: ``max_i |ln(x_i f_i^x) - ln(y_i f_i^y)|`` with
      ``f = phi`` (phi-phi, key ``fugacity_residual``) or ``f = gamma``
      (gamma-gamma, key ``equilibrium_residual``). The pressure (phi-phi) or the
      pure-liquid reference (gamma-gamma) cancels, so this is the equal-fugacity
      condition in log form.
    - ``delta_g_split_rt``: ``beta sum_i y_i ln(y_i f_i^y)
      + (1 - beta) sum_i x_i ln(x_i f_i^x) - sum_i z_i ln(z_i f_i^feed)``,
      which must be negative for the split to be an improvement on the
      single-phase feed. For phi-phi the feed term uses the *minimum-Gibbs*
      branch (the stability module's convention, i.e. the lowest single-phase
      Gibbs energy available), and since ADR-0019 each split phase uses its own
      minimum-Gibbs root as well, so all three terms are on the same footing
      and the quantity is exact rather than the upper bound it was while the
      two phases were pinned to fixed branches. An activity model has one
      branch, so for gamma-gamma the quantity was always exact.
    """
    beta = float(beta)
    balance = float(np.max(np.abs(z - (beta * y + (1.0 - beta) * x))))
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)

    delta_g = (
        beta * _reduced_g(y, ln_f_y)
        + (1.0 - beta) * _reduced_g(x, ln_f_x)
        - _reduced_g(z, ln_f_feed)
    )

    return {
        "mass_balance_residual": balance,
        residual_key: residual,
        "delta_g_split_rt": delta_g,
    }


def _post_split_report(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    phases: Sequence[tuple[str, np.ndarray]],
    settings: FlashSettings,
    vapor: Literal["none", "ideal"] = "none",
) -> _PostSplitReport:
    """Test every converged phase for stability and adjudicate the phase set.

    Coexisting phases share one tangent plane, so a stability test run on any
    of them finds its *partners* as stationary points with ``tpd = 0``
    (validation Case S-3). Up to the split's own convergence tolerance that
    zero is a small negative number - measured worst case -7.0e-09 over the
    in-repo Peng-Robinson grid - so it can dip below ``-tpd_tol`` for a loosely
    converged split. A minimizer that *is* one of the other converged phases is
    therefore reported as ``"marginal"``, not as an instability. "Is a partner"
    uses the same measure as the stability module's trivial-solution test,
    ``sum_i ln(w_i / x_i^partner)^2 < trivial_tol``, rather than a new
    tolerance.

    Anything else negative means the phase set is not an answer: a further
    phase exists at a lower Gibbs energy, and the minimizer found here is the
    seed for it (:mod:`chemthermo.flash._multiphase`).

    ``vapor`` is forwarded to :func:`chemthermo.stability_tp`, so a
    modified-Raoult result has **every** phase re-tested against **both**
    candidates (liquid and ideal vapor). That is what turns a vapor-liquid
    answer whose liquid is inside a miscibility gap - or a liquid-liquid answer
    that should be boiling - into a three-phase state instead of a
    plausible-looking wrong one.

    Returns:
        A :class:`_PostSplitReport`. The diagnostics keys are
        ``post_split_checked``, ``post_split_stable``, ``post_split_status``,
        ``post_split_tpd_min``, ``phase_stability_<name>`` and
        ``phase_stability_tpd_min_<name>``.
    """
    from ..stability import StabilitySettings, stability_tp

    stability_settings = settings.stability_settings or StabilitySettings()
    diagnostics: dict[str, float | int | str | bool] = {"post_split_checked": True}
    worst = math.inf
    instabilities: list[_PhaseInstability] = []
    inconclusive: list[str] = []

    for position, (name, composition) in enumerate(phases):
        partners = [other for index, (_, other) in enumerate(phases) if index != position]
        phase_mixture = Mixture(
            components=mixture.components,
            composition=Composition(
                fractions=tuple(float(value) for value in composition),
                basis=mixture.basis,
                normalize=True,
                tol=COMPOSITION_SUM_TOL,
            ),
        )
        result = stability_tp(
            phase_mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            eos=eos,
            activity_model=activity_model,
            vapor=vapor,
            settings=stability_settings,
        )

        verdict = result.status
        if verdict == "unstable" and result.trial_composition is not None:
            minimizer = np.array(result.trial_composition, dtype=float)
            if any(
                _is_same_phase(minimizer, partner, stability_settings.trivial_tol)
                for partner in partners
            ):
                verdict = "marginal"

        diagnostics[f"phase_stability_{name}"] = verdict
        diagnostics[f"phase_stability_tpd_min_{name}"] = float(result.tpd_min)
        worst = min(worst, float(result.tpd_min))
        if verdict == "unstable":
            assert result.trial_composition is not None
            instabilities.append(
                _PhaseInstability(
                    phase_name=name,
                    composition=np.array(result.trial_composition, dtype=float),
                    tpd=float(result.tpd_min),
                    branch=result.phase_branch,
                )
            )
        elif verdict == "inconclusive":
            inconclusive.append(name)

    diagnostics["post_split_tpd_min"] = worst if math.isfinite(worst) else float("nan")
    status = "unstable" if instabilities else ("inconclusive" if inconclusive else "stable")
    diagnostics["post_split_status"] = status
    diagnostics["post_split_stable"] = status == "stable"

    instabilities.sort(key=lambda failure: failure.tpd)
    return _PostSplitReport(
        diagnostics=diagnostics,
        status=status,
        tpd_min=worst if math.isfinite(worst) else float("nan"),
        instabilities=instabilities,
        inconclusive=inconclusive,
    )


def _post_split_stability(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    phases: Sequence[tuple[str, np.ndarray]],
    settings: FlashSettings,
    vapor: Literal["none", "ideal"] = "none",
) -> dict[str, float | int | str | bool]:
    """:func:`_post_split_report` with the historical raise-on-failure behavior.

    Used by the paths that cannot add a phase: the phi-phi and gamma-gamma
    splits (ADR-0011 leaves both at two phases because no state in this
    repository exercises a third one there). The check always runs;
    ``settings.post_split_stability`` decides whether a *failure* raises or is
    only reported, so the diagnostics show the failure either way.
    """
    report = _post_split_report(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=activity_model,
        phases=phases,
        settings=settings,
        vapor=vapor,
    )

    if report.status != "stable" and settings.post_split_stability:
        detail = ", ".join(
            [failure.phase_name for failure in report.instabilities] or report.inconclusive
        )
        raise ConvergenceError(
            "The converged two-phase solution is not a stable phase set: the post-split "
            f"stability test reports '{report.status}' for phase(s) {detail} "
            f"(most negative post-split tpd = {report.tpd_min:.6e}). A third phase is "
            "required, and flash_tp returns at most two phases in this release "
            "(multiphase flash is the next slice). Pass "
            "FlashSettings(post_split_stability=False) to receive the two-phase result "
            "anyway, with this failure recorded in diagnostics."
        )

    return report.diagnostics


def _is_same_phase(w: np.ndarray, other: np.ndarray, trivial_tol: float) -> bool:
    """True when ``w`` is the composition ``other`` in the trivial-solution metric."""
    mask = (w > 0.0) & (other > 0.0)
    if not np.any(mask):
        return False
    if np.any((w > 0.0) != (other > 0.0)):
        return False
    ln_ratio = np.log(w[mask] / other[mask])
    return bool(float(np.sum(ln_ratio * ln_ratio)) < trivial_tol)
