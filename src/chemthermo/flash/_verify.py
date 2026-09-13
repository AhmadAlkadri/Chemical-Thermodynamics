"""A phase split must be verified, not just converged.

Mass balance, the equal-fugacity/equal-activity residual, the Gibbs-energy
reduction against the single-phase feed (ADR-0008 decision 4), and the
post-split stability re-test of every converged phase (ADR-0009 decision 4).
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..core import Composition, Mixture
from ..exceptions import ConvergenceError
from ..models import ActivityModel, EquationOfState
from ..validation import COMPOSITION_SUM_TOL
from .settings import FlashSettings


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
      Gibbs energy available), while each split phase uses the branch the solver
      actually converged on. Since the min-Gibbs branch is never higher in Gibbs
      energy, this combination is a conservative (upper-bound) estimate of the
      true reduction. An activity model has one branch, so for gamma-gamma the
      quantity is exact.
    """
    beta = float(beta)
    balance = float(np.max(np.abs(z - (beta * y + (1.0 - beta) * x))))
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)

    def _reduced_g(fractions: np.ndarray, ln_f: np.ndarray) -> float:
        mask = fractions > 0.0
        return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_f[mask])))

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


def _post_split_stability(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    phases: Sequence[tuple[str, np.ndarray]],
    settings: FlashSettings,
) -> dict[str, float | int | str | bool]:
    """Test every converged phase for stability and adjudicate the phase set.

    Two coexisting phases share one tangent plane, so a stability test run on
    either of them finds its *partner* as a stationary point with ``tpd = 0``
    (validation Case S-3). Up to the split's own convergence tolerance that
    zero is a small negative number - measured worst case -7.0e-09 over the
    in-repo Peng-Robinson grid - so it can dip below ``-tpd_tol`` for a loosely
    converged split. A minimizer that *is* the partner phase is therefore
    reported as ``"marginal"``, not as an instability. "Is the partner" uses the
    same measure as the stability module's trivial-solution test,
    ``sum_i ln(w_i / x_i^partner)^2 < trivial_tol``, rather than a new
    tolerance.

    Anything else negative means the two-phase answer is not a stable phase set
    and a third phase is required, which this release cannot produce.

    The check always runs on the paths that can run it;
    ``settings.post_split_stability`` decides whether a *failure* raises or is
    only reported, so the diagnostics show the failure either way.

    Returns:
        Diagnostics keys ``post_split_checked``, ``post_split_stable``,
        ``post_split_status``, ``post_split_tpd_min``,
        ``phase_stability_<name>`` and ``phase_stability_tpd_min_<name>``.
    """
    from ..stability import StabilitySettings, stability_tp

    stability_settings = settings.stability_settings or StabilitySettings()
    diagnostics: dict[str, float | int | str | bool] = {"post_split_checked": True}
    worst = math.inf
    failures: list[str] = []
    inconclusive: list[str] = []

    for position, (name, composition) in enumerate(phases):
        partner = phases[1 - position][1] if len(phases) == 2 else None
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
            settings=stability_settings,
        )

        verdict = result.status
        if verdict == "unstable" and partner is not None and result.trial_composition is not None:
            if _is_same_phase(
                np.array(result.trial_composition, dtype=float),
                partner,
                stability_settings.trivial_tol,
            ):
                verdict = "marginal"

        diagnostics[f"phase_stability_{name}"] = verdict
        diagnostics[f"phase_stability_tpd_min_{name}"] = float(result.tpd_min)
        worst = min(worst, float(result.tpd_min))
        if verdict == "unstable":
            failures.append(name)
        elif verdict == "inconclusive":
            inconclusive.append(name)

    diagnostics["post_split_tpd_min"] = worst if math.isfinite(worst) else float("nan")
    status = "unstable" if failures else ("inconclusive" if inconclusive else "stable")
    diagnostics["post_split_status"] = status
    diagnostics["post_split_stable"] = status == "stable"

    if status != "stable" and settings.post_split_stability:
        detail = ", ".join(failures or inconclusive)
        raise ConvergenceError(
            "The converged two-phase solution is not a stable phase set: the post-split "
            f"stability test reports '{status}' for phase(s) {detail} "
            f"(most negative post-split tpd = {worst:.6e}). A third phase is required, and "
            "flash_tp returns at most two phases in this release (multiphase flash is the "
            "next slice). Pass FlashSettings(post_split_stability=False) to receive the "
            "two-phase result anyway, with this failure recorded in diagnostics."
        )

    return diagnostics


def _is_same_phase(w: np.ndarray, other: np.ndarray, trivial_tol: float) -> bool:
    """True when ``w`` is the composition ``other`` in the trivial-solution metric."""
    mask = (w > 0.0) & (other > 0.0)
    if not np.any(mask):
        return False
    if np.any((w > 0.0) != (other > 0.0)):
        return False
    ln_ratio = np.log(w[mask] / other[mask])
    return bool(float(np.sum(ln_ratio * ln_ratio)) < trivial_tol)
