"""The legacy Wilson-K-bound / Rachford-Rice-root phase-detection heuristic.

Kept unchanged (ADR-0008, ADR-0009) so the pre-``flash-auto-phase-detection``
behavior stays reachable and testable via
``FlashSettings(phase_detection="wilson-heuristic")``, and because gamma-phi
has no stability test to fall back on (ADR-0007). It is deliberately *not*
extended with the post-split stability check: reproducing the old behavior is
the whole point of this path.
"""

from __future__ import annotations

import numpy as np

from ..core import Mixture
from ..exceptions import ModelError
from ..models import ActivityModel, EquationOfState
from ._assemble import _single_phase_result, _two_phase_result
from ._common import wilson_k
from ._split import _rachford_rice, _solve_k_loop
from .results import FlashResult
from .settings import FlashSettings


def _flash_tp_wilson_heuristic(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Legacy TP VLE solver: Wilson K-bounds and a Rachford-Rice root check.

    Kept unchanged (ADR-0008, ADR-0009) so the pre-``flash-auto-phase-detection``
    behavior stays reachable and testable, and because gamma-phi has no
    stability test to fall back on (ADR-0007). It is deliberately *not*
    extended with the post-split stability check: reproducing the old behavior
    is the whole point of this path.
    """
    K = wilson_k(mixture, temperature, pressure)
    if np.any(K <= 0.0):
        raise ModelError("Non-positive K-values encountered in Wilson estimate.")

    k_min = float(np.min(K))
    k_max = float(np.max(K))
    base: dict[str, float | int | str | bool] = {
        "flash_mode": mode,
        "phase_detection": "wilson-heuristic",
        "k_seed": "wilson",
    }

    if np.all(K <= 1.0) or np.all(K >= 1.0):
        phase_name = "liquid" if np.all(K <= 1.0) else "vapor"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=0.0 if phase_name == "liquid" else 1.0,
            diagnostics={
                **base,
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "single_phase_k_bounds",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
            },
        )

    vapor_fraction, f0, f1 = _rachford_rice(z, K)
    if vapor_fraction is None:
        phase_name = "vapor" if f0 > 0.0 else "liquid"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=1.0 if phase_name == "vapor" else 0.0,
            diagnostics={
                **base,
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "rr_no_root",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
                "rr_f0": float(f0),
                "rr_f1": float(f1),
                "rr_status": "no_root",
            },
        )

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=activity_model,
        mode=mode,
        settings=settings,
        z=z,
        K=K,
        vapor_fraction=vapor_fraction,
    )
    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        split.x,
        split.y,
        split.vapor_fraction,
        names=("liquid", "vapor"),
        vapor_fraction=split.vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "max_delta_k": split.max_delta,
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": "VLE",
            "post_split_checked": False,
            "post_split_skipped_reason": (
                "gamma_phi_stability_unsupported"
                if mode == "gamma-phi"
                else "legacy_wilson_heuristic_path"
            ),
        },
    )
