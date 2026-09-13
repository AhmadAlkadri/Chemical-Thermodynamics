"""Antoine pure-component vapor pressure: the modified-Raoult reference fugacity.

This module is private. Nothing here is exported from ``chemthermo`` or from
``chemthermo.models``.

Equation form
-------------
The packaged databank (``src/chemthermo/data/components.json``) stores Antoine
coefficients in the **base-e, pressure-in-bar** form of M. D. Koretsky,
"Engineering and Chemical Thermodynamics" (2nd ed., 2012), Appendix A.1:

    ln( P^sat / [units] ) = A - B / (T/K + C)                             (1)

with ``units`` taken from the record's ``units`` field (``"bar"`` for every
packaged record; ``"kPa"`` and ``"Pa"`` are accepted by the schema and converted
here). The base is **e, not 10**, and ``C`` is added to ``T`` in kelvin, so
``C`` is negative for the packaged records. Sanity check used in the test suite:
water (A = 11.6834, B = 3816.44, C = -46.13) at 373.15 K gives
``exp(11.6834 - 3816.44 / 327.02) = 1.01309 bar``, i.e. one atmosphere at the
normal boiling point.

Why it lives here
-----------------
``P^sat_i(T)`` is the pure-liquid reference fugacity of the modified-Raoult
model (``f_i^0 = P^sat_i``, with ``phi_i^sat = 1`` and no Poynting correction),
so both :mod:`chemthermo.stability` and :mod:`chemthermo.flash` need it and
neither owns it.

Range policy
------------
An Antoine correlation is a fit over the temperature interval its source
states. Evaluating it outside that interval is extrapolation, and the error is
unbounded and silent. :func:`antoine_saturation_pressures` therefore raises
:class:`chemthermo.InputRangeError` when ``T`` is outside ``[Tmin_K, Tmax_K]``
of **any** component, rather than returning a number nobody can trust.
"""

from __future__ import annotations

import math

import numpy as np

from ..core import Mixture
from ..exceptions import InputRangeError, PropertyNotFoundError
from ..schemas import AntoineCoefficients

#: Multiplier converting a record's pressure unit to Pa.
_UNIT_TO_PA = {"bar": 1.0e5, "kPa": 1.0e3, "Pa": 1.0}


def antoine_temperature_range(mixture: Mixture) -> tuple[float, float]:
    """Return the ``(Tmin_K, Tmax_K)`` interval valid for every component.

    The intersection of the per-component Antoine validity intervals. When the
    intersection is empty the returned lower bound exceeds the upper bound,
    which :func:`antoine_saturation_pressures` turns into an error at the first
    evaluation.

    Raises:
        PropertyNotFoundError: If a component carries no Antoine record.
    """
    lower = -math.inf
    upper = math.inf
    for component in mixture.components:
        antoine = _require_antoine(component.name, component.antoine)
        lower = max(lower, float(antoine.Tmin_K))
        upper = min(upper, float(antoine.Tmax_K))
    return lower, upper


def antoine_saturation_pressures(mixture: Mixture, temperature_K: float) -> np.ndarray:
    """Pure-component saturation pressures in Pa, one per component.

    Args:
        mixture: Mixture supplying the component records.
        temperature_K: Temperature in K.

    Returns:
        ``P^sat_i(T)`` in Pa, from equation (1) of the module docstring.

    Raises:
        PropertyNotFoundError: If a component carries no Antoine record.
        InputRangeError: If ``temperature_K`` lies outside the component's
            ``[Tmin_K, Tmax_K]`` validity interval, or if the correlation
            returns a non-finite or non-positive pressure.
    """
    temperature = float(temperature_K)
    values: list[float] = []
    for component in mixture.components:
        antoine = _require_antoine(component.name, component.antoine)
        if temperature < float(antoine.Tmin_K) or temperature > float(antoine.Tmax_K):
            raise InputRangeError(
                f"Antoine correlation for '{component.name}' is valid over "
                f"[{antoine.Tmin_K:.2f}, {antoine.Tmax_K:.2f}] K; got {temperature:.4f} K. "
                "Extrapolating a vapor-pressure fit outside its stated range is not done "
                "silently."
            )
        denominator = temperature + float(antoine.C)
        if denominator == 0.0:
            raise InputRangeError(
                f"Antoine correlation for '{component.name}' is singular at {temperature:.4f} K."
            )
        scale = _UNIT_TO_PA.get(str(antoine.units))
        if scale is None:  # pragma: no cover - the schema restricts the unit set
            raise InputRangeError(
                f"Unsupported Antoine pressure unit '{antoine.units}' for '{component.name}'."
            )
        pressure = math.exp(float(antoine.A) - float(antoine.B) / denominator) * scale
        if not math.isfinite(pressure) or pressure <= 0.0:
            raise InputRangeError(
                f"Antoine correlation for '{component.name}' returned a non-physical "
                f"saturation pressure at {temperature:.4f} K."
            )
        values.append(pressure)
    return np.array(values, dtype=float)


def _require_antoine(name: str, antoine: AntoineCoefficients | None) -> AntoineCoefficients:
    if antoine is None:
        raise PropertyNotFoundError(
            f"Component '{name}' has no Antoine vapor-pressure record, so the modified-Raoult "
            "reference fugacity f_i^0 = Psat_i(T) cannot be evaluated."
        )
    return antoine
