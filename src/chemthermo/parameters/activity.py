"""Loaders for activity-model parameters."""

from __future__ import annotations

from dataclasses import dataclass

from ..exceptions import ModelError
from .nrtl import NRTLParameters


@dataclass(frozen=True)
class ActivityParameters:
    """Factory for activity-model parameter stores."""

    @staticmethod
    def load(model: str) -> NRTLParameters:
        key = model.strip().casefold()
        if key == "nrtl":
            return NRTLParameters.load()
        raise ModelError(f"Unsupported activity model '{model}'.")
