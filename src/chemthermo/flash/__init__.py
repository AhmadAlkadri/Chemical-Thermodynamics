"""Flash calculation interfaces and result types."""

from .api import flash_tp
from .results import FlashResult, PhaseResult
from .settings import FlashSettings

__all__ = ["FlashResult", "FlashSettings", "PhaseResult", "flash_tp"]
