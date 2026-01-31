"""Errors for the VLLE plugin boundary."""

from __future__ import annotations


class VLLEError(RuntimeError):
    """Base error for VLLE boundary issues."""


class VLLEPluginNotInstalledError(ModuleNotFoundError):
    """Raised when the VLLE plugin package is missing."""


class VLLEPluginError(VLLEError):
    """Raised when the VLLE plugin package is misconfigured."""
