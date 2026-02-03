"""Registry helpers for residual-Helmholtz EOS models."""

from __future__ import annotations

from typing import Callable

from ..exceptions import ModelError
from .api import EOSProtocol

EOSFactory = Callable[..., EOSProtocol]

_EOS_REGISTRY: dict[str, EOSFactory] = {}


def _normalize_key(name: str) -> str:
    return "".join(ch for ch in name.casefold() if ch.isalnum())


def register_eos(name: str, factory: EOSFactory, *, replace: bool = False) -> None:
    """Register an EOS factory under a normalized name."""

    key = _normalize_key(name)
    if not key:
        raise ModelError("EOS name must be non-empty.")
    if not replace and key in _EOS_REGISTRY:
        raise ModelError(f"EOS '{name}' is already registered.")
    _EOS_REGISTRY[key] = factory


def get_eos(name: str, **kwargs: object) -> EOSProtocol:
    """Return an EOS instance by name."""

    key = _normalize_key(name)
    try:
        factory = _EOS_REGISTRY[key]
    except KeyError as exc:
        raise ModelError(f"Unsupported EOS '{name}'.") from exc
    return factory(**kwargs)


def list_eos() -> list[str]:
    """Return registered EOS names (normalized)."""

    return sorted(_EOS_REGISTRY.keys())
