"""NRTL binary interaction parameter storage."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Iterable, Mapping, Sequence

import numpy as np

from ..core import Mixture
from ..data import normalize_name
from ..exceptions import ModelError

_NRTL_RESOURCE = ("data", "activity", "nrtl.json")
_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class NRTLParameters:
    """Binary NRTL parameters keyed by ordered component pairs.

    Payload schema (nrtl.json) requires:
        schema_version: 1
        model: "NRTL"
        pairs: list of objects with components, tau_12, tau_21, alpha_12, alpha_21

    Component names are normalized for lookup. Missing pairs raise ModelError
    when building arrays for a mixture.
    """

    tau: Mapping[tuple[str, str], float]
    alpha: Mapping[tuple[str, str], float]

    @classmethod
    def load(cls) -> "NRTLParameters":
        payload = _load_payload()
        return cls._from_payload(payload)

    @classmethod
    def _from_payload(cls, payload: Mapping[str, object]) -> "NRTLParameters":
        if payload.get("schema_version") != _SCHEMA_VERSION:
            raise ModelError(
                "Unsupported NRTL parameter schema "
                f"{payload.get('schema_version')!r}; expected {_SCHEMA_VERSION}."
            )
        if str(payload.get("model", "")).strip().casefold() != "nrtl":
            raise ModelError("NRTL parameter payload has an unexpected model label.")

        pairs = payload.get("pairs")
        if not isinstance(pairs, list):
            raise ModelError("NRTL parameter payload must define a 'pairs' list.")

        pair_entries: list[tuple[str, str, float, float, float, float]] = []
        for entry in pairs:
            if not isinstance(entry, dict):
                raise ModelError("NRTL pair entries must be objects.")
            components = entry.get("components")
            if not isinstance(components, list) or len(components) != 2:
                raise ModelError("NRTL pair entries must include exactly two components.")
            name_a, name_b = components
            try:
                tau_ab = float(entry["tau_12"])
                tau_ba = float(entry["tau_21"])
                alpha_ab = float(entry["alpha_12"])
                alpha_ba = float(entry["alpha_21"])
            except KeyError as exc:
                raise ModelError(f"NRTL pair entry missing key {exc}.") from exc
            except (TypeError, ValueError) as exc:
                raise ModelError("NRTL pair entry contains non-numeric parameters.") from exc
            pair_entries.append((str(name_a), str(name_b), tau_ab, tau_ba, alpha_ab, alpha_ba))

        return cls.from_pairs(pair_entries)

    @classmethod
    def from_pairs(
        cls,
        pairs: Iterable[tuple[str, str, float, float, float, float]],
    ) -> "NRTLParameters":
        tau: dict[tuple[str, str], float] = {}
        alpha: dict[tuple[str, str], float] = {}

        for name_a, name_b, tau_ab, tau_ba, alpha_ab, alpha_ba in pairs:
            canonical_a = normalize_name(name_a)
            canonical_b = normalize_name(name_b)
            if not canonical_a or not canonical_b:
                raise ModelError("NRTL component names must be non-empty.")
            if canonical_a == canonical_b:
                raise ModelError("NRTL pair components must be distinct.")

            _insert_pair(tau, canonical_a, canonical_b, tau_ab, tau_ba, "tau")
            _insert_pair(alpha, canonical_a, canonical_b, alpha_ab, alpha_ba, "alpha")

        return cls(tau=tau, alpha=alpha)

    def for_mixture(self, mixture: Mixture) -> tuple[np.ndarray, np.ndarray]:
        return self.for_components([component.name for component in mixture.components])

    def for_components(self, names: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
        canonical = [normalize_name(name) for name in names]
        count = len(canonical)

        tau = np.zeros((count, count), dtype=float)
        alpha = np.zeros((count, count), dtype=float)

        for i in range(count):
            for j in range(count):
                if i == j:
                    continue
                key = (canonical[i], canonical[j])
                if key not in self.tau or key not in self.alpha:
                    raise ModelError(
                        f"Missing NRTL parameters for pair {names[i]!r} -> {names[j]!r}."
                    )
                tau[i, j] = float(self.tau[key])
                alpha[i, j] = float(self.alpha[key])

        return tau, alpha


def _insert_pair(
    store: dict[tuple[str, str], float],
    name_a: str,
    name_b: str,
    value_ab: float,
    value_ba: float,
    label: str,
) -> None:
    """Insert ordered pair parameters into a symmetric lookup store."""
    key_ab = (name_a, name_b)
    key_ba = (name_b, name_a)
    if key_ab in store or key_ba in store:
        raise ModelError(f"Duplicate NRTL {label} parameters for pair '{name_a}'/'{name_b}'.")
    store[key_ab] = float(value_ab)
    store[key_ba] = float(value_ba)


@lru_cache(maxsize=1)
def _load_payload() -> dict[str, object]:
    """Load the raw NRTL parameter payload from package resources."""
    data_path = resources.files("chemthermo.parameters")
    for segment in _NRTL_RESOURCE:
        data_path = data_path / segment
    with data_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)
