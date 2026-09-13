"""PC-SAFT pure-component parameter storage.

Holds the three pure-component parameters every non-associating PC-SAFT
component needs (Gross & Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244):

- ``m``            - number of segments per chain (dimensionless),
- ``sigma_A``      - segment diameter, in Angstrom,
- ``epsilon_k_K``  - segment dispersion energy divided by the Boltzmann
  constant, ``epsilon / k_B``, in K.

The packaged defaults live in ``data/eos/pcsaft.json`` and carry a
``provenance`` block naming their source and how it was verified. Binary
interaction parameters are **not** stored here: ``k_ij`` is a property of the
model instance (``PCSAFTEOS(kij=...)``, ADR-0006), not of a component.

Follows the ``NRTLParameters`` pattern: a packaged JSON payload with a
``schema_version``, a ``model`` label and a list of records, plus a
``from_records`` constructor for user-supplied values. Unknown keys (top level
or inside a record) are ignored, so documentation such as ``provenance`` or a
per-record ``source`` needs no schema bump.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Iterable, Mapping, Sequence

import numpy as np

from ..data import normalize_name
from ..exceptions import ModelError

_PCSAFT_RESOURCE = ("data", "eos", "pcsaft.json")
_SCHEMA_VERSION = 1


class PCSAFTParameterError(ModelError):
    """Raised when PC-SAFT parameters are missing or invalid."""


@dataclass(frozen=True)
class PCSAFTRecord:
    """Pure-component PC-SAFT parameters for one non-associating compound.

    ``name`` is the canonical (normalized) component name. ``MW_g_mol`` is
    optional and unused by the equation of state itself; the packaged records
    omit it because the component databank is the source of truth for molar
    mass.
    """

    name: str
    m: float
    sigma_A: float
    epsilon_k_K: float
    MW_g_mol: float | None = None
    source: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise PCSAFTParameterError("PC-SAFT component names must be non-empty.")
        for label, value in (
            ("m", self.m),
            ("sigma_A", self.sigma_A),
            ("epsilon_k_K", self.epsilon_k_K),
        ):
            if not np.isfinite(value) or value <= 0.0:
                raise PCSAFTParameterError(
                    f"PC-SAFT parameter {label!r} for {self.name!r} must be finite and "
                    f"positive (got {value!r})."
                )
        if self.MW_g_mol is not None and (not np.isfinite(self.MW_g_mol) or self.MW_g_mol <= 0.0):
            raise PCSAFTParameterError(
                f"PC-SAFT parameter 'MW_g_mol' for {self.name!r} must be finite and "
                f"positive when given (got {self.MW_g_mol!r})."
            )


@dataclass(frozen=True)
class PCSAFTParameters:
    """Pure-component PC-SAFT parameters keyed by normalized component name."""

    records: Mapping[str, PCSAFTRecord]

    @classmethod
    def load(cls) -> "PCSAFTParameters":
        """Return the packaged parameter set (Gross & Sadowski 2001, Table 1)."""
        return cls._from_payload(_load_payload())

    @classmethod
    def _from_payload(cls, payload: Mapping[str, object]) -> "PCSAFTParameters":
        if payload.get("schema_version") != _SCHEMA_VERSION:
            raise PCSAFTParameterError(
                "Unsupported PC-SAFT parameter schema "
                f"{payload.get('schema_version')!r}; expected {_SCHEMA_VERSION}."
            )
        if str(payload.get("model", "")).strip().casefold() != "pc-saft":
            raise PCSAFTParameterError("PC-SAFT parameter payload has an unexpected model label.")

        components = payload.get("components")
        if not isinstance(components, list):
            raise PCSAFTParameterError("PC-SAFT parameter payload must define a 'components' list.")
        return cls.from_records(components)

    @classmethod
    def from_records(
        cls, records: Iterable[PCSAFTRecord | Mapping[str, object]]
    ) -> "PCSAFTParameters":
        """Build a parameter set from user-supplied records.

        Each entry is either a :class:`PCSAFTRecord` or a mapping with the keys
        ``name``, ``m``, ``sigma_A``, ``epsilon_k_K`` and the optional
        ``MW_g_mol`` / ``source``. Names are normalized with
        ``chemthermo.data.normalize_name``; duplicates are rejected.
        """
        store: dict[str, PCSAFTRecord] = {}
        for entry in records:
            record = _coerce_record(entry)
            if record.name in store:
                raise PCSAFTParameterError(
                    f"Duplicate PC-SAFT parameters for component {record.name!r}."
                )
            store[record.name] = record
        if not store:
            raise PCSAFTParameterError("PC-SAFT parameter set must contain at least one record.")
        return cls(records=store)

    def names(self) -> tuple[str, ...]:
        """Return the canonical component names this set covers, sorted."""
        return tuple(sorted(self.records))

    def record(self, name: str) -> PCSAFTRecord:
        """Return the record for one component, raising if it is missing."""
        canonical = normalize_name(name)
        try:
            return self.records[canonical]
        except KeyError as exc:
            raise PCSAFTParameterError(f"Missing PC-SAFT parameters for: {name}.") from exc

    def for_components(self, names: Sequence[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(m, sigma_A, epsilon_k_K)`` arrays in the requested order.

        Raises:
            PCSAFTParameterError: If ``names`` is empty or any component has no
                record.
        """
        if not names:
            raise PCSAFTParameterError("PC-SAFT requires at least one component.")

        canonical = [normalize_name(name) for name in names]
        missing = [original for original, key in zip(names, canonical) if key not in self.records]
        if missing:
            missing_list = ", ".join(str(name) for name in missing)
            raise PCSAFTParameterError(f"Missing PC-SAFT parameters for: {missing_list}.")

        picked = [self.records[key] for key in canonical]
        return (
            np.array([r.m for r in picked], dtype=float),
            np.array([r.sigma_A for r in picked], dtype=float),
            np.array([r.epsilon_k_K for r in picked], dtype=float),
        )


def _coerce_record(entry: PCSAFTRecord | Mapping[str, object]) -> PCSAFTRecord:
    if isinstance(entry, PCSAFTRecord):
        return PCSAFTRecord(
            name=normalize_name(entry.name),
            m=entry.m,
            sigma_A=entry.sigma_A,
            epsilon_k_K=entry.epsilon_k_K,
            MW_g_mol=entry.MW_g_mol,
            source=entry.source,
        )
    if not isinstance(entry, Mapping):
        raise PCSAFTParameterError("PC-SAFT component entries must be objects.")
    try:
        name = normalize_name(str(entry["name"]))
        m = float(entry["m"])  # type: ignore[arg-type]
        sigma_A = float(entry["sigma_A"])  # type: ignore[arg-type]
        epsilon_k_K = float(entry["epsilon_k_K"])  # type: ignore[arg-type]
    except KeyError as exc:
        raise PCSAFTParameterError(f"PC-SAFT component entry missing key {exc}.") from exc
    except (TypeError, ValueError) as exc:
        raise PCSAFTParameterError(
            "PC-SAFT component entry contains non-numeric parameters."
        ) from exc

    raw_mw = entry.get("MW_g_mol")
    try:
        mw = None if raw_mw is None else float(raw_mw)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise PCSAFTParameterError(
            "PC-SAFT component entry contains a non-numeric 'MW_g_mol'."
        ) from exc

    return PCSAFTRecord(
        name=name,
        m=m,
        sigma_A=sigma_A,
        epsilon_k_K=epsilon_k_K,
        MW_g_mol=mw,
        source=str(entry.get("source", "")),
    )


@lru_cache(maxsize=1)
def _load_payload() -> dict[str, object]:
    """Load the raw PC-SAFT parameter payload from package resources."""
    data_path = resources.files("chemthermo.parameters")
    for segment in _PCSAFT_RESOURCE:
        data_path = data_path / segment
    with data_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


@lru_cache(maxsize=1)
def get_pcsaft_parameters() -> PCSAFTParameters:
    """Return the default (packaged) PC-SAFT parameter set."""
    return PCSAFTParameters.load()
