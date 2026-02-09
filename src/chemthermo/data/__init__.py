"""Runtime access to the packaged component databank."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from typing import Any

from ..exceptions import PropertyNotFoundError
from .schema import ANTOINE_KEYS, COMPONENT_KEYS, PARAMETER_KEYS, SCHEMA_VERSION, TOP_LEVEL_KEYS


def normalize_name(name: str) -> str:
    """Normalize a component name for case-insensitive matching."""

    return " ".join(name.casefold().split())


@lru_cache(maxsize=1)
def load_component_database() -> dict[str, dict[str, Any]]:
    """Load the JSON databank and return a mapping by canonical name."""

    data_path = resources.files("chemthermo.data") / "components.json"
    with data_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    for key in TOP_LEVEL_KEYS:
        if key not in payload:
            raise ValueError(f"Missing required key '{key}' in component databank.")

    if payload["schema_version"] != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported schema version {payload['schema_version']}; expected {SCHEMA_VERSION}."
        )

    components = payload["components"]
    mapping: dict[str, dict[str, Any]] = {}
    
    # We validate against COMPONENT_KEYS which now includes MW, Tc, etc.
    # And we check that those keys are objects with PARAMETER_KEYS.
    # We do NOT use Pydantic here to avoid runtime overhead on import if possible,
    # or just keep it simple validation.
    
    for record in components:
        for key in COMPONENT_KEYS:
            if key not in record:
                raise ValueError(f"Component record missing '{key}'.")

        # Check parameter objects
        for param_name in ("MW", "Tc", "Pc", "omega"):
            param = record[param_name]
            for key in PARAMETER_KEYS:  # value, units, source_key
                 if key not in param:
                      raise ValueError(f"Component '{record.get('name')}' property '{param_name}' missing '{key}'.")

        if "antoine" in record and record["antoine"] is not None:
            for key in ANTOINE_KEYS:
                if key not in record["antoine"]:
                    raise ValueError(
                        f"Component '{record.get('name')}' missing Antoine key '{key}'."
                    )
        
        canonical = normalize_name(str(record["name"]))
        if canonical in mapping:
            raise ValueError(f"Duplicate component name after normalization: {record['name']!r}")
        mapping[canonical] = record

    return mapping


def get_component_record(name: str) -> dict[str, Any]:
    """Return a component record by name, raising if missing."""

    canonical = normalize_name(name)
    mapping = load_component_database()
    try:
        return mapping[canonical]
    except KeyError as exc:
        raise PropertyNotFoundError(f"Component '{name}' not found in databank.") from exc


def list_component_names() -> list[str]:
    """Return a sorted list of component names in the databank."""

    return [record["name"] for record in load_component_database().values()]
