"""Shared pytest fixtures.

Currently only the published NRTL validation fixture (Tessier, Brennecke and
Stadtherr 2000, Problem 1). It is deliberately kept out of the packaged
parameter data (`src/chemthermo/parameters/data/activity/nrtl.json`): it is a
citation-backed validation anchor, not a default the library should silently
apply to user mixtures.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest

import chemthermo as ct

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
TESSIER_FIXTURE_PATH = FIXTURES_DIR / "nrtl" / "tessier2000_problem1.json"


def load_tessier2000_problem1() -> dict[str, Any]:
    """Return the raw Tessier (2000) Problem 1 fixture payload."""
    with TESSIER_FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        payload: dict[str, Any] = json.load(handle)
    return payload


def tessier_component_names(payload: dict[str, Any]) -> list[str]:
    """Return the chemthermo databank names in the paper's component order."""
    return [str(entry["chemthermo_name"]) for entry in payload["components"]]


def tessier_nrtl_parameters(payload: dict[str, Any]) -> ct.NRTLParameters:
    """Build `NRTLParameters` from the fixture's tau/alpha matrices."""
    names = tessier_component_names(payload)
    tau = payload["tau"]
    alpha = payload["alpha"]
    pairs = [
        (
            names[i],
            names[j],
            float(tau[i][j]),
            float(tau[j][i]),
            float(alpha[i][j]),
            float(alpha[j][i]),
        )
        for i in range(len(names))
        for j in range(i + 1, len(names))
    ]
    return ct.NRTLParameters.from_pairs(pairs)


@pytest.fixture(scope="session")
def tessier2000_payload() -> dict[str, Any]:
    return load_tessier2000_problem1()


@pytest.fixture(scope="session")
def tessier2000_names(tessier2000_payload: dict[str, Any]) -> list[str]:
    return tessier_component_names(tessier2000_payload)


@pytest.fixture(scope="session")
def tessier2000_parameters(tessier2000_payload: dict[str, Any]) -> ct.NRTLParameters:
    return tessier_nrtl_parameters(tessier2000_payload)


@pytest.fixture(scope="session")
def tessier2000_model(tessier2000_parameters: ct.NRTLParameters) -> ct.NRTL:
    return ct.NRTL(parameters=tessier2000_parameters)


@pytest.fixture(scope="session")
def tessier2000_ln_gamma(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> Callable[[Any], Any]:
    """Return ln gamma(x) for the Tessier Problem 1 system as a plain callable.

    The temperature is immaterial: Table 1 gives dimensionless tau directly and
    the NRTL equation at fixed tau does not depend on T.

    The argument is rescaled to sum to one before the call. NRTL's ln gamma is
    homogeneous of degree zero in the mole numbers, so this is an exact
    identity, not an approximation; it only lets callers (finite-difference
    Jacobians, for instance) probe off the simplex without tripping the
    composition validator.
    """
    mixture = ct.Mixture.from_database(tessier2000_names, [1 / 3, 1 / 3, 1 / 3], normalize=True)

    def ln_gamma(x: Any) -> Any:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = tessier2000_model.activity_coefficients(
            mixture=mixture,
            temperature_K=298.15,
            composition=[float(value) for value in values],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma
