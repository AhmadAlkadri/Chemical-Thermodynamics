from __future__ import annotations

from typing import TypedDict

import pytest

import chemthermo as ct


class ExpectedStatus(TypedDict):
    status: str


class PCSAFTSanityCase(TypedDict):
    components: list[str]
    temperature_K: float
    pressure_Pa: float
    expected: ExpectedStatus


@pytest.fixture
def pcsaft_pure_sanity() -> PCSAFTSanityCase:
    return {
        "components": ["Methane"],
        "temperature_K": 300.0,
        "pressure_Pa": 101325.0,
        "expected": {"status": "missing_parameters"},
    }


@pytest.fixture
def pcsaft_binary_sanity() -> PCSAFTSanityCase:
    return {
        "components": ["Methane", "Ethane"],
        "temperature_K": 280.0,
        "pressure_Pa": 2.0e6,
        "expected": {"status": "missing_parameters"},
    }


def _pcsaft_status(case: PCSAFTSanityCase) -> str:
    eos = ct.get_eos("pcsaft", components=case["components"])
    assert eos.num_components() == len(case["components"])

    composition = [1.0 / len(case["components"])] * len(case["components"])
    try:
        eos.residual_helmholtz(
            temperature_K=case["temperature_K"],
            volume_m3=1.0,
            composition=composition,
        )
    except ct.PCSAFTParameterError:
        return "missing_parameters"
    except NotImplementedError:
        return "not_implemented"

    return "ok"


def test_pcsaft_pure_sanity(pcsaft_pure_sanity: PCSAFTSanityCase) -> None:
    assert "pcsaft" in ct.list_eos()
    status = _pcsaft_status(pcsaft_pure_sanity)
    assert status == pcsaft_pure_sanity["expected"]["status"]


def test_pcsaft_binary_sanity(pcsaft_binary_sanity: PCSAFTSanityCase) -> None:
    status = _pcsaft_status(pcsaft_binary_sanity)
    assert status == pcsaft_binary_sanity["expected"]["status"]
