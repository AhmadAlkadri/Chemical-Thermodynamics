from importlib import resources

from chemthermo.data import get_component_record
from chemthermo.parameters import NRTLParameters


def test_packaged_resource_files_exist() -> None:
    assert (resources.files("chemthermo.data") / "components.json").is_file()
    assert (resources.files("chemthermo.data") / "references.bib").is_file()
    assert (resources.files("chemthermo.parameters") / "data" / "activity" / "nrtl.json").is_file()


def test_packaged_runtime_data_loads() -> None:
    record = get_component_record("Methane")
    assert record["name"] == "Methane"

    tau, alpha = NRTLParameters.load().for_components(["Methane", "Ethane"])
    assert tau.shape == (2, 2)
    assert alpha.shape == (2, 2)
