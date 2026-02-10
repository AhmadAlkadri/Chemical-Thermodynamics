from importlib.metadata import version


def test_import() -> None:
    import chemthermo as ct

    assert hasattr(ct, "__version__")


def test_version_matches_installed_metadata() -> None:
    import chemthermo as ct

    assert ct.__version__ == version("chemthermo")
