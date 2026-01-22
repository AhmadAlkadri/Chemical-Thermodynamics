def test_import() -> None:
    import chemthermo as ct

    assert hasattr(ct, "__version__")
