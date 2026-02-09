"""Smoke tests for all examples."""

import runpy
import sys
from pathlib import Path
import pytest

# Locate examples directory relative to this test file
EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"

def get_example_scripts():
    """Return a list of all .py files in the examples directory."""
    return [
        p for p in EXAMPLES_DIR.glob("*.py")
        if p.name != "__init__.py" and p.name != "README.md"
    ]

@pytest.mark.parametrize("script_path", get_example_scripts(), ids=lambda p: p.name)
def test_example_runs_successfully(script_path):
    """Run the example script and ensure it doesn't raise an exception."""
    print(f"Running example: {script_path.name}")
    
    # We use runpy to execute the script in the same process, assuming they are clean.
    # Alternatively, subprocess could be used for isolation.
    # runpy is faster but risks side effects on global state (though mostly safe for these scripts).
    
    try:
        runpy.run_path(str(script_path), run_name="__main__")
    except Exception as e:
        pytest.fail(f"Example {script_path.name} failed with error: {e}")
