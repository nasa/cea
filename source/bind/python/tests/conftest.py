import sys
from pathlib import Path

import pytest

PROJECT_PYTHON_DIR = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def cea_module():
    try:
        import cea
    except Exception as exc:
        pytest.skip(f"cea import failed: {exc}")
    try:
        if hasattr(cea, "is_initialized") and not cea.is_initialized():
            pytest.skip("cea database not initialized")
    except Exception as exc:
        pytest.skip(f"cea initialization check failed: {exc}")
    return cea
