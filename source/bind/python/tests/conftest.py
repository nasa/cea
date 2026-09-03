import os
import sys
from pathlib import Path

import pytest

PROJECT_PYTHON_DIR = Path(__file__).resolve().parents[1]

# In CI, this fixture is exercising a build CI just produced, so an import or
# init failure is a real failure, not an optional/environment-dependent skip.
_IN_CI = bool(os.environ.get("CI") or os.environ.get("GITHUB_ACTIONS"))


def _skip_or_fail(message):
    if _IN_CI:
        pytest.fail(message)
    pytest.skip(message)


@pytest.fixture(scope="session")
def cea_module():
    try:
        import cea
    except Exception as exc:
        _skip_or_fail(f"cea import failed: {exc}")
    try:
        if hasattr(cea, "is_initialized") and not cea.is_initialized():
            _skip_or_fail("cea database not initialized")
    except Exception as exc:
        _skip_or_fail(f"cea initialization check failed: {exc}")
    return cea
