import sys
import os
from pathlib import Path
import pytest

# Repo root is 4 levels up from this file:
# tests/peet/conftest.py -> tests -> TEMPy_extensions_py3 -> pex (repo root)
_REPO = Path(__file__).parents[3]

_TEMPY_SRC = _REPO / 'TEMPy_py3' / 'src' / 'TEMPy'
_EXT_SRC   = _REPO / 'TEMPy_extensions_py3'

for _p in [str(_TEMPY_SRC), str(_EXT_SRC)]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Private test data — symlinked in the repo but not committed.
DATA_DIR = _REPO / 'not_to_be_distributed_test_data'

# The PEET run directory symlinks into flexo; resolve run1 to find the parent
# peet/ dir that holds the r05_1_* files.
def _peet_dir():
    run1 = DATA_DIR / 'run1'
    if run1.exists():
        return run1.resolve().parent
    return None

PEET_DIR = _peet_dir()


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "needs_data: test requires not_to_be_distributed_test_data"
    )


@pytest.fixture(scope='session')
def data_dir():
    if not DATA_DIR.exists():
        pytest.skip("test data directory not present")
    return DATA_DIR


@pytest.fixture(scope='session')
def peet_dir():
    if PEET_DIR is None or not PEET_DIR.exists():
        pytest.skip("peet data directory not resolvable")
    return PEET_DIR


@pytest.fixture
def cwd_run1(data_dir):
    """Change CWD to run1/ so prm relative paths (initMOTL, fnModParticle) resolve correctly."""
    old = os.getcwd()
    os.chdir(str(data_dir / 'run1'))
    yield
    os.chdir(old)
