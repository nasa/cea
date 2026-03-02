__version__ = "3.1.0"

# initialize libcea, loading in the default data files
from cea.lib.libcea import init as libcea_init
import os as _os

if not _os.environ.get("CEA_SKIP_INIT"):
    libcea_init()

# cleanup the root namespace
del libcea_init
del _os

# expose all public methods in the library to the root package.
from cea.lib.libcea import *
from cea.lib.libcea import __all__ as _libcea_all
from cea.lib.libcea import _version as lib_version
from cea.lib.libcea import _version_major as lib_version_major
from cea.lib.libcea import _version_minor as lib_version_minor
from cea.lib.libcea import _version_patch as lib_version_patch

from .constants import R
# Allow attribute-style access (e.g., `cea.units.atm_to_bar`)
from . import units as units

__all__ = list(_libcea_all)
__all__.extend([
    "__version__",
    "lib_version",
    "lib_version_major",
    "lib_version_minor",
    "lib_version_patch",
    "R",
    "units",
])

del _libcea_all
