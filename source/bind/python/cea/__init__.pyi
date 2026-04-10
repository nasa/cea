from typing import Any

from . import units as units
from .constants import R as R
from .lib.libcea import *
from .lib.libcea import _version as lib_version
from .lib.libcea import _version_major as lib_version_major
from .lib.libcea import _version_minor as lib_version_minor
from .lib.libcea import _version_patch as lib_version_patch

__version__: str
__all__: list[str]

def __getattr__(name: str) -> Any: ...
