from typing import Any, SupportsFloat

import numpy as np
import numpy.typing as npt

from . import units as units
from .constants import R as R
from .lib.libcea import *
from .lib.libcea import _version as lib_version
from .lib.libcea import _version_major as lib_version_major
from .lib.libcea import _version_minor as lib_version_minor
from .lib.libcea import _version_patch as lib_version_patch

__version__: str
__all__: list[str]

FloatArray = npt.NDArray[np.float64]
VectorLike = list[SupportsFloat] | tuple[SupportsFloat, ...] | FloatArray

def eq_solve(
    eq_type: int,
    reactants: list[str],
    T: float | None = None,
    H: float | None = None,
    S: float | None = None,
    U: float | None = None,
    P: float | None = None,
    V: float | None = None,
    T_reac: SupportsFloat | VectorLike | None = None,
    fuel_amounts: VectorLike | None = None,
    oxid_amounts: VectorLike | None = None,
    moles: bool = False,
    of_ratio: SupportsFloat | VectorLike | None = None,
    phi: SupportsFloat | VectorLike | None = None,
    r_eq: SupportsFloat | VectorLike | None = None,
    pct_fuel: SupportsFloat | VectorLike | None = None,
    only: list[str] | None = None,
    omit: list[str] | None = None,
    insert: list[str] | None = None,
    trace: float | None = None,
    transport: bool = False,
    ions: bool = False,
) -> EqSolution: ...

def __getattr__(name: str) -> Any: ...
