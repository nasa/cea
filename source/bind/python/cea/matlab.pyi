from __future__ import annotations

from types import SimpleNamespace
from typing import Sequence, SupportsFloat

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = SupportsFloat | VectorLike

__all__: list[str]


class MatlabEqSolution(SimpleNamespace):
    last_error: int
    converged: bool
    T: float
    P: float
    volume: float
    density: float
    M: float
    MW: float
    enthalpy: float
    energy: float
    entropy: float
    gibbs_energy: float
    gamma_s: float
    cp_fr: float
    cp_eq: float
    cp: float
    cv_fr: float
    cv_eq: float
    cv: float
    viscosity: float
    conductivity_fr: float
    conductivity_eq: float
    Pr_fr: float
    Pr_eq: float
    nj: FloatArray
    ln_nj: FloatArray
    n: float
    mass_fractions: dict[str, float]
    mole_fractions: dict[str, float]

def eq_solve(
    eq_type: int,
    reactants: list[str],
    T: float | None = None,
    H: float | None = None,
    S: float | None = None,
    U: float | None = None,
    P: float | None = None,
    V: float | None = None,
    T_reac: ScalarOrArray | None = None,
    fuel_amounts: VectorLike | None = None,
    oxid_amounts: VectorLike | None = None,
    moles: bool = False,
    of_ratio: ScalarOrArray | None = None,
    phi: ScalarOrArray | None = None,
    r_eq: ScalarOrArray | None = None,
    pct_fuel: ScalarOrArray | None = None,
    only: list[str] | None = None,
    omit: list[str] | None = None,
    insert: list[str] | None = None,
    trace: float | None = None,
    transport: bool = False,
    ions: bool = False,
) -> MatlabEqSolution: ...
