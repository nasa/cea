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


class MatlabRocketSolution(SimpleNamespace):
    last_error: int
    converged: bool
    num_pts: int
    T: FloatArray
    P: FloatArray
    volume: FloatArray
    density: FloatArray
    M: FloatArray
    MW: FloatArray
    enthalpy: FloatArray
    energy: FloatArray
    entropy: FloatArray
    gibbs_energy: FloatArray
    gamma_s: FloatArray
    cp_fr: FloatArray
    cp_eq: FloatArray
    cp: FloatArray
    cv_fr: FloatArray
    cv_eq: FloatArray
    cv: FloatArray
    Mach: FloatArray
    sonic_velocity: FloatArray
    ae_at: FloatArray
    c_star: FloatArray
    coefficient_of_thrust: FloatArray
    Isp: FloatArray
    Isp_vacuum: FloatArray
    viscosity: FloatArray
    conductivity_fr: FloatArray
    conductivity_eq: FloatArray
    Pr_fr: FloatArray
    Pr_eq: FloatArray
    nj: FloatArray
    ln_nj: FloatArray
    n: FloatArray
    mass_fractions: dict[str, FloatArray]
    mole_fractions: dict[str, FloatArray]


class MatlabShockSolution(SimpleNamespace):
    last_error: int
    converged: bool
    T: ScalarOrArray
    P: ScalarOrArray
    velocity: ScalarOrArray
    Mach: ScalarOrArray
    sonic_velocity: ScalarOrArray
    rho12: ScalarOrArray
    rho52: ScalarOrArray
    P21: ScalarOrArray
    P52: ScalarOrArray
    T21: ScalarOrArray
    T52: ScalarOrArray
    M21: ScalarOrArray
    M52: ScalarOrArray
    v2: ScalarOrArray
    u5_p_v2: ScalarOrArray
    volume: ScalarOrArray
    density: ScalarOrArray
    M: ScalarOrArray
    MW: ScalarOrArray
    enthalpy: ScalarOrArray
    energy: ScalarOrArray
    entropy: ScalarOrArray
    gibbs_energy: ScalarOrArray
    gamma_s: ScalarOrArray
    cp_fr: ScalarOrArray
    cp_eq: ScalarOrArray
    cp: ScalarOrArray
    cv_fr: ScalarOrArray
    cv_eq: ScalarOrArray
    cv: ScalarOrArray
    viscosity: ScalarOrArray
    conductivity_fr: ScalarOrArray
    conductivity_eq: ScalarOrArray
    Pr_fr: ScalarOrArray
    Pr_eq: ScalarOrArray
    nj: FloatArray
    ln_nj: FloatArray
    n: ScalarOrArray
    mass_fractions: dict[str, ScalarOrArray]
    mole_fractions: dict[str, ScalarOrArray]


class MatlabDetonationSolution(SimpleNamespace):
    last_error: int
    converged: bool
    P1: float
    T1: float
    H1: float
    M1: float
    gamma1: float
    sonic_velocity1: float
    P: float
    T: float
    density: float
    enthalpy: float
    energy: float
    gibbs_energy: float
    entropy: float
    Mach: float
    velocity: float
    sonic_velocity: float
    gamma_s: float
    P_P1: float
    T_T1: float
    M_M1: float
    rho_rho1: float
    cp_fr: float
    cv_fr: float
    cp_eq: float
    cv_eq: float
    M: float
    MW: float
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


def rocket_solve(
    reactants: list[str],
    pc: SupportsFloat,
    pi_p: ScalarOrArray | None = None,
    subar: ScalarOrArray | None = None,
    supar: ScalarOrArray | None = None,
    T_reac: ScalarOrArray | None = None,
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
    iac: bool = True,
    n_frz: int | None = None,
    hc: SupportsFloat | None = None,
    tc: SupportsFloat | None = None,
    mdot: SupportsFloat | None = None,
    ac_at: SupportsFloat | None = None,
    tc_est: SupportsFloat | None = None,
) -> MatlabRocketSolution: ...


def shock_solve(
    reactants: list[str],
    T0: SupportsFloat,
    p0: SupportsFloat,
    *,
    u1: SupportsFloat | VectorLike | None = None,
    Mach1: SupportsFloat | VectorLike | None = None,
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
    reflected: bool = True,
    incident_frozen: bool = False,
    reflected_frozen: bool = False,
) -> MatlabShockSolution: ...


def detonation_solve(
    reactants: list[str],
    T1: SupportsFloat,
    p1: SupportsFloat,
    *,
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
    frozen: bool = False,
) -> MatlabDetonationSolution: ...
