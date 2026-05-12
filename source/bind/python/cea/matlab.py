from __future__ import annotations

from types import SimpleNamespace
from typing import Optional, Sequence, SupportsFloat

import numpy as np
import numpy.typing as npt

from .lib.libcea import EqSolution
from .lib.libcea import eq_solve as _compiled_eq_solve

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = SupportsFloat | VectorLike

__all__ = ["eq_solve"]


def _flatten_eq_solution(soln: EqSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        T=float(soln.T),
        P=float(soln.P),
        volume=float(soln.volume),
        density=float(soln.density),
        M=float(soln.M),
        MW=float(soln.MW),
        enthalpy=float(soln.enthalpy),
        energy=float(soln.energy),
        entropy=float(soln.entropy),
        gibbs_energy=float(soln.gibbs_energy),
        gamma_s=float(soln.gamma_s),
        cp_fr=float(soln.cp_fr),
        cp_eq=float(soln.cp_eq),
        cp=float(soln.cp),
        cv_fr=float(soln.cv_fr),
        cv_eq=float(soln.cv_eq),
        cv=float(soln.cv),
        viscosity=float(soln.viscosity),
        conductivity_fr=float(soln.conductivity_fr),
        conductivity_eq=float(soln.conductivity_eq),
        Pr_fr=float(soln.Pr_fr),
        Pr_eq=float(soln.Pr_eq),
        nj=np.array(soln.nj, copy=True),
        ln_nj=np.array(soln.ln_nj, copy=True),
        n=float(soln.n),
        mass_fractions=dict(soln.mass_fractions),
        mole_fractions=dict(soln.mole_fractions),
    )


def eq_solve(
    eq_type: int,
    reactants: list[str],
    T: Optional[float] = None,
    H: Optional[float] = None,
    S: Optional[float] = None,
    U: Optional[float] = None,
    P: Optional[float] = None,
    V: Optional[float] = None,
    T_reac: Optional[ScalarOrArray] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[ScalarOrArray] = None,
    phi: Optional[ScalarOrArray] = None,
    r_eq: Optional[ScalarOrArray] = None,
    pct_fuel: Optional[ScalarOrArray] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
) -> SimpleNamespace:
    """
    Solve equilibrium problem using CEA methodology (Matlab-style interface).
    """
    soln = _compiled_eq_solve(
        eq_type,
        reactants,
        T=T,
        H=H,
        S=S,
        U=U,
        P=P,
        V=V,
        T_reac=T_reac,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=moles,
        of_ratio=of_ratio,
        phi=phi,
        r_eq=r_eq,
        pct_fuel=pct_fuel,
        only=only,
        omit=omit,
        insert=insert,
        trace=trace,
        transport=transport,
        ions=ions,
    )
    return _flatten_eq_solution(soln)
