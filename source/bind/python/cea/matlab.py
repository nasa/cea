from __future__ import annotations

import warnings
from typing import Optional, Sequence, SupportsFloat

import numpy as np
import numpy.typing as npt

from .lib.libcea import (
    ENERGY,
    ENTHALPY,
    ENTROPY,
    HP,
    SP,
    SV,
    TP,
    TV,
    UV,
    EqSolution,
    EqSolver,
    Mixture,
)

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = SupportsFloat | VectorLike

__all__ = ["eq_solve"]


def _build_mixtures(
    reactants: list[str],
    only: Optional[list[str]],
    omit: Optional[list[str]],
    ions: bool,
) -> tuple[Mixture, Mixture]:
    omit_list = [] if omit is None else omit
    only_list = [] if only is None else only

    reactants_mix = Mixture(reactants, ions=ions)
    if len(only_list) > 0:
        products_mix = Mixture(only_list, ions=ions)
    else:
        products_mix = Mixture(
            reactants,
            products_from_reactants=True,
            omit=omit_list,
            ions=ions,
        )

    return reactants_mix, products_mix


def _build_weights(
    reactants_mix: Mixture,
    reactants: list[str],
    fuel_amounts: Optional[VectorLike],
    oxid_amounts: Optional[VectorLike],
    moles: bool,
    of_ratio: Optional[ScalarOrArray],
    phi: Optional[ScalarOrArray],
    r_eq: Optional[ScalarOrArray],
    pct_fuel: Optional[ScalarOrArray],
) -> FloatArray:
    weights = np.zeros(len(reactants))
    if (fuel_amounts is None) or (oxid_amounts is None):
        raise TypeError("Fuel and oxidizer amounts not defined")

    fuel_amounts = np.asarray(fuel_amounts, dtype=np.double)
    oxid_amounts = np.asarray(oxid_amounts, dtype=np.double)
    if moles:
        fuel_weights = reactants_mix.moles_to_weights(fuel_amounts)
        oxid_weights = reactants_mix.moles_to_weights(oxid_amounts)
    else:
        fuel_weights = fuel_amounts
        oxid_weights = oxid_amounts

    if of_ratio is not None:
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio,
        )
    elif phi is not None:
        of_ratio_calc = reactants_mix.weight_eq_ratio_to_of_ratio(
            oxid_weights,
            fuel_weights,
            phi,
        )
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio_calc,
        )
    elif r_eq is not None:
        of_ratio_calc = reactants_mix.chem_eq_ratio_to_of_ratio(
            oxid_weights,
            fuel_weights,
            r_eq,
        )
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio_calc,
        )
    elif pct_fuel is not None:
        of_ratio_calc = (100.0 - pct_fuel) / pct_fuel
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio_calc,
        )
    else:
        weights = fuel_weights + oxid_weights

    return weights


def _resolve_state1(
    eq_type: int,
    reactants: list[str],
    reactants_mix: Mixture,
    weights: FloatArray,
    T: Optional[float],
    H: Optional[float],
    S: Optional[float],
    U: Optional[float],
    T_reac: Optional[ScalarOrArray],
) -> float:
    state1 = 0.0
    needs_state1 = (T is None) and (H is None) and (S is None) and (U is None)
    if T is not None:
        state1 = T
    elif H is not None:
        state1 = H
    elif S is not None:
        state1 = S
    elif U is not None:
        state1 = U

    if not needs_state1:
        return state1

    if T_reac is None:
        raise TypeError("Reactant temperature not defined")

    if type(T_reac) in [list, np.ndarray]:
        if len(T_reac) == 0:
            raise ValueError("T_reac cannot be empty")
        if len(T_reac) != len(reactants):
            raise ValueError("T_reac must have the same length as reactants")
    elif not isinstance(T_reac, (float, int)):
        raise ValueError("T_reac must be a float, list, or np.ndarray")

    if eq_type == TP or eq_type == TV:
        state1 = T_reac[0] if type(T_reac) in [list, np.ndarray] else T_reac
        warnings.warn(
            "Problem temperature not defined; using first reactant temperature"
        )
    elif eq_type == HP:
        state1 = reactants_mix.calc_property(ENTHALPY, weights, T_reac)
    elif eq_type == SP or eq_type == SV:
        state1 = reactants_mix.calc_property(ENTROPY, weights, T_reac)
    elif eq_type == UV:
        state1 = reactants_mix.calc_property(ENERGY, weights, T_reac)

    return state1


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
) -> EqSolution:
    """
    Solve equilibrium problem using CEA methodology (Matlab-style interface).
    """
    insert_list = [] if insert is None else insert
    reactants_mix, products_mix = _build_mixtures(reactants, only, omit, ions)

    if P is not None:
        state2 = P
    elif V is not None:
        state2 = V
    else:
        raise TypeError("Pressure or volume state not defined")

    weights = _build_weights(
        reactants_mix,
        reactants,
        fuel_amounts,
        oxid_amounts,
        moles,
        of_ratio,
        phi,
        r_eq,
        pct_fuel,
    )
    state1 = _resolve_state1(
        eq_type,
        reactants,
        reactants_mix,
        weights,
        T,
        H,
        S,
        U,
        T_reac,
    )

    solver_kwargs = {
        "reactants": reactants_mix,
        "transport": transport,
        "ions": ions,
        "insert": insert_list,
    }
    if trace is not None:
        solver_kwargs["trace"] = trace
    solver = EqSolver(products_mix, **solver_kwargs)

    soln = EqSolution(solver)
    solver.solve(soln, eq_type, state1, state2, weights)
    return soln
