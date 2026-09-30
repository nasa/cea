"""MATLAB-oriented compatibility wrappers built on top of the Python binding."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Callable, Optional, Sequence, SupportsFloat

import numpy as np
import numpy.typing as npt

from .constants import R
from .lib.libcea import ENTHALPY
from .lib.libcea import DetonationSolution
from .lib.libcea import DetonationSolver
from .lib.libcea import EqSolution
from .lib.libcea import Mixture
from .lib.libcea import RocketSolution
from .lib.libcea import RocketSolver
from .lib.libcea import ShockSolution
from .lib.libcea import ShockSolver
from .lib.libcea import eq_solve as _compiled_eq_solve

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = SupportsFloat | VectorLike

__all__ = ["eq_solve", "rocket_solve", "shock_solve", "detonation_solve"]


def _copy_array(value: object) -> FloatArray:
    return np.array(value, copy=True)


def _copy_scalar_or_array(value: object) -> float | bool | int | FloatArray:
    if isinstance(value, np.ndarray):
        return np.array(value, copy=True)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return np.array(value, dtype=np.float64, copy=True)
    if isinstance(value, (bool, int, float)):
        return value
    return float(value)  # type: ignore[arg-type]


def _copy_fraction_dict(
    fractions: dict[str, object],
) -> dict[str, float | bool | int | FloatArray]:
    return {name: _copy_scalar_or_array(value) for name, value in fractions.items()}


def _coerce_float(name: str, value: SupportsFloat | VectorLike) -> float:
    if isinstance(value, np.ndarray):
        if value.size != 1:
            raise ValueError(f"{name} must be a scalar")
        return float(value.reshape(-1)[0])
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 1:
            raise ValueError(f"{name} must be a scalar")
        return float(value[0])
    return float(value)


def _normalize_float_or_array(
    value: SupportsFloat | VectorLike | None,
) -> float | FloatArray | None:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return float(value.item())
        return np.array(value, dtype=np.float64, copy=True)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return np.asarray(value, dtype=np.float64)
    return float(value)


def _normalize_strings(values: Optional[Sequence[str]]) -> list[str]:
    return [] if values is None else list(values)


def _build_reactant_and_product_mixtures(
    reactants: list[str],
    *,
    only: Optional[Sequence[str]] = None,
    omit: Optional[Sequence[str]] = None,
    ions: bool = False,
) -> tuple[Mixture, Mixture]:
    only_list = _normalize_strings(only)
    omit_list = _normalize_strings(omit)

    reactants_mix = Mixture(reactants, ions=ions)
    if only_list:
        products_mix = Mixture(only_list, ions=ions)
    else:
        products_mix = Mixture(
            reactants,
            products_from_reactants=True,
            omit=omit_list,
            ions=ions,
        )
    return reactants_mix, products_mix


def _resolve_weights(
    reactants_mix: Mixture,
    reactants: Sequence[str],
    *,
    fuel_amounts: Optional[VectorLike],
    oxid_amounts: Optional[VectorLike],
    moles: bool,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
) -> FloatArray:
    weights = np.zeros(len(reactants), dtype=np.float64)

    if fuel_amounts is None or oxid_amounts is None:
        raise TypeError("Fuel and oxidizer amounts not defined")

    fuel_array = np.asarray(fuel_amounts, dtype=np.float64)
    oxid_array = np.asarray(oxid_amounts, dtype=np.float64)
    if moles:
        fuel_weights = reactants_mix.moles_to_weights(fuel_array)
        oxid_weights = reactants_mix.moles_to_weights(oxid_array)
    else:
        fuel_weights = fuel_array
        oxid_weights = oxid_array

    if of_ratio is not None:
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            _coerce_float("of_ratio", of_ratio),
        )
    elif phi is not None:
        of_ratio_calc = reactants_mix.weight_eq_ratio_to_of_ratio(
            oxid_weights,
            fuel_weights,
            _coerce_float("phi", phi),
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
            _coerce_float("r_eq", r_eq),
        )
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio_calc,
        )
    elif pct_fuel is not None:
        pct_fuel_value = _coerce_float("pct_fuel", pct_fuel)
        of_ratio_calc = (100.0 - pct_fuel_value) / pct_fuel_value
        weights = reactants_mix.of_ratio_to_weights(
            oxid_weights,
            fuel_weights,
            of_ratio_calc,
        )
    else:
        weights = fuel_weights + oxid_weights

    return np.asarray(weights, dtype=np.float64)


def _build_rocket_solver(
    products_mix: Mixture,
    reactants_mix: Mixture,
    *,
    trace: Optional[float],
    transport: bool,
    ions: bool,
    insert: Optional[Sequence[str]],
) -> RocketSolver:
    kwargs = {
        "reactants": reactants_mix,
        "transport": transport,
        "ions": ions,
        "insert": _normalize_strings(insert),
    }
    if trace is not None:
        kwargs["trace"] = trace
    return RocketSolver(products_mix, **kwargs)


def _build_shock_solver(
    products_mix: Mixture,
    reactants_mix: Mixture,
    *,
    trace: Optional[float],
    transport: bool,
    ions: bool,
    insert: Optional[Sequence[str]],
) -> ShockSolver:
    kwargs = {
        "reactants": reactants_mix,
        "transport": transport,
        "ions": ions,
        "insert": _normalize_strings(insert),
    }
    if trace is not None:
        kwargs["trace"] = trace
    return ShockSolver(products_mix, **kwargs)


def _build_detonation_solver(
    products_mix: Mixture,
    reactants_mix: Mixture,
    *,
    trace: Optional[float],
    transport: bool,
    ions: bool,
    insert: Optional[Sequence[str]],
) -> DetonationSolver:
    kwargs = {
        "reactants": reactants_mix,
        "transport": transport,
        "ions": ions,
        "insert": _normalize_strings(insert),
    }
    if trace is not None:
        kwargs["trace"] = trace
    return DetonationSolver(products_mix, **kwargs)


_FieldSpec = str | tuple[str, Callable[[object], object]]


def _flatten_fields(
    soln: object,
    fields: Sequence[_FieldSpec],
    wrapper: Callable[[object], object],
    *,
    overrides: Optional[dict[str, object]] = None,
) -> dict[str, object]:
    values: dict[str, object] = {}
    for field in fields:
        name, field_wrapper = field if isinstance(field, tuple) else (field, wrapper)
        if overrides is not None and name in overrides:
            values[name] = overrides[name]
        else:
            values[name] = field_wrapper(getattr(soln, name))
    return values


_EQ_FIELDS: tuple[_FieldSpec, ...] = (
    "T", "P", "volume", "density", "M", "MW", "enthalpy", "energy", "entropy",
    "gibbs_energy", "gamma_s", "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv",
    "viscosity", "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq",
    ("nj", _copy_array), ("ln_nj", _copy_array), "n",
)

_ROCKET_FIELDS: tuple[_FieldSpec, ...] = (
    "T", "P", "volume", "density", "M", "MW", "enthalpy", "energy", "entropy",
    "gibbs_energy", "gamma_s", "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv",
    "Mach", "sonic_velocity", "ae_at", "c_star", "coefficient_of_thrust",
    "Isp", "Isp_vacuum", "viscosity", "conductivity_fr", "conductivity_eq",
    "Pr_fr", "Pr_eq", "nj", "ln_nj", "n",
)


def _flatten_eq_solution(soln: EqSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        **_flatten_fields(soln, _EQ_FIELDS, float),
        mass_fractions=dict(soln.mass_fractions),
        mole_fractions=dict(soln.mole_fractions),
    )


def _flatten_rocket_solution(soln: RocketSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        num_pts=int(soln.num_pts),
        **_flatten_fields(soln, _ROCKET_FIELDS, _copy_array),
        mass_fractions=_copy_fraction_dict(soln.mass_fractions),
        mole_fractions=_copy_fraction_dict(soln.mole_fractions),
    )


def _reconstruct_shock_amounts(soln: ShockSolution) -> tuple[FloatArray, FloatArray, FloatArray]:
    mole_fractions = soln.mole_fractions
    species_names = list(mole_fractions)
    total_moles = 1.0 / np.asarray(soln.M, dtype=np.float64)
    nj = np.column_stack(
        [
            np.asarray(mole_fractions[name], dtype=np.float64) * total_moles
            for name in species_names
        ]
    )
    with np.errstate(divide="ignore"):
        ln_nj = np.log(nj)
    return total_moles, nj, ln_nj


_SHOCK_FIELDS: tuple[_FieldSpec, ...] = (
    "T", "P", "velocity", "Mach", "sonic_velocity", "rho12", "rho52", "P21",
    "P52", "T21", "T52", "M21", "M52", "v2", "u5_p_v2", "volume", "density",
    "M", "MW", "enthalpy", "energy", "entropy", "gibbs_energy", "gamma_s",
    "cp_fr", "cp_eq", "cp", "cv_fr", "cv_eq", "cv", "viscosity",
    "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq", "nj", "ln_nj", "n",
)

_DETONATION_FIELDS: tuple[_FieldSpec, ...] = (
    "P1", "T1", "H1", "M1", "gamma1", "sonic_velocity1", "P", "T", "density",
    "enthalpy", "energy", "gibbs_energy", "entropy", "Mach", "velocity",
    "sonic_velocity", "gamma_s", "P_P1", "T_T1", "M_M1", "rho_rho1", "cp_fr",
    "cv_fr", "cp_eq", "cv_eq", "M", "MW", "viscosity", "conductivity_fr",
    "conductivity_eq", "Pr_fr", "Pr_eq",
    ("nj", _copy_array), ("ln_nj", _copy_array), "n",
)


def _flatten_shock_solution(soln: ShockSolution) -> SimpleNamespace:
    total_moles, species_amounts, log_species_amounts = _reconstruct_shock_amounts(soln)
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        **_flatten_fields(
            soln,
            _SHOCK_FIELDS,
            _copy_scalar_or_array,
            overrides={"nj": species_amounts, "ln_nj": log_species_amounts, "n": total_moles},
        ),
        mass_fractions=_copy_fraction_dict(soln.mass_fractions),
        mole_fractions=_copy_fraction_dict(soln.mole_fractions),
    )


def _flatten_detonation_solution(soln: DetonationSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        **_flatten_fields(soln, _DETONATION_FIELDS, float),
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
    Solve an equilibrium problem and return a flat MATLAB-friendly namespace.
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


def rocket_solve(
    reactants: list[str],
    pc: SupportsFloat,
    pi_p: Optional[ScalarOrArray] = None,
    subar: Optional[ScalarOrArray] = None,
    supar: Optional[ScalarOrArray] = None,
    T_reac: Optional[ScalarOrArray] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    iac: bool = True,
    n_frz: Optional[int] = None,
    hc: Optional[SupportsFloat] = None,
    tc: Optional[SupportsFloat] = None,
    mdot: Optional[SupportsFloat] = None,
    ac_at: Optional[SupportsFloat] = None,
    tc_est: Optional[SupportsFloat] = None,
) -> SimpleNamespace:
    """
    Solve a rocket-performance problem and return flat array/scalar outputs.
    """
    reactants_mix, products_mix = _build_reactant_and_product_mixtures(
        reactants,
        only=only,
        omit=omit,
        ions=ions,
    )
    weights = _resolve_weights(
        reactants_mix,
        reactants,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=moles,
        of_ratio=of_ratio,
        phi=phi,
        r_eq=r_eq,
        pct_fuel=pct_fuel,
    )
    if hc is None and tc is None:
        if T_reac is None:
            raise TypeError("Reactant temperature not defined")
        hc = reactants_mix.calc_property(ENTHALPY, weights, T_reac) / R

    solver = _build_rocket_solver(
        products_mix,
        reactants_mix,
        trace=trace,
        transport=transport,
        ions=ions,
        insert=insert,
    )
    soln = RocketSolution(solver)
    solver.solve(
        soln,
        weights,
        float(pc),
        pi_p=_normalize_float_or_array(pi_p),
        subar=_normalize_float_or_array(subar),
        supar=_normalize_float_or_array(supar),
        iac=iac,
        n_frz=n_frz,
        hc=None if hc is None else float(hc),
        tc=None if tc is None else float(tc),
        mdot=None if mdot is None else float(mdot),
        ac_at=None if ac_at is None else float(ac_at),
        tc_est=None if tc_est is None else float(tc_est),
    )
    return _flatten_rocket_solution(soln)


def shock_solve(
    reactants: list[str],
    T0: SupportsFloat,
    p0: SupportsFloat,
    *,
    u1: Optional[SupportsFloat | VectorLike] = None,
    Mach1: Optional[SupportsFloat | VectorLike] = None,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    reflected: bool = True,
    incident_frozen: bool = False,
    reflected_frozen: bool = False,
) -> SimpleNamespace:
    """
    Solve an incident/reflected shock problem with MATLAB-friendly outputs.
    """
    reactants_mix, products_mix = _build_reactant_and_product_mixtures(
        reactants,
        only=only,
        omit=omit,
        ions=ions,
    )
    weights = _resolve_weights(
        reactants_mix,
        reactants,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=moles,
        of_ratio=of_ratio,
        phi=phi,
        r_eq=r_eq,
        pct_fuel=pct_fuel,
    )
    solver = _build_shock_solver(
        products_mix,
        reactants_mix,
        trace=trace,
        transport=transport,
        ions=ions,
        insert=insert,
    )
    soln = ShockSolution(solver, reflected=reflected)
    solver.solve(
        soln,
        weights,
        float(T0),
        float(p0),
        u1=None if u1 is None else _coerce_float("u1", u1),
        Mach1=None if Mach1 is None else _coerce_float("Mach1", Mach1),
        reflected=reflected,
        incident_frozen=incident_frozen,
        reflected_frozen=reflected_frozen,
    )
    return _flatten_shock_solution(soln)


def detonation_solve(
    reactants: list[str],
    T1: SupportsFloat,
    p1: SupportsFloat,
    *,
    fuel_amounts: Optional[VectorLike] = None,
    oxid_amounts: Optional[VectorLike] = None,
    moles: bool = False,
    of_ratio: Optional[SupportsFloat | VectorLike] = None,
    phi: Optional[SupportsFloat | VectorLike] = None,
    r_eq: Optional[SupportsFloat | VectorLike] = None,
    pct_fuel: Optional[SupportsFloat | VectorLike] = None,
    only: Optional[list[str]] = None,
    omit: Optional[list[str]] = None,
    insert: Optional[list[str]] = None,
    trace: Optional[float] = None,
    transport: bool = False,
    ions: bool = False,
    frozen: bool = False,
) -> SimpleNamespace:
    """
    Solve a detonation problem and return a flat MATLAB-friendly namespace.
    """
    reactants_mix, products_mix = _build_reactant_and_product_mixtures(
        reactants,
        only=only,
        omit=omit,
        ions=ions,
    )
    weights = _resolve_weights(
        reactants_mix,
        reactants,
        fuel_amounts=fuel_amounts,
        oxid_amounts=oxid_amounts,
        moles=moles,
        of_ratio=of_ratio,
        phi=phi,
        r_eq=r_eq,
        pct_fuel=pct_fuel,
    )
    solver = _build_detonation_solver(
        products_mix,
        reactants_mix,
        trace=trace,
        transport=transport,
        ions=ions,
        insert=insert,
    )
    soln = DetonationSolution(solver)
    solver.solve(soln, weights, float(T1), float(p1), frozen=frozen)
    return _flatten_detonation_solution(soln)
