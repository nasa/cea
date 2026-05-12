from __future__ import annotations

from types import SimpleNamespace
from typing import Optional, Sequence, SupportsFloat

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


def _flatten_rocket_solution(soln: RocketSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        num_pts=int(soln.num_pts),
        T=np.array(soln.T, copy=True),
        P=np.array(soln.P, copy=True),
        volume=np.array(soln.volume, copy=True),
        density=np.array(soln.density, copy=True),
        M=np.array(soln.M, copy=True),
        MW=np.array(soln.MW, copy=True),
        enthalpy=np.array(soln.enthalpy, copy=True),
        energy=np.array(soln.energy, copy=True),
        entropy=np.array(soln.entropy, copy=True),
        gibbs_energy=np.array(soln.gibbs_energy, copy=True),
        gamma_s=np.array(soln.gamma_s, copy=True),
        cp_fr=np.array(soln.cp_fr, copy=True),
        cp_eq=np.array(soln.cp_eq, copy=True),
        cp=np.array(soln.cp, copy=True),
        cv_fr=np.array(soln.cv_fr, copy=True),
        cv_eq=np.array(soln.cv_eq, copy=True),
        cv=np.array(soln.cv, copy=True),
        Mach=np.array(soln.Mach, copy=True),
        sonic_velocity=np.array(soln.sonic_velocity, copy=True),
        ae_at=np.array(soln.ae_at, copy=True),
        c_star=np.array(soln.c_star, copy=True),
        coefficient_of_thrust=np.array(soln.coefficient_of_thrust, copy=True),
        Isp=np.array(soln.Isp, copy=True),
        Isp_vacuum=np.array(soln.Isp_vacuum, copy=True),
        viscosity=np.array(soln.viscosity, copy=True),
        conductivity_fr=np.array(soln.conductivity_fr, copy=True),
        conductivity_eq=np.array(soln.conductivity_eq, copy=True),
        Pr_fr=np.array(soln.Pr_fr, copy=True),
        Pr_eq=np.array(soln.Pr_eq, copy=True),
        nj=np.array(soln.nj, copy=True),
        ln_nj=np.array(soln.ln_nj, copy=True),
        n=np.array(soln.n, copy=True),
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


def _flatten_shock_solution(soln: ShockSolution) -> SimpleNamespace:
    total_moles, species_amounts, log_species_amounts = _reconstruct_shock_amounts(soln)
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        T=_copy_scalar_or_array(soln.T),
        P=_copy_scalar_or_array(soln.P),
        velocity=_copy_scalar_or_array(soln.velocity),
        Mach=_copy_scalar_or_array(soln.Mach),
        sonic_velocity=_copy_scalar_or_array(soln.sonic_velocity),
        rho12=_copy_scalar_or_array(soln.rho12),
        rho52=_copy_scalar_or_array(soln.rho52),
        P21=_copy_scalar_or_array(soln.P21),
        P52=_copy_scalar_or_array(soln.P52),
        T21=_copy_scalar_or_array(soln.T21),
        T52=_copy_scalar_or_array(soln.T52),
        M21=_copy_scalar_or_array(soln.M21),
        M52=_copy_scalar_or_array(soln.M52),
        v2=_copy_scalar_or_array(soln.v2),
        u5_p_v2=_copy_scalar_or_array(soln.u5_p_v2),
        volume=_copy_scalar_or_array(soln.volume),
        density=_copy_scalar_or_array(soln.density),
        M=_copy_scalar_or_array(soln.M),
        MW=_copy_scalar_or_array(soln.MW),
        enthalpy=_copy_scalar_or_array(soln.enthalpy),
        energy=_copy_scalar_or_array(soln.energy),
        entropy=_copy_scalar_or_array(soln.entropy),
        gibbs_energy=_copy_scalar_or_array(soln.gibbs_energy),
        gamma_s=_copy_scalar_or_array(soln.gamma_s),
        cp_fr=_copy_scalar_or_array(soln.cp_fr),
        cp_eq=_copy_scalar_or_array(soln.cp_eq),
        cp=_copy_scalar_or_array(soln.cp),
        cv_fr=_copy_scalar_or_array(soln.cv_fr),
        cv_eq=_copy_scalar_or_array(soln.cv_eq),
        cv=_copy_scalar_or_array(soln.cv),
        viscosity=_copy_scalar_or_array(soln.viscosity),
        conductivity_fr=_copy_scalar_or_array(soln.conductivity_fr),
        conductivity_eq=_copy_scalar_or_array(soln.conductivity_eq),
        Pr_fr=_copy_scalar_or_array(soln.Pr_fr),
        Pr_eq=_copy_scalar_or_array(soln.Pr_eq),
        nj=species_amounts,
        ln_nj=log_species_amounts,
        n=total_moles,
        mass_fractions=_copy_fraction_dict(soln.mass_fractions),
        mole_fractions=_copy_fraction_dict(soln.mole_fractions),
    )


def _flatten_detonation_solution(soln: DetonationSolution) -> SimpleNamespace:
    return SimpleNamespace(
        last_error=int(soln.last_error),
        converged=bool(soln.converged),
        P1=float(soln.P1),
        T1=float(soln.T1),
        H1=float(soln.H1),
        M1=float(soln.M1),
        gamma1=float(soln.gamma1),
        sonic_velocity1=float(soln.sonic_velocity1),
        P=float(soln.P),
        T=float(soln.T),
        density=float(soln.density),
        enthalpy=float(soln.enthalpy),
        energy=float(soln.energy),
        gibbs_energy=float(soln.gibbs_energy),
        entropy=float(soln.entropy),
        Mach=float(soln.Mach),
        velocity=float(soln.velocity),
        sonic_velocity=float(soln.sonic_velocity),
        gamma_s=float(soln.gamma_s),
        P_P1=float(soln.P_P1),
        T_T1=float(soln.T_T1),
        M_M1=float(soln.M_M1),
        rho_rho1=float(soln.rho_rho1),
        cp_fr=float(soln.cp_fr),
        cv_fr=float(soln.cv_fr),
        cp_eq=float(soln.cp_eq),
        cv_eq=float(soln.cv_eq),
        M=float(soln.M),
        MW=float(soln.MW),
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
