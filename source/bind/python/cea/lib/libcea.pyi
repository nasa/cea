from __future__ import annotations

from typing import Any, Literal, Mapping, Sequence, SupportsFloat, overload

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
VectorLike = Sequence[SupportsFloat] | FloatArray
ScalarOrArray = float | FloatArray

WeightEnthalpyUnits = Literal[
    "J/kg",
    "kJ/kg",
    "cal/kg",
    "kcal/kg",
    "j/kg",
    "kj/kg",
]
MolarEnthalpyUnits = Literal[
    "J/mol",
    "kJ/mol",
    "cal/mol",
    "kcal/mol",
    "J/mole",
    "kJ/mole",
    "cal/mole",
    "kcal/mole",
    "j/mol",
    "kj/mol",
    "j/mole",
    "kj/mole",
]

SUCCESS: int
INVALID_FILENAME: int
INVALID_PROPERTY_TYPE: int
INVALID_EQUILIBRIUM_TYPE: int

LOG_CRITICAL: int
LOG_ERROR: int
LOG_WARNING: int
LOG_INFO: int
LOG_DEBUG: int
LOG_NONE: int

TP: int
HP: int
SP: int
TV: int
UV: int
SV: int

TEMPERATURE: int
PRESSURE: int
VOLUME: int
DENSITY: int
ENTHALPY: int
ENERGY: int
ENTROPY: int
GIBBS_ENERGY: int
GAMMA_S: int
FROZEN_CP: int
FROZEN_CV: int
EQUILIBRIUM_CP: int
EQUILIBRIUM_CV: int
VISCOSITY: int
FROZEN_CONDUCTIVITY: int
EQUILIBRIUM_CONDUCTIVITY: int
FROZEN_PRANDTL: int
EQUILIBRIUM_PRANDTL: int

def _version() -> str: ...
def _version_major() -> int: ...
def _version_minor() -> int: ...
def _version_patch() -> int: ...
def set_log_level(level: int) -> None: ...
def init(path: str | None = None, trans: bool = True) -> None: ...
def init_thermo(thermofile: str) -> None: ...
def init_trans(transfile: str) -> None: ...
def is_initialized() -> bool: ...


class Reactant:
    name: str
    formula: Mapping[str, float] | None
    molecular_weight: float | None
    enthalpy: float | None
    enthalpy_units: str | None
    temperature: float | None

    @overload
    def __init__(
        self,
        name: str,
        formula: Mapping[str, SupportsFloat] | None = None,
        molecular_weight: SupportsFloat = ...,
        enthalpy: SupportsFloat = ...,
        enthalpy_units: WeightEnthalpyUnits = ...,
        temperature: SupportsFloat | None = None,
    ) -> None: ...

    @overload
    def __init__(
        self,
        name: str,
        formula: Mapping[str, SupportsFloat] | None = None,
        molecular_weight: SupportsFloat | None = None,
        enthalpy: SupportsFloat | None = None,
        enthalpy_units: MolarEnthalpyUnits | None = None,
        temperature: SupportsFloat | None = None,
    ) -> None: ...


class Mixture:
    num_species: int

    def __init__(
        self,
        species: list[str | Reactant] | None = None,
        products_from_reactants: bool = False,
        omit: list[str] = ...,
        ions: bool = False,
    ) -> None: ...
    @property
    def species_names(self) -> list[str]: ...
    def moles_to_weights(self, moles: VectorLike) -> FloatArray: ...
    def weights_to_moles(self, weights: VectorLike) -> FloatArray: ...
    def chem_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: VectorLike,
        fuel_weights: VectorLike,
        chem_eq_ratio: float,
    ) -> float: ...
    def weight_eq_ratio_to_of_ratio(
        self,
        oxidant_weights: VectorLike,
        fuel_weights: VectorLike,
        weight_eq_ratio: float,
    ) -> float: ...
    def of_ratio_to_weights(
        self,
        oxidant_weights: VectorLike,
        fuel_weights: VectorLike,
        of_ratio: float,
    ) -> FloatArray: ...
    def calc_property(
        self,
        prop_type: int,
        weights: VectorLike,
        temperature: SupportsFloat | VectorLike,
        pressure: float | None = None,
    ) -> float: ...


class EqSolver:
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        reactants: Mixture,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    @property
    def num_reactants(self) -> int: ...
    @property
    def num_products(self) -> int: ...
    @property
    def num_gas(self) -> int: ...
    @property
    def num_condensed(self) -> int: ...
    @property
    def num_elements(self) -> int: ...
    @property
    def max_equations(self) -> int: ...
    def solve(
        self,
        soln: EqSolution,
        eq_type: int,
        state1: float,
        state2: float,
        amounts: VectorLike,
    ) -> None: ...


class EqSolution:
    @overload
    def __init__(self, solver: EqSolver) -> None: ...
    @overload
    def __init__(
        self,
        solver: EqSolver,
        T_init: float | None = None,
        nj_init: VectorLike | None = None,
    ) -> None: ...
    @property
    def T(self) -> float: ...
    @property
    def P(self) -> float: ...
    @property
    def volume(self) -> float: ...
    @property
    def density(self) -> float: ...
    @property
    def M(self) -> float: ...
    @property
    def MW(self) -> float: ...
    @property
    def enthalpy(self) -> float: ...
    @property
    def energy(self) -> float: ...
    @property
    def entropy(self) -> float: ...
    @property
    def gibbs_energy(self) -> float: ...
    @property
    def gamma_s(self) -> float: ...
    @property
    def cp_fr(self) -> float: ...
    @property
    def cp_eq(self) -> float: ...
    @property
    def cp(self) -> float: ...
    @property
    def cv_fr(self) -> float: ...
    @property
    def cv_eq(self) -> float: ...
    @property
    def cv(self) -> float: ...
    @property
    def viscosity(self) -> float: ...
    @property
    def conductivity_fr(self) -> float: ...
    @property
    def conductivity_eq(self) -> float: ...
    @property
    def Pr_fr(self) -> float: ...
    @property
    def Pr_eq(self) -> float: ...
    @property
    def nj(self) -> FloatArray: ...
    @property
    def ln_nj(self) -> FloatArray: ...
    @property
    def n(self) -> float: ...
    @property
    def mass_fractions(self) -> dict[str, float]: ...
    @property
    def mole_fractions(self) -> dict[str, float]: ...
    @property
    def converged(self) -> bool: ...


class EqDerivatives:
    def __init__(self, solver: EqSolver, solution: EqSolution) -> None: ...
    def compute_derivatives(self, check_closure_defect: bool = False) -> None: ...
    def compute_fd(self, h: float, verbose: bool = False, central: bool = False) -> None: ...


class RocketSolver:
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        reactants: Mixture,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    def solve(
        self,
        soln: RocketSolution,
        weights: VectorLike,
        pc: float,
        pi_p: SupportsFloat | VectorLike | None = None,
        subar: SupportsFloat | VectorLike | None = None,
        supar: SupportsFloat | VectorLike | None = None,
        iac: bool = True,
        n_frz: int | None = None,
        hc: float | None = None,
        tc: float | None = None,
        mdot: float | None = None,
        ac_at: float | None = None,
        tc_est: float | None = None,
    ) -> None: ...


class RocketSolution:
    def __init__(self, solver: RocketSolver) -> None: ...
    @property
    def num_pts(self) -> int: ...
    @property
    def T(self) -> FloatArray: ...
    @property
    def P(self) -> FloatArray: ...
    @property
    def volume(self) -> FloatArray: ...
    @property
    def density(self) -> FloatArray: ...
    @property
    def M(self) -> FloatArray: ...
    @property
    def MW(self) -> FloatArray: ...
    @property
    def enthalpy(self) -> FloatArray: ...
    @property
    def energy(self) -> FloatArray: ...
    @property
    def entropy(self) -> FloatArray: ...
    @property
    def gibbs_energy(self) -> FloatArray: ...
    @property
    def gamma_s(self) -> FloatArray: ...
    @property
    def cp_fr(self) -> FloatArray: ...
    @property
    def cp_eq(self) -> FloatArray: ...
    @property
    def cp(self) -> FloatArray: ...
    @property
    def cv_fr(self) -> FloatArray: ...
    @property
    def cv_eq(self) -> FloatArray: ...
    @property
    def cv(self) -> FloatArray: ...
    @property
    def Mach(self) -> FloatArray: ...
    @property
    def sonic_velocity(self) -> FloatArray: ...
    @property
    def ae_at(self) -> FloatArray: ...
    @property
    def c_star(self) -> FloatArray: ...
    @property
    def coefficient_of_thrust(self) -> FloatArray: ...
    @property
    def Isp(self) -> FloatArray: ...
    @property
    def Isp_vacuum(self) -> FloatArray: ...
    @property
    def viscosity(self) -> FloatArray: ...
    @property
    def conductivity_fr(self) -> FloatArray: ...
    @property
    def conductivity_eq(self) -> FloatArray: ...
    @property
    def Pr_fr(self) -> FloatArray: ...
    @property
    def Pr_eq(self) -> FloatArray: ...
    @property
    def nj(self) -> FloatArray: ...
    @property
    def ln_nj(self) -> FloatArray: ...
    @property
    def n(self) -> FloatArray: ...
    @property
    def mass_fractions(self) -> dict[str, FloatArray]: ...
    @property
    def mole_fractions(self) -> dict[str, FloatArray]: ...
    @property
    def converged(self) -> bool: ...


class ShockSolver:
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        reactants: Mixture,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    def solve(
        self,
        soln: ShockSolution,
        weights: VectorLike,
        T0: float,
        p0: float,
        u1: float | None = None,
        Mach1: float | None = None,
        reflected: bool = True,
        incident_frozen: bool = False,
        reflected_frozen: bool = False,
    ) -> None: ...


class ShockSolution:
    def __init__(self, solver: ShockSolver, reflected: bool = True) -> None: ...
    @property
    def T(self) -> ScalarOrArray: ...
    @property
    def P(self) -> ScalarOrArray: ...
    @property
    def velocity(self) -> ScalarOrArray: ...
    @property
    def Mach(self) -> ScalarOrArray: ...
    @property
    def sonic_velocity(self) -> ScalarOrArray: ...
    @property
    def rho12(self) -> ScalarOrArray: ...
    @property
    def rho52(self) -> ScalarOrArray: ...
    @property
    def P21(self) -> ScalarOrArray: ...
    @property
    def P52(self) -> ScalarOrArray: ...
    @property
    def T21(self) -> ScalarOrArray: ...
    @property
    def T52(self) -> ScalarOrArray: ...
    @property
    def M21(self) -> ScalarOrArray: ...
    @property
    def M52(self) -> ScalarOrArray: ...
    @property
    def v2(self) -> ScalarOrArray: ...
    @property
    def u5_p_v2(self) -> ScalarOrArray: ...
    @property
    def volume(self) -> ScalarOrArray: ...
    @property
    def density(self) -> ScalarOrArray: ...
    @property
    def M(self) -> ScalarOrArray: ...
    @property
    def MW(self) -> ScalarOrArray: ...
    @property
    def enthalpy(self) -> ScalarOrArray: ...
    @property
    def energy(self) -> ScalarOrArray: ...
    @property
    def entropy(self) -> ScalarOrArray: ...
    @property
    def gibbs_energy(self) -> ScalarOrArray: ...
    @property
    def gamma_s(self) -> ScalarOrArray: ...
    @property
    def cp_fr(self) -> ScalarOrArray: ...
    @property
    def cp_eq(self) -> ScalarOrArray: ...
    @property
    def cp(self) -> ScalarOrArray: ...
    @property
    def cv_fr(self) -> ScalarOrArray: ...
    @property
    def cv_eq(self) -> ScalarOrArray: ...
    @property
    def cv(self) -> ScalarOrArray: ...
    @property
    def viscosity(self) -> ScalarOrArray: ...
    @property
    def conductivity_fr(self) -> ScalarOrArray: ...
    @property
    def conductivity_eq(self) -> ScalarOrArray: ...
    @property
    def Pr_fr(self) -> ScalarOrArray: ...
    @property
    def Pr_eq(self) -> ScalarOrArray: ...
    @property
    def nj(self) -> FloatArray: ...
    @property
    def ln_nj(self) -> FloatArray: ...
    @property
    def n(self) -> ScalarOrArray: ...
    @property
    def mass_fractions(self) -> dict[str, ScalarOrArray]: ...
    @property
    def mole_fractions(self) -> dict[str, ScalarOrArray]: ...
    @property
    def converged(self) -> bool: ...


class DetonationSolver:
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        reactants: Mixture,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    @overload
    def __init__(
        self,
        products: Mixture,
        *,
        transport: bool = False,
        ions: bool = False,
        trace: float = -1.0,
        insert: Sequence[str] = (),
        smooth_truncation: bool = False,
        truncation_width: float = -1.0,
    ) -> None: ...
    def solve(
        self,
        soln: DetonationSolution,
        weights: VectorLike,
        T1: float,
        p1: float,
        frozen: bool = False,
    ) -> None: ...


class DetonationSolution:
    def __init__(self, solver: DetonationSolver) -> None: ...
    @property
    def P1(self) -> float: ...
    @property
    def T1(self) -> float: ...
    @property
    def H1(self) -> float: ...
    @property
    def M1(self) -> float: ...
    @property
    def gamma1(self) -> float: ...
    @property
    def sonic_velocity1(self) -> float: ...
    @property
    def P(self) -> float: ...
    @property
    def T(self) -> float: ...
    @property
    def density(self) -> float: ...
    @property
    def enthalpy(self) -> float: ...
    @property
    def energy(self) -> float: ...
    @property
    def gibbs_energy(self) -> float: ...
    @property
    def entropy(self) -> float: ...
    @property
    def Mach(self) -> float: ...
    @property
    def velocity(self) -> float: ...
    @property
    def sonic_velocity(self) -> float: ...
    @property
    def gamma_s(self) -> float: ...
    @property
    def P_P1(self) -> float: ...
    @property
    def T_T1(self) -> float: ...
    @property
    def M_M1(self) -> float: ...
    @property
    def rho_rho1(self) -> float: ...
    @property
    def cp_fr(self) -> float: ...
    @property
    def cv_fr(self) -> float: ...
    @property
    def cp_eq(self) -> float: ...
    @property
    def cv_eq(self) -> float: ...
    @property
    def M(self) -> float: ...
    @property
    def MW(self) -> float: ...
    @property
    def viscosity(self) -> float: ...
    @property
    def conductivity_fr(self) -> float: ...
    @property
    def conductivity_eq(self) -> float: ...
    @property
    def Pr_fr(self) -> float: ...
    @property
    def Pr_eq(self) -> float: ...
    @property
    def nj(self) -> FloatArray: ...
    @property
    def ln_nj(self) -> FloatArray: ...
    @property
    def n(self) -> float: ...
    @property
    def mass_fractions(self) -> dict[str, float]: ...
    @property
    def mole_fractions(self) -> dict[str, float]: ...
    @property
    def converged(self) -> bool: ...


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


__all__: list[str]
def __getattr__(name: str) -> Any: ...
