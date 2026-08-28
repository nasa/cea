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
NOT_CONVERGED: int
LAST_VALID_SOLUTION: int
FORTRAN_ABORT: int

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

NUM_REACTANTS: int
NUM_PRODUCTS: int
NUM_GAS: int
NUM_CONDENSED: int
NUM_ELEMENTS: int
MAX_EQUATIONS: int

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

ROCKET_TEMPERATURE: int
ROCKET_PRESSURE: int
ROCKET_VOLUME: int
ROCKET_DENSITY: int
ROCKET_M: int
ROCKET_MW: int
ROCKET_ENTHALPY: int
ROCKET_ENERGY: int
ROCKET_ENTROPY: int
ROCKET_GIBBS_ENERGY: int
ROCKET_GAMMA_S: int
ROCKET_FROZEN_CP: int
ROCKET_FROZEN_CV: int
ROCKET_EQUILIBRIUM_CP: int
ROCKET_EQUILIBRIUM_CV: int
MACH: int
SONIC_VELOCITY: int
AE_AT: int
C_STAR: int
COEFFICIENT_OF_THRUST: int
ISP: int
ISP_VACUUM: int
ROCKET_VISCOSITY: int
ROCKET_FROZEN_CONDUCTIVITY: int
ROCKET_EQUILIBRIUM_CONDUCTIVITY: int
ROCKET_FROZEN_PRANDTL: int
ROCKET_EQUILIBRIUM_PRANDTL: int

SHOCK_TEMPERATURE: int
SHOCK_PRESSURE: int
SHOCK_VELOCITY: int
SHOCK_MACH: int
SHOCK_SONIC_VELOCITY: int
SHOCK_RHO12: int
SHOCK_RHO52: int
SHOCK_P21: int
SHOCK_P52: int
SHOCK_T21: int
SHOCK_T52: int
SHOCK_M21: int
SHOCK_M52: int
SHOCK_V2: int
SHOCK_U5_P_V2: int
SHOCK_VOLUME: int
SHOCK_DENSITY: int
SHOCK_M: int
SHOCK_MW: int
SHOCK_ENTHALPY: int
SHOCK_ENERGY: int
SHOCK_ENTROPY: int
SHOCK_GIBBS_ENERGY: int
SHOCK_GAMMA_S: int
SHOCK_FROZEN_CP: int
SHOCK_FROZEN_CV: int
SHOCK_EQUILIBRIUM_CP: int
SHOCK_EQUILIBRIUM_CV: int
SHOCK_VISCOSITY: int
SHOCK_FROZEN_CONDUCTIVITY: int
SHOCK_EQUILIBRIUM_CONDUCTIVITY: int
SHOCK_FROZEN_PRANDTL: int
SHOCK_EQUILIBRIUM_PRANDTL: int

DETONATION_P1: int
DETONATION_T1: int
DETONATION_H1: int
DETONATION_M1: int
DETONATION_GAMMA1: int
DETONATION_V_SONIC1: int
DETONATION_PRESSURE: int
DETONATION_TEMPERATURE: int
DETONATION_DENSITY: int
DETONATION_ENTHALPY: int
DETONATION_ENERGY: int
DETONATION_GIBBS_ENERGY: int
DETONATION_ENTROPY: int
DETONATION_MACH: int
DETONATION_VELOCITY: int
DETONATION_SONIC_VELOCITY: int
DETONATION_GAMMA: int
DETONATION_P_P1: int
DETONATION_T_T1: int
DETONATION_M_M1: int
DETONATION_RHO_RHO1: int
DETONATION_FROZEN_CP: int
DETONATION_FROZEN_CV: int
DETONATION_EQUILIBRIUM_CP: int
DETONATION_EQUILIBRIUM_CV: int
DETONATION_M: int
DETONATION_MW: int
DETONATION_VISCOSITY: int
DETONATION_FROZEN_CONDUCTIVITY: int
DETONATION_EQUILIBRIUM_CONDUCTIVITY: int
DETONATION_FROZEN_PRANDTL: int
DETONATION_EQUILIBRIUM_PRANDTL: int

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

    def get_valid_temperature_range(self) -> tuple[float, float]: ...

    @overload
    def __init__(
        self,
        name: str,
        formula: Mapping[str, SupportsFloat] | None = None,
        molecular_weight: SupportsFloat | None = None,
        enthalpy: None = None,
        enthalpy_units: None = None,
        temperature: SupportsFloat | None = None,
    ) -> None: ...

    @overload
    def __init__(
        self,
        name: str,
        formula: Mapping[str, SupportsFloat] | None = None,
        molecular_weight: SupportsFloat,
        enthalpy: SupportsFloat,
        enthalpy_units: WeightEnthalpyUnits,
        temperature: SupportsFloat | None = None,
    ) -> None: ...

    @overload
    def __init__(
        self,
        name: str,
        formula: Mapping[str, SupportsFloat] | None = None,
        molecular_weight: SupportsFloat | None = None,
        enthalpy: SupportsFloat,
        enthalpy_units: MolarEnthalpyUnits,
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
    ) -> None:
        """Solve with both thermodynamic state inputs preserved in double precision."""
        ...


class EqSolution:
    last_error: int
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


class EqDerivatives:
    last_error: int
    def __init__(self, solver: EqSolver, solution: EqSolution) -> None: ...
    def compute_derivatives(self, check_closure_defect: bool = False) -> None: ...
    def compute_fd(self, h: float, verbose: bool = False, central: bool = False) -> None: ...
    @property
    def dT_dstate1(self) -> float: ...
    @property
    def dT_dstate1_fd(self) -> float: ...
    @property
    def dT_dstate2(self) -> float: ...
    @property
    def dT_dstate2_fd(self) -> float: ...
    @property
    def dn_dstate1(self) -> float: ...
    @property
    def dn_dstate1_fd(self) -> float: ...
    @property
    def dn_dstate2(self) -> float: ...
    @property
    def dn_dstate2_fd(self) -> float: ...
    @property
    def dH_dstate1(self) -> float: ...
    @property
    def dH_dstate1_fd(self) -> float: ...
    @property
    def dH_dstate2(self) -> float: ...
    @property
    def dH_dstate2_fd(self) -> float: ...
    @property
    def dU_dstate1(self) -> float: ...
    @property
    def dU_dstate1_fd(self) -> float: ...
    @property
    def dU_dstate2(self) -> float: ...
    @property
    def dU_dstate2_fd(self) -> float: ...
    @property
    def dG_dstate1(self) -> float: ...
    @property
    def dG_dstate1_fd(self) -> float: ...
    @property
    def dG_dstate2(self) -> float: ...
    @property
    def dG_dstate2_fd(self) -> float: ...
    @property
    def dS_dstate1(self) -> float: ...
    @property
    def dS_dstate1_fd(self) -> float: ...
    @property
    def dS_dstate2(self) -> float: ...
    @property
    def dS_dstate2_fd(self) -> float: ...
    @property
    def dT_dw0(self) -> FloatArray: ...
    @property
    def dT_dw0_fd(self) -> FloatArray: ...
    @property
    def dn_dw0(self) -> FloatArray: ...
    @property
    def dn_dw0_fd(self) -> FloatArray: ...
    @property
    def dnj_dstate1(self) -> FloatArray: ...
    @property
    def dnj_dstate1_fd(self) -> FloatArray: ...
    @property
    def dnj_dstate2(self) -> FloatArray: ...
    @property
    def dnj_dstate2_fd(self) -> FloatArray: ...
    @property
    def dH_dw0(self) -> FloatArray: ...
    @property
    def dH_dw0_fd(self) -> FloatArray: ...
    @property
    def dU_dw0(self) -> FloatArray: ...
    @property
    def dU_dw0_fd(self) -> FloatArray: ...
    @property
    def dG_dw0(self) -> FloatArray: ...
    @property
    def dG_dw0_fd(self) -> FloatArray: ...
    @property
    def dS_dw0(self) -> FloatArray: ...
    @property
    def dS_dw0_fd(self) -> FloatArray: ...
    @property
    def dnj_dw0(self) -> FloatArray: ...
    @property
    def dnj_dw0_fd(self) -> FloatArray: ...


class RocketSolver:
    """Rocket solver; omitted reactants defaults to the products mixture."""

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
    ) -> None:
        """Require 1D weights of length num_reactants in reactant species order.

        If reactants was omitted, weights must follow products.species_names.
        Invalid weight shapes or lengths raise ValueError.
        """
        ...


class RocketSolution:
    last_error: int
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
    """Stations [1, 2, 5]; speeds/Mach use incident, incident, reflected shock frames."""
    last_error: int
    def __init__(self, solver: ShockSolver, reflected: bool = True) -> None: ...
    @property
    def T(self) -> ScalarOrArray: ...
    @property
    def P(self) -> ScalarOrArray: ...
    @property
    def velocity(self) -> ScalarOrArray:
        """Gas speeds [u1, u2, u5] in the shock frames; u5 is also the wall-frame wave speed magnitude."""
        ...
    @property
    def Mach(self) -> ScalarOrArray:
        """velocity / sonic_velocity, including the legacy reflected-frozen sound-speed basis."""
        ...
    @property
    def sonic_velocity(self) -> ScalarOrArray:
        """sqrt(R*T*gamma_s/M); reflected-frozen uses incident frozen Cp at reflected T."""
        ...
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
    def M21(self) -> ScalarOrArray:
        """Incident molecular-weight ratio M2/M1, not a Mach ratio."""
        ...
    @property
    def M52(self) -> ScalarOrArray:
        """Reflected molecular-weight ratio M5/M2, not a Mach ratio."""
        ...
    @property
    def v2(self) -> ScalarOrArray:
        """Station-2 gas speed in the wall frame, u1-u2 [m/s]."""
        ...
    @property
    def u5_p_v2(self) -> ScalarOrArray:
        """Upstream gas speed in the reflected-shock frame [m/s]."""
        ...
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
    def solve(
        self,
        soln: DetonationSolution,
        weights: VectorLike,
        T1: float,
        p1: float,
        frozen: bool = False,
    ) -> None: ...


class DetonationSolution:
    last_error: int
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

__all__: list[str]
def __getattr__(name: str) -> Any: ...
