# Import string stuff
from libc.string cimport const_char

# Import the python version information
from cpython.version cimport PY_MAJOR_VERSION

cdef extern from "cea.h":
    ctypedef enum cea_err:
        CEA_SUCCESS
        CEA_INVALID_FILENAME
        CEA_INVALID_PROPERTY_TYPE
        CEA_INVALID_EQUILIBRIUM_TYPE
        CEA_INVALID_INDEX
        CEA_INVALID_SIZE
        CEA_NOT_CONVERGED

    ctypedef enum cea_log_level:
        CEA_LOG_CRITICAL
        CEA_LOG_ERROR
        CEA_LOG_WARNING
        CEA_LOG_INFO
        CEA_LOG_DEBUG
        CEA_LOG_NONE

    ctypedef enum cea_equilibrium_type:
        CEA_TP
        CEA_HP
        CEA_SP
        CEA_TV
        CEA_UV
        CEA_SV

    ctypedef enum cea_equilibrium_size:
        CEA_NUM_REACTANTS
        CEA_NUM_PRODUCTS
        CEA_NUM_GAS
        CEA_NUM_CONDENSED
        CEA_NUM_ELEMENTS
        CEA_MAX_EQUATIONS

    ctypedef enum cea_property_type:
        CEA_TEMPERATURE
        CEA_PRESSURE
        CEA_VOLUME
        CEA_DENSITY
        CEA_M
        CEA_MW
        CEA_ENTHALPY
        CEA_ENERGY
        CEA_ENTROPY
        CEA_GIBBS_ENERGY
        CEA_GAMMA_S
        CEA_FROZEN_CP
        CEA_FROZEN_CV
        CEA_EQUILIBRIUM_CP
        CEA_EQUILIBRIUM_CV
        CEA_VISCOSITY
        CEA_FROZEN_CONDUCTIVITY
        CEA_EQUILIBRIUM_CONDUCTIVITY
        CEA_FROZEN_PRANDTL
        CEA_EQUILIBRIUM_PRANDTL

    ctypedef enum cea_rocket_property_type:
        CEA_ROCKET_TEMPERATURE
        CEA_ROCKET_PRESSURE
        CEA_ROCKET_VOLUME
        CEA_ROCKET_DENSITY
        CEA_ROCKET_M
        CEA_ROCKET_MW
        CEA_ROCKET_ENTHALPY
        CEA_ROCKET_ENERGY
        CEA_ROCKET_ENTROPY
        CEA_ROCKET_GIBBS_ENERGY
        CEA_ROCKET_GAMMA_S
        CEA_ROCKET_FROZEN_CP
        CEA_ROCKET_FROZEN_CV
        CEA_ROCKET_EQUILIBRIUM_CP
        CEA_ROCKET_EQUILIBRIUM_CV
        CEA_MACH
        CEA_SONIC_VELOCITY
        CEA_AE_AT
        CEA_C_STAR
        CEA_COEFFICIENT_OF_THRUST
        CEA_ISP
        CEA_ISP_VACUUM
        CEA_ROCKET_VISCOSITY
        CEA_ROCKET_FROZEN_CONDUCTIVITY
        CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY
        CEA_ROCKET_FROZEN_PRANDTL
        CEA_ROCKET_EQUILIBRIUM_PRANDTL

    ctypedef enum cea_shock_property_type:
        CEA_SHOCK_TEMPERATURE
        CEA_SHOCK_PRESSURE
        CEA_SHOCK_VELOCITY
        CEA_SHOCK_MACH
        CEA_SHOCK_SONIC_VELOCITY
        CEA_SHOCK_RHO12
        CEA_SHOCK_RHO52
        CEA_SHOCK_P21
        CEA_SHOCK_P52
        CEA_SHOCK_T21
        CEA_SHOCK_T52
        CEA_SHOCK_M21
        CEA_SHOCK_M52
        CEA_SHOCK_V2
        CEA_SHOCK_U5_P_V2
        CEA_SHOCK_VOLUME
        CEA_SHOCK_DENSITY
        CEA_SHOCK_M
        CEA_SHOCK_MW
        CEA_SHOCK_ENTHALPY
        CEA_SHOCK_ENERGY
        CEA_SHOCK_ENTROPY
        CEA_SHOCK_GIBBS_ENERGY
        CEA_SHOCK_GAMMA_S
        CEA_SHOCK_FROZEN_CP
        CEA_SHOCK_FROZEN_CV
        CEA_SHOCK_EQUILIBRIUM_CP
        CEA_SHOCK_EQUILIBRIUM_CV
        CEA_SHOCK_VISCOSITY
        CEA_SHOCK_FROZEN_CONDUCTIVITY
        CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY
        CEA_SHOCK_FROZEN_PRANDTL
        CEA_SHOCK_EQUILIBRIUM_PRANDTL

    ctypedef enum cea_detonation_property_type:
        CEA_DETONATION_P1
        CEA_DETONATION_T1
        CEA_DETONATION_H1
        CEA_DETONATION_M1
        CEA_DETONATION_GAMMA1
        CEA_DETONATION_V_SONIC1
        CEA_DETONATION_PRESSURE
        CEA_DETONATION_TEMPERATURE
        CEA_DETONATION_DENSITY
        CEA_DETONATION_ENTHALPY
        CEA_DETONATION_ENERGY
        CEA_DETONATION_GIBBS_ENERGY
        CEA_DETONATION_ENTROPY
        CEA_DETONATION_MACH
        CEA_DETONATION_VELOCITY
        CEA_DETONATION_SONIC_VELOCITY
        CEA_DETONATION_GAMMA
        CEA_DETONATION_SONIC_VELOCITY
        CEA_DETONATION_P_P1
        CEA_DETONATION_T_T1
        CEA_DETONATION_M_M1
        CEA_DETONATION_RHO_RHO1
        CEA_DETONATION_FROZEN_CP
        CEA_DETONATION_FROZEN_CV
        CEA_DETONATION_EQUILIBRIUM_CP
        CEA_DETONATION_EQUILIBRIUM_CV
        CEA_DETONATION_M
        CEA_DETONATION_MW
        CEA_DETONATION_VISCOSITY
        CEA_DETONATION_FROZEN_CONDUCTIVITY
        CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY
        CEA_DETONATION_FROZEN_PRANDTL
        CEA_DETONATION_EQUILIBRIUM_PRANDTL

    # Typdefs
    ctypedef int cea_int
    ctypedef double cea_real
    ctypedef const char* cea_string
    ctypedef double* cea_array
    ctypedef unsigned char cea_bool

    ctypedef struct cea_mixture_t
    ctypedef cea_mixture_t* cea_mixture

    ctypedef struct cea_eqsolver_t
    ctypedef cea_eqsolver_t* cea_eqsolver

    ctypedef struct cea_eqsolution_t
    ctypedef cea_eqsolution_t* cea_eqsolution

    ctypedef struct cea_eqpartials_t
    ctypedef cea_eqpartials_t* cea_eqpartials

    ctypedef struct cea_rocket_solver_t
    ctypedef cea_rocket_solver_t* cea_rocket_solver

    ctypedef struct cea_rocket_solution_t
    ctypedef cea_rocket_solution_t* cea_rocket_solution

    ctypedef struct cea_shock_solver_t
    ctypedef cea_shock_solver_t* cea_shock_solver

    ctypedef struct cea_shock_solution_t
    ctypedef cea_shock_solution_t* cea_shock_solution

    ctypedef struct cea_detonation_solver_t
    ctypedef cea_detonation_solver_t* cea_detonation_solver

    ctypedef struct cea_detonation_solution_t
    ctypedef cea_detonation_solution_t* cea_detonation_solution

    ctypedef struct cea_solver_opts:
        cea_real trace
        cea_bool ions
        cea_bool transport
        cea_mixture reactants
        cea_int ninsert
        const cea_string* insert

    cpdef cea_err cea_solver_opts_init(cea_solver_opts *opts)
    cpdef cea_err cea_species_name_len(cea_int *name_len)

    # Version
    cpdef cea_err cea_version_major(cea_int *major)
    cpdef cea_err cea_version_minor(cea_int *minor)
    cpdef cea_err cea_version_patch(cea_int *patch)
    cpdef cea_err cea_set_log_level(const cea_log_level level)

    # Initialization
    cpdef cea_err cea_init()
    cpdef cea_err cea_init_thermo(const cea_string thermofile)
    cpdef cea_err cea_init_trans(const cea_string transfile)
    cpdef cea_err cea_is_initialized(cea_int *initialized)

    # Mixture
    cpdef cea_err cea_mixture_create(cea_mixture *mix, const cea_int nspecies, const cea_string species[])
    cpdef cea_err cea_mixture_create_w_ions(cea_mixture *mix, const cea_int nspecies, const cea_string species[])
    cpdef cea_err cea_mixture_create_from_reactants(cea_mixture *mix, const cea_int nreac, const cea_string reactants[],
                                                    const cea_int nomit, const cea_string omit[])
    cpdef cea_err cea_mixture_create_from_reactants_w_ions(cea_mixture *mix, const cea_int nreac, const cea_string reactants[],
                                                           const cea_int nomit, const cea_string omit[])
    cpdef cea_err cea_mixture_destroy(cea_mixture *mix)
    cpdef cea_err cea_mixture_get_num_species(const cea_mixture mix, cea_int *num_species)
    cpdef cea_err cea_mixture_get_species_name(const cea_mixture *mix, const cea_int i_species, cea_string *species)
    cpdef cea_err cea_mixture_get_species_names(const cea_mixture *mix, const cea_int nspecies, cea_string *species[])
    cpdef cea_err cea_mixture_get_species_name_buf(const cea_mixture *mix, const cea_int i_species, char *species,
                                                   const cea_int buf_len)
    cpdef cea_err cea_mixture_get_species_names_buf(const cea_mixture *mix, const cea_int nspecies, char *species,
                                                    const cea_int stride)
    cpdef cea_err cea_string_free(cea_string species)
    cpdef cea_err cea_string_array_free(cea_string species[], const cea_int nspecies)
    cpdef cea_err cea_mixture_moles_to_weights(const cea_mixture mix, const cea_int len, const cea_real moles[],
                                               cea_real weights[])
    cpdef cea_err cea_mixture_weights_to_moles(const cea_mixture mix, const cea_int len, const cea_real weights[],
                                               cea_real moles[])
    cpdef cea_err cea_mixture_chem_eq_ratio_to_of_ratio(const cea_mixture mix, const cea_int len,
        const cea_real oxidant_weights[], const cea_real fuel_weights[], const cea_real chem_eq_ratio, cea_real *of_ratio)
    cpdef cea_err cea_mixture_weight_eq_ratio_to_of_ratio(const cea_mixture mix, const cea_int len,
        const cea_real oxidant_weights[], const cea_real fuel_weights[], const cea_real weight_eq_ratio, cea_real *of_ratio)
    cpdef cea_err cea_mixture_of_ratio_to_weights(const cea_mixture mix, const cea_int len,
        const cea_real oxidant_weights[], const cea_real fuel_weights[], const cea_real of_ratio, cea_real reactant_weights[])
    cpdef cea_err cea_mixture_calc_property(const cea_mixture mix, const cea_property_type prop_type,
        const cea_int len_weights, const cea_real weights[], const cea_real temperature, cea_real *value)
    cpdef cea_err cea_mixture_calc_property_multitemp(const cea_mixture mix, const cea_property_type prop_type,
        const cea_int len_weights, const cea_real weights[], const cea_int len_temperatures,
        const cea_real temperatures[], cea_real *value)
    cpdef cea_err cea_mixture_calc_property_tp(const cea_mixture mix, const cea_property_type prop_type,
        const cea_int len_weights, const cea_real weights[], const cea_real temperature,
        const cea_real pressure, cea_real *value)
    cpdef cea_err cea_mixture_calc_property_tp_multitemp(const cea_mixture mix, const cea_property_type prop_type,
        const cea_int len_weights, const cea_real weights[], const cea_int len_temperatures,
        const cea_real temperatures[], const cea_real pressure, cea_real *value)

    # Equilibrium Solver
    cpdef cea_err cea_eqsolver_create(cea_eqsolver *solver, const cea_mixture products)
    cpdef cea_err cea_eqsolver_create_with_reactants(cea_eqsolver *solver, const cea_mixture products,
                                                     const cea_mixture reactants)
    cpdef cea_err cea_eqsolver_create_with_options(cea_eqsolver *solver, const cea_mixture products,
                                                   const cea_solver_opts opts)
    cpdef cea_err cea_eqsolver_destroy(cea_eqsolver *solver)
    cpdef cea_err cea_eqsolver_solve(const cea_eqsolver solver, const cea_equilibrium_type eq_type,
        const cea_real state1, const cea_real state2, const cea_array amounts, cea_eqsolution soln)
    cpdef cea_err cea_eqsolver_solve_with_partials(const cea_eqsolver solver, const cea_equilibrium_type eq_type,
        const cea_real state1, const cea_real state2, const cea_array amounts, cea_eqsolution soln, cea_eqpartials partials)
    cpdef cea_err cea_eqsolver_get_size(const cea_eqsolver solver, const cea_equilibrium_size eq_variable,
        cea_int *value)

    # Equilibrium Solution
    cpdef cea_err cea_eqsolution_create(cea_eqsolution *solution, const cea_eqsolver solver)
    cpdef cea_err cea_eqsolution_destroy(cea_eqsolution *solution)
    cpdef cea_err cea_eqsolution_get_property(const cea_eqsolution solution, const cea_property_type type,
                                              cea_real *value)
    cpdef cea_err cea_eqsolution_get_weights(const cea_eqsolution solution, const cea_int np,
                                             cea_real weights[], const cea_bool log)
    cpdef cea_err cea_eqsolution_set_T(const cea_eqsolution solution, const cea_real T)
    cpdef cea_err cea_eqsolution_set_nj(const cea_eqsolution solution, const cea_eqsolver solver,
                                        const cea_int np, const cea_real nj[])
    cpdef cea_err cea_eqsolution_get_species_amounts(const cea_eqsolution solution, const cea_int np,
                                                     cea_real amounts[], const cea_bool mass)
    cpdef cea_err cea_eqsolution_get_moles(const cea_eqsolution solution, cea_real *value)
    cpdef cea_err cea_eqsolution_get_converged(const cea_eqsolution solution, bint *converged)

    # Equilibrium Partials
    cpdef cea_err cea_eqpartials_create(cea_eqpartials *partials, const cea_eqsolver solver)
    cpdef cea_err cea_eqpartials_destroy(cea_eqpartials *partials)

    # Rocket Solver
    cpdef cea_err cea_rocket_solver_create(cea_rocket_solver *solver, const cea_mixture products)
    cpdef cea_err cea_rocket_solver_create_with_reactants(cea_rocket_solver *solver, const cea_mixture products,
                                                          const cea_mixture reactants)
    cpdef cea_err cea_rocket_solver_create_with_options(cea_rocket_solver *solver, const cea_mixture products,
                                                        const cea_solver_opts opts)
    cpdef cea_err cea_rocket_solver_destroy(cea_rocket_solver *solver)
    # pi_p is optional when n_pi_p == 0.
    cpdef cea_err cea_rocket_solver_solve_iac(const cea_rocket_solver solver, cea_rocket_solution soln,
                                              const cea_array weights, const cea_real pc, const cea_array pi_p,
                                              const cea_int n_pi_p, const cea_array subar, const cea_int nsubar,
                                              const cea_array supar, const cea_int nsupar, const cea_int n_frz,
                                              const cea_real hc_or_tc, const cea_bool use_hc,
                                              const cea_real tc_est, const cea_bool use_tc_est)
    # pi_p is optional when n_pi_p == 0.
    cpdef cea_err cea_rocket_solver_solve_fac(const cea_rocket_solver solver, cea_rocket_solution soln,
                                              const cea_array weights, const cea_real pc, const cea_array pi_p,
                                              const cea_int n_pi_p, const cea_array subar, const cea_int nsubar,
                                              const cea_array supar, const cea_int nsupar, const cea_int n_frz,
                                              const cea_real hc_or_tc, const cea_bool use_hc,
                                              const cea_real mdot_or_acat, const cea_bool use_mdot,
                                              const cea_real tc_est, const cea_bool use_tc_est)
    cpdef cea_err cea_rocket_solver_get_size(const cea_rocket_solver solver, const cea_equilibrium_size eq_variable,
                                             cea_int *value)

    # Rocket Solution
    cpdef cea_err cea_rocket_solution_create(cea_rocket_solution *solution, const cea_rocket_solver solver)
    cpdef cea_err cea_rocket_solution_destroy(cea_rocket_solution *solution)
    cpdef cea_err cea_rocket_solution_get_size(cea_rocket_solution solution, cea_int *num_pts)
    cpdef cea_err cea_rocket_solution_get_property(const cea_rocket_solution solution,
                                                   const cea_rocket_property_type type, const cea_int len,
                                                   cea_real value[])
    cpdef cea_err cea_rocket_solution_get_weights(const cea_rocket_solution solution, const cea_int np,
                                                  const cea_int station, cea_real weights[], const cea_bool log)
    cpdef cea_err cea_rocket_solution_get_species_amounts(const cea_rocket_solution solution, const cea_int np,
                                                          const cea_int station, cea_real amounts[],
                                                          const cea_bool mass)
    cpdef cea_err cea_rocket_solution_get_moles(const cea_rocket_solution solution, cea_real *value)
    cpdef cea_err cea_rocket_solution_get_converged(const cea_rocket_solution solution, bint *converged)
    # cpdef cea_err cea_rocket_solution_get_eq_solutions(const cea_rocket_solution solution, cea_int npts,
    #                                                    cea_eqsolution *eq_solution[])
    # cpdef cea_err cea_rocket_solution_destroy_eq_solutions(const cea_rocket_solution solution, cea_int npts,
    #                                                        cea_eqsolution *eq_solution[])

    # Shock Solver
    cpdef cea_err cea_shock_solver_create(cea_shock_solver *solver, const cea_mixture products)
    cpdef cea_err cea_shock_solver_create_with_reactants(cea_shock_solver *solver, const cea_mixture products,
                                                         const cea_mixture reactants)
    cpdef cea_err cea_shock_solver_create_with_options(cea_shock_solver *solver, const cea_mixture products,
                                                       const cea_solver_opts options)
    cpdef cea_err cea_shock_solver_destroy(cea_shock_solver *solver)
    cpdef cea_err cea_shock_solver_get_size(const cea_shock_solver solver, const cea_equilibrium_size eq_variable,
                                            cea_int *value)
    cpdef cea_err cea_shock_solver_solve(const cea_shock_solver solver, cea_shock_solution soln,
                                         const cea_array weights, const cea_real T0, const cea_real p0,
                                         const cea_real mach_or_u1, const cea_bool use_mach, const cea_bool refl,
                                         const cea_bool incd_froz, const cea_bool refl_froz)

    # Shock Solution
    cpdef cea_err cea_shock_solution_create(cea_shock_solution *solution, const cea_int num_pts)
    cpdef cea_err cea_shock_solution_destroy(cea_shock_solution *solution)
    cpdef cea_err cea_shock_solution_get_scalar_property(const cea_shock_solution solution,
                                                         const cea_shock_property_type type, cea_real *value)
    cpdef cea_err cea_shock_solution_get_property(const cea_shock_solution solution,
                                                  const cea_shock_property_type type, const cea_int len,
                                                  cea_real value[])
    cpdef cea_err cea_shock_solution_get_weights(const cea_shock_solution solution, const cea_int np,
                                                 const cea_int station, cea_real weights[], const cea_bool log)
    cpdef cea_err cea_shock_solution_get_species_amounts(const cea_shock_solution solution, const cea_int np,
                                                         const cea_int station, cea_real amounts[],
                                                         const cea_bool mass)
    cpdef cea_err cea_shock_solution_get_moles(const cea_shock_solution solution, cea_real *value)
    cpdef cea_err cea_shock_solution_get_converged(const cea_shock_solution solution, bint *converged)

    # Detonation Solver
    cpdef cea_err cea_detonation_solver_create(cea_detonation_solver *solver, const cea_mixture products)
    cpdef cea_err cea_detonation_solver_create_with_reactants(cea_detonation_solver *solver, const cea_mixture products,
                                                              const cea_mixture reactants)
    cpdef cea_err cea_detonation_solver_create_with_options(cea_detonation_solver *solver, const cea_mixture products,
                                                            const cea_solver_opts options)
    cpdef cea_err cea_detonation_solver_destroy(cea_detonation_solver *solver)
    cpdef cea_err cea_detonation_solver_get_size(const cea_detonation_solver solver, const cea_equilibrium_size eq_variable,
                                                 cea_int *value)
    cpdef cea_err cea_detonation_solver_solve(const cea_detonation_solver solver, cea_detonation_solution soln,
                                              const cea_array weights, const cea_real T1, const cea_real p1,
                                              const cea_bool frozen)
    cpdef cea_err cea_detonation_solution_get_weights(const cea_detonation_solution solution, const cea_int np,
                                                      cea_real weights[], const cea_bool log)
    cpdef cea_err cea_detonation_solution_get_species_amounts(const cea_detonation_solution solution, const cea_int np,
                                                              cea_real amounts[], const cea_bool mass)
    cpdef cea_err cea_detonation_solution_get_moles(const cea_detonation_solution solution, cea_real *value)
    cpdef cea_err cea_detonation_solution_get_converged(const cea_detonation_solution solution, bint *converged)

    # Detonation Solution
    cpdef cea_err cea_detonation_solution_create(cea_detonation_solution *solution)
    cpdef cea_err cea_detonation_solution_destroy(cea_detonation_solution *solution)
    cpdef cea_err cea_detonation_solution_get_property(const cea_detonation_solution solution,
                                                       const cea_detonation_property_type type, const cea_int len,
                                                       cea_real *value)
