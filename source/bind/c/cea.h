#include "stdbool.h"
#include "cea_enum.h"

#ifdef __cplusplus
extern "C"
{
#endif

  // Enumerations
  typedef enum
  {
    CEA_LOG_LEVEL_ENUM
  } cea_log_level;
  typedef enum
  {
    CEA_EQUILIBRIUM_TYPE_ENUM
  } cea_equilibrium_type;
  typedef enum
  {
    CEA_EQUILIBRIUM_SIZE_ENUM
  } cea_equilibrium_size;
  typedef enum
  {
    CEA_PROPERTY_TYPE_ENUM
  } cea_property_type;
  typedef enum
  {
    CEA_ROCKET_PROPERTY_TYPE_ENUM
  } cea_rocket_property_type;
  typedef enum
  {
    CEA_SHOCK_PROPERTY_TYPE_ENUM
  } cea_shock_property_type;
  typedef enum
  {
    CEA_DETONATION_PROPERTY_TYPE_ENUM
  } cea_detonation_property_type;
  typedef enum
  {
    CEA_ERROR_CODE_ENUM
  } cea_error_code;

  // Type Shorthand, Incomplete types
  typedef struct cea_mixture_t *cea_reactant;
  typedef struct cea_mixture_t *cea_mixture;
  typedef struct cea_eqsolver_t *cea_eqsolver;
  typedef struct cea_eqsolution *cea_eqsolution;
  typedef struct cea_eqpartials *cea_eqpartials;
  typedef struct cea_rocket_solver *cea_rocket_solver;
  typedef struct cea_rocket_solution *cea_rocket_solution;
  typedef struct cea_shock_solver *cea_shock_solver;
  typedef struct cea_shock_solution *cea_shock_solution;
  typedef struct cea_detonation_solver *cea_detonation_solver;
  typedef struct cea_detonation_solution *cea_detonation_solution;
  typedef cea_error_code cea_err;
  typedef const char *cea_string;
  typedef int cea_int;
  typedef double cea_real;
  typedef double *cea_array;

  typedef struct
  {
    cea_string name;
    cea_int num_elements;
    const cea_string *elements;
    const cea_real *coefficients;
    bool has_molecular_weight;
    cea_real molecular_weight;
    bool has_enthalpy;
    cea_real enthalpy;
    cea_string enthalpy_units;
    bool has_temperature;
    cea_real temperature;
    cea_string temperature_units;
  } cea_reactant_input;

  // Struct types for optional arguments
  typedef struct
  {
    cea_real trace;
    bool ions;
    bool transport;
    cea_mixture reactants;
    cea_int ninsert;
    const cea_string *insert;
  } cea_solver_opts;

  // Initialize optional arguments
  cea_err cea_solver_opts_init(cea_solver_opts *opts);

  // Input string ownership
  // All input strings/arrays (species, reactants, omit, insert, reactant_input fields) are copied by the
  // C/Fortran layer during the call. Callers retain ownership and may free or reuse their buffers after
  // the call returns.

  // Species name utilities
  cea_err cea_species_name_len(cea_int *name_len);

  cea_int cea_test_add(cea_int a, cea_int b);

  // Version Information
  cea_err cea_version_major(cea_int *major);
  cea_err cea_version_minor(cea_int *minor);
  cea_err cea_version_patch(cea_int *patch);

  // Logging Control
  cea_err cea_set_log_level(const cea_log_level level);

  // Initialization (not thread safe)
  cea_err cea_init();
  cea_err cea_init_thermo(const cea_string thermofile);
  cea_err cea_init_trans(const cea_string transfile);
  cea_err cea_is_initialized(cea_int *initialized);

  //----------------------------------------------------------------------
  // Mixture API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_mixture_create(
      cea_mixture *mix,
      const cea_int nspecies,
      const cea_string species[]);

  cea_err cea_mixture_create_w_ions(
      cea_mixture *mix,
      const cea_int nspecies,
      const cea_string species[]);

  cea_err cea_mixture_create_from_reactants(
      cea_mixture *mix,
      const cea_int nreactants,
      const cea_string reactants[],
      const cea_int nomit,
      const cea_string omit[]);

  cea_err cea_mixture_create_from_reactants_w_ions(
      cea_mixture *mix,
      const cea_int nreactants,
      const cea_string reactants[],
      const cea_int nomit,
      const cea_string omit[]);

  cea_err cea_mixture_create_from_input_reactants(
      cea_mixture *mix,
      const cea_int nreactants,
      const cea_reactant_input reactants[]);

  cea_err cea_mixture_create_from_input_reactants_w_ions(
      cea_mixture *mix,
      const cea_int nreactants,
      const cea_reactant_input reactants[]);

  cea_err cea_mixture_destroy(
      cea_mixture *mix);

  // Get values
  cea_err cea_mixture_get_num_species(
      const cea_mixture mix,
      cea_int *num_species);

  // Deprecated: returns a heap-allocated string that must be freed with cea_string_free.
  cea_err cea_mixture_get_species_name(
      const cea_mixture *mix,
      const cea_int i_species,
      cea_string *species);

  // Deprecated: returns heap-allocated strings that must be freed with cea_string_array_free.
  cea_err cea_mixture_get_species_names(
      const cea_mixture *mix,
      const cea_int nspecies,
      cea_string *species[]);

  // Buffer-based species name retrieval
  cea_err cea_mixture_get_species_name_buf(
      const cea_mixture *mix,
      const cea_int i_species,
      char *species,
      const cea_int buf_len);

  cea_err cea_mixture_get_species_names_buf(
      const cea_mixture *mix,
      const cea_int nspecies,
      char *species,
      const cea_int stride);

  // Free functions for deprecated allocating getters
  cea_err cea_string_free(cea_string species);
  cea_err cea_string_array_free(cea_string species[], const cea_int nspecies);

  // Unit Conversion
  cea_err cea_mixture_moles_to_weights(
      const cea_mixture mix,
      const cea_int len,
      const cea_real moles[],
      cea_real weights[]);

  cea_err cea_mixture_weights_to_moles(
      const cea_mixture mix,
      const cea_int len,
      const cea_real weights[],
      cea_real moles[]);

  cea_err cea_mixture_per_mole_to_per_weight(
      const cea_mixture mix,
      const cea_int len,
      const cea_real per_mole[],
      cea_real per_weight[]);

  cea_err cea_mixture_per_weight_to_per_mole(
      const cea_mixture mix,
      const cea_int len,
      const cea_real per_weight[],
      cea_real per_mole[]);

  // Fuel/Oxidant Mixture Conversion
  cea_err cea_mixture_chem_eq_ratio_to_of_ratio(
      const cea_mixture mix,
      const cea_int len,
      const cea_real oxidant_weights[],
      const cea_real fuel_weights[],
      const cea_real chem_eq_ratio,
      cea_real *of_ratio);

  cea_err cea_mixture_weight_eq_ratio_to_of_ratio(
      const cea_mixture mix,
      const cea_int len,
      const cea_real oxidant_weights[],
      const cea_real fuel_weights[],
      const cea_real weight_eq_ratio,
      cea_real *of_ratio);

  cea_err cea_mixture_of_ratio_to_weights(
      const cea_mixture mix,
      const cea_int len,
      const cea_real oxidant_weights[],
      const cea_real fuel_weights[],
      const cea_real of_ratio,
      cea_real reactant_weights[]);

  // Property Calculation
  cea_err cea_mixture_calc_property(
      const cea_mixture mix,
      const cea_property_type type,
      const cea_int len_weights,
      const cea_real weights[],
      const cea_real temperature,
      cea_real *value);

  cea_err cea_mixture_calc_property_multitemp(
      const cea_mixture mix,
      const cea_property_type type,
      const cea_int len_weights,
      const cea_real weights[],
      const cea_int len_temperatures,
      const cea_real temperatures[],
      cea_real *value);

  cea_err cea_mixture_calc_property_tp(
      const cea_mixture mix,
      const cea_property_type type,
      const cea_int len_weights,
      const cea_real weights[],
      const cea_real temperature,
      const cea_real pressure,
      cea_real *value);

  cea_err cea_mixture_calc_property_tp_multitemp(
      const cea_mixture mix,
      const cea_property_type type,
      const cea_int len_weights,
      const cea_real weights[],
      const cea_int len_temperatures,
      const cea_real temperatures[],
      const cea_real pressure,
      cea_real *value);

  //----------------------------------------------------------------------
  // Equilibrium Solver API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_eqsolver_create(
      cea_eqsolver *solver,
      const cea_mixture products);

  cea_err cea_eqsolver_create_with_reactants(
      cea_eqsolver *solver,
      const cea_mixture products,
      const cea_mixture reactants);

  cea_err cea_eqsolver_create_with_options(
      cea_eqsolver *solver,
      const cea_mixture products,
      const cea_solver_opts options);

  cea_err cea_eqsolver_destroy(
      cea_eqsolver *solver);

  // Solve
  cea_err cea_eqsolver_solve(
      const cea_eqsolver solver,
      const cea_equilibrium_type type,
      const cea_real state1,
      const cea_real state2,
      const cea_array amounts,
      cea_eqsolution soln);

  cea_err cea_eqsolver_solve_with_partials(
      const cea_eqsolver solver,
      const cea_equilibrium_type type,
      const cea_real state1,
      const cea_real state2,
      const cea_array amounts,
      cea_eqsolution soln,
      cea_eqpartials eqpartials);

  // Querry functions
  cea_err cea_eqsolver_get_size(
      const cea_eqsolver solver,
      const cea_equilibrium_size eq_variable,
      cea_int *value);

  // // Function interface
  // cea_eqsolution cea_solve_eq(
  //     const cea_equilibrium_type problem_type,
  //     const cea_real state1,
  //     const cea_real state2,
  //     const cea_int nreactants,
  //     const cea_string creactants[],
  //     const cea_array amounts,
  //     const cea_int nproducts,
  //     const cea_string cproducts[],
  //     const cea_int ninsert,
  //     const cea_string insert[],
  //     const bool set_trace,
  //     const cea_real trace,
  //     const bool transport,
  //     const bool ions
  // );

  //----------------------------------------------------------------------
  // Equilibrium Solution API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_eqsolution_create(
      cea_eqsolution *soln,
      const cea_eqsolver solver);
  cea_err cea_eqsolution_destroy(
      cea_eqsolution *soln);

  // Property Queries
  cea_err cea_eqsolution_get_property(
      const cea_eqsolution soln,
      const cea_property_type type,
      cea_real *value);

  cea_err cea_eqsolution_get_weights(
      const cea_eqsolution soln,
      const cea_int np,
      cea_real weights[],
      const bool log);

  cea_err cea_eqsolution_set_T(
      const cea_eqsolution soln,
      const cea_real T);

  cea_err cea_eqsolution_set_nj(
      const cea_eqsolution soln,
      const cea_eqsolver solver,
      const cea_int np,
      const cea_real nj[]);

  cea_err cea_eqsolution_get_species_amounts(
      const cea_eqsolution soln,
      const cea_int np,
      cea_real amounts[],
      const bool mass);

  cea_err cea_eqsolution_get_moles(
      const cea_eqsolution soln,
      cea_real *moles);

  cea_err cea_eqsolution_get_converged(
      const cea_eqsolution soln,
      int *converged);

  //----------------------------------------------------------------------
  // Equilibrium Partials API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_eqpartials_create(
      cea_eqpartials *eqpartials,
      const cea_eqsolver solver);

  cea_err cea_eqpartials_destroy(
      cea_eqpartials *eqpartials);

  //----------------------------------------------------------------------
  // Rocket Solver API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_rocket_solver_create(
      cea_rocket_solver *solver,
      const cea_mixture products);

  cea_err cea_rocket_solver_create_with_reactants(
      cea_rocket_solver *solver,
      const cea_mixture products,
      const cea_mixture reactants);

  cea_err cea_rocket_solver_create_with_options(
      cea_rocket_solver *solver,
      const cea_mixture products,
      const cea_solver_opts options);

  cea_err cea_rocket_solver_destroy(
      cea_rocket_solver *solver);

  // Get
  cea_err cea_rocket_solver_get_size(
      const cea_rocket_solver solver,
      const cea_equilibrium_size eq_variable,
      cea_int *value);

  // Solve
  cea_err cea_rocket_solver_solve_iac(
      const cea_rocket_solver solver,
      cea_rocket_solution soln,
      const cea_array weights,
      const cea_real pc,
      // Optional: ignored when n_pi_p == 0 (pi_p may be NULL)
      const cea_array pi_p,
      const cea_int n_pi_p,
      const cea_array subar,
      const cea_int nsubar,
      const cea_array supar,
      const cea_int nsupar,
      const cea_int n_frz,
      const cea_real hc_or_tc,
      const bool use_hc,
      const cea_real tc_est,
      const bool use_tc_est);

  cea_err cea_rocket_solver_solve_fac(
      const cea_rocket_solver solver,
      cea_rocket_solution soln,
      const cea_array weights,
      const cea_real pc,
      // Optional: ignored when n_pi_p == 0 (pi_p may be NULL)
      const cea_array pi_p,
      const cea_int n_pi_p,
      const cea_array subar,
      const cea_int nsubar,
      const cea_array supar,
      const cea_int nsupar,
      const cea_int n_frz,
      const cea_real hc_or_tc,
      const bool use_hc,
      const cea_real mdot_or_acat,
      const bool use_mdot,
      const cea_real tc_est,
      const bool use_tc_est);

  //----------------------------------------------------------------------
  // Rocket Solution API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_rocket_solution_create(
      cea_rocket_solution *soln,
      const cea_rocket_solver solver);

  cea_err cea_rocket_solution_destroy(
      cea_rocket_solution *soln);

  // Property Queries
  cea_err cea_rocket_solution_get_size(
      const cea_rocket_solution soln,
      cea_int *num_pts);

  cea_err cea_rocket_solution_get_property(
      const cea_rocket_solution soln,
      const cea_rocket_property_type type,
      const cea_int len,
      cea_real *value);

  cea_err cea_rocket_solution_get_weights(
      const cea_rocket_solution soln,
      const cea_int np,
      const cea_int station,
      cea_real weights[],
      const bool log);

  cea_err cea_rocket_solution_get_species_amounts(
      const cea_rocket_solution soln,
      const cea_int np,
      const cea_int station,
      cea_real amounts[],
      const bool mass);

  cea_err cea_rocket_solution_get_moles(
      const cea_rocket_solution soln,
      cea_real *moles);

  cea_err cea_rocket_solution_get_converged(
      const cea_rocket_solution soln,
      int *converged);

  // cea_err cea_rocket_solution_get_eq_solutions(
  //     const cea_rocket_solution soln,
  //     cea_int npts,
  //     cea_eqsolution *eq_solns[]
  // );

  // cea_err cea_rocket_solution_destroy_eq_solutions(
  //     const cea_rocket_solution soln,
  //     cea_int npts,
  //     cea_eqsolution *eq_solns[]
  // );

  //----------------------------------------------------------------------
  // Shock Solver API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_shock_solver_create(
      cea_shock_solver *solver,
      const cea_mixture products);

  cea_err cea_shock_solver_create_with_reactants(
      cea_shock_solver *solver,
      const cea_mixture products,
      const cea_mixture reactants);

  cea_err cea_shock_solver_create_with_options(
      cea_shock_solver *solver,
      const cea_mixture products,
      const cea_solver_opts options);

  cea_err cea_shock_solver_destroy(
      cea_shock_solver *solver);

  // Get
  cea_err cea_shock_solver_get_size(
      const cea_shock_solver solver,
      const cea_equilibrium_size eq_variable,
      cea_int *value);

  // Solve
  cea_err cea_shock_solver_solve(
      const cea_shock_solver solver,
      cea_shock_solution soln,
      const cea_array weights,
      const cea_real T0,
      const cea_real p0,
      const cea_real mach1_or_u1,
      const bool use_mach,
      const bool refl,
      const bool incd_froz,
      const bool refl_froz);

  //----------------------------------------------------------------------
  // Shock Solution API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_shock_solution_create(
      cea_shock_solution *soln,
      cea_int num_pts);
  cea_err cea_shock_solution_destroy(
      cea_shock_solution *soln);

  // Property Queries
  cea_err cea_shock_solution_get_property(
      const cea_shock_solution soln,
      const cea_shock_property_type type,
      const cea_int len,
      cea_real *value);

  cea_err cea_shock_solution_get_scalar_property(
      const cea_shock_solution soln,
      const cea_shock_property_type type,
      cea_real *value);

  cea_err cea_shock_solution_get_weights(
      const cea_shock_solution soln,
      const cea_int np,
      const cea_int station,
      cea_real weights[],
      const bool log);

  cea_err cea_shock_solution_get_species_amounts(
      const cea_shock_solution soln,
      const cea_int np,
      const cea_int station,
      cea_real amounts[],
      const bool mass);

  cea_err cea_shock_solution_get_moles(
      const cea_shock_solution soln,
      cea_real *moles);

  cea_err cea_shock_solution_get_converged(
      const cea_shock_solution soln,
      int *converged);

  //----------------------------------------------------------------------
  // Detonation Solver API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_detonation_solver_create(
      cea_detonation_solver *solver,
      const cea_mixture products);

  cea_err cea_detonation_solver_create_with_reactants(
      cea_detonation_solver *solver,
      const cea_mixture products,
      const cea_mixture reactants);

  cea_err cea_detonation_solver_create_with_options(
      cea_detonation_solver *solver,
      const cea_mixture products,
      const cea_solver_opts options);

  cea_err cea_detonation_solver_destroy(
      cea_detonation_solver *solver);

  // Get
  cea_err cea_detonation_solver_get_size(
      const cea_detonation_solver solver,
      const cea_equilibrium_size eq_variable,
      cea_int *value);

  // Solve
  cea_err cea_detonation_solver_solve(
      const cea_detonation_solver solver,
      cea_detonation_solution soln,
      const cea_array weights,
      const cea_real T1,
      const cea_real p1,
      const bool frozen);

  //----------------------------------------------------------------------
  // Detonation Solution API
  //----------------------------------------------------------------------

  // Create/Destroy
  cea_err cea_detonation_solution_create(
      cea_detonation_solution *soln);

  cea_err cea_detonation_solution_destroy(
      cea_detonation_solution *soln);

  // Property Queries
  cea_err cea_detonation_solution_get_property(
      const cea_detonation_solution soln,
      const cea_detonation_property_type type,
      const cea_int len,
      cea_real *value);

  cea_err cea_detonation_solution_get_weights(
      const cea_detonation_solution soln,
      const cea_int np,
      cea_real weights[],
      const bool log);

  cea_err cea_detonation_solution_get_species_amounts(
      const cea_detonation_solution soln,
      const cea_int np,
      cea_real amounts[],
      const bool mass);

  cea_err cea_detonation_solution_get_moles(
      const cea_detonation_solution soln,
      cea_real *moles);

  cea_err cea_detonation_solution_get_converged(
      const cea_detonation_solution soln,
      int *converged);

#ifdef __cplusplus
}
#endif
