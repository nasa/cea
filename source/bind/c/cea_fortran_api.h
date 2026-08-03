#ifndef CEA_FORTRAN_API_H
#define CEA_FORTRAN_API_H

#include "cea.h"
#include <setjmp.h>

#ifdef __cplusplus
extern "C" {
#endif

int cea_bindc_abort_recovery_is_active(void);
void cea_bindc_abort_recover(const char *message, cea_int message_len);
void cea_bindc_abort_clear_message(void);
void cea_bindc_abort_set_active(int active);
jmp_buf *cea_bindc_abort_jmp_buf(void);

cea_err cea_init_fortran(void);
cea_err cea_init_thermo_fortran(const cea_string thermofile);
cea_err cea_init_trans_fortran(const cea_string transfile);
cea_err cea_reactant_get_valid_temperature_range_fortran(
    const cea_string name,
    cea_real *temperature_min,
    cea_real *temperature_max);

cea_err cea_mixture_create_fortran(
    cea_mixture *mix,
    const cea_int nspecies,
    const cea_string species[]);
cea_err cea_mixture_create_w_ions_fortran(
    cea_mixture *mix,
    const cea_int nspecies,
    const cea_string species[]);
cea_err cea_mixture_create_from_reactants_fortran(
    cea_mixture *mix,
    const cea_int nreac,
    const cea_string reactants[],
    const cea_int nomit,
    const cea_string omit[]);
cea_err cea_mixture_create_from_reactants_w_ions_fortran(
    cea_mixture *mix,
    const cea_int nreac,
    const cea_string reactants[],
    const cea_int nomit,
    const cea_string omit[]);
cea_err cea_mixture_create_from_input_reactants_fortran(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[]);
cea_err cea_mixture_create_from_input_reactants_w_ions_fortran(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[]);
cea_err cea_mixture_create_products_from_input_reactants_fortran(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[],
    const cea_int nomit,
    const cea_string omit[]);
cea_err cea_mixture_create_products_from_input_reactants_w_ions_fortran(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[],
    const cea_int nomit,
    const cea_string omit[]);

cea_err cea_mixture_of_ratio_to_chem_eq_ratio(
    const cea_mixture mix,
    const cea_int len,
    const cea_real oxidant_weights[],
    const cea_real fuel_weights[],
    const cea_real of_ratio,
    cea_real *chem_eq_ratio);
cea_err cea_mixture_of_ratio_to_weight_eq_ratio(
    const cea_mixture mix,
    const cea_int len,
    const cea_real oxidant_weights[],
    const cea_real fuel_weights[],
    const cea_real of_ratio,
    cea_real *weight_eq_ratio);

cea_err cea_eqsolver_create_fortran(cea_eqsolver *solver, const cea_mixture products);
cea_err cea_eqsolver_create_with_reactants_fortran(
    cea_eqsolver *solver,
    const cea_mixture products,
    const cea_mixture reactants);
cea_err cea_eqsolver_create_with_options_fortran(
    cea_eqsolver *solver,
    const cea_mixture products,
    const cea_solver_opts opts);
cea_err cea_eqsolver_solve_fortran(
    const cea_eqsolver solver,
    const cea_equilibrium_type eq_type,
    const cea_real state1,
    const cea_real state2,
    const cea_real amounts[],
    cea_eqsolution solution);
cea_err cea_eqsolver_solve_with_partials_fortran(
    const cea_eqsolver solver,
    const cea_equilibrium_type eq_type,
    const cea_real state1,
    const cea_real state2,
    const cea_real amounts[],
    cea_eqsolution solution,
    cea_eqpartials partials);

cea_err cea_eqderivatives_compute_derivatives_fortran(
    const cea_eqderivatives derivs,
    const cea_eqsolver solver,
    const cea_eqsolution solution,
    const bool check_closure_defect);
cea_err cea_eqderivatives_compute_fd_fortran(
    const cea_eqderivatives derivs,
    const cea_eqsolver solver,
    const cea_eqsolution solution,
    const cea_real h,
    const bool verbose,
    const bool central);

cea_err cea_rocket_solver_create_fortran(cea_rocket_solver *solver, const cea_mixture products);
cea_err cea_rocket_solver_create_with_reactants_fortran(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants);
cea_err cea_rocket_solver_create_with_options_fortran(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts);
cea_err cea_rocket_solver_solve_iac_fortran(
    const cea_rocket_solver solver,
    cea_rocket_solution solution,
    const cea_real weights[],
    const cea_real pc,
    const cea_real pi_p[],
    const cea_int n_pi_p,
    const cea_real subar[],
    const cea_int nsubar,
    const cea_real supar[],
    const cea_int nsupar,
    const cea_int n_frz,
    const cea_real hc_or_tc,
    const bool use_hc,
    const cea_real tc_est,
    const bool use_tc_est);
cea_err cea_rocket_solver_solve_fac_fortran(
    const cea_rocket_solver solver,
    cea_rocket_solution solution,
    const cea_real weights[],
    const cea_real pc,
    const cea_real pi_p[],
    const cea_int n_pi_p,
    const cea_real subar[],
    const cea_int nsubar,
    const cea_real supar[],
    const cea_int nsupar,
    const cea_int n_frz,
    const cea_real hc_or_tc,
    const bool use_hc,
    const cea_real mdot_or_acat,
    const bool use_mdot,
    const cea_real tc_est,
    const bool use_tc_est);

cea_err cea_shock_solver_create_fortran(cea_shock_solver *solver, const cea_mixture products);
cea_err cea_shock_solver_create_with_reactants_fortran(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants);
cea_err cea_shock_solver_create_with_options_fortran(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts);
cea_err cea_shock_solver_solve_fortran(
    const cea_shock_solver solver,
    cea_shock_solution solution,
    const cea_real weights[],
    const cea_real T0,
    const cea_real p0,
    const cea_real mach1_or_u1,
    const bool use_mach,
    const bool reflected,
    const bool incident_frozen,
    const bool reflected_frozen);

cea_err cea_detonation_solver_create_fortran(cea_detonation_solver *solver, const cea_mixture products);
cea_err cea_detonation_solver_create_with_reactants_fortran(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants);
cea_err cea_detonation_solver_create_with_options_fortran(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts);
cea_err cea_detonation_solver_solve_fortran(
    const cea_detonation_solver solver,
    cea_detonation_solution solution,
    const cea_real weights[],
    const cea_real T1,
    const cea_real p1,
    const bool frozen);

#ifdef __cplusplus
}
#endif

#endif
