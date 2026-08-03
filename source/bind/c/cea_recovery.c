#include "cea_fortran_api.h"

#include <stddef.h>
#include <setjmp.h>
#define CEA_GUARD_BEGIN()                  \
  int prev_active = cea_bindc_abort_recovery_is_active(); \
  if (prev_active == 0)                    \
  {                                        \
    cea_bindc_abort_clear_message();       \
    cea_bindc_abort_set_active(1);         \
    if (setjmp(*cea_bindc_abort_jmp_buf()) != 0) \
    {                                      \
      cea_bindc_abort_set_active(0);       \
      return CEA_FORTRAN_ABORT;            \
    }                                      \
  }

#define CEA_GUARD_END(ierr_expr)           \
  do                                       \
  {                                        \
    cea_err ierr__ = (ierr_expr);          \
    if (prev_active == 0)                  \
    {                                      \
      cea_bindc_abort_set_active(0);       \
    }                                      \
    return ierr__;                         \
  } while (0)

cea_err cea_init(void)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_init_fortran());
}

cea_err cea_init_thermo(const cea_string thermofile)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_init_thermo_fortran(thermofile));
}

cea_err cea_init_trans(const cea_string transfile)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_init_trans_fortran(transfile));
}

cea_err cea_reactant_get_valid_temperature_range(
    const cea_string name,
    cea_real *temperature_min,
    cea_real *temperature_max)
{
  if (name == NULL || temperature_min == NULL || temperature_max == NULL)
  {
    return CEA_INVALID_SIZE;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_reactant_get_valid_temperature_range_fortran(name, temperature_min, temperature_max));
}

cea_err cea_mixture_create(cea_mixture *mix, const cea_int nspecies, const cea_string species[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_fortran(mix, nspecies, species));
}

cea_err cea_mixture_create_w_ions(cea_mixture *mix, const cea_int nspecies, const cea_string species[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_w_ions_fortran(mix, nspecies, species));
}

cea_err cea_mixture_create_from_reactants(
    cea_mixture *mix,
    const cea_int nreac,
    const cea_string reactants[],
    const cea_int nomit,
    const cea_string omit[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_from_reactants_fortran(mix, nreac, reactants, nomit, omit));
}

cea_err cea_mixture_create_from_reactants_w_ions(
    cea_mixture *mix,
    const cea_int nreac,
    const cea_string reactants[],
    const cea_int nomit,
    const cea_string omit[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_from_reactants_w_ions_fortran(mix, nreac, reactants, nomit, omit));
}

cea_err cea_mixture_create_from_input_reactants(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_from_input_reactants_fortran(mix, nreactants, reactants));
}

cea_err cea_mixture_create_from_input_reactants_w_ions(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_from_input_reactants_w_ions_fortran(mix, nreactants, reactants));
}

cea_err cea_mixture_create_products_from_input_reactants(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[],
    const cea_int nomit,
    const cea_string omit[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_products_from_input_reactants_fortran(mix, nreactants, reactants, nomit, omit));
}

cea_err cea_mixture_create_products_from_input_reactants_w_ions(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[],
    const cea_int nomit,
    const cea_string omit[])
{
  if (mix != NULL)
  {
    *mix = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_mixture_create_products_from_input_reactants_w_ions_fortran(mix, nreactants, reactants, nomit, omit));
}

cea_err cea_eqsolver_create(cea_eqsolver *solver, const cea_mixture products)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqsolver_create_fortran(solver, products));
}

cea_err cea_eqsolver_create_with_reactants(cea_eqsolver *solver, const cea_mixture products, const cea_mixture reactants)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqsolver_create_with_reactants_fortran(solver, products, reactants));
}

cea_err cea_eqsolver_create_with_options(cea_eqsolver *solver, const cea_mixture products, const cea_solver_opts opts)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqsolver_create_with_options_fortran(solver, products, opts));
}

cea_err cea_eqsolver_solve(
    const cea_eqsolver solver,
    const cea_equilibrium_type eq_type,
    const cea_real state1,
    const cea_real state2,
    const cea_real amounts[],
    cea_eqsolution solution)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqsolver_solve_fortran(solver, eq_type, state1, state2, amounts, solution));
}

cea_err cea_eqsolver_solve_with_partials(
    const cea_eqsolver solver,
    const cea_equilibrium_type eq_type,
    const cea_real state1,
    const cea_real state2,
    const cea_real amounts[],
    cea_eqsolution solution,
    cea_eqpartials partials)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqsolver_solve_with_partials_fortran(solver, eq_type, state1, state2, amounts, solution, partials));
}

cea_err cea_eqderivatives_compute_derivatives(
    const cea_eqderivatives derivs,
    const cea_eqsolver solver,
    const cea_eqsolution solution,
    const bool check_closure_defect)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqderivatives_compute_derivatives_fortran(derivs, solver, solution, check_closure_defect));
}

cea_err cea_eqderivatives_compute_fd(
    const cea_eqderivatives derivs,
    const cea_eqsolver solver,
    const cea_eqsolution solution,
    const cea_real h,
    const bool verbose,
    const bool central)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_eqderivatives_compute_fd_fortran(derivs, solver, solution, h, verbose, central));
}

cea_err cea_rocket_solver_create(cea_rocket_solver *solver, const cea_mixture products)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_rocket_solver_create_fortran(solver, products));
}

cea_err cea_rocket_solver_create_with_reactants(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_rocket_solver_create_with_reactants_fortran(solver, products, reactants));
}

cea_err cea_rocket_solver_create_with_options(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_rocket_solver_create_with_options_fortran(solver, products, opts));
}

cea_err cea_rocket_solver_solve_iac(
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
    const bool use_tc_est)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_rocket_solver_solve_iac_fortran(
      solver, solution, weights, pc, pi_p, n_pi_p, subar, nsubar, supar, nsupar, n_frz, hc_or_tc, use_hc, tc_est,
      use_tc_est));
}

cea_err cea_rocket_solver_solve_fac(
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
    const bool use_tc_est)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_rocket_solver_solve_fac_fortran(
      solver, solution, weights, pc, pi_p, n_pi_p, subar, nsubar, supar, nsupar, n_frz, hc_or_tc, use_hc,
      mdot_or_acat, use_mdot, tc_est, use_tc_est));
}

cea_err cea_shock_solver_create(cea_shock_solver *solver, const cea_mixture products)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_shock_solver_create_fortran(solver, products));
}

cea_err cea_shock_solver_create_with_reactants(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_shock_solver_create_with_reactants_fortran(solver, products, reactants));
}

cea_err cea_shock_solver_create_with_options(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_shock_solver_create_with_options_fortran(solver, products, opts));
}

cea_err cea_shock_solver_solve(
    const cea_shock_solver solver,
    cea_shock_solution solution,
    const cea_real weights[],
    const cea_real T0,
    const cea_real p0,
    const cea_real mach1_or_u1,
    const bool use_mach,
    const bool reflected,
    const bool incident_frozen,
    const bool reflected_frozen)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_shock_solver_solve_fortran(
      solver, solution, weights, T0, p0, mach1_or_u1, use_mach, reflected, incident_frozen, reflected_frozen));
}

cea_err cea_detonation_solver_create(cea_detonation_solver *solver, const cea_mixture products)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_detonation_solver_create_fortran(solver, products));
}

cea_err cea_detonation_solver_create_with_reactants(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_detonation_solver_create_with_reactants_fortran(solver, products, reactants));
}

cea_err cea_detonation_solver_create_with_options(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_solver_opts opts)
{
  if (solver != NULL)
  {
    *solver = NULL;
  }
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_detonation_solver_create_with_options_fortran(solver, products, opts));
}

cea_err cea_detonation_solver_solve(
    const cea_detonation_solver solver,
    cea_detonation_solution solution,
    const cea_real weights[],
    const cea_real T1,
    const cea_real p1,
    const bool frozen)
{
  CEA_GUARD_BEGIN();
  CEA_GUARD_END(cea_detonation_solver_solve_fortran(solver, solution, weights, T1, p1, frozen));
}
