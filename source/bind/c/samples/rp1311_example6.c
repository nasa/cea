#include "stdlib.h"
#include "stdio.h"
#include "cea.h"

#define TRUE 1
#define FALSE 0
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))

int main(void) {

    //------------------------------------------------------------------
    // Problem Specification
    //------------------------------------------------------------------

    // Species
    const cea_string reactants[]     = {"H2(L)", "O2(L)"};
    const cea_real fuel_weights[]    = {    1.0,     0.0};
    const cea_real oxidant_weights[] = {    0.0,     1.0};
    const int nr = LEN(reactants);

    // Products
    const cea_string omitted_products[] = {""};

    // Thermo States
    const cea_real p0 = 1.0;
    const cea_real T0 = 298.15;
    const cea_real eq_ratio = 1.0;


    //----------------------------------------------------------------
    // CEA Setup
    //----------------------------------------------------------------

    cea_set_log_level(CEA_LOG_DEBUG);
    cea_init();

    // Mixtures
    cea_mixture reac, prod;
    cea_mixture_create(&reac, nr, reactants);
    cea_mixture_create_from_reactants(&prod,
        nr, reactants,
        0, omitted_products
    );

    // Solver
    cea_detonation_solver solver;
    cea_detonation_solver_create_with_reactants(&solver, prod, reac);

    // Solution
    cea_detonation_solution soln;
    cea_detonation_solution_create(&soln);


    //----------------------------------------------------------------
    // Solve the detonation problem
    //----------------------------------------------------------------

    cea_real of_ratio;
    cea_real* weights = (cea_real*) calloc(nr, sizeof(cea_real));
    cea_mixture_chem_eq_ratio_to_of_ratio(reac, LEN(reactants), oxidant_weights, fuel_weights, eq_ratio, &of_ratio);
    cea_mixture_of_ratio_to_weights(reac, LEN(reactants), oxidant_weights, fuel_weights, of_ratio, weights);

    cea_detonation_solver_solve(solver, soln, weights, T0, p0, FALSE);

    cea_real temperature, pressure, velocity, mach, sonic_velocity, gamma, enthalpy;
    cea_detonation_solution_get_property(soln, CEA_DETONATION_TEMPERATURE, 1, &temperature);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_PRESSURE, 1, &pressure);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_VELOCITY, 1, &velocity);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_MACH, 1, &mach);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_SONIC_VELOCITY, 1, &sonic_velocity);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_GAMMA, 1, &gamma);
    cea_detonation_solution_get_property(soln, CEA_DETONATION_ENTHALPY, 1, &enthalpy);

    printf(
        "%10s %10s %12s %8s %12s %10s %8s\n",
        "T (K)", "P (bar)", "Velocity (m/s)", "Mach", "Sonic Velocity (m/s)", "Gamma", "H/R"
    );

    printf(
        "%10.3f  %9.3f  %13.3f  %7.4f  %19.3f  %9.4f  %6.3f\n",
        temperature, pressure, velocity, mach, sonic_velocity, gamma, enthalpy
    );

    //----------------------------------------------------------------
    // CEA Cleanup
    //----------------------------------------------------------------
    cea_detonation_solution_destroy(&soln);
    cea_detonation_solver_destroy(&solver);
    cea_mixture_destroy(&prod);
    cea_mixture_destroy(&reac);
    free(weights);

    return 0;

}
