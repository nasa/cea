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
    const cea_string reactants[]     = {"H2",  "O2",  "Ar" };
    const cea_real moles[]           = {0.050, 0.050, 0.900};
    const int nr = LEN(reactants);

    // Products
    const cea_string omitted_products[] = {""};

    // Thermo States
    const cea_real p0 = 0.1;
    const cea_real T0 = 300.0;

    // Shock Parameters
    const cea_real u1 = 1400.0;


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
    cea_shock_solver solver;
    cea_shock_solver_create_with_reactants(&solver, prod, reac);

    // Solution
    const cea_int num_pts = 3;
    cea_shock_solution soln;
    cea_shock_solution_create(&soln, num_pts);


    //----------------------------------------------------------------
    // Solve the shock problem
    //----------------------------------------------------------------

    cea_real* weights = (cea_real*) calloc(nr, sizeof(cea_real));
    cea_mixture_moles_to_weights(reac, LEN(reactants), moles, weights);

    cea_shock_solver_solve(solver, soln, weights, T0, p0, u1, FALSE, TRUE, FALSE, FALSE);

    cea_real* temperature = (cea_real*) calloc(num_pts, sizeof(cea_real));
    cea_real* pressure = (cea_real*) calloc(num_pts, sizeof(cea_real));
    cea_real* velocity = (cea_real*) calloc(num_pts, sizeof(cea_real));
    cea_real* mach = (cea_real*) calloc(num_pts, sizeof(cea_real));
    cea_real* v_sonic = (cea_real*) calloc(num_pts, sizeof(cea_real));
    cea_shock_solution_get_property(soln, CEA_SHOCK_TEMPERATURE, num_pts, temperature);
    cea_shock_solution_get_property(soln, CEA_SHOCK_PRESSURE, num_pts, pressure);
    cea_shock_solution_get_property(soln, CEA_SHOCK_VELOCITY, num_pts, velocity);
    cea_shock_solution_get_property(soln, CEA_SHOCK_MACH, num_pts, mach);
    cea_shock_solution_get_property(soln, CEA_SHOCK_SONIC_VELOCITY, num_pts, v_sonic);

    printf(
        "%14s %14s %14s %14s %14s \n",
        "T (K)", "P (bar)", "u (m/s)", "Mach", "Son. Vel. (m/s)"
    );

    for (int i = 0; i < num_pts; ++i) {
        printf(
            "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n",
            temperature[i], pressure[i], velocity[i], mach[i], v_sonic[i]
        );
    }

    //----------------------------------------------------------------
    // CEA Cleanup
    //----------------------------------------------------------------
    cea_shock_solution_destroy(&soln);
    cea_shock_solver_destroy(&solver);
    cea_mixture_destroy(&prod);
    cea_mixture_destroy(&reac);
    free(weights);
    free(temperature);
    free(pressure);
    free(velocity);
    free(mach);
    free(v_sonic);

    return 0;

}
