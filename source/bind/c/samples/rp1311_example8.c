#include "stdlib.h"
#include "stdio.h"
#include "cea.h"

#define TRUE 1
#define FALSE 0
#define LEN(x)  (sizeof(x) / sizeof((x)[0]))
#define BAR     10000.00
#define R      8314.51

int main(void) {

    //------------------------------------------------------------------
    // Problem Specification
    //------------------------------------------------------------------

    // Species
    const cea_string reactants[]     = {"H2(L)", "O2(L)"};
    const cea_real reactant_temps[]  = {  20.27,   90.17};
    const cea_real fuel_weights[]    = {    1.0,     0.0};
    const cea_real oxidant_weights[] = {    0.0,     1.0};
    const int nr = LEN(reactants);

    // Products
    const cea_string omitted_products[] = {""};

    // Thermo States
    const cea_real pressures = 53.3172*BAR;
    const cea_real of_ratio  = 5.55157;


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
    cea_rocket_solver solver;
    cea_rocket_solver_create_with_reactants(&solver, prod, reac);

    // Solution
    const cea_int num_pts = 5;
    cea_rocket_solution soln;
    cea_rocket_solution_create(&soln, solver);


    //----------------------------------------------------------------
    // Rocket Solve
    //----------------------------------------------------------------

    cea_real* weights = calloc(nr, sizeof(cea_real));
    cea_mixture_of_ratio_to_weights(reac, LEN(reactants), oxidant_weights, fuel_weights, of_ratio, weights);

    const cea_real pip[] = { 10.0 };
    const cea_real subar[] = { 1.58 };
    const cea_real supar[] = { 25.0 };
    const cea_real pc = 53.3172;
    cea_real hc;

    cea_mixture_calc_property_multitemp(reac, CEA_ENTHALPY, LEN(reactants), weights, LEN(reactants), reactant_temps, &hc);
    hc = hc/R;

    cea_rocket_solver_solve_iac(solver, soln, weights, pc, pip, 1, subar, 1, supar, 1, 0, hc, TRUE, 0.0, FALSE);

    cea_real* temperature = calloc(num_pts, sizeof(cea_real));
    cea_real* pressure = calloc(num_pts, sizeof(cea_real));
    cea_real* gamma = calloc(num_pts, sizeof(cea_real));
    cea_real* mw = calloc(num_pts, sizeof(cea_real));
    cea_real* mach = calloc(num_pts, sizeof(cea_real));
    cea_real* area_ratio = calloc(num_pts, sizeof(cea_real));
    cea_real* isp = calloc(num_pts, sizeof(cea_real));
    cea_real* cstar = calloc(num_pts, sizeof(cea_real));
    cea_real* cf = calloc(num_pts, sizeof(cea_real));
    cea_rocket_solution_get_property(soln, CEA_ROCKET_TEMPERATURE, num_pts, temperature);
    cea_rocket_solution_get_property(soln, CEA_ROCKET_PRESSURE, num_pts, pressure);
    cea_rocket_solution_get_property(soln, CEA_ROCKET_GAMMA_S, num_pts, gamma);
    cea_rocket_solution_get_property(soln, CEA_ROCKET_M, num_pts, mw);
    cea_rocket_solution_get_property(soln, CEA_MACH, num_pts, mach);
    cea_rocket_solution_get_property(soln, CEA_AE_AT, num_pts, area_ratio);
    cea_rocket_solution_get_property(soln, CEA_ISP, num_pts, isp);
    cea_rocket_solution_get_property(soln, CEA_C_STAR, num_pts, cstar);
    cea_rocket_solution_get_property(soln, CEA_COEFFICIENT_OF_THRUST, num_pts, cf);

    printf(
        "%10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
        "T (K)", "P (bar)", "gamma_s", "M (1/n)", "Mach", "Ae/At", "ISP", "C* (m/s)", "C_f"
    );

    for (int i = 0; i < num_pts; ++i) {
        printf(
            "%10.4f  %10.4f  %10.4f  %10.4f %10.4f  %10.4e  %10.4e  %10.4e  %10.4e\n",
            temperature[i], pressure[i], gamma[i], mw[i], mach[i],
            area_ratio[i], isp[i], cstar[i], cf[i]
        );
    }

    //----------------------------------------------------------------
    // CEA Cleanup
    //----------------------------------------------------------------
    cea_rocket_solution_destroy(&soln);
    cea_rocket_solver_destroy(&solver);
    cea_mixture_destroy(&prod);
    cea_mixture_destroy(&reac);
    free(weights);
    free(temperature);
    free(pressure);
    free(gamma);
    free(mw);
    free(mach);
    free(area_ratio);
    free(isp);
    free(cstar);
    free(cf);

    return 0;

}
