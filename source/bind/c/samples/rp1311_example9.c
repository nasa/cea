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

    cea_set_log_level(CEA_LOG_NONE);
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
    cea_rocket_solution soln;
    cea_rocket_solution_create(&soln, solver);

    const cea_int num_pts;
    cea_rocket_solution_get_size(soln, &num_pts);

    //----------------------------------------------------------------
    // Rocket Solve
    //----------------------------------------------------------------

    cea_real* weights = calloc(nr, sizeof(cea_real));
    cea_mixture_of_ratio_to_weights(reac, LEN(reactants), oxidant_weights, fuel_weights, of_ratio, weights);

    const cea_real pip[] = { 10.0, 100.0, 1000.0 };
    const cea_real subar[] = {0.0};
    const cea_real supar[] = { 25.0, 50.0, 75.0 };
    const cea_real pc = 53.3172;
    const cea_real ac_at = 1.58;
    cea_real hc;

    cea_mixture_calc_property_multitemp(reac, CEA_ENTHALPY, LEN(reactants), weights, LEN(reactants), reactant_temps, &hc);
    hc = hc/R;

    cea_rocket_solver_solve_fac(solver, soln, weights, pc, pip, 3, subar, 0, supar, 3, 0, hc, TRUE, ac_at, FALSE, 0.0, FALSE);

    cea_real* temperature = calloc(num_pts, sizeof(cea_real));
    cea_real* pressure = calloc(num_pts, sizeof(cea_real));
    cea_real* gamma = calloc(num_pts, sizeof(cea_real));
    cea_real* mach = calloc(num_pts, sizeof(cea_real));
    cea_real* area_ratio = calloc(num_pts, sizeof(cea_real));
    cea_real* isp = calloc(num_pts, sizeof(cea_real));
    cea_real* isp_vac = calloc(num_pts, sizeof(cea_real));
    cea_real* cstar = calloc(num_pts, sizeof(cea_real));
    cea_real* cf = calloc(num_pts, sizeof(cea_real));
    cea_rocket_solution_get_property(soln, CEA_ROCKET_TEMPERATURE, num_pts, temperature);
    cea_rocket_solution_get_property(soln, CEA_ROCKET_PRESSURE, num_pts, pressure);
    cea_rocket_solution_get_property(soln, CEA_ROCKET_GAMMA_S, num_pts, gamma);
    cea_rocket_solution_get_property(soln, CEA_MACH, num_pts, mach);
    cea_rocket_solution_get_property(soln, CEA_AE_AT, num_pts, area_ratio);
    cea_rocket_solution_get_property(soln, CEA_ISP, num_pts, isp);
    cea_rocket_solution_get_property(soln, CEA_ISP_VACUUM, num_pts, isp_vac);
    cea_rocket_solution_get_property(soln, CEA_C_STAR, num_pts, cstar);
    cea_rocket_solution_get_property(soln, CEA_COEFFICIENT_OF_THRUST, num_pts, cf);

    printf("P, bar           ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.3f ", pressure[i]);
    }
    printf("\n");

    printf("T, K             ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", temperature[i]);
    }
    printf("\n");

    printf("gamma_s          ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", gamma[i]);
    }
    printf("\n");

    printf("Mach             ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", mach[i]);
    }
    printf("\n");

    printf("Ae/At            ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", area_ratio[i]);
    }
    printf("\n");

    printf("C*               ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", cstar[i]);
    }
    printf("\n");

    printf("Cf               ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", cf[i]);
    }
    printf("\n");

    printf("Isp              ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", isp[i]);
    }
    printf("\n");

    printf("Isp, vac.        ");
    for (int i = 0; i < num_pts; i++) {
        printf(" %9.2f ", isp_vac[i]);
    }
    printf("\n");

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
    free(mach);
    free(area_ratio);
    free(isp);
    free(isp_vac);
    free(cstar);
    free(cf);

    return 0;
}
