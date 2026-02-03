#include "stdio.h"
#include "cea.h"

#define LEN(x)  (sizeof(x) / sizeof((x)[0]))
#define ATM     1.01325

int main(void) {

    //------------------------------------------------------------------
    // TP Problem Specification
    //------------------------------------------------------------------

    // Reactants
    const cea_string reactants[]   = { "H2", "Air" };
    const cea_real fuel_weights[]    = {  1.0,   0.0 };
    const cea_real oxidant_weights[] = {  0.0,   1.0 };

    // Products
    const cea_string products[] = {
        "Ar",    "C",    "CO",   "CO2",  "H",
        "H2",    "H2O",  "HNO",  "HO2",  "HNO2",
        "HNO3",  "N",    "NH",   "NO",   "N2",
        "N2O3",  "O",    "O2",   "OH",   "O3"
    };

    // Mixture States
    const cea_real densities[] = { 1.0e3*9.1864e-5, 1.0e3*8.0877e-6, 1.0e3*6.6054e-7 };
    const cea_real temperatures[] = { 3000.0 };
    const cea_real phi_eq_ratios[] = { 1.0 };


    //------------------------------------------------------------------
    // CEA Setup
    //------------------------------------------------------------------

    cea_set_log_level(CEA_LOG_DEBUG);
    cea_init();

    // Mixtures
    cea_mixture reac, prod;
    cea_mixture_create(&reac, LEN(reactants), reactants);
    cea_mixture_create(&prod, LEN(products),  products);

    // EqSolver
    cea_eqsolver solver;
    cea_solver_opts opts;
    cea_solver_opts_init(&opts);
    opts.reactants = reac;
    opts.transport = 1;

    cea_eqsolver_create_with_options(&solver, prod, opts);

    // EqSolution
    cea_eqsolution soln;
    cea_eqsolution_create(&soln, solver);

    // EqPartials
    cea_eqpartials partials;
    cea_eqpartials_create(&partials, solver);

    //------------------------------------------------------------------
    // Unit Conversions
    //------------------------------------------------------------------

    cea_real of_ratios[LEN(phi_eq_ratios)];
    for (int ir=0; ir < LEN(phi_eq_ratios); ++ir) {
        cea_mixture_chem_eq_ratio_to_of_ratio(reac, LEN(reactants),
            oxidant_weights,
            fuel_weights,
            phi_eq_ratios[ir],
            &of_ratios[ir]
        );
    }

    //------------------------------------------------------------------
    // Equilibrium Solve
    //------------------------------------------------------------------

    printf(
        "%10s  %10s  %10s  %12s  %12s  %12s\n",
        "O/F Ratio", "P (atm)", "T (K)",
        "H (cal/g)", "Cp (cal/g-K)", "Viscosity (mP)"
    );

    for (int ir=0; ir < LEN(of_ratios);    ++ir) {

        cea_real weights[LEN(reactants)];
        cea_mixture_of_ratio_to_weights(reac, LEN(reactants), oxidant_weights, fuel_weights, of_ratios[ir], weights);

        for (int ip=0; ip < LEN(densities);    ++ip) {
        for (int it=0; it < LEN(temperatures); ++it) {

            cea_eqsolver_solve_with_partials(solver, CEA_TV, temperatures[it], 1.0/densities[ip], weights, soln, partials);

            cea_real pressure, enthalpy, heat_capacity, viscosity;
            cea_eqsolution_get_property(soln, CEA_PRESSURE, &pressure);
            cea_eqsolution_get_property(soln, CEA_ENTHALPY, &enthalpy);
            cea_eqsolution_get_property(soln, CEA_EQUILIBRIUM_CP, &heat_capacity);
            cea_eqsolution_get_property(soln, CEA_VISCOSITY, &viscosity);
            pressure /= ATM;
            enthalpy /= 4.184;
            heat_capacity /= 4.184;

            printf(
                "%10.2f  %10.2f  %10.2f  %12.5e  %10.5f  %10.5f\n",
                of_ratios[ir], pressure, temperatures[it],
                enthalpy, heat_capacity, viscosity
            );

        }}
    }


    //----------------------------------------------------------------
    // CEA Cleanup
    //----------------------------------------------------------------
    cea_eqpartials_destroy(&partials);
    cea_eqsolution_destroy(&soln);
    cea_eqsolver_destroy(&solver);
    cea_mixture_destroy(&prod);
    cea_mixture_destroy(&reac);

    return 0;

}
