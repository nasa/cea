#include "stdlib.h"
#include "stdio.h"
#include "cea.h"

#define LEN(x) (sizeof(x) / sizeof((x)[0]))
#define NUM_PTS 2

int main(void) {
    const cea_string reactants[] = {"N2 ", "CH4"};
    const cea_real moles[] = {0.985, 0.015};
    const cea_string omitted_products[] = {""};
    const cea_int nr = LEN(reactants);
    const cea_real p0 = 0.0002986421052631579;
    const cea_real T0 = 291.0;
    const cea_real u1 = 1900.0;

    int status = EXIT_FAILURE;
    cea_err ierr;
    cea_mixture reac = NULL;
    cea_mixture prod = NULL;
    cea_shock_solver solver = NULL;
    cea_shock_solution soln = NULL;
    cea_real *weights = NULL;
    cea_real pressure[NUM_PTS];
    cea_real temperature[NUM_PTS];
    int converged = 1;
    cea_solver_opts opts;

    cea_set_log_level(CEA_LOG_NONE);
    ierr = cea_init();
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "cea_init failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_mixture_create_w_ions(&reac, nr, reactants);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "reactant create failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_mixture_create_from_reactants_w_ions(&prod, nr, reactants, 0, omitted_products);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "product create failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_solver_opts_init(&opts);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "opts init failed: %d\n", ierr);
        goto cleanup;
    }
    opts.trace = 1.0e-10;
    opts.ions = true;
    opts.transport = true;
    opts.reactants = reac;

    ierr = cea_shock_solver_create_with_options(&solver, prod, opts);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "solver create failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_shock_solution_create(&soln, NUM_PTS);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "solution create failed: %d\n", ierr);
        goto cleanup;
    }

    weights = (cea_real *) calloc(nr, sizeof(cea_real));
    if (weights == NULL) {
        fprintf(stderr, "weights allocation failed\n");
        goto cleanup;
    }

    ierr = cea_mixture_moles_to_weights(reac, nr, moles, weights);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "weights conversion failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_shock_solver_solve(solver, soln, weights, T0, p0, u1, false, false, false, false);
    if (ierr != CEA_LAST_VALID_SOLUTION) {
        fprintf(stderr, "unexpected solve status: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_shock_solution_get_converged(soln, &converged);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "get_converged failed: %d\n", ierr);
        goto cleanup;
    }
    if (converged != 0) {
        fprintf(stderr, "expected non-converged shock status\n");
        goto cleanup;
    }

    ierr = cea_shock_solution_get_property(soln, CEA_SHOCK_PRESSURE, NUM_PTS, pressure);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "pressure query failed: %d\n", ierr);
        goto cleanup;
    }

    ierr = cea_shock_solution_get_property(soln, CEA_SHOCK_TEMPERATURE, NUM_PTS, temperature);
    if (ierr != CEA_SUCCESS) {
        fprintf(stderr, "temperature query failed: %d\n", ierr);
        goto cleanup;
    }

    if (pressure[1] <= 0.0 || temperature[1] <= 0.0) {
        fprintf(stderr, "expected retained incident state to remain populated\n");
        goto cleanup;
    }

    status = EXIT_SUCCESS;

cleanup:
    if (soln != NULL) cea_shock_solution_destroy(&soln);
    if (solver != NULL) cea_shock_solver_destroy(&solver);
    if (prod != NULL) cea_mixture_destroy(&prod);
    if (reac != NULL) cea_mixture_destroy(&reac);
    free(weights);
    return status;
}
