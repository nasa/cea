#include "cea.h"

#include <stddef.h>

int main(void)
{
    const cea_string reactants_names[] = {"N2 ", "CH4"};
    const cea_string omitted[] = {""};
    cea_mixture reac = NULL;
    cea_mixture prod = NULL;
    cea_eqsolver solver = NULL;
    cea_solver_opts opts;

    if (cea_init() != CEA_SUCCESS) return 1;
    if (cea_mixture_create_w_ions(&reac, 2, reactants_names) != CEA_SUCCESS) return 2;
    if (cea_mixture_create_from_reactants_w_ions(&prod, 2, reactants_names, 0, omitted) != CEA_SUCCESS) return 3;

    // No reactants supplied: opts.reactants must be nullified rather than dereferenced.
    cea_solver_opts_init(&opts);
    if (cea_eqsolver_create_with_options(&solver, prod, opts) != CEA_SUCCESS) return 4;
    cea_eqsolver_destroy(&solver);

    // Negative ninsert is rejected.
    cea_solver_opts_init(&opts);
    opts.reactants = reac;
    opts.ninsert = -1;
    if (cea_eqsolver_create_with_options(&solver, prod, opts) != CEA_INVALID_SIZE) return 5;

    // ninsert > 0 requires a non-null insert array.
    cea_solver_opts_init(&opts);
    opts.reactants = reac;
    opts.ninsert = 2;
    opts.insert = NULL;
    if (cea_eqsolver_create_with_options(&solver, prod, opts) != CEA_INVALID_SIZE) return 6;

    // Insert species names longer than species_name_len are rejected.
    {
        const cea_string insert[] = {"THIS_NAME_IS_WAY_TOO_LONG_FOR_A_SPECIES"};
        cea_solver_opts_init(&opts);
        opts.reactants = reac;
        opts.ninsert = 1;
        opts.insert = insert;
        if (cea_eqsolver_create_with_options(&solver, prod, opts) != CEA_INVALID_SIZE) return 7;
    }

    cea_mixture_destroy(&prod);
    cea_mixture_destroy(&reac);

    return 0;
}
