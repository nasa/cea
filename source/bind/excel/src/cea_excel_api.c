#include "cea_excel.h"

#include "cea.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <strings.h>
#endif

#ifdef _WIN32
#define excel_strcasecmp _stricmp
#else
#define excel_strcasecmp strcasecmp
#endif

#define CEA_EXCEL_GAS_CONSTANT 8314.5100

typedef struct {
    char **items;
    int count;
} excel_string_list;

typedef struct {
    const char *name;
    int type;
    int transport;
} excel_property;

static const excel_property eq_properties[] = {
    {"T", CEA_TEMPERATURE, 0},
    {"P", CEA_PRESSURE, 0},
    {"volume", CEA_VOLUME, 0},
    {"density", CEA_DENSITY, 0},
    {"M", CEA_M, 0},
    {"MW", CEA_MW, 0},
    {"enthalpy", CEA_ENTHALPY, 0},
    {"energy", CEA_ENERGY, 0},
    {"entropy", CEA_ENTROPY, 0},
    {"gibbs_energy", CEA_GIBBS_ENERGY, 0},
    {"gamma_s", CEA_GAMMA_S, 0},
    {"cp_fr", CEA_FROZEN_CP, 0},
    {"cv_fr", CEA_FROZEN_CV, 0},
    {"cp_eq", CEA_EQUILIBRIUM_CP, 0},
    {"cv_eq", CEA_EQUILIBRIUM_CV, 0},
    {"viscosity", CEA_VISCOSITY, 1},
    {"conductivity_fr", CEA_FROZEN_CONDUCTIVITY, 1},
    {"conductivity_eq", CEA_EQUILIBRIUM_CONDUCTIVITY, 1},
    {"Pr_fr", CEA_FROZEN_PRANDTL, 1},
    {"Pr_eq", CEA_EQUILIBRIUM_PRANDTL, 1},
};

static const excel_property rocket_properties[] = {
    {"T", CEA_ROCKET_TEMPERATURE, 0},
    {"P", CEA_ROCKET_PRESSURE, 0},
    {"volume", CEA_ROCKET_VOLUME, 0},
    {"density", CEA_ROCKET_DENSITY, 0},
    {"M", CEA_ROCKET_M, 0},
    {"MW", CEA_ROCKET_MW, 0},
    {"enthalpy", CEA_ROCKET_ENTHALPY, 0},
    {"energy", CEA_ROCKET_ENERGY, 0},
    {"entropy", CEA_ROCKET_ENTROPY, 0},
    {"gibbs_energy", CEA_ROCKET_GIBBS_ENERGY, 0},
    {"gamma_s", CEA_ROCKET_GAMMA_S, 0},
    {"cp_fr", CEA_ROCKET_FROZEN_CP, 0},
    {"cv_fr", CEA_ROCKET_FROZEN_CV, 0},
    {"cp_eq", CEA_ROCKET_EQUILIBRIUM_CP, 0},
    {"cv_eq", CEA_ROCKET_EQUILIBRIUM_CV, 0},
    {"Mach", CEA_MACH, 0},
    {"sonic_velocity", CEA_SONIC_VELOCITY, 0},
    {"ae_at", CEA_AE_AT, 0},
    {"c_star", CEA_C_STAR, 0},
    {"coefficient_of_thrust", CEA_COEFFICIENT_OF_THRUST, 0},
    {"Isp", CEA_ISP, 0},
    {"Isp_vacuum", CEA_ISP_VACUUM, 0},
    {"viscosity", CEA_ROCKET_VISCOSITY, 1},
    {"conductivity_fr", CEA_ROCKET_FROZEN_CONDUCTIVITY, 1},
    {"conductivity_eq", CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY, 1},
    {"Pr_fr", CEA_ROCKET_FROZEN_PRANDTL, 1},
    {"Pr_eq", CEA_ROCKET_EQUILIBRIUM_PRANDTL, 1},
};

static const excel_property shock_properties[] = {
    {"T", CEA_SHOCK_TEMPERATURE, 0},
    {"P", CEA_SHOCK_PRESSURE, 0},
    {"velocity", CEA_SHOCK_VELOCITY, 0},
    {"Mach", CEA_SHOCK_MACH, 0},
    {"sonic_velocity", CEA_SHOCK_SONIC_VELOCITY, 0},
    {"rho12", CEA_SHOCK_RHO12, 0},
    {"rho52", CEA_SHOCK_RHO52, 0},
    {"P21", CEA_SHOCK_P21, 0},
    {"P52", CEA_SHOCK_P52, 0},
    {"T21", CEA_SHOCK_T21, 0},
    {"T52", CEA_SHOCK_T52, 0},
    {"M21", CEA_SHOCK_M21, 0},
    {"M52", CEA_SHOCK_M52, 0},
    {"v2", CEA_SHOCK_V2, 0},
    {"u5_p_v2", CEA_SHOCK_U5_P_V2, 0},
    {"volume", CEA_SHOCK_VOLUME, 0},
    {"density", CEA_SHOCK_DENSITY, 0},
    {"M", CEA_SHOCK_M, 0},
    {"MW", CEA_SHOCK_MW, 0},
    {"enthalpy", CEA_SHOCK_ENTHALPY, 0},
    {"energy", CEA_SHOCK_ENERGY, 0},
    {"entropy", CEA_SHOCK_ENTROPY, 0},
    {"gibbs_energy", CEA_SHOCK_GIBBS_ENERGY, 0},
    {"gamma_s", CEA_SHOCK_GAMMA_S, 0},
    {"cp_fr", CEA_SHOCK_FROZEN_CP, 0},
    {"cv_fr", CEA_SHOCK_FROZEN_CV, 0},
    {"cp_eq", CEA_SHOCK_EQUILIBRIUM_CP, 0},
    {"cv_eq", CEA_SHOCK_EQUILIBRIUM_CV, 0},
    {"viscosity", CEA_SHOCK_VISCOSITY, 1},
    {"conductivity_fr", CEA_SHOCK_FROZEN_CONDUCTIVITY, 1},
    {"conductivity_eq", CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY, 1},
    {"Pr_fr", CEA_SHOCK_FROZEN_PRANDTL, 1},
    {"Pr_eq", CEA_SHOCK_EQUILIBRIUM_PRANDTL, 1},
};

static const excel_property detonation_properties[] = {
    {"P1", CEA_DETONATION_P1, 0},
    {"T1", CEA_DETONATION_T1, 0},
    {"H1", CEA_DETONATION_H1, 0},
    {"M1", CEA_DETONATION_M1, 0},
    {"gamma1", CEA_DETONATION_GAMMA1, 0},
    {"sonic_velocity1", CEA_DETONATION_V_SONIC1, 0},
    {"P", CEA_DETONATION_PRESSURE, 0},
    {"T", CEA_DETONATION_TEMPERATURE, 0},
    {"density", CEA_DETONATION_DENSITY, 0},
    {"enthalpy", CEA_DETONATION_ENTHALPY, 0},
    {"energy", CEA_DETONATION_ENERGY, 0},
    {"gibbs_energy", CEA_DETONATION_GIBBS_ENERGY, 0},
    {"entropy", CEA_DETONATION_ENTROPY, 0},
    {"Mach", CEA_DETONATION_MACH, 0},
    {"velocity", CEA_DETONATION_VELOCITY, 0},
    {"sonic_velocity", CEA_DETONATION_SONIC_VELOCITY, 0},
    {"gamma_s", CEA_DETONATION_GAMMA, 0},
    {"P_P1", CEA_DETONATION_P_P1, 0},
    {"T_T1", CEA_DETONATION_T_T1, 0},
    {"M_M1", CEA_DETONATION_M_M1, 0},
    {"rho_rho1", CEA_DETONATION_RHO_RHO1, 0},
    {"cp_fr", CEA_DETONATION_FROZEN_CP, 0},
    {"cv_fr", CEA_DETONATION_FROZEN_CV, 0},
    {"cp_eq", CEA_DETONATION_EQUILIBRIUM_CP, 0},
    {"cv_eq", CEA_DETONATION_EQUILIBRIUM_CV, 0},
    {"M", CEA_DETONATION_M, 0},
    {"MW", CEA_DETONATION_MW, 0},
    {"viscosity", CEA_DETONATION_VISCOSITY, 1},
    {"conductivity_fr", CEA_DETONATION_FROZEN_CONDUCTIVITY, 1},
    {"conductivity_eq", CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY, 1},
    {"Pr_fr", CEA_DETONATION_FROZEN_PRANDTL, 1},
    {"Pr_eq", CEA_DETONATION_EQUILIBRIUM_PRANDTL, 1},
};

static void excel_set_message(char *message, int message_len, const char *fmt, ...)
{
    va_list args;

    if (message == NULL || message_len <= 0) {
        return;
    }

    va_start(args, fmt);
    vsnprintf(message, (size_t)message_len, fmt, args);
    va_end(args);
    message[message_len - 1] = '\0';
}

static int excel_catch_cea_error(cea_err ierr, char *message, int message_len)
{
    if (ierr == CEA_SUCCESS) {
        excel_set_message(message, message_len, "OK");
        return CEA_EXCEL_SUCCESS;
    }

    if (ierr == CEA_FORTRAN_ABORT) {
        cea_excel_last_error(message, message_len);
    } else {
        excel_set_message(message, message_len, "CEA error %d", (int)ierr);
    }
    return CEA_EXCEL_ERROR_CEA;
}

static int excel_ensure_initialized(char *message, int message_len)
{
    cea_int initialized = 0;
    cea_err ierr = cea_is_initialized(&initialized);

    if (ierr != CEA_SUCCESS) {
        return excel_catch_cea_error(ierr, message, message_len);
    }
    if (initialized) {
        return CEA_EXCEL_SUCCESS;
    }
    return excel_catch_cea_error(cea_init(), message, message_len);
}

static int excel_is_separator(char c)
{
    return c == '\n' || c == '\r' || c == '\t' || c == ',' || c == ';';
}

static char *excel_strdup_range(const char *start, size_t len)
{
    char *out = (char *)malloc(len + 1U);
    if (out == NULL) {
        return NULL;
    }
    memcpy(out, start, len);
    out[len] = '\0';
    return out;
}

static void excel_trim_span(const char **start, const char **end)
{
    while (*start < *end && isspace((unsigned char)**start)) {
        (*start)++;
    }
    while (*end > *start && isspace((unsigned char)*((*end) - 1))) {
        (*end)--;
    }
}

static void excel_free_string_list(excel_string_list *list)
{
    int i;

    if (list == NULL) {
        return;
    }
    for (i = 0; i < list->count; i++) {
        free(list->items[i]);
    }
    free(list->items);
    list->items = NULL;
    list->count = 0;
}

static int excel_parse_string_list(
    const char *text,
    excel_string_list *list,
    char *message,
    int message_len)
{
    const char *p;
    const char *token_start;
    const char *token_end;
    int count = 0;
    int idx = 0;

    list->items = NULL;
    list->count = 0;
    if (text == NULL || text[0] == '\0') {
        return CEA_EXCEL_SUCCESS;
    }

    p = text;
    while (*p != '\0') {
        while (*p != '\0' && excel_is_separator(*p)) {
            p++;
        }
        token_start = p;
        while (*p != '\0' && !excel_is_separator(*p)) {
            p++;
        }
        token_end = p;
        excel_trim_span(&token_start, &token_end);
        if (token_end > token_start) {
            count++;
        }
    }

    if (count == 0) {
        return CEA_EXCEL_SUCCESS;
    }

    list->items = (char **)calloc((size_t)count, sizeof(char *));
    if (list->items == NULL) {
        excel_set_message(message, message_len, "Failed to allocate string list");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }

    p = text;
    while (*p != '\0') {
        while (*p != '\0' && excel_is_separator(*p)) {
            p++;
        }
        token_start = p;
        while (*p != '\0' && !excel_is_separator(*p)) {
            p++;
        }
        token_end = p;
        excel_trim_span(&token_start, &token_end);
        if (token_end > token_start) {
            list->items[idx] = excel_strdup_range(token_start, (size_t)(token_end - token_start));
            if (list->items[idx] == NULL) {
                excel_free_string_list(list);
                excel_set_message(message, message_len, "Failed to allocate string token");
                return CEA_EXCEL_ERROR_ALLOCATION;
            }
            idx++;
        }
    }
    list->count = count;
    return CEA_EXCEL_SUCCESS;
}

static const char **excel_const_items(excel_string_list *list)
{
    return (const char **)list->items;
}

static int excel_create_mixture(
    const char *reactant_names,
    int ions,
    cea_mixture *mix,
    excel_string_list *names,
    char *message,
    int message_len)
{
    int status;
    cea_err ierr;

    status = excel_parse_string_list(reactant_names, names, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    if (names->count <= 0) {
        excel_set_message(message, message_len, "No reactant names were provided");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }
    ierr = ions ? cea_mixture_create_w_ions(mix, names->count, excel_const_items(names)) :
                  cea_mixture_create(mix, names->count, excel_const_items(names));
    return excel_catch_cea_error(ierr, message, message_len);
}

static int excel_create_products(
    excel_string_list *reactants,
    const char *only_names,
    const char *omit_names,
    int ions,
    cea_mixture *products,
    excel_string_list *only,
    excel_string_list *omit,
    char *message,
    int message_len)
{
    int status;
    cea_err ierr;

    status = excel_parse_string_list(only_names, only, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_parse_string_list(omit_names, omit, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    if (only->count > 0) {
        ierr = ions ? cea_mixture_create_w_ions(products, only->count, excel_const_items(only)) :
                      cea_mixture_create(products, only->count, excel_const_items(only));
    } else {
        ierr = ions ? cea_mixture_create_from_reactants_w_ions(
                          products, reactants->count, excel_const_items(reactants),
                          omit->count, excel_const_items(omit)) :
                      cea_mixture_create_from_reactants(
                          products, reactants->count, excel_const_items(reactants),
                          omit->count, excel_const_items(omit));
    }
    return excel_catch_cea_error(ierr, message, message_len);
}

static int excel_convert_base_weights(
    cea_mixture mix,
    const double base_amounts[],
    const int bases[],
    int nreactants,
    double weights[],
    char *message,
    int message_len)
{
    int i;
    double *moles;
    double *mole_weights;
    cea_err ierr;

    if (base_amounts == NULL || bases == NULL || weights == NULL || nreactants <= 0) {
        excel_set_message(message, message_len, "Invalid reactant amount inputs");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    moles = (double *)calloc((size_t)nreactants, sizeof(double));
    mole_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (moles == NULL || mole_weights == NULL) {
        free(moles);
        free(mole_weights);
        excel_set_message(message, message_len, "Failed to allocate reactant amount buffers");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }

    for (i = 0; i < nreactants; i++) {
        weights[i] = 0.0;
        if (bases[i] == CEA_EXCEL_BASIS_MOLE) {
            moles[i] = base_amounts[i];
        } else if (bases[i] == CEA_EXCEL_BASIS_WEIGHT) {
            weights[i] = base_amounts[i];
        } else {
            free(moles);
            free(mole_weights);
            excel_set_message(message, message_len, "Invalid reactant basis at row %d", i + 1);
            return CEA_EXCEL_ERROR_INVALID_INPUT;
        }
    }

    ierr = cea_mixture_moles_to_weights(mix, nreactants, moles, mole_weights);
    if (ierr != CEA_SUCCESS) {
        free(moles);
        free(mole_weights);
        return excel_catch_cea_error(ierr, message, message_len);
    }
    for (i = 0; i < nreactants; i++) {
        if (bases[i] == CEA_EXCEL_BASIS_MOLE) {
            weights[i] = mole_weights[i];
        }
    }

    free(moles);
    free(mole_weights);
    return CEA_EXCEL_SUCCESS;
}

static int excel_partition_base_weights(
    cea_mixture mix,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double fuel_weights[],
    double oxidant_weights[],
    double inert_weights[],
    char *message,
    int message_len)
{
    int i;
    int has_fuel = 0;
    int has_oxidizer = 0;
    double *base_weights;
    int status;

    if (roles == NULL) {
        excel_set_message(message, message_len, "Reactant roles are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    base_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (base_weights == NULL) {
        excel_set_message(message, message_len, "Failed to allocate base weights");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }
    status = excel_convert_base_weights(
        mix, base_amounts, bases, nreactants, base_weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        free(base_weights);
        return status;
    }

    for (i = 0; i < nreactants; i++) {
        fuel_weights[i] = 0.0;
        oxidant_weights[i] = 0.0;
        inert_weights[i] = 0.0;
        if (roles[i] == CEA_EXCEL_ROLE_FUEL) {
            fuel_weights[i] = base_weights[i];
            has_fuel = 1;
        } else if (roles[i] == CEA_EXCEL_ROLE_OXIDIZER) {
            oxidant_weights[i] = base_weights[i];
            has_oxidizer = 1;
        } else if (roles[i] == CEA_EXCEL_ROLE_INERT) {
            inert_weights[i] = base_weights[i];
        } else {
            free(base_weights);
            excel_set_message(message, message_len, "Invalid reactant role at row %d", i + 1);
            return CEA_EXCEL_ERROR_INVALID_INPUT;
        }
    }
    free(base_weights);

    if (!has_fuel || !has_oxidizer) {
        excel_set_message(message, message_len, "At least one fuel and one oxidizer are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }
    return CEA_EXCEL_SUCCESS;
}

static void excel_add_inerts(double weights[], const double inert_weights[], int nreactants)
{
    int i;

    if (inert_weights == NULL) {
        return;
    }
    for (i = 0; i < nreactants; i++) {
        weights[i] += inert_weights[i];
    }
}

static int excel_resolve_weights(
    cea_mixture mix,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    const double explicit_weights[],
    int explicit_weights_len,
    int amount_mode,
    double amount_value,
    double weights[],
    char *message,
    int message_len)
{
    int i;
    int status;
    double of_ratio = amount_value;
    double *fuel_weights;
    double *oxidant_weights;
    double *inert_weights;
    cea_err ierr;

    if (amount_mode == CEA_EXCEL_AMOUNT_WEIGHTS) {
        if (explicit_weights == NULL || explicit_weights_len < nreactants) {
            excel_set_message(message, message_len, "Explicit weights are missing or too short");
            return CEA_EXCEL_ERROR_INVALID_INPUT;
        }
        for (i = 0; i < nreactants; i++) {
            weights[i] = explicit_weights[i];
        }
        return CEA_EXCEL_SUCCESS;
    }

    fuel_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    oxidant_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    inert_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (fuel_weights == NULL || oxidant_weights == NULL || inert_weights == NULL) {
        free(fuel_weights);
        free(oxidant_weights);
        free(inert_weights);
        excel_set_message(message, message_len, "Failed to allocate split reactant weights");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }

    status = excel_partition_base_weights(
        mix, base_amounts, roles, bases, nreactants,
        fuel_weights, oxidant_weights, inert_weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        free(fuel_weights);
        free(oxidant_weights);
        free(inert_weights);
        return status;
    }

    if (amount_mode == CEA_EXCEL_AMOUNT_PHI) {
        ierr = cea_mixture_weight_eq_ratio_to_of_ratio(
            mix, nreactants, oxidant_weights, fuel_weights, amount_value, &of_ratio);
    } else if (amount_mode == CEA_EXCEL_AMOUNT_R_EQ) {
        ierr = cea_mixture_chem_eq_ratio_to_of_ratio(
            mix, nreactants, oxidant_weights, fuel_weights, amount_value, &of_ratio);
    } else if (amount_mode == CEA_EXCEL_AMOUNT_PCT_FUEL) {
        if (amount_value <= 0.0 || amount_value >= 100.0) {
            free(fuel_weights);
            free(oxidant_weights);
            free(inert_weights);
            excel_set_message(message, message_len, "pct_fuel must be between 0 and 100");
            return CEA_EXCEL_ERROR_INVALID_INPUT;
        }
        of_ratio = (100.0 - amount_value) / amount_value;
        ierr = CEA_SUCCESS;
    } else if (amount_mode == CEA_EXCEL_AMOUNT_OF) {
        ierr = CEA_SUCCESS;
    } else {
        free(fuel_weights);
        free(oxidant_weights);
        free(inert_weights);
        excel_set_message(message, message_len, "Invalid amount mode %d", amount_mode);
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }
    if (ierr != CEA_SUCCESS) {
        free(fuel_weights);
        free(oxidant_weights);
        free(inert_weights);
        return excel_catch_cea_error(ierr, message, message_len);
    }

    ierr = cea_mixture_of_ratio_to_weights(
        mix, nreactants, oxidant_weights, fuel_weights, of_ratio, weights);
    if (ierr == CEA_SUCCESS) {
        excel_add_inerts(weights, inert_weights, nreactants);
    }

    free(fuel_weights);
    free(oxidant_weights);
    free(inert_weights);
    return excel_catch_cea_error(ierr, message, message_len);
}

static int excel_append_output(
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    const char *header,
    double value)
{
    int used;
    int written;

    if (*nvalues >= values_cap) {
        return CEA_EXCEL_ERROR_OUTPUT_TOO_SMALL;
    }
    values[*nvalues] = value;

    if (headers != NULL && headers_len > 0) {
        used = (int)strlen(headers);
        written = snprintf(
            headers + used,
            (size_t)(headers_len - used),
            "%s%s",
            used == 0 ? "" : "\t",
            header);
        if (written < 0 || written >= headers_len - used) {
            return CEA_EXCEL_ERROR_OUTPUT_TOO_SMALL;
        }
    }

    (*nvalues)++;
    return CEA_EXCEL_SUCCESS;
}

static const excel_property *excel_find_property(
    const excel_property properties[],
    int nproperties,
    const char *name)
{
    int i;

    for (i = 0; i < nproperties; i++) {
        if (excel_strcasecmp(properties[i].name, name) == 0) {
            return &properties[i];
        }
    }
    return NULL;
}

static int excel_selected_properties(
    const excel_property *defaults[],
    int ndefaults,
    const excel_property all_properties[],
    int nall_properties,
    const char *property_names,
    int transport,
    const excel_property ***selected,
    int *nselected,
    char *message,
    int message_len)
{
    int status;
    int i;
    excel_string_list names = {0};
    const excel_property **out;

    *selected = NULL;
    *nselected = 0;
    status = excel_parse_string_list(property_names, &names, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }

    if (names.count > 0) {
        out = (const excel_property **)calloc((size_t)names.count, sizeof(excel_property *));
        if (out == NULL) {
            excel_free_string_list(&names);
            excel_set_message(message, message_len, "Failed to allocate property list");
            return CEA_EXCEL_ERROR_ALLOCATION;
        }
        for (i = 0; i < names.count; i++) {
            out[i] = excel_find_property(all_properties, nall_properties, names.items[i]);
            if (out[i] == NULL) {
                free(out);
                excel_set_message(message, message_len, "Unknown property '%s'", names.items[i]);
                excel_free_string_list(&names);
                return CEA_EXCEL_ERROR_INVALID_INPUT;
            }
            if (out[i]->transport && !transport) {
                out[i] = NULL;
            }
        }
        *nselected = names.count;
        *selected = out;
    } else {
        out = (const excel_property **)calloc((size_t)nall_properties, sizeof(excel_property *));
        if (out == NULL) {
            excel_free_string_list(&names);
            excel_set_message(message, message_len, "Failed to allocate default property list");
            return CEA_EXCEL_ERROR_ALLOCATION;
        }
        for (i = 0; i < ndefaults; i++) {
            out[*nselected] = defaults[i];
            (*nselected)++;
        }
        if (transport) {
            for (i = 0; i < nall_properties; i++) {
                if (all_properties[i].transport) {
                    out[*nselected] = &all_properties[i];
                    (*nselected)++;
                }
            }
        }
        *selected = out;
    }

    excel_free_string_list(&names);
    return CEA_EXCEL_SUCCESS;
}

static const char *excel_rocket_station_name(int station)
{
    if (station == 0) {
        return "chamber";
    }
    if (station == 1) {
        return "throat";
    }
    return "exit";
}

static const char *excel_shock_station_name(int station)
{
    if (station == 0) {
        return "state1";
    }
    if (station == 1) {
        return "state2";
    }
    if (station == 2) {
        return "state5";
    }
    return "station";
}

static int excel_append_vector_outputs(
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    const char *property,
    const double vector[],
    int nvector,
    const char *(*station_name)(int),
    char *message,
    int message_len)
{
    int i;
    int status;
    char label[128];

    for (i = 0; i < nvector; i++) {
        if (station_name == excel_rocket_station_name && i >= 2) {
            snprintf(label, sizeof(label), "%s_%s_%d", property, station_name(i), i - 1);
        } else {
            snprintf(label, sizeof(label), "%s_%s", property, station_name(i));
        }
        status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, vector[i]);
        if (status != CEA_EXCEL_SUCCESS) {
            excel_set_message(message, message_len, "Output buffer is too small");
            return status;
        }
    }
    return CEA_EXCEL_SUCCESS;
}

static int excel_find_species_index(cea_mixture products, const char *name, int *index)
{
    int i;
    cea_int nspecies = 0;
    cea_int name_len = 0;
    char *buffer;

    *index = -1;
    if (cea_mixture_get_num_species(products, &nspecies) != CEA_SUCCESS ||
        cea_species_name_len(&name_len) != CEA_SUCCESS) {
        return CEA_EXCEL_ERROR_CEA;
    }

    buffer = (char *)malloc((size_t)name_len + 1U);
    if (buffer == NULL) {
        return CEA_EXCEL_ERROR_ALLOCATION;
    }
    for (i = 0; i < nspecies; i++) {
        if (cea_mixture_get_species_name_buf(&products, i, buffer, name_len + 1) != CEA_SUCCESS) {
            free(buffer);
            return CEA_EXCEL_ERROR_CEA;
        }
        if (excel_strcasecmp(buffer, name) == 0) {
            *index = i;
            free(buffer);
            return CEA_EXCEL_SUCCESS;
        }
    }
    free(buffer);
    return CEA_EXCEL_SUCCESS;
}

static double excel_state1_to_internal(int eq_type, double state1)
{
    if (eq_type == CEA_HP || eq_type == CEA_SP || eq_type == CEA_UV || eq_type == CEA_SV) {
        return state1 * 1000.0 / CEA_EXCEL_GAS_CONSTANT;
    }
    return state1;
}

int cea_excel_last_error(char *message, int message_len)
{
    cea_int last_len = 0;
    cea_err ierr;

    if (message == NULL) {
        return CEA_EXCEL_ERROR_NULL_BUFFER;
    }
    if (message_len <= 0) {
        return CEA_EXCEL_ERROR_INVALID_LENGTH;
    }

    message[0] = '\0';
    ierr = cea_last_error_message_len(&last_len);
    if (ierr != CEA_SUCCESS || last_len <= 0) {
        excel_set_message(message, message_len, "No CEA error message is available");
        return CEA_EXCEL_SUCCESS;
    }
    ierr = cea_last_error_message_buf(message, message_len);
    if (ierr != CEA_SUCCESS) {
        excel_set_message(message, message_len, "Failed to read CEA error message");
        return CEA_EXCEL_ERROR_CEA;
    }
    return CEA_EXCEL_SUCCESS;
}

int cea_excel_to_si(
    double value,
    const char *units,
    double *si_value,
    char *message,
    int message_len)
{
    char normalized[32];
    size_t i;
    size_t n;

    if (units == NULL || si_value == NULL) {
        excel_set_message(message, message_len, "Units and output value are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    n = strlen(units);
    if (n >= sizeof(normalized)) {
        excel_set_message(message, message_len, "Unit string is too long");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }
    for (i = 0; i < n; i++) {
        normalized[i] = (char)tolower((unsigned char)units[i]);
    }
    normalized[n] = '\0';

    if (strcmp(normalized, "k") == 0 || strcmp(normalized, "bar") == 0 ||
        strcmp(normalized, "j/mole") == 0 || strcmp(normalized, "kg/m**3") == 0 ||
        strcmp(normalized, "m**3/kg") == 0 || strcmp(normalized, "kj/kg") == 0 ||
        strcmp(normalized, "kj/kg-k") == 0) {
        *si_value = value;
    } else if (strcmp(normalized, "c") == 0) {
        *si_value = value + 273.15;
    } else if (strcmp(normalized, "f") == 0) {
        *si_value = (value - 32.0) / 1.8 + 273.15;
    } else if (strcmp(normalized, "r") == 0) {
        *si_value = value / 1.8;
    } else if (strcmp(normalized, "atm") == 0) {
        *si_value = value * 1.01325;
    } else if (strcmp(normalized, "psi") == 0 || strcmp(normalized, "psia") == 0) {
        *si_value = value * 0.0689475729;
    } else if (strcmp(normalized, "mmhg") == 0) {
        *si_value = value * 0.001333;
    } else if (strcmp(normalized, "pa") == 0) {
        *si_value = value * 1.0e-5;
    } else if (strcmp(normalized, "kpa") == 0) {
        *si_value = value * 0.01;
    } else if (strcmp(normalized, "mpa") == 0) {
        *si_value = value * 10.0;
    } else if (strcmp(normalized, "kj/mole") == 0) {
        *si_value = value * 1.0e3;
    } else if (strcmp(normalized, "kcal/mole") == 0) {
        *si_value = value * 4.184e3;
    } else if (strcmp(normalized, "cal/mole") == 0) {
        *si_value = value * 4.184;
    } else if (strcmp(normalized, "g/cm**3") == 0 || strcmp(normalized, "g/cc") == 0) {
        *si_value = value * 1.0e3;
    } else if (strcmp(normalized, "cm**3/g") == 0 || strcmp(normalized, "cc/g") == 0) {
        *si_value = value * 1.0e-3;
    } else {
        excel_set_message(message, message_len, "Unknown units '%s'", units);
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    excel_set_message(message, message_len, "OK");
    return CEA_EXCEL_SUCCESS;
}

static int excel_ratio_helper(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double input_value,
    double *output_value,
    int mode,
    char *message,
    int message_len)
{
    int status;
    double *fuel_weights;
    double *oxidant_weights;
    double *inert_weights;
    cea_mixture mix = NULL;
    excel_string_list names = {0};
    cea_err ierr = CEA_SUCCESS;

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, 0, &mix, &names, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        excel_free_string_list(&names);
        return status;
    }
    if (nreactants != names.count) {
        cea_mixture_destroy(&mix);
        excel_free_string_list(&names);
        excel_set_message(message, message_len, "Reactant count does not match name list");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    fuel_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    oxidant_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    inert_weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (fuel_weights == NULL || oxidant_weights == NULL || inert_weights == NULL) {
        free(fuel_weights);
        free(oxidant_weights);
        free(inert_weights);
        cea_mixture_destroy(&mix);
        excel_free_string_list(&names);
        excel_set_message(message, message_len, "Failed to allocate ratio helper weights");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }
    status = excel_partition_base_weights(
        mix, base_amounts, roles, bases, nreactants,
        fuel_weights, oxidant_weights, inert_weights, message, message_len);
    free(inert_weights);
    if (status == CEA_EXCEL_SUCCESS) {
        if (mode == CEA_EXCEL_AMOUNT_OF) {
            ierr = cea_mixture_of_ratio_to_weights(
                mix, nreactants, oxidant_weights, fuel_weights, input_value, output_value);
        } else if (mode == CEA_EXCEL_AMOUNT_R_EQ) {
            ierr = cea_mixture_chem_eq_ratio_to_of_ratio(
                mix, nreactants, oxidant_weights, fuel_weights, input_value, output_value);
        } else if (mode == CEA_EXCEL_AMOUNT_PHI) {
            ierr = cea_mixture_weight_eq_ratio_to_of_ratio(
                mix, nreactants, oxidant_weights, fuel_weights, input_value, output_value);
        } else if (mode == -CEA_EXCEL_AMOUNT_R_EQ) {
            ierr = cea_mixture_of_ratio_to_chem_eq_ratio(
                mix, nreactants, oxidant_weights, fuel_weights, input_value, output_value);
        } else if (mode == -CEA_EXCEL_AMOUNT_PHI) {
            ierr = cea_mixture_of_ratio_to_weight_eq_ratio(
                mix, nreactants, oxidant_weights, fuel_weights, input_value, output_value);
        }
        status = excel_catch_cea_error(ierr, message, message_len);
    }

    free(fuel_weights);
    free(oxidant_weights);
    cea_mixture_destroy(&mix);
    excel_free_string_list(&names);
    return status;
}

int cea_excel_weights_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double weights[],
    char *message,
    int message_len)
{
    return excel_ratio_helper(
        reactant_names, base_amounts, roles, bases, nreactants,
        of_ratio, weights, CEA_EXCEL_AMOUNT_OF, message, message_len);
}

int cea_excel_of_from_equivalence(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double equivalence,
    double *of_ratio,
    char *message,
    int message_len)
{
    return excel_ratio_helper(
        reactant_names, base_amounts, roles, bases, nreactants,
        equivalence, of_ratio, CEA_EXCEL_AMOUNT_R_EQ, message, message_len);
}

int cea_excel_of_from_phi(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double phi,
    double *of_ratio,
    char *message,
    int message_len)
{
    return excel_ratio_helper(
        reactant_names, base_amounts, roles, bases, nreactants,
        phi, of_ratio, CEA_EXCEL_AMOUNT_PHI, message, message_len);
}

int cea_excel_equivalence_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double *equivalence,
    char *message,
    int message_len)
{
    return excel_ratio_helper(
        reactant_names, base_amounts, roles, bases, nreactants,
        of_ratio, equivalence, -CEA_EXCEL_AMOUNT_R_EQ, message, message_len);
}

int cea_excel_phi_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double *phi,
    char *message,
    int message_len)
{
    return excel_ratio_helper(
        reactant_names, base_amounts, roles, bases, nreactants,
        of_ratio, phi, -CEA_EXCEL_AMOUNT_PHI, message, message_len);
}

static int excel_array_conversion(
    const char *reactant_names,
    const double input[],
    int nreactants,
    double output[],
    int conversion,
    char *message,
    int message_len)
{
    int status;
    cea_err ierr = CEA_SUCCESS;
    cea_mixture mix = NULL;
    excel_string_list names = {0};

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, 0, &mix, &names, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        excel_free_string_list(&names);
        return status;
    }
    if (nreactants != names.count) {
        cea_mixture_destroy(&mix);
        excel_free_string_list(&names);
        excel_set_message(message, message_len, "Reactant count does not match name list");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    if (conversion == 0) {
        ierr = cea_mixture_moles_to_weights(mix, nreactants, input, output);
    } else if (conversion == 1) {
        ierr = cea_mixture_weights_to_moles(mix, nreactants, input, output);
    } else if (conversion == 2) {
        ierr = cea_mixture_per_weight_to_per_mole(mix, nreactants, input, output);
    } else if (conversion == 3) {
        ierr = cea_mixture_per_mole_to_per_weight(mix, nreactants, input, output);
    }

    cea_mixture_destroy(&mix);
    excel_free_string_list(&names);
    return excel_catch_cea_error(ierr, message, message_len);
}

int cea_excel_weights_from_moles(
    const char *reactant_names,
    const double moles[],
    int nreactants,
    double weights[],
    char *message,
    int message_len)
{
    return excel_array_conversion(reactant_names, moles, nreactants, weights, 0, message, message_len);
}

int cea_excel_moles_from_weights(
    const char *reactant_names,
    const double weights[],
    int nreactants,
    double moles[],
    char *message,
    int message_len)
{
    return excel_array_conversion(reactant_names, weights, nreactants, moles, 1, message, message_len);
}

int cea_excel_per_mole_from_per_weight(
    const char *reactant_names,
    const double per_weight[],
    int nreactants,
    double per_mole[],
    char *message,
    int message_len)
{
    return excel_array_conversion(reactant_names, per_weight, nreactants, per_mole, 2, message, message_len);
}

int cea_excel_per_weight_from_per_mole(
    const char *reactant_names,
    const double per_mole[],
    int nreactants,
    double per_weight[],
    char *message,
    int message_len)
{
    return excel_array_conversion(reactant_names, per_mole, nreactants, per_weight, 3, message, message_len);
}

int cea_excel_calc_thermo(
    const char *reactant_names,
    const double weights[],
    int nreactants,
    int property_type,
    const double temperatures[],
    int ntemperatures,
    double pressure,
    int use_pressure,
    double *value,
    char *message,
    int message_len)
{
    int status;
    cea_err ierr;
    cea_mixture mix = NULL;
    excel_string_list names = {0};

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, 0, &mix, &names, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        excel_free_string_list(&names);
        return status;
    }
    if (nreactants != names.count || weights == NULL || temperatures == NULL || value == NULL) {
        cea_mixture_destroy(&mix);
        excel_free_string_list(&names);
        excel_set_message(message, message_len, "Invalid calc_thermo inputs");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    if (use_pressure) {
        if (ntemperatures == nreactants) {
            ierr = cea_mixture_calc_property_tp_multitemp(
                mix, property_type, nreactants, weights, ntemperatures, temperatures, pressure, value);
        } else {
            ierr = cea_mixture_calc_property_tp(
                mix, property_type, nreactants, weights, temperatures[0], pressure, value);
        }
    } else {
        if (ntemperatures == nreactants) {
            ierr = cea_mixture_calc_property_multitemp(
                mix, property_type, nreactants, weights, ntemperatures, temperatures, value);
        } else {
            ierr = cea_mixture_calc_property(
                mix, property_type, nreactants, weights, temperatures[0], value);
        }
    }

    cea_mixture_destroy(&mix);
    excel_free_string_list(&names);
    return excel_catch_cea_error(ierr, message, message_len);
}

static int excel_append_eq_species(
    cea_mixture products,
    cea_eqsolution soln,
    const char *species_names,
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    char *message,
    int message_len)
{
    int status;
    int i;
    int idx;
    cea_int np = 0;
    double *mass = NULL;
    double *mole = NULL;
    char label[160];
    excel_string_list species = {0};

    status = excel_parse_string_list(species_names, &species, message, message_len);
    if (status != CEA_EXCEL_SUCCESS || species.count == 0) {
        excel_free_string_list(&species);
        return status;
    }
    if (cea_mixture_get_num_species(products, &np) != CEA_SUCCESS) {
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }
    mass = (double *)calloc((size_t)np, sizeof(double));
    mole = (double *)calloc((size_t)np, sizeof(double));
    if (mass == NULL || mole == NULL) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        excel_set_message(message, message_len, "Failed to allocate species buffers");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }
    if (cea_eqsolution_get_species_amounts(soln, np, mass, true) != CEA_SUCCESS ||
        cea_eqsolution_get_species_amounts(soln, np, mole, false) != CEA_SUCCESS) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }

    for (i = 0; i < species.count; i++) {
        status = excel_find_species_index(products, species.items[i], &idx);
        if (status != CEA_EXCEL_SUCCESS || idx < 0) {
            free(mass);
            free(mole);
            excel_set_message(message, message_len, "Species '%s' is not in the product set", species.items[i]);
            excel_free_string_list(&species);
            return status == CEA_EXCEL_SUCCESS ? CEA_EXCEL_ERROR_INVALID_INPUT : status;
        }
        snprintf(label, sizeof(label), "mass_%s", species.items[i]);
        status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mass[idx]);
        if (status == CEA_EXCEL_SUCCESS) {
            snprintf(label, sizeof(label), "mole_%s", species.items[i]);
            status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mole[idx]);
        }
        if (status != CEA_EXCEL_SUCCESS) {
            free(mass);
            free(mole);
            excel_free_string_list(&species);
            excel_set_message(message, message_len, "Output buffer is too small");
            return status;
        }
    }

    free(mass);
    free(mole);
    excel_free_string_list(&species);
    return CEA_EXCEL_SUCCESS;
}

int cea_excel_eq_solve(
    int eq_type,
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    const double explicit_weights[],
    int explicit_weights_len,
    int amount_mode,
    double amount_value,
    const double reactant_temperatures[],
    int nreactant_temperatures,
    double state1,
    double state2,
    const char *only_names,
    const char *omit_names,
    const char *insert_names,
    const char *property_names,
    const char *species_names,
    int transport,
    int ions,
    double trace,
    double values[],
    int values_cap,
    char *headers,
    int headers_len,
    int *nvalues,
    int *converged,
    char *message,
    int message_len)
{
    const excel_property *defaults[] = {
        &eq_properties[0], &eq_properties[1], &eq_properties[5], &eq_properties[6],
        &eq_properties[8], &eq_properties[10], &eq_properties[13]
    };
    const excel_property **selected = NULL;
    int nselected = 0;
    int status;
    int i;
    double prop_value;
    double *weights = NULL;
    cea_err ierr;
    cea_err solve_ierr;
    cea_mixture reactants = NULL;
    cea_mixture products = NULL;
    cea_eqsolver solver = NULL;
    cea_eqsolution soln = NULL;
    cea_solver_opts opts;
    excel_string_list reactant_list = {0};
    excel_string_list only_list = {0};
    excel_string_list omit_list = {0};
    excel_string_list insert_list = {0};

    (void)reactant_temperatures;
    (void)nreactant_temperatures;

    if (headers != NULL && headers_len > 0) {
        headers[0] = '\0';
    }
    if (nvalues != NULL) {
        *nvalues = 0;
    }
    if (converged != NULL) {
        *converged = 0;
    }
    if (values == NULL || nvalues == NULL || converged == NULL) {
        excel_set_message(message, message_len, "Output buffers are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, ions, &reactants, &reactant_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    if (nreactants != reactant_list.count) {
        excel_set_message(message, message_len, "Reactant count does not match name list");
        status = CEA_EXCEL_ERROR_INVALID_INPUT;
        goto cleanup;
    }
    status = excel_create_products(
        &reactant_list, only_names, omit_names, ions, &products, &only_list, &omit_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (weights == NULL) {
        excel_set_message(message, message_len, "Failed to allocate solve weights");
        status = CEA_EXCEL_ERROR_ALLOCATION;
        goto cleanup;
    }
    status = excel_resolve_weights(
        reactants, base_amounts, roles, bases, nreactants, explicit_weights, explicit_weights_len,
        amount_mode, amount_value, weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    status = excel_parse_string_list(insert_names, &insert_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    cea_solver_opts_init(&opts);
    opts.reactants = reactants;
    opts.transport = transport ? true : false;
    opts.ions = ions ? true : false;
    opts.trace = trace;
    opts.ninsert = insert_list.count;
    opts.insert = excel_const_items(&insert_list);
    ierr = cea_eqsolver_create_with_options(&solver, products, opts);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }
    ierr = cea_eqsolution_create(&soln, solver);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }

    ierr = cea_eqsolver_solve(
        solver, (cea_equilibrium_type)eq_type, excel_state1_to_internal(eq_type, state1), state2, weights, soln);
    solve_ierr = ierr;
    *converged = 0;
    cea_eqsolution_get_converged(soln, converged);
    status = (int)ierr;
    if (ierr != CEA_SUCCESS && ierr != CEA_NOT_CONVERGED && ierr != CEA_LAST_VALID_SOLUTION) {
        excel_catch_cea_error(ierr, message, message_len);
        goto cleanup;
    }

    status = excel_selected_properties(
        defaults, (int)(sizeof(defaults) / sizeof(defaults[0])),
        eq_properties, (int)(sizeof(eq_properties) / sizeof(eq_properties[0])),
        property_names, transport, &selected, &nselected, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    for (i = 0; i < nselected; i++) {
        if (selected[i] == NULL) {
            continue;
        }
        ierr = cea_eqsolution_get_property(soln, (cea_property_type)selected[i]->type, &prop_value);
        if (ierr != CEA_SUCCESS) {
            status = excel_catch_cea_error(ierr, message, message_len);
            goto cleanup;
        }
        status = excel_append_output(
            values, values_cap, nvalues, headers, headers_len, selected[i]->name, prop_value);
        if (status != CEA_EXCEL_SUCCESS) {
            excel_set_message(message, message_len, "Output buffer is too small");
            goto cleanup;
        }
    }
    status = excel_append_eq_species(
        products, soln, species_names, values, values_cap, nvalues, headers, headers_len, message, message_len);
    if (status == CEA_EXCEL_SUCCESS) {
        excel_set_message(message, message_len,
            solve_ierr == CEA_SUCCESS ? "OK" : "CEA solve warning %d", (int)solve_ierr);
        status = (int)solve_ierr;
    }

cleanup:
    free(selected);
    free(weights);
    if (soln != NULL) {
        cea_eqsolution_destroy(&soln);
    }
    if (solver != NULL) {
        cea_eqsolver_destroy(&solver);
    }
    if (products != NULL) {
        cea_mixture_destroy(&products);
    }
    if (reactants != NULL) {
        cea_mixture_destroy(&reactants);
    }
    excel_free_string_list(&insert_list);
    excel_free_string_list(&omit_list);
    excel_free_string_list(&only_list);
    excel_free_string_list(&reactant_list);
    return status;
}

static int excel_append_rocket_species(
    cea_mixture products,
    cea_rocket_solution soln,
    int num_pts,
    const char *species_names,
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    char *message,
    int message_len)
{
    int status;
    int i;
    int station;
    int idx;
    cea_int np = 0;
    double *mass = NULL;
    double *mole = NULL;
    char label[192];
    excel_string_list species = {0};

    status = excel_parse_string_list(species_names, &species, message, message_len);
    if (status != CEA_EXCEL_SUCCESS || species.count == 0) {
        excel_free_string_list(&species);
        return status;
    }
    if (cea_mixture_get_num_species(products, &np) != CEA_SUCCESS) {
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }
    mass = (double *)calloc((size_t)np, sizeof(double));
    mole = (double *)calloc((size_t)np, sizeof(double));
    if (mass == NULL || mole == NULL) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        excel_set_message(message, message_len, "Failed to allocate species buffers");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }

    for (station = 0; station < num_pts; station++) {
        if (cea_rocket_solution_get_species_amounts(soln, np, station, mass, true) != CEA_SUCCESS ||
            cea_rocket_solution_get_species_amounts(soln, np, station, mole, false) != CEA_SUCCESS) {
            free(mass);
            free(mole);
            excel_free_string_list(&species);
            return CEA_EXCEL_ERROR_CEA;
        }
        for (i = 0; i < species.count; i++) {
            status = excel_find_species_index(products, species.items[i], &idx);
            if (status != CEA_EXCEL_SUCCESS || idx < 0) {
                free(mass);
                free(mole);
                excel_set_message(message, message_len, "Species '%s' is not in the product set", species.items[i]);
                excel_free_string_list(&species);
                return status == CEA_EXCEL_SUCCESS ? CEA_EXCEL_ERROR_INVALID_INPUT : status;
            }
            snprintf(label, sizeof(label), "mass_%s_%s", species.items[i], excel_rocket_station_name(station));
            if (station >= 2) {
                snprintf(label, sizeof(label), "mass_%s_exit_%d", species.items[i], station - 1);
            }
            status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mass[idx]);
            if (status == CEA_EXCEL_SUCCESS) {
                snprintf(label, sizeof(label), "mole_%s_%s", species.items[i], excel_rocket_station_name(station));
                if (station >= 2) {
                    snprintf(label, sizeof(label), "mole_%s_exit_%d", species.items[i], station - 1);
                }
                status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mole[idx]);
            }
            if (status != CEA_EXCEL_SUCCESS) {
                free(mass);
                free(mole);
                excel_free_string_list(&species);
                excel_set_message(message, message_len, "Output buffer is too small");
                return status;
            }
        }
    }

    free(mass);
    free(mole);
    excel_free_string_list(&species);
    return CEA_EXCEL_SUCCESS;
}

int cea_excel_rocket_solve(
    int finite_area,
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    const double explicit_weights[],
    int explicit_weights_len,
    int amount_mode,
    double amount_value,
    const double reactant_temperatures[],
    int nreactant_temperatures,
    double pc,
    const double pi_p[],
    int n_pi_p,
    const double subar[],
    int nsubar,
    const double supar[],
    int nsupar,
    int n_frz,
    double hc_or_tc,
    int use_hc,
    double mdot_or_acat,
    int use_mdot,
    double tc_est,
    int use_tc_est,
    const char *omit_names,
    const char *insert_names,
    const char *property_names,
    const char *species_names,
    int transport,
    int ions,
    double trace,
    double values[],
    int values_cap,
    char *headers,
    int headers_len,
    int *nvalues,
    int *converged,
    char *message,
    int message_len)
{
    const excel_property *defaults[] = {
        &rocket_properties[0], &rocket_properties[1], &rocket_properties[4],
        &rocket_properties[10], &rocket_properties[15], &rocket_properties[17],
        &rocket_properties[18], &rocket_properties[19], &rocket_properties[20],
        &rocket_properties[21]
    };
    const excel_property **selected = NULL;
    int nselected = 0;
    int status;
    int i;
    int num_pts = 0;
    double *weights = NULL;
    double *vector = NULL;
    double hc_or_tc_internal = hc_or_tc;
    cea_err ierr;
    cea_err solve_ierr;
    cea_mixture reactants = NULL;
    cea_mixture products = NULL;
    cea_rocket_solver solver = NULL;
    cea_rocket_solution soln = NULL;
    cea_solver_opts opts;
    excel_string_list reactant_list = {0};
    excel_string_list only_list = {0};
    excel_string_list omit_list = {0};
    excel_string_list insert_list = {0};

    if (headers != NULL && headers_len > 0) {
        headers[0] = '\0';
    }
    if (nvalues != NULL) {
        *nvalues = 0;
    }
    if (converged != NULL) {
        *converged = 0;
    }
    if (values == NULL || nvalues == NULL || converged == NULL) {
        excel_set_message(message, message_len, "Output buffers are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, ions, &reactants, &reactant_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    if (nreactants != reactant_list.count) {
        excel_set_message(message, message_len, "Reactant count does not match name list");
        status = CEA_EXCEL_ERROR_INVALID_INPUT;
        goto cleanup;
    }
    status = excel_create_products(
        &reactant_list, NULL, omit_names, ions, &products, &only_list, &omit_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (weights == NULL) {
        excel_set_message(message, message_len, "Failed to allocate solve weights");
        status = CEA_EXCEL_ERROR_ALLOCATION;
        goto cleanup;
    }
    status = excel_resolve_weights(
        reactants, base_amounts, roles, bases, nreactants, explicit_weights, explicit_weights_len,
        amount_mode, amount_value, weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    status = excel_parse_string_list(insert_names, &insert_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    cea_solver_opts_init(&opts);
    opts.reactants = reactants;
    opts.transport = transport ? true : false;
    opts.ions = ions ? true : false;
    opts.trace = trace;
    opts.ninsert = insert_list.count;
    opts.insert = excel_const_items(&insert_list);
    ierr = cea_rocket_solver_create_with_options(&solver, products, opts);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }
    ierr = cea_rocket_solution_create(&soln, solver);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }

    if (use_hc) {
        hc_or_tc_internal = hc_or_tc * 1000.0 / CEA_EXCEL_GAS_CONSTANT;
    } else if (nreactant_temperatures == nreactants && reactant_temperatures != NULL) {
        ierr = cea_mixture_calc_property_multitemp(
            reactants, CEA_ENTHALPY, nreactants, weights,
            nreactant_temperatures, reactant_temperatures, &hc_or_tc_internal);
        if (ierr != CEA_SUCCESS) {
            status = excel_catch_cea_error(ierr, message, message_len);
            goto cleanup;
        }
        hc_or_tc_internal = hc_or_tc_internal / CEA_EXCEL_GAS_CONSTANT;
        use_hc = 1;
    }

    ierr = finite_area ?
        cea_rocket_solver_solve_fac(
            solver, soln, weights, pc, pi_p, n_pi_p, subar, nsubar, supar, nsupar,
            n_frz, hc_or_tc_internal, use_hc ? true : false,
            mdot_or_acat, use_mdot ? true : false, tc_est, use_tc_est ? true : false) :
        cea_rocket_solver_solve_iac(
            solver, soln, weights, pc, pi_p, n_pi_p, subar, nsubar, supar, nsupar,
            n_frz, hc_or_tc_internal, use_hc ? true : false, tc_est, use_tc_est ? true : false);
    solve_ierr = ierr;
    *converged = 0;
    cea_rocket_solution_get_converged(soln, converged);
    status = (int)ierr;
    if (ierr != CEA_SUCCESS && ierr != CEA_NOT_CONVERGED) {
        excel_catch_cea_error(ierr, message, message_len);
        goto cleanup;
    }
    if (cea_rocket_solution_get_size(soln, &num_pts) != CEA_SUCCESS || num_pts <= 0) {
        excel_set_message(message, message_len, "Failed to query rocket solution size");
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }
    vector = (double *)calloc((size_t)num_pts, sizeof(double));
    if (vector == NULL) {
        excel_set_message(message, message_len, "Failed to allocate rocket property vector");
        status = CEA_EXCEL_ERROR_ALLOCATION;
        goto cleanup;
    }

    status = excel_selected_properties(
        defaults, (int)(sizeof(defaults) / sizeof(defaults[0])),
        rocket_properties, (int)(sizeof(rocket_properties) / sizeof(rocket_properties[0])),
        property_names, transport, &selected, &nselected, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    for (i = 0; i < nselected; i++) {
        if (selected[i] == NULL) {
            continue;
        }
        ierr = cea_rocket_solution_get_property(
            soln, (cea_rocket_property_type)selected[i]->type, num_pts, vector);
        if (ierr != CEA_SUCCESS) {
            status = excel_catch_cea_error(ierr, message, message_len);
            goto cleanup;
        }
        status = excel_append_vector_outputs(
            values, values_cap, nvalues, headers, headers_len,
            selected[i]->name, vector, num_pts, excel_rocket_station_name, message, message_len);
        if (status != CEA_EXCEL_SUCCESS) {
            goto cleanup;
        }
    }
    status = excel_append_rocket_species(
        products, soln, num_pts, species_names, values, values_cap, nvalues,
        headers, headers_len, message, message_len);
    if (status == CEA_EXCEL_SUCCESS) {
        excel_set_message(message, message_len,
            solve_ierr == CEA_SUCCESS ? "OK" : "CEA solve warning %d", (int)solve_ierr);
        status = (int)solve_ierr;
    }

cleanup:
    free(selected);
    free(vector);
    free(weights);
    if (soln != NULL) {
        cea_rocket_solution_destroy(&soln);
    }
    if (solver != NULL) {
        cea_rocket_solver_destroy(&solver);
    }
    if (products != NULL) {
        cea_mixture_destroy(&products);
    }
    if (reactants != NULL) {
        cea_mixture_destroy(&reactants);
    }
    excel_free_string_list(&insert_list);
    excel_free_string_list(&omit_list);
    excel_free_string_list(&only_list);
    excel_free_string_list(&reactant_list);
    return status;
}

static int excel_append_shock_species(
    cea_mixture products,
    cea_shock_solution soln,
    int num_pts,
    const char *species_names,
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    char *message,
    int message_len)
{
    int status;
    int i;
    int station;
    int idx;
    cea_int np = 0;
    double *mass = NULL;
    double *mole = NULL;
    char label[192];
    excel_string_list species = {0};

    status = excel_parse_string_list(species_names, &species, message, message_len);
    if (status != CEA_EXCEL_SUCCESS || species.count == 0) {
        excel_free_string_list(&species);
        return status;
    }
    if (cea_mixture_get_num_species(products, &np) != CEA_SUCCESS) {
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }
    mass = (double *)calloc((size_t)np, sizeof(double));
    mole = (double *)calloc((size_t)np, sizeof(double));
    if (mass == NULL || mole == NULL) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        excel_set_message(message, message_len, "Failed to allocate species buffers");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }

    for (station = 0; station < num_pts; station++) {
        if (cea_shock_solution_get_species_amounts(soln, np, station, mass, true) != CEA_SUCCESS ||
            cea_shock_solution_get_species_amounts(soln, np, station, mole, false) != CEA_SUCCESS) {
            free(mass);
            free(mole);
            excel_free_string_list(&species);
            return CEA_EXCEL_ERROR_CEA;
        }
        for (i = 0; i < species.count; i++) {
            status = excel_find_species_index(products, species.items[i], &idx);
            if (status != CEA_EXCEL_SUCCESS || idx < 0) {
                free(mass);
                free(mole);
                excel_set_message(message, message_len, "Species '%s' is not in the product set", species.items[i]);
                excel_free_string_list(&species);
                return status == CEA_EXCEL_SUCCESS ? CEA_EXCEL_ERROR_INVALID_INPUT : status;
            }
            snprintf(label, sizeof(label), "mass_%s_%s", species.items[i], excel_shock_station_name(station));
            status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mass[idx]);
            if (status == CEA_EXCEL_SUCCESS) {
                snprintf(label, sizeof(label), "mole_%s_%s", species.items[i], excel_shock_station_name(station));
                status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mole[idx]);
            }
            if (status != CEA_EXCEL_SUCCESS) {
                free(mass);
                free(mole);
                excel_free_string_list(&species);
                excel_set_message(message, message_len, "Output buffer is too small");
                return status;
            }
        }
    }

    free(mass);
    free(mole);
    excel_free_string_list(&species);
    return CEA_EXCEL_SUCCESS;
}

int cea_excel_shock_solve(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    const double explicit_weights[],
    int explicit_weights_len,
    int amount_mode,
    double amount_value,
    double T0,
    double p0,
    double mach1_or_u1,
    int use_mach,
    int reflected,
    int incident_frozen,
    int reflected_frozen,
    const char *omit_names,
    const char *insert_names,
    const char *property_names,
    const char *species_names,
    int transport,
    int ions,
    double trace,
    double values[],
    int values_cap,
    char *headers,
    int headers_len,
    int *nvalues,
    int *converged,
    char *message,
    int message_len)
{
    const excel_property *defaults[] = {
        &shock_properties[0], &shock_properties[1], &shock_properties[2],
        &shock_properties[3], &shock_properties[4], &shock_properties[5],
        &shock_properties[7], &shock_properties[9], &shock_properties[17],
        &shock_properties[23]
    };
    const excel_property **selected = NULL;
    int nselected = 0;
    int status;
    int i;
    int num_pts = reflected ? 3 : 1;
    double *weights = NULL;
    double *vector = NULL;
    cea_err ierr;
    cea_err solve_ierr;
    cea_mixture reactants = NULL;
    cea_mixture products = NULL;
    cea_shock_solver solver = NULL;
    cea_shock_solution soln = NULL;
    cea_solver_opts opts;
    excel_string_list reactant_list = {0};
    excel_string_list only_list = {0};
    excel_string_list omit_list = {0};
    excel_string_list insert_list = {0};

    if (headers != NULL && headers_len > 0) {
        headers[0] = '\0';
    }
    if (nvalues != NULL) {
        *nvalues = 0;
    }
    if (converged != NULL) {
        *converged = 0;
    }
    if (values == NULL || nvalues == NULL || converged == NULL) {
        excel_set_message(message, message_len, "Output buffers are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, ions, &reactants, &reactant_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    if (nreactants != reactant_list.count) {
        excel_set_message(message, message_len, "Reactant count does not match name list");
        status = CEA_EXCEL_ERROR_INVALID_INPUT;
        goto cleanup;
    }
    status = excel_create_products(
        &reactant_list, NULL, omit_names, ions, &products, &only_list, &omit_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    weights = (double *)calloc((size_t)nreactants, sizeof(double));
    vector = (double *)calloc((size_t)num_pts, sizeof(double));
    if (weights == NULL || vector == NULL) {
        excel_set_message(message, message_len, "Failed to allocate shock buffers");
        status = CEA_EXCEL_ERROR_ALLOCATION;
        goto cleanup;
    }
    status = excel_resolve_weights(
        reactants, base_amounts, roles, bases, nreactants, explicit_weights, explicit_weights_len,
        amount_mode, amount_value, weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    status = excel_parse_string_list(insert_names, &insert_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    cea_solver_opts_init(&opts);
    opts.reactants = reactants;
    opts.transport = transport ? true : false;
    opts.ions = ions ? true : false;
    opts.trace = trace;
    opts.ninsert = insert_list.count;
    opts.insert = excel_const_items(&insert_list);
    ierr = cea_shock_solver_create_with_options(&solver, products, opts);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }
    ierr = cea_shock_solution_create(&soln, num_pts);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }

    ierr = cea_shock_solver_solve(
        solver, soln, weights, T0, p0, mach1_or_u1, use_mach ? true : false,
        reflected ? true : false, incident_frozen ? true : false, reflected_frozen ? true : false);
    solve_ierr = ierr;
    *converged = 0;
    cea_shock_solution_get_converged(soln, converged);
    status = (int)ierr;
    if (ierr != CEA_SUCCESS && ierr != CEA_NOT_CONVERGED && ierr != CEA_LAST_VALID_SOLUTION) {
        excel_catch_cea_error(ierr, message, message_len);
        goto cleanup;
    }

    status = excel_selected_properties(
        defaults, (int)(sizeof(defaults) / sizeof(defaults[0])),
        shock_properties, (int)(sizeof(shock_properties) / sizeof(shock_properties[0])),
        property_names, transport, &selected, &nselected, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    for (i = 0; i < nselected; i++) {
        if (selected[i] == NULL) {
            continue;
        }
        ierr = cea_shock_solution_get_property(
            soln, (cea_shock_property_type)selected[i]->type, num_pts, vector);
        if (ierr != CEA_SUCCESS) {
            status = excel_catch_cea_error(ierr, message, message_len);
            goto cleanup;
        }
        status = excel_append_vector_outputs(
            values, values_cap, nvalues, headers, headers_len,
            selected[i]->name, vector, num_pts, excel_shock_station_name, message, message_len);
        if (status != CEA_EXCEL_SUCCESS) {
            goto cleanup;
        }
    }
    status = excel_append_shock_species(
        products, soln, num_pts, species_names, values, values_cap, nvalues,
        headers, headers_len, message, message_len);
    if (status == CEA_EXCEL_SUCCESS) {
        excel_set_message(message, message_len,
            solve_ierr == CEA_SUCCESS ? "OK" : "CEA solve warning %d", (int)solve_ierr);
        status = (int)solve_ierr;
    }

cleanup:
    free(selected);
    free(vector);
    free(weights);
    if (soln != NULL) {
        cea_shock_solution_destroy(&soln);
    }
    if (solver != NULL) {
        cea_shock_solver_destroy(&solver);
    }
    if (products != NULL) {
        cea_mixture_destroy(&products);
    }
    if (reactants != NULL) {
        cea_mixture_destroy(&reactants);
    }
    excel_free_string_list(&insert_list);
    excel_free_string_list(&omit_list);
    excel_free_string_list(&only_list);
    excel_free_string_list(&reactant_list);
    return status;
}

static int excel_append_detonation_species(
    cea_mixture products,
    cea_detonation_solution soln,
    const char *species_names,
    double values[],
    int values_cap,
    int *nvalues,
    char *headers,
    int headers_len,
    char *message,
    int message_len)
{
    int status;
    int i;
    int idx;
    cea_int np = 0;
    double *mass = NULL;
    double *mole = NULL;
    char label[160];
    excel_string_list species = {0};

    status = excel_parse_string_list(species_names, &species, message, message_len);
    if (status != CEA_EXCEL_SUCCESS || species.count == 0) {
        excel_free_string_list(&species);
        return status;
    }
    if (cea_mixture_get_num_species(products, &np) != CEA_SUCCESS) {
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }
    mass = (double *)calloc((size_t)np, sizeof(double));
    mole = (double *)calloc((size_t)np, sizeof(double));
    if (mass == NULL || mole == NULL) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        excel_set_message(message, message_len, "Failed to allocate species buffers");
        return CEA_EXCEL_ERROR_ALLOCATION;
    }
    if (cea_detonation_solution_get_species_amounts(soln, np, mass, true) != CEA_SUCCESS ||
        cea_detonation_solution_get_species_amounts(soln, np, mole, false) != CEA_SUCCESS) {
        free(mass);
        free(mole);
        excel_free_string_list(&species);
        return CEA_EXCEL_ERROR_CEA;
    }

    for (i = 0; i < species.count; i++) {
        status = excel_find_species_index(products, species.items[i], &idx);
        if (status != CEA_EXCEL_SUCCESS || idx < 0) {
            free(mass);
            free(mole);
            excel_set_message(message, message_len, "Species '%s' is not in the product set", species.items[i]);
            excel_free_string_list(&species);
            return status == CEA_EXCEL_SUCCESS ? CEA_EXCEL_ERROR_INVALID_INPUT : status;
        }
        snprintf(label, sizeof(label), "mass_%s", species.items[i]);
        status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mass[idx]);
        if (status == CEA_EXCEL_SUCCESS) {
            snprintf(label, sizeof(label), "mole_%s", species.items[i]);
            status = excel_append_output(values, values_cap, nvalues, headers, headers_len, label, mole[idx]);
        }
        if (status != CEA_EXCEL_SUCCESS) {
            free(mass);
            free(mole);
            excel_free_string_list(&species);
            excel_set_message(message, message_len, "Output buffer is too small");
            return status;
        }
    }

    free(mass);
    free(mole);
    excel_free_string_list(&species);
    return CEA_EXCEL_SUCCESS;
}

int cea_excel_detonation_solve(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    const double explicit_weights[],
    int explicit_weights_len,
    int amount_mode,
    double amount_value,
    double T1,
    double p1,
    int frozen,
    const char *omit_names,
    const char *insert_names,
    const char *property_names,
    const char *species_names,
    int transport,
    int ions,
    double trace,
    double values[],
    int values_cap,
    char *headers,
    int headers_len,
    int *nvalues,
    int *converged,
    char *message,
    int message_len)
{
    const excel_property *defaults[] = {
        &detonation_properties[6], &detonation_properties[7], &detonation_properties[8],
        &detonation_properties[13], &detonation_properties[14], &detonation_properties[15],
        &detonation_properties[16], &detonation_properties[17], &detonation_properties[18],
        &detonation_properties[20], &detonation_properties[25]
    };
    const excel_property **selected = NULL;
    int nselected = 0;
    int status;
    int i;
    double prop_value;
    double *weights = NULL;
    cea_err ierr;
    cea_err solve_ierr;
    cea_mixture reactants = NULL;
    cea_mixture products = NULL;
    cea_detonation_solver solver = NULL;
    cea_detonation_solution soln = NULL;
    cea_solver_opts opts;
    excel_string_list reactant_list = {0};
    excel_string_list only_list = {0};
    excel_string_list omit_list = {0};
    excel_string_list insert_list = {0};

    if (headers != NULL && headers_len > 0) {
        headers[0] = '\0';
    }
    if (nvalues != NULL) {
        *nvalues = 0;
    }
    if (converged != NULL) {
        *converged = 0;
    }
    if (values == NULL || nvalues == NULL || converged == NULL) {
        excel_set_message(message, message_len, "Output buffers are required");
        return CEA_EXCEL_ERROR_INVALID_INPUT;
    }

    status = excel_ensure_initialized(message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        return status;
    }
    status = excel_create_mixture(reactant_names, ions, &reactants, &reactant_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    if (nreactants != reactant_list.count) {
        excel_set_message(message, message_len, "Reactant count does not match name list");
        status = CEA_EXCEL_ERROR_INVALID_INPUT;
        goto cleanup;
    }
    status = excel_create_products(
        &reactant_list, NULL, omit_names, ions, &products, &only_list, &omit_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    weights = (double *)calloc((size_t)nreactants, sizeof(double));
    if (weights == NULL) {
        excel_set_message(message, message_len, "Failed to allocate solve weights");
        status = CEA_EXCEL_ERROR_ALLOCATION;
        goto cleanup;
    }
    status = excel_resolve_weights(
        reactants, base_amounts, roles, bases, nreactants, explicit_weights, explicit_weights_len,
        amount_mode, amount_value, weights, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }

    status = excel_parse_string_list(insert_names, &insert_list, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    cea_solver_opts_init(&opts);
    opts.reactants = reactants;
    opts.transport = transport ? true : false;
    opts.ions = ions ? true : false;
    opts.trace = trace;
    opts.ninsert = insert_list.count;
    opts.insert = excel_const_items(&insert_list);
    ierr = cea_detonation_solver_create_with_options(&solver, products, opts);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }
    ierr = cea_detonation_solution_create(&soln);
    if (excel_catch_cea_error(ierr, message, message_len) != CEA_EXCEL_SUCCESS) {
        status = CEA_EXCEL_ERROR_CEA;
        goto cleanup;
    }

    ierr = cea_detonation_solver_solve(solver, soln, weights, T1, p1, frozen ? true : false);
    solve_ierr = ierr;
    *converged = 0;
    cea_detonation_solution_get_converged(soln, converged);
    status = (int)ierr;
    if (ierr != CEA_SUCCESS && ierr != CEA_NOT_CONVERGED) {
        excel_catch_cea_error(ierr, message, message_len);
        goto cleanup;
    }

    status = excel_selected_properties(
        defaults, (int)(sizeof(defaults) / sizeof(defaults[0])),
        detonation_properties, (int)(sizeof(detonation_properties) / sizeof(detonation_properties[0])),
        property_names, transport, &selected, &nselected, message, message_len);
    if (status != CEA_EXCEL_SUCCESS) {
        goto cleanup;
    }
    for (i = 0; i < nselected; i++) {
        if (selected[i] == NULL) {
            continue;
        }
        ierr = cea_detonation_solution_get_property(
            soln, (cea_detonation_property_type)selected[i]->type, 1, &prop_value);
        if (ierr != CEA_SUCCESS) {
            status = excel_catch_cea_error(ierr, message, message_len);
            goto cleanup;
        }
        status = excel_append_output(
            values, values_cap, nvalues, headers, headers_len, selected[i]->name, prop_value);
        if (status != CEA_EXCEL_SUCCESS) {
            excel_set_message(message, message_len, "Output buffer is too small");
            goto cleanup;
        }
    }
    status = excel_append_detonation_species(
        products, soln, species_names, values, values_cap, nvalues, headers, headers_len, message, message_len);
    if (status == CEA_EXCEL_SUCCESS) {
        excel_set_message(message, message_len,
            solve_ierr == CEA_SUCCESS ? "OK" : "CEA solve warning %d", (int)solve_ierr);
        status = (int)solve_ierr;
    }

cleanup:
    free(selected);
    free(weights);
    if (soln != NULL) {
        cea_detonation_solution_destroy(&soln);
    }
    if (solver != NULL) {
        cea_detonation_solver_destroy(&solver);
    }
    if (products != NULL) {
        cea_mixture_destroy(&products);
    }
    if (reactants != NULL) {
        cea_mixture_destroy(&reactants);
    }
    excel_free_string_list(&insert_list);
    excel_free_string_list(&omit_list);
    excel_free_string_list(&only_list);
    excel_free_string_list(&reactant_list);
    return status;
}
