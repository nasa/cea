#ifndef CEA_EXCEL_H
#define CEA_EXCEL_H

#ifdef _WIN32
    #ifdef CEA_EXCEL_BUILD
        #define CEA_EXCEL_API __declspec(dllexport)
    #else
        #define CEA_EXCEL_API __declspec(dllimport)
    #endif
#else
    #define CEA_EXCEL_API __attribute__((visibility("default")))
#endif

#define CEA_EXCEL_SUCCESS 0
#define CEA_EXCEL_ERROR_NULL_BUFFER 1
#define CEA_EXCEL_ERROR_INVALID_LENGTH 2
#define CEA_EXCEL_ERROR_TRUNCATED 3
#define CEA_EXCEL_ERROR_INVALID_INPUT 4
#define CEA_EXCEL_ERROR_ALLOCATION 5
#define CEA_EXCEL_ERROR_CEA 6
#define CEA_EXCEL_ERROR_OUTPUT_TOO_SMALL 7

#define CEA_EXCEL_ROLE_FUEL 1
#define CEA_EXCEL_ROLE_OXIDIZER 2
#define CEA_EXCEL_ROLE_INERT 3

#define CEA_EXCEL_BASIS_WEIGHT 1
#define CEA_EXCEL_BASIS_MOLE 2

#define CEA_EXCEL_AMOUNT_WEIGHTS 0
#define CEA_EXCEL_AMOUNT_OF 1
#define CEA_EXCEL_AMOUNT_PHI 2
#define CEA_EXCEL_AMOUNT_R_EQ 3
#define CEA_EXCEL_AMOUNT_PCT_FUEL 4

#ifdef __cplusplus
extern "C" {
#endif

CEA_EXCEL_API int cea_excel_version(char *version, int version_len);
CEA_EXCEL_API int cea_excel_test_add(int a, int b);
CEA_EXCEL_API int cea_excel_last_error(char *message, int message_len);
CEA_EXCEL_API int cea_excel_to_si(
    double value,
    const char *units,
    double *si_value,
    char *message,
    int message_len);

CEA_EXCEL_API int cea_excel_weights_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double weights[],
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_of_from_equivalence(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double equivalence,
    double *of_ratio,
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_of_from_phi(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double phi,
    double *of_ratio,
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_equivalence_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double *equivalence,
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_phi_from_of(
    const char *reactant_names,
    const double base_amounts[],
    const int roles[],
    const int bases[],
    int nreactants,
    double of_ratio,
    double *phi,
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_weights_from_moles(
    const char *reactant_names,
    const double moles[],
    int nreactants,
    double weights[],
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_moles_from_weights(
    const char *reactant_names,
    const double weights[],
    int nreactants,
    double moles[],
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_per_mole_from_per_weight(
    const char *reactant_names,
    const double per_weight[],
    int nreactants,
    double per_mole[],
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_per_weight_from_per_mole(
    const char *reactant_names,
    const double per_mole[],
    int nreactants,
    double per_weight[],
    char *message,
    int message_len);
CEA_EXCEL_API int cea_excel_calc_thermo(
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
    int message_len);

CEA_EXCEL_API int cea_excel_eq_solve(
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
    int message_len);
CEA_EXCEL_API int cea_excel_rocket_solve(
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
    int message_len);
CEA_EXCEL_API int cea_excel_shock_solve(
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
    int message_len);
CEA_EXCEL_API int cea_excel_detonation_solve(
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
    int message_len);

#ifdef __cplusplus
}
#endif

#endif
