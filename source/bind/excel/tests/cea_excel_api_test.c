#include "cea_excel.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define LEN(x) ((int)(sizeof(x) / sizeof((x)[0])))

static int expect_true(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        return 1;
    }
    return 0;
}

static int expect_status_ok(int status, const char *message, const char *detail)
{
    if (status != CEA_EXCEL_SUCCESS) {
        fprintf(stderr, "%s: status=%d detail=%s\n", message, status, detail);
        return 1;
    }
    return 0;
}

static int test_helpers(void)
{
    const char *reactants = "H2\nAir";
    const double amounts[] = {1.0, 1.0};
    const int roles[] = {CEA_EXCEL_ROLE_FUEL, CEA_EXCEL_ROLE_OXIDIZER};
    const int bases[] = {CEA_EXCEL_BASIS_MOLE, CEA_EXCEL_BASIS_MOLE};
    double weights[2] = {0.0, 0.0};
    double of_ratio = 0.0;
    double phi = 0.0;
    double enthalpy = 0.0;
    double temperatures[] = {298.15, 298.15};
    char message[512];
    int status;

    status = cea_excel_of_from_equivalence(
        reactants, amounts, roles, bases, LEN(amounts), 1.0, &of_ratio, message, LEN(message));
    if (expect_status_ok(status, "CEA_OF_FROM_EQUIVALENCE", message)) return 1;
    if (expect_true(isfinite(of_ratio) && of_ratio > 0.0, "Expected finite positive O/F ratio")) return 1;

    status = cea_excel_weights_from_of(
        reactants, amounts, roles, bases, LEN(amounts), of_ratio, weights, message, LEN(message));
    if (expect_status_ok(status, "CEA_WEIGHTS_FROM_OF", message)) return 1;
    if (expect_true(isfinite(weights[0]) && isfinite(weights[1]), "Expected finite reactant weights")) return 1;

    status = cea_excel_phi_from_of(
        reactants, amounts, roles, bases, LEN(amounts), of_ratio, &phi, message, LEN(message));
    if (expect_status_ok(status, "CEA_PHI_FROM_OF", message)) return 1;
    if (expect_true(fabs(phi - 1.0) < 1.0e-10, "Expected phi to round-trip near 1")) return 1;

    status = cea_excel_calc_thermo(
        reactants, weights, LEN(weights), 6, temperatures, LEN(temperatures), 0.0, 0,
        &enthalpy, message, LEN(message));
    if (expect_status_ok(status, "CEA_CALC_THERMO", message)) return 1;
    if (expect_true(isfinite(enthalpy), "Expected finite reactant enthalpy")) return 1;

    return 0;
}

static int test_eq_solve(void)
{
    const char *reactants = "H2\nAir";
    const double amounts[] = {1.0, 1.0};
    const int roles[] = {CEA_EXCEL_ROLE_FUEL, CEA_EXCEL_ROLE_OXIDIZER};
    const int bases[] = {CEA_EXCEL_BASIS_MOLE, CEA_EXCEL_BASIS_MOLE};
    const double dummy[] = {0.0};
    double values[128];
    char headers[2048];
    char message[512];
    int nvalues = 0;
    int converged = 0;
    int status;

    status = cea_excel_eq_solve(
        0, reactants, amounts, roles, bases, LEN(amounts), dummy, 0,
        CEA_EXCEL_AMOUNT_R_EQ, 1.0, dummy, 0, 3000.0, 1.01325, "",
        "", "", "T\nP\nenthalpy\ncp_eq", "H2O\nOH", 0, 0, -1.0,
        values, LEN(values), headers, LEN(headers), &nvalues, &converged, message, LEN(message));
    if (expect_status_ok(status, "CEA_TP_SOLVE", message)) return 1;
    if (expect_true(converged != 0, "Expected TP solve to converge")) return 1;
    if (expect_true(nvalues >= 8, "Expected selected properties and species outputs")) return 1;
    if (expect_true(strstr(headers, "enthalpy") != NULL, "Expected enthalpy header")) return 1;
    if (expect_true(isfinite(values[0]), "Expected finite first solve output")) return 1;

    return 0;
}

static int test_rocket_solve(void)
{
    const char *reactants = "H2(L)\nO2(L)";
    const double amounts[] = {1.0, 1.0};
    const int roles[] = {CEA_EXCEL_ROLE_FUEL, CEA_EXCEL_ROLE_OXIDIZER};
    const int bases[] = {CEA_EXCEL_BASIS_WEIGHT, CEA_EXCEL_BASIS_WEIGHT};
    const double temps[] = {20.27, 90.17};
    const double pi_p[] = {10.0};
    const double subar[] = {1.58};
    const double supar[] = {25.0};
    const double dummy[] = {0.0};
    double values[256];
    char headers[4096];
    char message[512];
    int nvalues = 0;
    int converged = 0;
    int status;

    status = cea_excel_rocket_solve(
        0, reactants, amounts, roles, bases, LEN(amounts), dummy, 0,
        CEA_EXCEL_AMOUNT_OF, 5.55157, temps, LEN(temps), 53.3172,
        pi_p, LEN(pi_p), subar, LEN(subar), supar, LEN(supar), 0,
        0.0, 0, 0.0, 0, 0.0, 0, "", "", "T\nP\nIsp", "H2O", 0, 0,
        -1.0, values, LEN(values), headers, LEN(headers), &nvalues, &converged,
        message, LEN(message));
    if (expect_status_ok(status, "CEA_ROCKET_IAC_SOLVE", message)) return 1;
    if (expect_true(converged != 0, "Expected rocket solve to converge")) return 1;
    if (expect_true(nvalues > 3, "Expected rocket vector outputs")) return 1;
    if (expect_true(strstr(headers, "mass_H2O_chamber") != NULL, "Expected rocket mass fraction header")) return 1;
    if (expect_true(strstr(headers, "mole_H2O_chamber") != NULL, "Expected rocket mole fraction header")) return 1;

    return 0;
}

static int test_shock_solve(void)
{
    const char *reactants = "H2\nO2\nAr";
    const double amounts[] = {0.050, 0.050, 0.900};
    const int roles[] = {CEA_EXCEL_ROLE_FUEL, CEA_EXCEL_ROLE_OXIDIZER, CEA_EXCEL_ROLE_INERT};
    const int bases[] = {CEA_EXCEL_BASIS_MOLE, CEA_EXCEL_BASIS_MOLE, CEA_EXCEL_BASIS_MOLE};
    const double dummy[] = {0.0};
    double weights[3];
    double values[256];
    char headers[4096];
    char message[512];
    int nvalues = 0;
    int converged = 0;
    int status;

    status = cea_excel_weights_from_moles(reactants, amounts, LEN(amounts), weights, message, LEN(message));
    if (expect_status_ok(status, "CEA_WEIGHTS_FROM_MOLES", message)) return 1;

    status = cea_excel_shock_solve(
        reactants, amounts, roles, bases, LEN(amounts), weights, LEN(weights),
        CEA_EXCEL_AMOUNT_WEIGHTS, 0.0, 300.0, 0.1, 1400.0, 0, 1, 0, 0,
        "", "", "T\nP\nMach\nP21\nT21\nM21\nrho12\nv2\nu5_p_v2", "OH", 0, 0, -1.0, values, LEN(values),
        headers, LEN(headers), &nvalues, &converged, message, LEN(message));
    if (expect_status_ok(status, "CEA_SHOCK_SOLVE", message)) return 1;
    if (expect_true(converged != 0, "Expected shock solve to converge")) return 1;
    if (expect_true(nvalues == 21, "Expected vector, scalar, and species shock outputs")) return 1;
    if (expect_true(strstr(headers, "P21") != NULL, "Expected scalar shock property header")) return 1;
    if (expect_true(strstr(headers, "P21_state1") == NULL, "Scalar shock property must not have a station suffix")) return 1;
    if (expect_true(strstr(headers, "mass_OH_state1") != NULL, "Expected shock mass fraction header")) return 1;
    if (expect_true(strstr(headers, "mole_OH_state1") != NULL, "Expected shock mole fraction header")) return 1;

    (void)dummy;
    return 0;
}

static int test_detonation_solve(void)
{
    const char *reactants = "H2(L)\nO2(L)";
    const double amounts[] = {1.0, 1.0};
    const int roles[] = {CEA_EXCEL_ROLE_FUEL, CEA_EXCEL_ROLE_OXIDIZER};
    const int bases[] = {CEA_EXCEL_BASIS_WEIGHT, CEA_EXCEL_BASIS_WEIGHT};
    const double dummy[] = {0.0};
    double values[128];
    char headers[2048];
    char message[512];
    int nvalues = 0;
    int converged = 0;
    int status;

    status = cea_excel_detonation_solve(
        reactants, amounts, roles, bases, LEN(amounts), dummy, 0,
        CEA_EXCEL_AMOUNT_R_EQ, 1.0, 298.15, 1.0, 0, "", "",
        "T\nP\nvelocity\nMach", "H2O", 0, 0, -1.0, values, LEN(values),
        headers, LEN(headers), &nvalues, &converged, message, LEN(message));
    if (expect_status_ok(status, "CEA_DETONATION_SOLVE", message)) return 1;
    if (expect_true(converged != 0, "Expected detonation solve to converge")) return 1;
    if (expect_true(nvalues >= 6, "Expected detonation property and species outputs")) return 1;

    return 0;
}

int main(void)
{
    if (test_helpers()) return 1;
    if (test_eq_solve()) return 1;
    if (test_rocket_solve()) return 1;
    if (test_shock_solve()) return 1;
    if (test_detonation_solve()) return 1;
    return 0;
}
