#include "stdio.h"
#include "cea.h"

#define LEN(x) (sizeof(x) / sizeof((x)[0]))

int main(void)
{
  const cea_real pressures[] = {34.473652, 17.236826, 8.618413, 3.447365, 0.344737};
  const cea_real reactant_weights[] = {0.7206, 0.1858, 0.09, 0.002, 0.0016};
  const cea_real reactant_temps[] = {298.15, 298.15, 298.15, 298.15, 298.15};
  const cea_string omitted_products[] = {
      "CCN", "CNC", "C2N2", "C2O", "C3H4,allene", "C3H4,propyne", "C3H4,cyclo-", "C3",
      "C3H5,allyl", "C3H6,propylene", "C3H6,cyclo-", "C3H3,propargyl", "C3H6O", "C3H7,n-propyl",
      "C3H7,i-propyl", "Jet-A(g)", "C3O2", "C4", "C4H2", "C3H8O,2propanol", "C4H4,1,3-cyclo-",
      "C4H6,butadiene", "C4H6,2-butyne", "C3H8O,1propanol", "C4H8,tr2-butene", "C4H8,isobutene",
      "C4H8,cyclo-", "C4H6,cyclo-", "(CH3COOH)2", "C4H9,n-butyl", "C4H9,i-butyl", "C4H8,1-butene",
      "C4H9,s-butyl", "C4H9,t-butyl", "C4H10,isobutane", "C4H8,cis2-buten", "C4H10,n-butane",
      "C4N2", "C5", "C3H8", "Jet-A(L)", "C6H6(L)", "H2O(s)", "H2O(L)"};

  const cea_string binder_elements[] = {"C", "H", "O", "S"};
  const cea_real binder_coeffs[] = {1.0, 1.86955, 0.031256, 0.008415};
  const cea_string unit_enthalpy = "j/mole";
  const cea_string unit_temp = "k";

  cea_reactant_input reactants[] = {
      {.name = "NH4CLO4(I)"},
      {.name = "CHOS-Binder",
       .num_elements = 4,
       .elements = binder_elements,
       .coefficients = binder_coeffs,
       .has_enthalpy = true,
       .enthalpy = -12548.159088,
       .enthalpy_units = unit_enthalpy,
       .has_temperature = true,
       .temperature = 298.15,
       .temperature_units = unit_temp},
      {.name = "AL(cr)"},
      {.name = "MgO(cr)"},
      {.name = "H2O(L)"}};

  cea_set_log_level(CEA_LOG_NONE);
  cea_init();

  cea_mixture reac, prod;
  cea_mixture_create_from_input_reactants(&reac, LEN(reactants), reactants);
  cea_mixture_create_products_from_input_reactants(
      &prod, LEN(reactants), reactants, LEN(omitted_products), omitted_products);

  cea_eqsolver solver;
  cea_eqsolver_create_with_reactants(&solver, prod, reac);

  cea_eqsolution soln;
  cea_eqsolution_create(&soln, solver);

  cea_real h0;
  cea_mixture_calc_property_multitemp(
      reac, CEA_ENTHALPY, LEN(reactants), reactant_weights, LEN(reactants), reactant_temps, &h0);

  printf("%12s %12s %12s %12s\n", "P(bar)", "T(K)", "Gamma_s", "Cp(cal/g-K)");
  for (int i = 0; i < LEN(pressures); ++i)
  {
    cea_eqsolver_solve(solver, CEA_HP, h0 / 8314.51, pressures[i], reactant_weights, soln);
    cea_real temperature, gamma_s, cp_eq;
    cea_eqsolution_get_property(soln, CEA_TEMPERATURE, &temperature);
    cea_eqsolution_get_property(soln, CEA_GAMMA_S, &gamma_s);
    cea_eqsolution_get_property(soln, CEA_EQUILIBRIUM_CP, &cp_eq);
    cp_eq = cp_eq / 4.184;
    printf("%12.6f %12.3f %12.5f %12.5f\n", pressures[i], temperature, gamma_s, cp_eq);
  }

  cea_eqsolution_destroy(&soln);
  cea_eqsolver_destroy(&solver);
  cea_mixture_destroy(&prod);
  cea_mixture_destroy(&reac);
  return 0;
}
