import numpy as np
import cea

"""
Example 1 from RP-1311
- TP equilibrium
"""

# Species
reac_names = ["H2", "Air"]
prod_names = ["Ar",   "C",   "CO",  "CO2", "H",
              "H2",   "H2O", "HNO", "HO2", "HNO2",
              "HNO3", "N",   "NH",  "NO",  "N2",
              "N2O3", "O",   "O2",  "OH",  "O3"]

# Thermo states
pressures = cea.units.atm_to_bar(np.array([1.0, 0.1, 0.01]))
temperatures = np.array([3000.0, 2000.0])

# Mixture states
fuel_moles = np.array([1.0, 0.0])
oxidant_moles = np.array([0.0, 1.0])
chem_eq_ratios = np.array([1.0, 1.5])

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(prod_names)

# Solver
solver = cea.EqSolver(prod, reactants=reac)
solution = cea.EqSolution(solver)

# Unit conversions
fuel_weights = reac.moles_to_weights(fuel_moles)
oxidant_weights = reac.moles_to_weights(oxidant_moles)
of_ratios = len(chem_eq_ratios)*[0.0]
for i, eqrat in enumerate(chem_eq_ratios):
    of_ratios[i] = reac.chem_eq_ratio_to_of_ratio(oxidant_weights,
                                                  fuel_weights,
                                                  eqrat)

n = len(chem_eq_ratios)*len(pressures)*len(temperatures)
of_ratio_out = np.zeros(n)
T_out = np.zeros(n)
P_out = np.zeros(n)
rho = np.zeros(n)
volume = np.zeros(n)
enthalpy = np.zeros(n)
energy = np.zeros(n)
gibbs = np.zeros(n)
entropy = np.zeros(n)
molecular_weight_M = np.zeros(n)
molecular_weight_MW = np.zeros(n)
gamma_s = np.zeros(n)
cp_eq = np.zeros(n)
cp_fr = np.zeros(n)
cv_eq = np.zeros(n)
cv_fr = np.zeros(n)
mole_fractions = {}
i = 0

# Equilibrium solve
for of_ratio in of_ratios:
    for p in pressures:
        for t in temperatures:
            weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
            solver.solve(solution, cea.TP, t, p, weights)

            # Store the output
            of_ratio_out[i] = of_ratio
            T_out[i] = t
            P_out[i] = cea.units.bar_to_atm(p)
            if solution.converged:
                rho[i] = solution.density*1.e-3
                volume[i] = solution.volume*1.e3
                enthalpy[i] = cea.units.joule_to_cal(solution.enthalpy)
                energy[i] = cea.units.joule_to_cal(solution.energy)
                gibbs[i] = cea.units.joule_to_cal(solution.gibbs_energy)
                entropy[i] = cea.units.joule_to_cal(solution.entropy)
                molecular_weight_M[i] = solution.M
                molecular_weight_MW[i] = solution.MW
                gamma_s[i] = solution.gamma_s
                cp_eq[i] = cea.units.joule_to_cal(solution.cp_eq)
                cp_fr[i] = cea.units.joule_to_cal(solution.cp_fr)
                cv_eq[i] = cea.units.joule_to_cal(solution.cv_eq)
                cv_fr[i] = cea.units.joule_to_cal(solution.cv_fr)

            if i == 0:
                for prod in solution.mole_fractions:
                    mole_fractions[prod] = np.array([solution.mole_fractions[prod]])
            else:
                for prod in mole_fractions:
                    mole_fractions[prod] = np.append(mole_fractions[prod], solution.mole_fractions[prod])

            i += 1

print("o/f             ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(of_ratio_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(of_ratio_out[i]))

print("P, atm          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(P_out[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(T_out[i]))

print("Density, g/cc   ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(rho[i]), end=" ")
    else:
        print("{0:10.3e}".format(rho[i]))

print("Volume, cc/g    ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(volume[i]), end=" ")
    else:
        print("{0:10.3e}".format(volume[i]))

print("H, cal/g        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(enthalpy[i]), end=" ")
    else:
        print("{0:10.3f}".format(enthalpy[i]))

print("U, cal/g        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(energy[i]), end=" ")
    else:
        print("{0:10.3f}".format(energy[i]))

print("G, cal/g        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.1f}".format(gibbs[i]), end=" ")
    else:
        print("{0:10.1f}".format(gibbs[i]))

print("S, cal/g-K      ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(entropy[i]), end=" ")
    else:
        print("{0:10.3f}".format(entropy[i]))

print("M, (1/n)        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(molecular_weight_M[i]), end=" ")
    else:
        print("{0:10.3f}".format(molecular_weight_M[i]))

print("MW              ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(molecular_weight_MW[i]), end=" ")
    else:
        print("{0:10.3f}".format(molecular_weight_MW[i]))

print("Gamma_s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(gamma_s[i]), end=" ")
    else:
        print("{0:10.4f}".format(gamma_s[i]))

print("Cp_eq, cal/g-K  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_eq[i]))

print("Cp_fr, cal/g-K  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_fr[i]))

print("Cv_eq, cal/g-K  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_eq[i]))

print("Cv_fr, cal/g-K  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_fr[i]))

print()
print("MOLE FRACTIONS")
print("")
trace_species = []
for prod in mole_fractions:
    if np.any(mole_fractions[prod] > 5e-6):
        print("{0:15s}".format(prod), end=" ")
        for j in range(n):
            if j < n-1:
                print("{0:10.5g}".format(mole_fractions[prod][j]), end=" ")
            else:
                print("{0:10.5g}".format(mole_fractions[prod][j]))
    else:
        trace_species.append(prod)

print()
print("TRACE SPECIES:")
max_cols = 10
nrows = (len(trace_species) + max_cols - 1) // max_cols
for i in range(nrows):
    print(" ".join("{0:10s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))
