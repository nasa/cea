import numpy as np
import cea

"""
Example 2 from RP-1311
- TV equilibrium
- transport properties
"""

# Options
transport = True

# Species
reac_names = ["H2", "Air"]
prod_names = ["Ar", "C", "CO", "CO2", "H", "H2", "H2O",
              "HNO", "HO2", "HNO2", "HNO3", "N", "NH",
              "NO", "N2", "N2O3", "O", "O2", "OH", "O3"]

# Thermo states
densities = 1.0e3*np.array([9.1864e-5, 8.0877e-6, 6.6054e-7])  # kg/m^3
temperatures = np.array([3000.0])
phis = np.array([1.0])
fuel_moles = np.array([1.0, 0.0])
oxidant_moles = np.array([0.0, 1.0])

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(prod_names)

# Solver
solver = cea.EqSolver(prod, reactants=reac, transport=transport)
solution = cea.EqSolution(solver)

# Unit conversions
fuel_weights = reac.moles_to_weights(fuel_moles)
oxidant_weights = reac.moles_to_weights(oxidant_moles)
of_ratios = len(phis)*[0.0]
for i, phi in enumerate(phis):
    of_ratios[i] = reac.weight_eq_ratio_to_of_ratio(oxidant_weights,
                                                    fuel_weights,
                                                    phi)

n = len(phis)*len(densities)*len(temperatures)
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
visc = np.zeros(n)
cond_fr = np.zeros(n)
cond_eq = np.zeros(n)
prandtl_fr = np.zeros(n)
prandtl_eq = np.zeros(n)
mole_fractions = {}
trace_species = []
i = 0

# Equilibrium solve
for of_ratio in of_ratios:
    for density in densities:
        for t in temperatures:
            weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
            solver.solve(solution, cea.TV, t, 1.0/density, weights)

            # Store the output
            of_ratio_out[i] = of_ratio
            T_out[i] = t
            if solution.converged:
                rho[i] = solution.density*1.e-3
                P_out[i] = cea.units.bar_to_atm(solution.P)
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
                visc[i] = solution.viscosity
                cond_fr[i] = cea.units.joule_to_cal(solution.conductivity_fr)
                cond_eq[i] = cea.units.joule_to_cal(solution.conductivity_eq)
                prandtl_fr[i] = solution.Pr_fr
                prandtl_eq[i] = solution.Pr_eq

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

print("")
print("Transport properties:")
print("")

print("Viscosity, mP   ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(visc[i]), end=" ")
    else:
        print("{0:10.4f}".format(visc[i]))

print("")
print("with equilibrium reaction:")
print("")

print("Cp, cal/g-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_eq[i]))

print("Cv, cal/g-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_eq[i]))

print("Conductivity    ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cond_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cond_eq[i]))

print("Prandtl number  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(prandtl_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(prandtl_eq[i]))

print("")
print("with frozen reaction:")
print("")

print("Cp, cal/g-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_fr[i]))

print("Cv, cal/g-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_fr[i]))

print("Conductivity    ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cond_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cond_fr[i]))

print("Prandtl number  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(prandtl_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(prandtl_fr[i]))

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
