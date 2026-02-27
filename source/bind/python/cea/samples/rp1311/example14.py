import numpy as np
import cea
import os

"""
Example 14 from RP-1311
- TP equilibrium
- condensed species
"""

# Species
reac_names = ["H2(L)", "O2(L)"]

# Thermo states
p = cea.units.atm_to_bar(0.05)
temperatures = np.array([1000.0, 500.0, 350.0, 305.0, 304.3, 304.2, 304.0, 300.0])

# Mixture states
fuel_moles = np.array([100.0, 0.0])
oxidant_moles = np.array([0.0, 60.0])

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

# Solver
solver = cea.EqSolver(prod, reactants=reac)
solution = cea.EqSolution(solver)

# Unit conversions
fuel_weights = reac.moles_to_weights(fuel_moles)
oxidant_weights = reac.moles_to_weights(oxidant_moles)
weights = fuel_weights + oxidant_weights

n = len(temperatures)
T_out = np.zeros(n)
P_out = np.zeros(n)
rho = np.zeros(n)
enthalpy = np.zeros(n)
energy = np.zeros(n)
gibbs = np.zeros(n)
entropy = np.zeros(n)
molecular_weight_M = np.zeros(n)
molecular_weight_MW = np.zeros(n)
gamma_s = np.zeros(n)
cp_eq = np.zeros(n)
mole_fractions = {}
i = 0

# Equilibrium solve
for t in temperatures:
    solver.solve(solution, cea.TP, t, p, weights)

    # Store the output
    T_out[i] = t
    P_out[i] = p
    if solution.converged:
        rho[i] = solution.density
        enthalpy[i] = solution.enthalpy
        energy[i] = solution.energy
        gibbs[i] = solution.gibbs_energy
        entropy[i] = solution.entropy
        molecular_weight_M[i] = solution.M
        molecular_weight_MW[i] = solution.MW
        gamma_s[i] = solution.gamma_s
        cp_eq[i] = solution.cp_eq

        if i == 0:
            for prod in solution.mole_fractions:
                mole_fractions[prod] = np.array([solution.mole_fractions[prod]])
        else:
            for prod in mole_fractions:
                mole_fractions[prod] = np.append(mole_fractions[prod], solution.mole_fractions[prod])

        i += 1

print("P, bar          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.5f}".format(P_out[i]), end=" ")
    else:
        print("{0:10.5f}".format(P_out[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(T_out[i]))

print("Density, kg/m^3 ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(rho[i]), end=" ")
    else:
        print("{0:10.3e}".format(rho[i]))

print("H, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.2f}".format(enthalpy[i]), end=" ")
    else:
        print("{0:10.2f}".format(enthalpy[i]))

print("U, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.2f}".format(energy[i]), end=" ")
    else:
        print("{0:10.2f}".format(energy[i]))

print("G, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.1f}".format(gibbs[i]), end=" ")
    else:
        print("{0:10.1f}".format(gibbs[i]))

print("S, kJ/kg-K      ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(entropy[i]), end=" ")
    else:
        print("{0:10.4f}".format(entropy[i]))

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

print("Cp, kJ/kg-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_eq[i]))

print("Gamma_s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(gamma_s[i]), end=" ")
    else:
        print("{0:10.4f}".format(gamma_s[i]))

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
