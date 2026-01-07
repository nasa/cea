import numpy as np
import cea

# Inputs
p0 = cea.units.atm_to_bar(1.0)  # Fixed-pressure state (bar)
T0 = 500.0   # Initial reactant temperature (K)
phi = 2.0       # Equivalence ratio

# Mixtures
reac = cea.Mixture(["H2", "O2"])
prod = cea.Mixture(["H", "H2", "H2O", "O", "O2", "OH"])
fuel_weights = np.array([100.0, 0.0])
oxidant_weights = np.array([0.0, 100.0])

# Solver
solver = cea.EqSolver(prod, reactants=reac)
solution = cea.EqSolution(solver)

# Get the weights for one mole of each species
of_ratio = reac.weight_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, phi)
weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

# Get the fixed-enthalpy
h0 = reac.calc_property(cea.ENTHALPY, weights, T0)/cea.R

# Equilibrium solve
solver.solve(solution, cea.HP, h0, p0, weights)

if solution.converged:
    rho = solution.density*1.e-3
    volume = solution.volume*1.e3
    enthalpy = solution.enthalpy
    energy = solution.energy
    gibbs = solution.gibbs_energy
    entropy = solution.entropy
    molecular_weight_M = solution.M
    molecular_weight_MW = solution.MW
    gamma_s = solution.gamma_s
    cp_eq = solution.cp_eq
    cp_fr = solution.cp_fr
    cv_eq = solution.cv_eq
    cv_fr = solution.cv_fr

mole_fractions = {}
for prod in solution.mole_fractions:
    mole_fractions[prod] = np.array([solution.mole_fractions[prod]])

print("o/f             ", end="")
print("{0:10.3f}".format(of_ratio))

print("P, bar          ", end="")
print("{0:10.3f}".format(p0))

print("T, K            ", end="")
print("{0:10.3f}".format(T0))

print("Density, kg/m^3  ", end="")
print("{0:10.3e}".format(rho))

print("Volume, m^3/kg   ", end="")
print("{0:10.3e}".format(volume))

print("H, kJ/kg        ", end="")
print("{0:10.3f}".format(enthalpy))

print("U, kJ/kg        ", end="")
print("{0:10.3f}".format(energy))

print("G, kJ/kg        ", end="")
print("{0:10.1f}".format(gibbs))

print("S, kJ/kg-K      ", end="")
print("{0:10.3f}".format(entropy))

print("M, (1/n)        ", end="")
print("{0:10.3f}".format(molecular_weight_M))

print("MW              ", end="")
print("{0:10.3f}".format(molecular_weight_MW))

print("Gamma_s         ", end="")
print("{0:10.4f}".format(gamma_s))

print("Cp_eq, kJ/kg-K  ", end="")
print("{0:10.4f}".format(cp_eq))

print("Cp_fr, kJ/kg-K  ", end="")
print("{0:10.4f}".format(cp_fr))

print("Cv_eq, kJ/kg-K  ", end="")
print("{0:10.4f}".format(cv_eq))

print("Cv_fr, kJ/kg-K  ", end="")
print("{0:10.4f}".format(cv_fr))

print()
print("MOLE FRACTIONS")
print("")
trace_species = []
for prod in mole_fractions:
    if np.any(mole_fractions[prod] > 5e-6):
        print("{0:15s}".format(prod), end=" ")
        print("{0:10.5g}".format(mole_fractions[prod][0]))
    else:
        trace_species.append(prod)

print()
print("TRACE SPECIES:")
max_cols = 10
nrows = (len(trace_species) + max_cols - 1) // max_cols
for i in range(nrows):
    print(" ".join("{0:10s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))
