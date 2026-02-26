import numpy as np
import cea

"""
Example 13 from RP-1311
- IAC rocket problem
- uses BeO(L) as an insert
- set trace threshold
- condensed species
"""


# Define the reactants
reac_names = ["N2H4(L)", "Be(a)", "H2O2(L)"]
fuel_weights = np.array([0.8, 0.2, 0.0])
oxidant_weights = np.array([0.0, 0.0, 1.0])
reac_T = np.array([298.15, 298.15, 298.15])  # Reactant temperatures (K)
pct_fuel = 67.0
trace = 1e-10

# Rocket states
pc = cea.units.psi_to_bar(3000)  # Chamber pressure (bar)
pi_p = [3.0, 10.0, 30.0, 300.0]   # Pressure ratio

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)
insert = ["BeO(L)"]

# Solver
solver = cea.RocketSolver(prod, reactants=reac, trace=trace, insert=insert)
solution = cea.RocketSolution(solver)

# Get the reactant weights
of_ratio = (100.0 - pct_fuel) / pct_fuel
weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

# Get the chamber enthalpy
hc = reac.calc_property(cea.ENTHALPY, weights, reac_T)/cea.R

# Solve the IAC rocket problem
solver.solve(solution, weights, pc, pi_p, iac=True, hc=hc)

# Print the results
num_pts = solution.num_pts
T = solution.T
P = cea.units.bar_to_atm(solution.P)
rho = 1e-3*solution.density
enthalpy = cea.units.joule_to_cal(solution.enthalpy)
energy = cea.units.joule_to_cal(solution.energy)
gibbs = cea.units.joule_to_cal(solution.gibbs_energy)
entropy = cea.units.joule_to_cal(solution.entropy)
n = solution.n
M_1n = solution.M
MW = solution.MW
cp_eq = cea.units.joule_to_cal(solution.cp_eq)
Mach = solution.Mach
gamma_s = solution.gamma_s
v_sonic = solution.sonic_velocity
ae_at = solution.ae_at
c_star = cea.units.m_per_s_to_ft_per_s(solution.c_star)
Cf = solution.coefficient_of_thrust
Isp = solution.Isp/9.80665
Isp_vac = solution.Isp_vacuum/9.80665

print("P, atm         ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(P[i]), end=" ")
    else:
        print("{0:10.3f}".format(P[i]))

print("T, K           ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(T[i]), end=" ")
    else:
        print("{0:10.3f}".format(T[i]))

print("Density, g/cm^3", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3g}".format(rho[i]), end=" ")
    else:
        print("{0:10.3g}".format(rho[i]))

print("H, cal/g       ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.2f}".format(enthalpy[i]), end=" ")
    else:
        print("{0:10.2f}".format(enthalpy[i]))

print("U, cal/g       ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.2f}".format(energy[i]), end=" ")
    else:
        print("{0:10.2f}".format(energy[i]))

print("G, cal/g       ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.1f}".format(gibbs[i]), end=" ")
    else:
        print("{0:10.1f}".format(gibbs[i]))

print("S, cal/g-K     ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(entropy[i]), end=" ")
    else:
        print("{0:10.3f}".format(entropy[i]))

print("M, (1/n)       ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(M_1n[i]), end=" ")
    else:
        print("{0:10.3f}".format(M_1n[i]))

print("MW             ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(MW[i]), end=" ")
    else:
        print("{0:10.3f}".format(MW[i]))

print("Cp, cal/g-K    ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(cp_eq[i]), end=" ")
    else:
        print("{0:10.3f}".format(cp_eq[i]))

print("Gamma_s        ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.4f}".format(gamma_s[i]), end=" ")
    else:
        print("{0:10.4f}".format(gamma_s[i]))

print("Son. vel., m/s ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.2f}".format(v_sonic[i]), end=" ")
    else:
        print("{0:10.2f}".format(v_sonic[i]))

print("Mach           ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(Mach[i]), end=" ")
    else:
        print("{0:10.3f}".format(Mach[i]))

print("")
print("PERFORMANCE PARAMETERS")
print("")

print("Ae/At          ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.4f}".format(ae_at[i]), end=" ")
    else:
        print("{0:10.4f}".format(ae_at[i]))

print("C*, ft/s       ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.2f}".format(c_star[i]), end=" ")
    else:
        print("{0:10.2f}".format(c_star[i]))

print("Cf             ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.4f}".format(Cf[i]), end=" ")
    else:
        print("{0:10.4f}".format(Cf[i]))

print("Ivac, lb-s/lb  ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(Isp_vac[i]), end=" ")
    else:
        print("{0:10.3f}".format(Isp_vac[i]))

print("Isp, lb-s/lb   ", end=" ")
for i in range(num_pts):
    if i < num_pts-1:
        print("{0:10.3f}".format(Isp[i]), end=" ")
    else:
        print("{0:10.3f}".format(Isp[i]))

print()
print("MOLE FRACTIONS")
print("")
trace_species = []
for prod in solution.mole_fractions:
    if np.any(solution.mole_fractions[prod] > trace):
        print("{0:15s}".format(prod), end=" ")
        for j in range(len(solution.mole_fractions[prod])):
            if j < len(solution.mole_fractions[prod])-1:
                print("{0:10.5g}".format(solution.mole_fractions[prod][j]), end=" ")
            else:
                print("{0:10.5g}".format(solution.mole_fractions[prod][j]))
    else:
        trace_species.append(prod)

print()
print("TRACE SPECIES:")
max_cols = 8
nrows = (len(trace_species) + max_cols - 1) // max_cols
for i in range(nrows):
    print(" ".join("{0:15s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))
