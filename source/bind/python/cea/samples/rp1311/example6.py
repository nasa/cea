import numpy as np
import cea

"""
Example 6 from RP-1311
- Chapman-Jouguet detonation
- transport properties
"""

# Define the reactants
reac_names = ["H2", "O2"]

# Detonation states
p1 = np.array([1.0, 20.0]) # Initial pressure (bar)
T1 = np.array([298.15, 500.0])  # Initial temperature (K)

# Mixture states
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])
eq_ratios = np.array([1.0])

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

# Solver
solver = cea.DetonationSolver(prod, reactants=reac, transport=True)
solution = cea.DetonationSolution(solver)

# Compute o/f ratios from equivalence ratios
of_ratios = len(eq_ratios)*[0.0]
for i, eqrat in enumerate(eq_ratios):
    of_ratios[i] = reac.chem_eq_ratio_to_of_ratio(oxidant_weights,
                                                  fuel_weights,
                                                  eqrat)

n = len(eq_ratios)*len(p1)*len(T1)
of_ratio_out = np.zeros(n)
T1_out = np.zeros(n)
P1_out = np.zeros(n)
H1 = np.zeros(n)
M1 = np.zeros(n)
gamma1 = np.zeros(n)
v_sonic1 = np.zeros(n)
P = np.zeros(n)
T = np.zeros(n)
rho = np.zeros(n)
enthalpy = np.zeros(n)
energy = np.zeros(n)
gibbs = np.zeros(n)
entropy = np.zeros(n)
mach = np.zeros(n)
velocity = np.zeros(n)
v_sonic = np.zeros(n)
gamma_s = np.zeros(n)
P_P1 = np.zeros(n)
T_T1 = np.zeros(n)
M_M1 = np.zeros(n)
rho_rho1 = np.zeros(n)
cp_eq = np.zeros(n)
cp_fr = np.zeros(n)
cv_eq = np.zeros(n)
cv_fr = np.zeros(n)
molecular_weight_M = np.zeros(n)
molecular_weight_MW = np.zeros(n)
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
    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
    for t in T1:
        for p in p1:
            # Solve the detonation problem
            solver.solve(solution, weights, t, p)
            # Store the output
            of_ratio_out[i] = of_ratio
            T1_out[i] = solution.T1
            P1_out[i] = cea.units.bar_to_atm(solution.P1)
            H1[i] = cea.units.joule_to_cal(solution.H1)
            M1[i] = solution.M1
            gamma1[i] = solution.gamma1
            v_sonic1[i] = solution.sonic_velocity1
            P[i] = cea.units.bar_to_atm(solution.P)
            T[i] = solution.T
            rho[i] = solution.density*1.e-3
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
            v_sonic[i] = solution.sonic_velocity
            P_P1[i] = solution.P_P1
            T_T1[i] = solution.T_T1
            M_M1[i] = solution.M_M1
            rho_rho1[i] = solution.rho_rho1
            mach[i] = solution.Mach
            velocity[i] = solution.velocity
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

print()
print("UNBURNED GAS")
print()
print("o/f             ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(of_ratio_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(of_ratio_out[i]))

print("P1, atm         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P1_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(P1_out[i]))

print("T1, K           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T1_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(T1_out[i]))

print("H1, cal/g       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(H1[i]), end=" ")
    else:
        print("{0:10.3f}".format(H1[i]))

print("M1 (1/n)        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M1[i]), end=" ")
    else:
        print("{0:10.3f}".format(M1[i]))

print("Gamma1          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(gamma1[i]), end=" ")
    else:
        print("{0:10.3f}".format(gamma1[i]))

print("Son. Vel.1, m/s ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(v_sonic1[i]), end=" ")
    else:
        print("{0:10.3f}".format(v_sonic1[i]))

print()
print("BURNED GAS")
print()
print("P, atm          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P[i]), end=" ")
    else:
        print("{0:10.3f}".format(P[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T[i]), end=" ")
    else:
        print("{0:10.3f}".format(T[i]))

print("Density, g/cc   ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(rho[i]), end=" ")
    else:
        print("{0:10.3e}".format(rho[i]))

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

print("Son. Vel., m/s  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(v_sonic[i]), end=" ")
    else:
        print("{0:10.4f}".format(v_sonic[i]))

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
print("DETONATION PARAMETERS")
print()

print("P/P1            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(P_P1[i]), end=" ")
    else:
        print("{0:10.4f}".format(P_P1[i]))

print("T/T1            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(T_T1[i]), end=" ")
    else:
        print("{0:10.4f}".format(T_T1[i]))

print("M/M1            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(M_M1[i]), end=" ")
    else:
        print("{0:10.4f}".format(M_M1[i]))

print("rho/rho1        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(rho_rho1[i]), end=" ")
    else:
        print("{0:10.4f}".format(rho_rho1[i]))

print("Mach            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(mach[i]), end=" ")
    else:
        print("{0:10.4f}".format(mach[i]))

print("Velocity        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(velocity[i]), end=" ")
    else:
        print("{0:10.4f}".format(velocity[i]))

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
max_cols = 8
nrows = (len(trace_species) + max_cols - 1) // max_cols
for i in range(nrows):
    print(" ".join("{0:15s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))
