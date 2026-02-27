import numpy as np
import cea

"""
Example 7 from RP-1311
- Shock tube problem
- Incident and reflected shocks with equilibrium analysis
"""

# Define the reactants
reac_names = ["H2", "O2", "Ar"]
moles = np.array([0.05, 0.05, 0.9])

# Shock states
p0 = cea.units.mmhg_to_bar(10.0)  # Unshocked pressure (bar)
T0 = 300.0   # Unshocked temperature (K)
u1 = np.array([1100, 1200, 1250, 1300, 1350, 1400])  # Initial velocity (m/s)

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

# Solver
solver = cea.ShockSolver(prod, reactants=reac)
solution = cea.ShockSolution(solver, reflected=True)

# Get the reactant weights
weights = reac.moles_to_weights(moles)

# Initialize output arrays
n = len(u1)
mach1 = np.zeros(n)
u1_out = np.zeros(n)
P1_out = np.zeros(n)
T1_out = np.zeros(n)
rho1 = np.zeros(n)
H1 = np.zeros(n)
U1 = np.zeros(n)
G1 = np.zeros(n)
S1 = np.zeros(n)
M1 = np.zeros(n)
cp1 = np.zeros(n)
gamma1 = np.zeros(n)
v_sonic1 = np.zeros(n)

u2 = np.zeros(n)
P2 = np.zeros(n)
T2 = np.zeros(n)
rho2 = np.zeros(n)
H2 = np.zeros(n)
U2 = np.zeros(n)
G2 = np.zeros(n)
S2 = np.zeros(n)
M2 = np.zeros(n)
cp2 = np.zeros(n)
gamma2 = np.zeros(n)
v_sonic2 = np.zeros(n)
P2_P1 = np.zeros(n)
T2_T1 = np.zeros(n)
M2_M1 = np.zeros(n)
rho2_rho1 = np.zeros(n)
v2 = np.zeros(n)

u5 = np.zeros(n)
P5 = np.zeros(n)
T5 = np.zeros(n)
rho5 = np.zeros(n)
H5 = np.zeros(n)
U5 = np.zeros(n)
G5 = np.zeros(n)
S5 = np.zeros(n)
M5 = np.zeros(n)
cp5 = np.zeros(n)
gamma5 = np.zeros(n)
v_sonic5 = np.zeros(n)
P5_P2 = np.zeros(n)
T5_T2 = np.zeros(n)
M5_M2 = np.zeros(n)
rho5_rho2 = np.zeros(n)
u5_p_v2 = np.zeros(n)
mole_fractions_incd = {}
mole_fractions_refl = {}

i = 0
for u in u1:
    # Solve the shock problem
    solver.solve(solution, weights, T0, p0, u1=u, reflected=True)

    # Store the output
    u1_out[i] = u
    P1_out[i] = p0
    T1_out[i] = T0
    if solution.converged:
        mach1[i] = solution.Mach[0]
        rho1[i] = solution.density[0]
        H1[i] = solution.enthalpy[0]
        U1[i] = solution.energy[0]
        G1[i] = solution.gibbs_energy[0]
        S1[i] = solution.entropy[0]
        M1[i] = solution.M[0]
        cp1[i] = solution.cp_eq[0]
        gamma1[i] = solution.gamma_s[0]
        v_sonic1[i] = solution.sonic_velocity[0]

        u2[i] = solution.velocity[1]
        P2[i] = solution.P[1]
        T2[i] = solution.T[1]
        rho2[i] = solution.density[1]
        H2[i] = solution.enthalpy[1]
        U2[i] = solution.energy[1]
        G2[i] = solution.gibbs_energy[1]
        S2[i] = solution.entropy[1]
        M2[i] = solution.M[1]
        cp2[i] = solution.cp_eq[1]
        gamma2[i] = solution.gamma_s[1]
        v_sonic2[i] = solution.sonic_velocity[1]
        P2_P1[i] = solution.P21
        T2_T1[i] = solution.T21
        M2_M1[i] = solution.M21
        rho2_rho1[i] = 1.0/solution.rho12
        v2[i] = solution.v2

        u5[i] = solution.velocity[2]
        P5[i] = solution.P[2]
        T5[i] = solution.T[2]
        rho5[i] = solution.density[2]
        H5[i] = solution.enthalpy[2]
        U5[i] = solution.energy[2]
        G5[i] = solution.gibbs_energy[2]
        S5[i] = solution.entropy[2]
        M5[i] = solution.M[2]
        cp5[i] = solution.cp_eq[2]
        gamma5[i] = solution.gamma_s[2]
        v_sonic5[i] = solution.sonic_velocity[2]
        P5_P2[i] = solution.P52
        T5_T2[i] = solution.T52
        M5_M2[i] = solution.M52
        rho5_rho2[i] = solution.rho52
        u5_p_v2[i] = solution.u5_p_v2

        if i == 0:
            for prod in solution.mole_fractions:
                    mole_fractions_incd[prod] = np.array([solution.mole_fractions[prod][1]])
                    mole_fractions_refl[prod] = np.array([solution.mole_fractions[prod][2]])
        else:
            for prod in mole_fractions_incd:
                mole_fractions_incd[prod] = np.append(mole_fractions_incd[prod], solution.mole_fractions[prod][1])
            for prod in mole_fractions_refl:
                mole_fractions_refl[prod] = np.append(mole_fractions_refl[prod], solution.mole_fractions[prod][2])

    i += 1

# Print the results
print()
print("INITIAL GAS (1)")
print()

print("Mach Number1    ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(mach1[i]), end=" ")
    else:
        print("{0:10.3f}".format(mach1[i]))

print("u1, m/s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(u1_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(u1_out[i]))

print("P, bar          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(P1_out[i]), end=" ")
    else:
        print("{0:10.4f}".format(P1_out[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T1_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(T1_out[i]))

print("Density, kg/m^3 ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4e}".format(rho1[i]), end=" ")
    else:
        print("{0:10.4e}".format(rho1[i]))

print("H, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(H1[i]), end=" ")
    else:
        print("{0:10.3f}".format(H1[i]))

print("U, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(U1[i]), end=" ")
    else:
        print("{0:10.3f}".format(U1[i]))

print("G, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(G1[i]), end=" ")
    else:
        print("{0:10.3f}".format(G1[i]))

print("S, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(S1[i]), end=" ")
    else:
        print("{0:10.3f}".format(S1[i]))

print("M, (1/n)        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M1[i]), end=" ")
    else:
        print("{0:10.3f}".format(M1[i]))

print("Cp, kJ/kg-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(cp1[i]), end=" ")
    else:
        print("{0:10.3f}".format(cp1[i]))

print("Gamma_s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(gamma1[i]), end=" ")
    else:
        print("{0:10.3f}".format(gamma1[i]))

print("Son. Vel., m/s  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(v_sonic1[i]), end=" ")
    else:
        print("{0:10.3f}".format(v_sonic1[i]))

print()
print("SHOCKED GAS (2) - Incident, Equilibrium")
print()
print("u2, m/s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(u2[i]), end=" ")
    else:
        print("{0:10.3f}".format(u2[i]))

print("P, bar          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(P2[i]), end=" ")
    else:
        print("{0:10.4f}".format(P2[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T2[i]), end=" ")
    else:
        print("{0:10.3f}".format(T2[i]))

print("Density, kg/m^3 ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4e}".format(rho2[i]), end=" ")
    else:
        print("{0:10.4e}".format(rho2[i]))

print("H, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(H2[i]), end=" ")
    else:
        print("{0:10.3f}".format(H2[i]))

print("U, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(U2[i]), end=" ")
    else:
        print("{0:10.3f}".format(U2[i]))

print("G, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(G2[i]), end=" ")
    else:
        print("{0:10.3f}".format(G2[i]))

print("S, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(S2[i]), end=" ")
    else:
        print("{0:10.3f}".format(S2[i]))

print("M, (1/n)        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M2[i]), end=" ")
    else:
        print("{0:10.3f}".format(M2[i]))

print("Cp, kJ/kg-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(cp2[i]), end=" ")
    else:
        print("{0:10.3f}".format(cp2[i]))

print("Gamma_s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(gamma2[i]), end=" ")
    else:
        print("{0:10.3f}".format(gamma2[i]))

print("Son. Vel., m/s  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(v_sonic2[i]), end=" ")
    else:
        print("{0:10.3f}".format(v_sonic2[i]))

print()
print("P2/P1           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P2_P1[i]), end=" ")
    else:
        print("{0:10.3f}".format(P2_P1[i]))

print("T2/T1           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T2_T1[i]), end=" ")
    else:
        print("{0:10.3f}".format(T2_T1[i]))

print("M2/M1           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M2_M1[i]), end=" ")
    else:
        print("{0:10.3f}".format(M2_M1[i]))

print("rho2/rho1       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(rho2_rho1[i]), end=" ")
    else:
        print("{0:10.3f}".format(rho2_rho1[i]))

print("v2, m/s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(v2[i]), end=" ")
    else:
        print("{0:10.3f}".format(v2[i]))

print()
print("MOLE FRACTIONS")
print()
trace_species = []
for prod in mole_fractions_incd:
    if np.any(mole_fractions_incd[prod] > 5e-6):
        print("{0:15s}".format(prod), end=" ")
        for j in range(n):
            if j < n-1:
                print("{0:10.5g}".format(mole_fractions_incd[prod][j]), end=" ")
            else:
                print("{0:10.5g}".format(mole_fractions_incd[prod][j]))
    else:
        trace_species.append(prod)

if len(trace_species) > 0:
    print()
    print("TRACE SPECIES:")
    max_cols = 8
    nrows = (len(trace_species) + max_cols - 1) // max_cols
    for i in range(nrows):
        print(" ".join("{0:15s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))

print()
print("SHOCKED GAS (5) - Reflected, Equilibrium")
print()
print("u5, m/s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(u5[i]), end=" ")
    else:
        print("{0:10.3f}".format(u5[i]))

print("P, bar          ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(P5[i]), end=" ")
    else:
        print("{0:10.4f}".format(P5[i]))

print("T, K            ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T5[i]), end=" ")
    else:
        print("{0:10.3f}".format(T5[i]))

print("Density, kg/m^3 ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4e}".format(rho5[i]), end=" ")
    else:
        print("{0:10.4e}".format(rho5[i]))

print("H, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(H5[i]), end=" ")
    else:
        print("{0:10.3f}".format(H5[i]))

print("U, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(U5[i]), end=" ")
    else:
        print("{0:10.3f}".format(U5[i]))

print("G, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(G5[i]), end=" ")
    else:
        print("{0:10.3f}".format(G5[i]))

print("S, kJ/kg        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(S5[i]), end=" ")
    else:
        print("{0:10.3f}".format(S5[i]))

print("M, (1/n)        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M5[i]), end=" ")
    else:
        print("{0:10.3f}".format(M5[i]))

print("Cp, kJ/kg-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(cp5[i]), end=" ")
    else:
        print("{0:10.3f}".format(cp5[i]))

print("Gamma_s         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(gamma5[i]), end=" ")
    else:
        print("{0:10.3f}".format(gamma5[i]))

print("Son. Vel., m/s  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(v_sonic5[i]), end=" ")
    else:
        print("{0:10.3f}".format(v_sonic5[i]))

print()
print("P5/P2           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P5_P2[i]), end=" ")
    else:
        print("{0:10.3f}".format(P5_P2[i]))

print("T5/T2           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T5_T2[i]), end=" ")
    else:
        print("{0:10.3f}".format(T5_T2[i]))

print("M5/M2           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(M5_M2[i]), end=" ")
    else:
        print("{0:10.3f}".format(M5_M2[i]))

print("rho5/rho2       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(rho5_rho2[i]), end=" ")
    else:
        print("{0:10.3f}".format(rho5_rho2[i]))

print("u5+v2, m/s      ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(u5_p_v2[i]), end=" ")
    else:
        print("{0:10.3f}".format(u5_p_v2[i]))

print()
print("MOLE FRACTIONS")
print()
trace_species = []
for prod in mole_fractions_refl:
    if np.any(mole_fractions_refl[prod] > 5e-6):
        print("{0:15s}".format(prod), end=" ")
        for j in range(n):
            if j < n-1:
                print("{0:10.5g}".format(mole_fractions_refl[prod][j]), end=" ")
            else:
                print("{0:10.5g}".format(mole_fractions_refl[prod][j]))
    else:
        trace_species.append(prod)

print()
print("TRACE SPECIES:")
max_cols = 8
nrows = (len(trace_species) + max_cols - 1) // max_cols
for i in range(nrows):
    print(" ".join("{0:15s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))
