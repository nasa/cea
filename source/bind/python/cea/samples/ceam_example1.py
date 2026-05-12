import numpy as np
import cea
import cea.matlab
import os

R = cea.R

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
i = 0

# Equilibrium solve
for r_eq in chem_eq_ratios:
    for p in pressures:
        for t in temperatures:
            solution = cea.matlab.eq_solve(
                cea.TP,
                reac_names,
                T=t,
                P=p,
                fuel_amounts=fuel_moles,
                oxid_amounts=oxidant_moles,
                moles=True,
                r_eq=r_eq,
                only=prod_names,
            )

            # Store the output
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

            i += 1

print("P, atm         ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(P_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(P_out[i]))

print("T, K           ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(T_out[i]), end=" ")
    else:
        print("{0:10.3f}".format(T_out[i]))

print("Density, g/cc  ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(rho[i]), end=" ")
    else:
        print("{0:10.3e}".format(rho[i]))

print("Volume, cc/g   ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3e}".format(volume[i]), end=" ")
    else:
        print("{0:10.3e}".format(volume[i]))

print("H, cal/g       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(enthalpy[i]), end=" ")
    else:
        print("{0:10.3f}".format(enthalpy[i]))

print("U, cal/g       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(energy[i]), end=" ")
    else:
        print("{0:10.3f}".format(energy[i]))

print("G, cal/g       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.1f}".format(gibbs[i]), end=" ")
    else:
        print("{0:10.1f}".format(gibbs[i]))

print("S, cal/g-K     ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(entropy[i]), end=" ")
    else:
        print("{0:10.3f}".format(entropy[i]))

print("M, (1/n)       ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(molecular_weight_M[i]), end=" ")
    else:
        print("{0:10.3f}".format(molecular_weight_M[i]))

print("MW             ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.3f}".format(molecular_weight_MW[i]), end=" ")
    else:
        print("{0:10.3f}".format(molecular_weight_MW[i]))

print("Gamma_s        ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(gamma_s[i]), end=" ")
    else:
        print("{0:10.4f}".format(gamma_s[i]))

print("Cp_eq, cal/g-K ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_eq[i]))

print("Cp_fr, cal/g-K ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cp_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cp_fr[i]))

print("Cv_eq, cal/g-K ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_eq[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_eq[i]))

print("Cv_fr, cal/g-K ", end="")
for i in range(n):
    if i < n-1:
        print("{0:10.4f}".format(cv_fr[i]), end=" ")
    else:
        print("{0:10.4f}".format(cv_fr[i]))
