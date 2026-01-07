import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

import cea

# Compute HP equilibrium of H2 and O2 at various equivalence ratios and initial temperatures
# Total points: 108,500

start = datetime.now()

# Thermo states
p0 = cea.units.atm_to_bar(1.0)  # Fixed-pressure state (bar)
n_T0 = 350
T0 = np.linspace(300.0, 2000.0, n_T0)   # Initial reactant temperature (K)

# Mixtures
reac = cea.Mixture(["H2", "O2"])
prod = cea.Mixture(["H", "H2", "H2O", "O", "O2", "OH"])

# Values of equivalence ratio, phi
n_phi = 310
phi = np.linspace(0.5, 2.0, n_phi)

# Solver
solver = cea.EqSolver(prod, reactants=reac)
soln = cea.EqSolution(solver)

# Initialize arrays to store some solution values
X, Y = np.meshgrid(T0, phi)
T_vals = np.zeros((n_phi, n_T0))
H_vals = np.zeros((n_phi, n_T0))
H2_vals = np.zeros((n_phi, n_T0))
H2O_vals = np.zeros((n_phi, n_T0))
O_vals = np.zeros((n_phi, n_T0))
O2_vals = np.zeros((n_phi, n_T0))
OH_vals = np.zeros((n_phi, n_T0))
convg_vals = np.zeros((n_phi, n_T0))

# Solve the HP equilibrium at each condition
for i in range(len(phi)):
    of_ratio = reac.weight_eq_ratio_to_of_ratio(np.array([0.0, 1.0]), np.array([1.0, 0.0]), phi[i])
    weights = reac.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), of_ratio)
    for j in range(len(T0)):
        # Get the fixed-enthalpy
        h0 = reac.calc_property(cea.ENTHALPY, weights, T0[j])/cea.R

        # Equilibrium solve
        solver.solve(soln, cea.HP, h0, p0, weights)

        convg_vals[i,j] = soln.converged
        if soln.converged:
            T_vals[i,j] = soln.T
            H_vals[i,j]   = soln.mole_fractions["H"]
            H2_vals[i,j]  = soln.mole_fractions["H2"]
            H2O_vals[i,j] = soln.mole_fractions["H2O"]
            O_vals[i,j]   = soln.mole_fractions["O"]
            O2_vals[i,j]  = soln.mole_fractions["O2"]
            OH_vals[i,j]  = soln.mole_fractions["OH"]
        else:
            print("Not converged for phi = ", phi[i], " and T0 = ", T0[j])

print("Time: ", datetime.now() - start)

def fmt(x):
    s = f"{100*x:.1f}"
    if s.endswith("0"):
        s = f"{100*x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

'''
Plot the reaction temperature values
'''
fig = plt.figure(figsize=(10, 10), dpi=100)
levels = 10

xticks = [300, 500, 1000, 1500, 2000]
yticks = np.linspace(0.5, 2.0, 4)

ax = fig.add_subplot(1, 1, 1)
CF = ax.contourf(X, Y, T_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax.contour(X, Y, T_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax.clabel(CS, CS.levels[4:], fmt="%.0f K")

# ax.axhline(1.0, linestyle='--', color="white", linewidth=2)

ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(xticks)
ax.set_yticks(yticks)

fs = 20
plt.title("Reaction Temperature", fontsize=fs)
fs = 16
plt.xlabel(r"Initial Temperature (K)", fontsize=fs)
plt.ylabel(r"Equivalence Ratio, $\phi$", fontsize=fs)

plt.show()
#plt.savefig("reaction_temperature.pdf", bbox_inches='tight')

'''
Plot the species concentrations
'''
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
levels = 10
fs = 16
params = {'axes.titlesize':'xx-large',
          'xtick.labelsize':'large',
          'ytick.labelsize':'large'}
plt.rcParams.update(params)

# -----------------------------------------
# H Concentration
# -----------------------------------------
CF = ax[0,0].contourf(X, Y, H_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[0,0].contour(X, Y, H_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[0,0].clabel(CS, CS.levels[:-5], fmt=fmt)#"%2.1e")

ax[0,0].spines['left'].set_visible(False)
ax[0,0].spines['right'].set_visible(False)
ax[0,0].spines['top'].set_visible(False)
ax[0,0].spines['bottom'].set_visible(False)
ax[0,0].yaxis.set_ticks_position('left')
ax[0,0].xaxis.set_ticks_position('bottom')
ax[0,0].set_xticks([])
ax[0,0].set_yticks(yticks)

ax[0,0].title.set_text(r'$n_{j}$: H')

# -----------------------------------------
# H2 Concentration
# -----------------------------------------
CF = ax[0,1].contourf(X, Y, H2_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[0,1].contour(X, Y, H2_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[0,1].clabel(CS, CS.levels, fmt=fmt)#"%2.1e")

ax[0,1].spines['left'].set_visible(False)
ax[0,1].spines['right'].set_visible(False)
ax[0,1].spines['top'].set_visible(False)
ax[0,1].spines['bottom'].set_visible(False)
ax[0,1].yaxis.set_ticks_position('left')
ax[0,1].xaxis.set_ticks_position('bottom')
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])

ax[0,1].title.set_text(r'$n_{j}$: H$_{2}$')

# -----------------------------------------
# H2O Concentration
# -----------------------------------------
CF = ax[0,2].contourf(X, Y, H2O_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[0,2].contour(X, Y, H2O_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[0,2].clabel(CS, CS.levels, fmt=fmt)#"%2.1e")

ax[0,2].spines['left'].set_visible(False)
ax[0,2].spines['right'].set_visible(False)
ax[0,2].spines['top'].set_visible(False)
ax[0,2].spines['bottom'].set_visible(False)
ax[0,2].yaxis.set_ticks_position('left')
ax[0,2].xaxis.set_ticks_position('bottom')
ax[0,2].set_xticks([])
ax[0,2].set_yticks([])

ax[0,2].title.set_text(r'$n_{j}$: H$_{2}$O')

# -----------------------------------------
# O Concentration
# -----------------------------------------
CF = ax[1,0].contourf(X, Y, O_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[1,0].contour(X, Y, O_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[1,0].clabel(CS, CS.levels[:-5], fmt=fmt)#="%2.1e")

ax[1,0].spines['left'].set_visible(False)
ax[1,0].spines['right'].set_visible(False)
ax[1,0].spines['top'].set_visible(False)
ax[1,0].spines['bottom'].set_visible(False)
ax[1,0].yaxis.set_ticks_position('left')
ax[1,0].xaxis.set_ticks_position('bottom')
ax[1,0].set_xticks(xticks)
ax[1,0].set_yticks(yticks)

ax[1,0].title.set_text(r'$n_{j}$: O')

# -----------------------------------------
# O2 Concentration
# -----------------------------------------
CF = ax[1,1].contourf(X, Y, O2_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[1,1].contour(X, Y, O2_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[1,1].clabel(CS, CS.levels[::2], fmt=fmt)#="%2.1e")

ax[1,1].spines['left'].set_visible(False)
ax[1,1].spines['right'].set_visible(False)
ax[1,1].spines['top'].set_visible(False)
ax[1,1].spines['bottom'].set_visible(False)
ax[1,1].yaxis.set_ticks_position('left')
ax[1,1].xaxis.set_ticks_position('bottom')
ax[1,1].set_xticks(xticks)
ax[1,1].set_yticks([])

ax[1,1].title.set_text(r'$n_{j}$: O$_{2}$')

# -----------------------------------------
# OH Concentration
# -----------------------------------------
CF = ax[1,2].contourf(X, Y, OH_vals, origin="lower", levels=levels, cmap="coolwarm")
CS = ax[1,2].contour(X, Y, OH_vals, origin="lower", levels=levels, colors="black", linewidths=0.5)
ax[1,2].clabel(CS, CS.levels, fmt=fmt)#"%2.1e")

ax[1,2].spines['left'].set_visible(False)
ax[1,2].spines['right'].set_visible(False)
ax[1,2].spines['top'].set_visible(False)
ax[1,2].spines['bottom'].set_visible(False)
ax[1,2].yaxis.set_ticks_position('left')
ax[1,2].xaxis.set_ticks_position('bottom')
ax[1,2].set_xticks(xticks)
ax[1,2].set_yticks([])

ax[1,2].title.set_text(r'$n_{j}$: OH')

fs = 20
fig.suptitle("Species Mole Fraction", fontsize=fs)
fs = 16
fig.text(0.5, 0.04, r"Initial Temperature (K)", fontsize=fs, ha='center')
fig.text(0.04, 0.5, r"Equivalence Ratio, $\phi$", fontsize=fs, va='center', rotation='vertical')

plt.show()
#plt.savefig("species_concentrations.pdf", bbox_inches='tight')
