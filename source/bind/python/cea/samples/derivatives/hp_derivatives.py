import numpy as np
import matplotlib.pyplot as plt
import cea

"""
Example demonstrating equilibrium derivatives for HP problem with finite-difference convergence study
- Single HP equilibrium calculation (based on Example 3 from RP-1311)
- Compute analytic derivatives
- Sweep finite-difference step sizes from 1e-1 to 1e-16
- Plot convergence behavior for representative outputs

Note: This example uses central differences (central=True) for second-order accuracy.
Set central=False to use forward differences, which are faster but less accurate.
"""

# Global flag to control error metric
USE_ABSOLUTE_ERROR = False  # Set to True for absolute error, False for relative error

def rel_error(analytical, fd):
    """
    Compute error between analytical and finite-difference derivatives.

    The error metric is controlled by the global USE_ABSOLUTE_ERROR flag:
    - If USE_ABSOLUTE_ERROR=True: Returns absolute error = |fd - analytical|
    - If USE_ABSOLUTE_ERROR=False: Returns relative error = |fd - analytical| / |analytical|

    Note: Relative error is NOT bounded by 1.0. It can be arbitrarily large when
    the finite-difference approximation is poor (either due to truncation error
    with large h, or round-off error with tiny h). This is expected and normal
    in convergence studies.
    """
    abs_error = abs(fd - analytical)

    if USE_ABSOLUTE_ERROR:
        # Return absolute error
        return abs_error
    else:
        # Return relative error
        if abs(analytical) < 1e-15:
            # For near-zero analytical values, return absolute error
            return abs_error if abs(fd) > 1e-15 else 0.0
        return abs_error / abs(analytical)

# Initialize CEA
cea.init()

# Species (from Example 3)
reac_names = ["Air", "C7H8(L)", "C8H18(L),n-octa"]
omit_names = ["CCN", "CNC", "C2N2", "C2O",
              "C3H4,allene", "C3H4,propyne", "C3H4,cyclo-", "C3",
              "C3H5,allyl", "C3H6,propylene", "C3H6,cyclo-", "C3H3,propargyl",
              "C3H6O", "C3H7,n-propyl", "C3H7,i-propyl", "Jet-A(g)",
              "C3O2", "C4", "C4H2", "C3H8O,2propanol",
              "C4H4,1,3-cyclo-", "C4H6,butadiene", "C4H6,2-butyne", "C3H8O,1propanol",
              "C4H8,tr2-butene", "C4H8,isobutene", "C4H8,cyclo-", "C4H6,cyclo-",
              "(CH3COOH)2", "C4H9,n-butyl", "C4H9,i-butyl", "C4H8,1-butene",
              "C4H9,s-butyl", "C4H9,t-butyl", "C4H10,isobutane", "C4H8,cis2-buten",
              "C4H10,n-butane", "C4N2", "C5", "C3H8",
              "C5H6,1,3cyclo-", "C5H8,cyclo-", "C5H10,1-pentene", "C10H21,n-decyl",
              "C5H10,cyclo-", "C5H11,pentyl", "C5H11,t-pentyl", "C12H10,biphenyl",
              "C5H12,n-pentane", "C5H12,i-pentane", "CH3C(CH3)2CH3", "C12H9,o-bipheny",
              "C6H6", "C6H5OH,phenol", "C6H10,cyclo-", "C6H2",
              "C6H12,1-hexene", "C6H12,cyclo-", "C6H13,n-hexyl", "C6H5,phenyl",
              "C7H7,benzyl", "C7H8", "C7H8O,cresol-mx", "C6H5O,phenoxy",
              "C7H14,1-heptene", "C7H15,n-heptyl", "C7H16,n-heptane", "C10H8,azulene",
              "C8H8,styrene", "C8H10,ethylbenz", "C8H16,1-octene", "C10H8,napthlene",
              "C8H17,n-octyl", "C8H18,isooctane", "C8H18,n-octane", "C9H19,n-nonyl",
              "Jet-A(L)", "C6H6(L)", "H2O(s)", "H2O(L)"]

# Single test condition
pressure = 100.0  # bar
trace = 1e-15

# Mixture composition (from Example 3)
fuel_weights = np.array([0.0, 0.4, 0.6])
oxidant_weights = np.array([1.0, 0.0, 0.0])
T_reac = np.array([700.0, 298.15, 298.15])
of_ratio = 17.0

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True, omit=omit_names)

# Solver
solver = cea.EqSolver(prod, reactants=reac, trace=trace, smooth_truncation=True)
solution = cea.EqSolution(solver)

# Compute reactant weights and enthalpy
weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
h0 = reac.calc_property(cea.ENTHALPY, weights, T_reac)

# Equilibrium solve (HP problem)
solver.solve(solution, cea.HP, h0/cea.R, pressure, weights)

if not solution.converged:
    print("ERROR: Solution did not converge")
    exit(1)

print("=" * 70)
print("EQUILIBRIUM SOLUTION (HP)")
print("=" * 70)
print(f"o/f ratio       : {of_ratio:.4f}")
print(f"P, bar          : {pressure:.4f}")
print(f"H, kJ/kg        : {solution.enthalpy:.4f}")
print(f"T, K            : {solution.T:.2f}")
print(f"S, kJ/(kg·K)    : {solution.entropy:.6f}")
print(f"M, (1/n)        : {solution.M:.4f}")
print(f"Gamma_s         : {solution.gamma_s:.6f}")

# Get species names
species_names = prod.species_names
print(f"\nNumber of product species: {len(species_names)}")

print("\n" + "=" * 70)
print("COMPUTING ANALYTIC DERIVATIVES")
print("=" * 70)

# Create derivatives object and compute analytic derivatives
derivs = cea.EqDerivatives(solver, solution)
derivs.compute_derivatives(check_closure_defect=False)

# Get analytic derivatives (for HP: state1=H, state2=P)
dT_dH_analytic = derivs.dT_dstate1
dT_dP_analytic = derivs.dT_dstate2
dT_dw0_analytic = derivs.dT_dw0

dS_dH_analytic = derivs.dS_dstate1
dS_dP_analytic = derivs.dS_dstate2
dS_dw0_analytic = derivs.dS_dw0

dnj_dH_analytic = derivs.dnj_dstate1
dnj_dP_analytic = derivs.dnj_dstate2
dnj_dw0_analytic = derivs.dnj_dw0

print(f"dT/dH (analytic)       : {dT_dH_analytic:.6e} K/(kJ/kg)")
print(f"dT/dP (analytic)       : {dT_dP_analytic:.6e} K/bar")
print(f"dS/dH (analytic)       : {dS_dH_analytic:.6e} (kg·K)/(kJ)")
print(f"dS/dP (analytic)       : {dS_dP_analytic:.6e} kJ/(kg·K·bar)")

print("\n" + "=" * 70)
print("FINITE-DIFFERENCE CONVERGENCE STUDY")
print("=" * 70)
print(f"Error metric: {'ABSOLUTE' if USE_ABSOLUTE_ERROR else 'RELATIVE'}")

# Define step sizes for convergence study
h_values = np.logspace(-1, -16, 50)
print(f"Testing {len(h_values)} step sizes from {h_values[0]:.1e} to {h_values[-1]:.1e}")

# Storage for errors
errors_dT_dH = []
errors_dT_dP = []
errors_dT_dw0 = []
errors_dS_dH = []
errors_dS_dP = []
errors_dS_dw0 = []
errors_dnj_dH = []
errors_dnj_dP = []
errors_dnj_dw0 = []

# Select a few representative species for plotting
# Choose species with significant mole fractions
mole_fracs = solution.mole_fractions
significant_species = [(name, mole_fracs[name]) for name in species_names if mole_fracs[name] > 1e-4]
significant_species.sort(key=lambda x: x[1], reverse=True)
plot_species_names = [name for name, _ in significant_species[:5]]  # Top 5 species
plot_species_indices = [species_names.index(name) for name in plot_species_names]

print(f"\nPlotting species: {', '.join(plot_species_names)}")

# Sweep through step sizes
for i, h in enumerate(h_values):
    if (i + 1) % 5 == 0:
        print(f"  Computing h = {h:.1e} ({i+1}/{len(h_values)})")

    # Compute finite-difference derivatives (central differences for 2nd-order accuracy)
    derivs.compute_fd(h=h, verbose=False, central=True)

    # Temperature derivatives
    dT_dH_fd = derivs.dT_dstate1_fd
    dT_dP_fd = derivs.dT_dstate2_fd
    dT_dw0_fd = derivs.dT_dw0_fd

    errors_dT_dH.append(rel_error(dT_dH_analytic, dT_dH_fd))
    errors_dT_dP.append(rel_error(dT_dP_analytic, dT_dP_fd))
    errors_dT_dw0.append(rel_error(dT_dw0_analytic[0], dT_dw0_fd[0]))  # First component

    # Entropy derivatives
    dS_dH_fd = derivs.dS_dstate1_fd
    dS_dP_fd = derivs.dS_dstate2_fd
    dS_dw0_fd = derivs.dS_dw0_fd

    errors_dS_dH.append(rel_error(dS_dH_analytic, dS_dH_fd))
    errors_dS_dP.append(rel_error(dS_dP_analytic, dS_dP_fd))
    errors_dS_dw0.append(rel_error(dS_dw0_analytic[0], dS_dw0_fd[0]))  # First component

    # Species derivatives
    dnj_dH_fd = derivs.dnj_dstate1_fd
    dnj_dP_fd = derivs.dnj_dstate2_fd
    dnj_dw0_fd = derivs.dnj_dw0_fd

    # Store errors for selected species
    errs_dH = [rel_error(dnj_dH_analytic[idx], dnj_dH_fd[idx]) for idx in plot_species_indices]
    errs_dP = [rel_error(dnj_dP_analytic[idx], dnj_dP_fd[idx]) for idx in plot_species_indices]
    errs_dw0 = [rel_error(dnj_dw0_analytic[idx, 0], dnj_dw0_fd[idx, 0]) for idx in plot_species_indices]

    errors_dnj_dH.append(errs_dH)
    errors_dnj_dP.append(errs_dP)
    errors_dnj_dw0.append(errs_dw0)

# Convert to numpy arrays
errors_dT_dH = np.array(errors_dT_dH)
errors_dT_dP = np.array(errors_dT_dP)
errors_dT_dw0 = np.array(errors_dT_dw0)
errors_dS_dH = np.array(errors_dS_dH)
errors_dS_dP = np.array(errors_dS_dP)
errors_dS_dw0 = np.array(errors_dS_dw0)
errors_dnj_dH = np.array(errors_dnj_dH)
errors_dnj_dP = np.array(errors_dnj_dP)
errors_dnj_dw0 = np.array(errors_dnj_dw0)

print("\nConvergence study complete!")

# Print error statistics
print("\n" + "=" * 70)
print("ERROR STATISTICS")
print("=" * 70)
print("Note: Relative errors can be >> 1.0 when FD approximation is poor.")
print("This is expected - the plots show the characteristic U-shaped curve.\n")

print("∂T/∂H errors:")
print(f"  Min: {np.min(errors_dT_dH):.2e}  Max: {np.max(errors_dT_dH):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dT_dH)]:.1e}")

print("∂T/∂P errors:")
print(f"  Min: {np.min(errors_dT_dP):.2e}  Max: {np.max(errors_dT_dP):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dT_dP)]:.1e}")

print("∂S/∂H errors:")
print(f"  Min: {np.min(errors_dS_dH):.2e}  Max: {np.max(errors_dS_dH):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dS_dH)]:.1e}")

print("∂S/∂P errors:")
print(f"  Min: {np.min(errors_dS_dP):.2e}  Max: {np.max(errors_dS_dP):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dS_dP)]:.1e}")

print("\nGenerating plots...")

# Create plots
plt.style.use('default')

# Figure 1: Temperature derivatives
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5))
fig1.suptitle('Temperature Derivative Convergence (HP)', fontsize=14, fontweight='bold')

# Find optimal points
opt_idx_dT_dH = np.argmin(errors_dT_dH)
opt_idx_dT_dP = np.argmin(errors_dT_dP)
opt_idx_dT_dw0 = np.argmin(errors_dT_dw0)

axes1[0].loglog(h_values, errors_dT_dH, 'o-', linewidth=2, markersize=4, label='Error')
axes1[0].axvline(h_values[opt_idx_dT_dH], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dH]:.1e}')
axes1[0].set_xlabel('Step size h', fontsize=11)
axes1[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[0].set_title('∂T/∂H', fontsize=12)
axes1[0].grid(True, which='both', alpha=0.3)
axes1[0].legend(fontsize=9)
axes1[0].invert_xaxis()

axes1[1].loglog(h_values, errors_dT_dP, 'o-', linewidth=2, markersize=4, label='Error')
axes1[1].axvline(h_values[opt_idx_dT_dP], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dP]:.1e}')
axes1[1].set_xlabel('Step size h', fontsize=11)
axes1[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[1].set_title('∂T/∂P', fontsize=12)
axes1[1].grid(True, which='both', alpha=0.3)
axes1[1].legend(fontsize=9)
axes1[1].invert_xaxis()

axes1[2].loglog(h_values, errors_dT_dw0, 'o-', linewidth=2, markersize=4, label='Error')
axes1[2].axvline(h_values[opt_idx_dT_dw0], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dw0]:.1e}')
axes1[2].set_xlabel('Step size h', fontsize=11)
axes1[2].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[2].set_title('∂T/∂w₀[0]', fontsize=12)
axes1[2].grid(True, which='both', alpha=0.3)
axes1[2].legend(fontsize=9)
axes1[2].invert_xaxis()

plt.tight_layout()

# Figure 2: Entropy derivatives
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
fig2.suptitle('Entropy Derivative Convergence (HP)', fontsize=14, fontweight='bold')

# Find optimal points
opt_idx_dS_dH = np.argmin(errors_dS_dH)
opt_idx_dS_dP = np.argmin(errors_dS_dP)
opt_idx_dS_dw0 = np.argmin(errors_dS_dw0)

axes2[0].loglog(h_values, errors_dS_dH, 'o-', linewidth=2, markersize=4, label='Error')
axes2[0].axvline(h_values[opt_idx_dS_dH], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dS_dH]:.1e}')
axes2[0].set_xlabel('Step size h', fontsize=11)
axes2[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[0].set_title('∂S/∂H', fontsize=12)
axes2[0].grid(True, which='both', alpha=0.3)
axes2[0].legend(fontsize=9)
axes2[0].invert_xaxis()

axes2[1].loglog(h_values, errors_dS_dP, 'o-', linewidth=2, markersize=4, label='Error')
axes2[1].axvline(h_values[opt_idx_dS_dP], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dS_dP]:.1e}')
axes2[1].set_xlabel('Step size h', fontsize=11)
axes2[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[1].set_title('∂S/∂P', fontsize=12)
axes2[1].grid(True, which='both', alpha=0.3)
axes2[1].legend(fontsize=9)
axes2[1].invert_xaxis()

axes2[2].loglog(h_values, errors_dS_dw0, 'o-', linewidth=2, markersize=4, label='Error')
axes2[2].axvline(h_values[opt_idx_dS_dw0], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dS_dw0]:.1e}')
axes2[2].set_xlabel('Step size h', fontsize=11)
axes2[2].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[2].set_title('∂S/∂w₀[0]', fontsize=12)
axes2[2].grid(True, which='both', alpha=0.3)
axes2[2].legend(fontsize=9)
axes2[2].invert_xaxis()

plt.tight_layout()

# Figure 3: Species derivatives
fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
fig3.suptitle('Species Concentration Derivative Convergence (HP)', fontsize=14, fontweight='bold')

for i, (species_name, idx) in enumerate(zip(plot_species_names, plot_species_indices)):
    axes3[0].loglog(h_values, errors_dnj_dH[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)
    axes3[1].loglog(h_values, errors_dnj_dP[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)
    axes3[2].loglog(h_values, errors_dnj_dw0[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)

axes3[0].set_xlabel('Step size h', fontsize=11)
axes3[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes3[0].set_title('∂nⱼ/∂H', fontsize=12)
axes3[0].grid(True, which='both', alpha=0.3)
axes3[0].legend(fontsize=9, loc='best')
axes3[0].invert_xaxis()

axes3[1].set_xlabel('Step size h', fontsize=11)
axes3[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes3[1].set_title('∂nⱼ/∂P', fontsize=12)
axes3[1].grid(True, which='both', alpha=0.3)
axes3[1].legend(fontsize=9, loc='best')
axes3[1].invert_xaxis()

axes3[2].set_xlabel('Step size h', fontsize=11)
axes3[2].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes3[2].set_title('∂nⱼ/∂w₀[0]', fontsize=12)
axes3[2].grid(True, which='both', alpha=0.3)
axes3[2].legend(fontsize=9, loc='best')
axes3[2].invert_xaxis()

plt.tight_layout()

# Save figures
fig1.savefig('hp_temperature_derivative_convergence.pdf', dpi=300, bbox_inches='tight')
fig2.savefig('hp_entropy_derivative_convergence.pdf', dpi=300, bbox_inches='tight')
fig3.savefig('hp_species_derivative_convergence.pdf', dpi=300, bbox_inches='tight')

print("\n" + "=" * 70)
print("FIGURES SAVED")
print("=" * 70)
print("  - hp_temperature_derivative_convergence.pdf")
print("  - hp_entropy_derivative_convergence.pdf")
print("  - hp_species_derivative_convergence.pdf")
print("\n" + "=" * 70)
print("HP derivatives demonstration complete!")
print("=" * 70)

plt.show()
