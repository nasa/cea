import numpy as np
import matplotlib.pyplot as plt
import cea

"""
Example demonstrating equilibrium derivatives with finite-difference convergence study
- Single SV equilibrium calculation
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

# Species
reac_names = ["H2", "Air"]
prod_names = ["Ar",   "C",   "CO",  "CO2", "H",
              "H2",   "H2O", "HNO", "HO2", "HNO2",
              "HNO3", "N",   "NH",  "NO",  "N2",
              "N2O3", "O",   "O2",  "OH",  "O3"]

# Mixture composition
fuel_moles = np.array([1.0, 0.0])
oxidant_moles = np.array([0.0, 1.0])
chem_eq_ratio = 1.0

# Mixtures
reac = cea.Mixture(reac_names)
prod = cea.Mixture(prod_names)

# Solver
solver = cea.EqSolver(prod, reactants=reac, smooth_truncation=True)
solution = cea.EqSolution(solver)

# Unit conversions
fuel_weights = reac.moles_to_weights(fuel_moles)
oxidant_weights = reac.moles_to_weights(oxidant_moles)
of_ratio = reac.chem_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, chem_eq_ratio)
weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

# Initial conditions
temperature_init = 2000.0
density = 1.0e3 * 9.1864e-5  # kg/m^3
volume = 1.0 / density  # m^3/kg
pressure_init = (1e-5)*0.1*cea.R*temperature_init/volume
entropy = reac.calc_property(cea.ENTROPY, weights, temperature_init, pressure_init)

# Now solve the SV problem using the extracted S and V
solver.solve(solution, cea.SV, entropy/cea.R, volume, weights)

if not solution.converged:
    print("ERROR: SV solution did not converge")
    exit(1)

print("\n" + "=" * 70)
print("SV EQUILIBRIUM SOLUTION")
print("=" * 70)
print(f"o/f ratio       : {of_ratio:.4f}")
print(f"S, kJ/(kg·K)    : {entropy:.6f}")
print(f"V, m³/kg        : {volume:.6e}")
print(f"T, K            : {solution.T:.2f}")
print(f"P, atm          : {cea.units.bar_to_atm(solution.P):.4f}")
print(f"H, kJ/kg        : {solution.enthalpy:.4f}")
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

# Get analytic derivatives (for SV: state1=S, state2=V)
dH_dS_analytic = derivs.dH_dstate1
dH_dV_analytic = derivs.dH_dstate2
dH_dw0_analytic = derivs.dH_dw0

dT_dS_analytic = derivs.dT_dstate1
dT_dV_analytic = derivs.dT_dstate2
dT_dw0_analytic = derivs.dT_dw0

dnj_dS_analytic = derivs.dnj_dstate1
dnj_dV_analytic = derivs.dnj_dstate2
dnj_dw0_analytic = derivs.dnj_dw0

print(f"dH/dS (analytic)       : {dH_dS_analytic:.6e} kJ/kg / [kJ/(kg·K)]")
print(f"dH/dV (analytic)       : {dH_dV_analytic:.6e} kJ/kg / [m³/kg]")
print(f"dT/dS (analytic)       : {dT_dS_analytic:.6e} K / [kJ/(kg·K)]")
print(f"dT/dV (analytic)       : {dT_dV_analytic:.6e} K / [m³/kg]")

print("\n" + "=" * 70)
print("FINITE-DIFFERENCE CONVERGENCE STUDY")
print("=" * 70)
print(f"Error metric: {'ABSOLUTE' if USE_ABSOLUTE_ERROR else 'RELATIVE'}")

# Define step sizes for convergence study
h_values = np.logspace(-1, -16, 50)
print(f"Testing {len(h_values)} step sizes from {h_values[0]:.1e} to {h_values[-1]:.1e}")

# Storage for errors
errors_dH_dS = []
errors_dH_dV = []
errors_dH_dw0 = []
errors_dT_dS = []
errors_dT_dV = []
errors_dT_dw0 = []
errors_dnj_dS = []
errors_dnj_dV = []
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

    # Enthalpy derivatives
    dH_dS_fd = derivs.dH_dstate1_fd
    dH_dV_fd = derivs.dH_dstate2_fd
    dH_dw0_fd = derivs.dH_dw0_fd

    errors_dH_dS.append(rel_error(dH_dS_analytic, dH_dS_fd))
    errors_dH_dV.append(rel_error(dH_dV_analytic, dH_dV_fd))
    errors_dH_dw0.append(rel_error(dH_dw0_analytic[0], dH_dw0_fd[0]))  # First component

    # Temperature derivatives
    dT_dS_fd = derivs.dT_dstate1_fd
    dT_dV_fd = derivs.dT_dstate2_fd
    dT_dw0_fd = derivs.dT_dw0_fd

    errors_dT_dS.append(rel_error(dT_dS_analytic, dT_dS_fd))
    errors_dT_dV.append(rel_error(dT_dV_analytic, dT_dV_fd))
    errors_dT_dw0.append(rel_error(dT_dw0_analytic[0], dT_dw0_fd[0]))  # First component

    # Species derivatives
    dnj_dS_fd = derivs.dnj_dstate1_fd
    dnj_dV_fd = derivs.dnj_dstate2_fd
    dnj_dw0_fd = derivs.dnj_dw0_fd

    # Store errors for selected species
    errs_dS = [rel_error(dnj_dS_analytic[idx], dnj_dS_fd[idx]) for idx in plot_species_indices]
    errs_dV = [rel_error(dnj_dV_analytic[idx], dnj_dV_fd[idx]) for idx in plot_species_indices]
    errs_dw0 = [rel_error(dnj_dw0_analytic[idx, 0], dnj_dw0_fd[idx, 0]) for idx in plot_species_indices]

    errors_dnj_dS.append(errs_dS)
    errors_dnj_dV.append(errs_dV)
    errors_dnj_dw0.append(errs_dw0)

# Convert to numpy arrays
errors_dH_dS = np.array(errors_dH_dS)
errors_dH_dV = np.array(errors_dH_dV)
errors_dH_dw0 = np.array(errors_dH_dw0)
errors_dT_dS = np.array(errors_dT_dS)
errors_dT_dV = np.array(errors_dT_dV)
errors_dT_dw0 = np.array(errors_dT_dw0)
errors_dnj_dS = np.array(errors_dnj_dS)
errors_dnj_dV = np.array(errors_dnj_dV)
errors_dnj_dw0 = np.array(errors_dnj_dw0)

print("\nConvergence study complete!")

# Print error statistics
print("\n" + "=" * 70)
print("ERROR STATISTICS")
print("=" * 70)
print("Note: Relative errors can be >> 1.0 when FD approximation is poor.")
print("This is expected - the plots show the characteristic U-shaped curve.\n")

print("∂H/∂S errors:")
print(f"  Min: {np.min(errors_dH_dS):.2e}  Max: {np.max(errors_dH_dS):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dH_dS)]:.1e}")

print("∂H/∂V errors:")
print(f"  Min: {np.min(errors_dH_dV):.2e}  Max: {np.max(errors_dH_dV):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dH_dV)]:.1e}")

print("∂T/∂S errors:")
print(f"  Min: {np.min(errors_dT_dS):.2e}  Max: {np.max(errors_dT_dS):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dT_dS)]:.1e}")

print("∂T/∂V errors:")
print(f"  Min: {np.min(errors_dT_dV):.2e}  Max: {np.max(errors_dT_dV):.2e}")
print(f"  Optimal h: {h_values[np.argmin(errors_dT_dV)]:.1e}")

print("\nGenerating plots...")

# Create plots
plt.style.use('default')

# Figure 1: Enthalpy derivatives
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5))
fig1.suptitle('Enthalpy Derivative Convergence (SV)', fontsize=14, fontweight='bold')

# Find optimal points
opt_idx_dH_dS = np.argmin(errors_dH_dS)
opt_idx_dH_dV = np.argmin(errors_dH_dV)
opt_idx_dH_dw0 = np.argmin(errors_dH_dw0)

axes1[0].loglog(h_values, errors_dH_dS, 'o-', linewidth=2, markersize=4, label='Error')
axes1[0].axvline(h_values[opt_idx_dH_dS], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dH_dS]:.1e}')
axes1[0].set_xlabel('Step size h', fontsize=11)
axes1[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[0].set_title('∂H/∂S', fontsize=12)
axes1[0].grid(True, which='both', alpha=0.3)
axes1[0].legend(fontsize=9)
axes1[0].invert_xaxis()

axes1[1].loglog(h_values, errors_dH_dV, 'o-', linewidth=2, markersize=4, label='Error')
axes1[1].axvline(h_values[opt_idx_dH_dV], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dH_dV]:.1e}')
axes1[1].set_xlabel('Step size h', fontsize=11)
axes1[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[1].set_title('∂H/∂V', fontsize=12)
axes1[1].grid(True, which='both', alpha=0.3)
axes1[1].legend(fontsize=9)
axes1[1].invert_xaxis()

axes1[2].loglog(h_values, errors_dH_dw0, 'o-', linewidth=2, markersize=4, label='Error')
axes1[2].axvline(h_values[opt_idx_dH_dw0], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dH_dw0]:.1e}')
axes1[2].set_xlabel('Step size h', fontsize=11)
axes1[2].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes1[2].set_title('∂H/∂w₀[0]', fontsize=12)
axes1[2].grid(True, which='both', alpha=0.3)
axes1[2].legend(fontsize=9)
axes1[2].invert_xaxis()

plt.tight_layout()

# Figure 2: Temperature derivatives
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
fig2.suptitle('Temperature Derivative Convergence (SV)', fontsize=14, fontweight='bold')

# Find optimal points
opt_idx_dT_dS = np.argmin(errors_dT_dS)
opt_idx_dT_dV = np.argmin(errors_dT_dV)
opt_idx_dT_dw0 = np.argmin(errors_dT_dw0)

axes2[0].loglog(h_values, errors_dT_dS, 'o-', linewidth=2, markersize=4, label='Error')
axes2[0].axvline(h_values[opt_idx_dT_dS], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dS]:.1e}')
axes2[0].set_xlabel('Step size h', fontsize=11)
axes2[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[0].set_title('∂T/∂S', fontsize=12)
axes2[0].grid(True, which='both', alpha=0.3)
axes2[0].legend(fontsize=9)
axes2[0].invert_xaxis()

axes2[1].loglog(h_values, errors_dT_dV, 'o-', linewidth=2, markersize=4, label='Error')
axes2[1].axvline(h_values[opt_idx_dT_dV], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dV]:.1e}')
axes2[1].set_xlabel('Step size h', fontsize=11)
axes2[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[1].set_title('∂T/∂V', fontsize=12)
axes2[1].grid(True, which='both', alpha=0.3)
axes2[1].legend(fontsize=9)
axes2[1].invert_xaxis()

axes2[2].loglog(h_values, errors_dT_dw0, 'o-', linewidth=2, markersize=4, label='Error')
axes2[2].axvline(h_values[opt_idx_dT_dw0], color='r', linestyle='--', alpha=0.7, linewidth=1.5, label=f'Optimal h={h_values[opt_idx_dT_dw0]:.1e}')
axes2[2].set_xlabel('Step size h', fontsize=11)
axes2[2].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes2[2].set_title('∂T/∂w₀[0]', fontsize=12)
axes2[2].grid(True, which='both', alpha=0.3)
axes2[2].legend(fontsize=9)
axes2[2].invert_xaxis()

plt.tight_layout()

# Figure 3: Species derivatives
fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
fig3.suptitle('Species Concentration Derivative Convergence (SV)', fontsize=14, fontweight='bold')

for i, (species_name, idx) in enumerate(zip(plot_species_names, plot_species_indices)):
    axes3[0].loglog(h_values, errors_dnj_dS[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)
    axes3[1].loglog(h_values, errors_dnj_dV[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)
    axes3[2].loglog(h_values, errors_dnj_dw0[:, i], 'o-', linewidth=2, markersize=4,
                    label=species_name, alpha=0.8)

axes3[0].set_xlabel('Step size h', fontsize=11)
axes3[0].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes3[0].set_title('∂nⱼ/∂S', fontsize=12)
axes3[0].grid(True, which='both', alpha=0.3)
axes3[0].legend(fontsize=9, loc='best')
axes3[0].invert_xaxis()

axes3[1].set_xlabel('Step size h', fontsize=11)
axes3[1].set_ylabel('Absolute error' if USE_ABSOLUTE_ERROR else 'Relative error', fontsize=11)
axes3[1].set_title('∂nⱼ/∂V', fontsize=12)
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
fig1.savefig('sv_enthalpy_derivative_convergence.pdf', dpi=300, bbox_inches='tight')
fig2.savefig('sv_temperature_derivative_convergence.pdf', dpi=300, bbox_inches='tight')
fig3.savefig('sv_species_derivative_convergence.pdf', dpi=300, bbox_inches='tight')

print("\n" + "=" * 70)
print("FIGURES SAVED")
print("=" * 70)
print("  - sv_enthalpy_derivative_convergence.pdf")
print("  - sv_temperature_derivative_convergence.pdf")
print("  - sv_species_derivative_convergence.pdf")
print("\n" + "=" * 70)
print("SV derivatives demonstration complete!")
print("=" * 70)

plt.show()
