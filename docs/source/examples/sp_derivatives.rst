Equilibrium Derivatives (SP)
============================

This example walks through ``source/bind/python/cea/samples/derivatives/sp_derivatives.py``.
It demonstrates how to compute analytic equilibrium derivatives for an SP
(entropy–pressure) problem and verify them against central finite differences.

Imports and initialisation
--------------------------

.. code-block:: python

    import numpy as np
    import cea

    cea.init()

Problem setup
-------------

Define reactant and product species, thermodynamic state, and mixture
composition. The setup is identical to the standard SP equilibrium workflow.

.. code-block:: python

    reac_names = ["H2", "Air"]
    prod_names = ["Ar",   "C",   "CO",  "CO2", "H",
                  "H2",   "H2O", "HNO", "HO2", "HNO2",
                  "HNO3", "N",   "NH",  "NO",  "N2",
                  "N2O3", "O",   "O2",  "OH",  "O3"]

    pressure = cea.units.atm_to_bar(1.0)
    fuel_moles    = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 1.0])

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(prod_names)

The :class:`~cea.EqSolver` is constructed with ``smooth_truncation=True``,
which replaces the hard species-threshold cutoff with a smooth ramp. This
produces continuously differentiable residuals and is required for accurate
analytic derivatives when trace species are present.

.. code-block:: python

    solver   = cea.EqSolver(prod, reactants=reac, smooth_truncation=True)
    solution = cea.EqSolution(solver)

Convert composition, compute the mixture entropy at a reference state (2000 K,
1 atm), then solve the SP problem at that entropy and pressure.

.. code-block:: python

    fuel_weights    = reac.moles_to_weights(fuel_moles)
    oxidant_weights = reac.moles_to_weights(oxidant_moles)
    of_ratio        = reac.chem_eq_ratio_to_of_ratio(
                          oxidant_weights, fuel_weights, 1.0)
    weights         = reac.of_ratio_to_weights(
                          oxidant_weights, fuel_weights, of_ratio)

    entropy = reac.calc_property(cea.ENTROPY, weights, 2000.0, pressure)

    # entropy is divided by cea.R to obtain the dimensionless form
    # expected by the SP solver
    solver.solve(solution, cea.SP, entropy / cea.R, pressure, weights)

    print(f"T, K      : {solution.T:.2f}")
    print(f"H, kJ/kg  : {solution.enthalpy:.4f}")

Computing analytic derivatives
-------------------------------

Construct an :class:`~cea.EqDerivatives` object from the converged solver and
solution, then call :meth:`~cea.EqDerivatives.compute_derivatives`.

.. code-block:: python

    derivs = cea.EqDerivatives(solver, solution)
    derivs.compute_derivatives(check_closure_defect=False)

For an SP problem the two thermodynamic state inputs are entropy (state1) and
pressure (state2). The derivative attributes follow this naming convention:

- ``dX_dstate1`` — derivative of output ``X`` with respect to S
- ``dX_dstate2`` — derivative of output ``X`` with respect to P
- ``dX_dw0``     — derivative of output ``X`` with respect to the reactant
  weight vector (one value per reactant species)

The outputs covered below are enthalpy (H), temperature (T), and the species
mole-number vector (nⱼ).

.. code-block:: python

    print("\n--- Analytic derivatives ---")
    print(f"dH/dS     : {derivs.dH_dstate1:.6e}")
    print(f"dH/dP     : {derivs.dH_dstate2:.6e}")
    print(f"dT/dS     : {derivs.dT_dstate1:.6e}")
    print(f"dT/dP     : {derivs.dT_dstate2:.6e}")
    print(f"dH/dw0[0] : {derivs.dH_dw0[0]:.6e}")
    print(f"dT/dw0[0] : {derivs.dT_dw0[0]:.6e}")

    species_names = prod.species_names
    for j, name in enumerate(species_names):
        print(f"dnj/dS [{name}] : {derivs.dnj_dstate1[j]:.6e}")
        print(f"dnj/dP [{name}] : {derivs.dnj_dstate2[j]:.6e}")

Finite-difference verification
--------------------------------

:meth:`~cea.EqDerivatives.compute_fd` perturbs each input and re-solves the
equilibrium problem to build finite-difference approximations. The key options
are:

- ``h``        — step size (dimensionless, applied as a relative perturbation
  to each input)
- ``central``  — use central differences (``True``) for second-order accuracy,
  or forward differences (``False``) for first-order accuracy
- ``verbose``  — print solver diagnostics for each perturbed solve

.. code-block:: python

    # Central differences at h = 1e-6 give second-order accuracy
    derivs.compute_fd(h=1e-6, central=True, verbose=False)

After the call the FD approximations are available under the same attribute
names with an ``_fd`` suffix (e.g. ``dH_dstate1_fd``).

.. code-block:: python

    def rel_err(a, b):
        return abs(b - a) / abs(a) if abs(a) > 1e-15 else abs(b - a)

    fmt = "{:<12} {:>14} {:>14} {:>12} {:>10}"
    print("\n--- Finite-difference verification (h=1e-6, central) ---")
    print(fmt.format("Derivative", "Analytic", "FD", "|Abs err|", "Rel err"))
    print("-" * 66)

    rows = [
        ("dH/dS",     derivs.dH_dstate1,  derivs.dH_dstate1_fd),
        ("dH/dP",     derivs.dH_dstate2,  derivs.dH_dstate2_fd),
        ("dT/dS",     derivs.dT_dstate1,  derivs.dT_dstate1_fd),
        ("dT/dP",     derivs.dT_dstate2,  derivs.dT_dstate2_fd),
        ("dH/dw0[0]", derivs.dH_dw0[0],   derivs.dH_dw0_fd[0]),
        ("dT/dw0[0]", derivs.dT_dw0[0],   derivs.dT_dw0_fd[0]),
    ]
    for name, a, fd in rows:
        abs_e = abs(fd - a)
        rel_e = rel_err(a, fd)
        print(fmt.format(name,
                         f"{a:.6e}", f"{fd:.6e}",
                         f"{abs_e:.2e}", f"{rel_e:.2e}"))

    print()
    print("dnj/dS and dnj/dP by species:")
    fmt2 = "  {:<6} {:>14} {:>14} {:>12} {:>10}"
    print(fmt2.format("Species", "dS analytic", "dS FD",
                      "|Abs err|", "Rel err"))
    print("  " + "-" * 58)
    for j, name in enumerate(species_names):
        a  = derivs.dnj_dstate1[j]
        fd = derivs.dnj_dstate1_fd[j]
        print(fmt2.format(name, f"{a:.6e}", f"{fd:.6e}",
                          f"{abs(fd-a):.2e}", f"{rel_err(a,fd):.2e}"))

At ``h=1e-6`` with central differences the relative error for well-resolved
derivatives (large mole fractions, large sensitivities) is typically in the
range ``1e-9`` to ``1e-11``. Species with very small mole fractions may show
larger relative errors because round-off in the perturbed solves dominates.
