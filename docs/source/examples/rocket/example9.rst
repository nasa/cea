Example 9 from RP-1311
======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example9.py` in the CEA repository.

Here we describe how to run example 9 from RP-1311 [1]_ using the Python API.
This is a rocket problem assuming a finite-area combustor (FAC), using the contraction ratio :math:`A_c /A_t` to define the combustor area.
The reactants are H\ :sub:`2`\ (L) and O\ :sub:`2`\ (L) , and aside from being an FAC rocket problem, the example is otherwise similar to :doc:`Example 8 <example8>`.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :data:`cea.R` inline when normalizing enthalpy.

Declare the reactants and set their amounts and initial temperatures.
The initial reactant temperatures, `T_reactant`, will be used later to compute the chamber enthalpy.
The amounts of each are specified through the `fuel_weights`, `oxidant_weights`, and `of_ratio` variables specified here.
Setting the `fuel_weights` array equal to `[1.0, 0.0]` means that `H2(L)` constitutes 100\% of the fuel, and similarly,
setting the `oxidant_weights` array equal to `[0.0, 1.0]` means that `O2(L)` constitutes 100\% of the oxidant.
These values will be used in conjunction later with `of_ratio` to compute the overall weight fraction array of the reactant mixture.

.. code-block:: python

    reac_names = ["H2(L)", "O2(L)"]
    T_reactant = np.array([20.27, 90.17])  # Reactant temperatures (K)
    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    of_ratio = 5.55157

Next, set some states for the rocket analysis. We will pass these values into the :class:`~cea.RocketSolver` later.

.. code-block:: python

    pc = 53.3172  # Chamber pressure (bar)
    pi_p = [10.0, 100.0, 1000.0]   # Pressure ratio
    supar = [25.0, 50.0, 75.0]  # Supersonic area ratio
    ac_at = 1.58  # Area ratio chamber to throat

Instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

Now instantiate the :class:`~cea.RocketSolver` and :class:`~cea.RocketSolution` objects.

.. code-block:: python

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)

Now we will use the reactant :class:`~cea.Mixture` object to compute the overall weight fraction array of the reactants:

.. code-block:: python

    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

And compute the chamber enthalpy value based on the reactant weights and temperatures. We will pass this in later when we call :meth:`~cea.RocketSolver.solve`.
Note that this value is normalized by `R` here.

.. code-block:: python

    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R

Now we can solve the :meth:`~cea.RocketSolver.solve` function:

.. code-block:: python

    solver.solve(solution, weights, pc, pi_p, supar=supar, ac_at=ac_at, iac=False, hc=hc)

Finally, query the solution variables and print them out:

.. code-block:: python

    num_pts = solution.num_pts
    T = solution.T
    P = solution.P
    rho = solution.density
    enthalpy = solution.enthalpy
    energy = solution.energy
    gibbs = solution.gibbs_energy
    entropy = solution.entropy
    n = solution.n
    M_1n = solution.M
    MW = solution.MW
    cp_eq = solution.cp_eq
    cp_fr = solution.cp_fr
    cv_eq = solution.cv_eq
    cv_fr = solution.cv_fr
    Mach = solution.Mach
    gamma_s = solution.gamma_s
    v_sonic = solution.sonic_velocity
    ae_at = solution.ae_at
    c_star = solution.c_star
    Cf = solution.coefficient_of_thrust
    Isp = solution.Isp
    Isp_vac = solution.Isp_vacuum

    print("P, bar         ", end=" ")
    for i in range(num_pts):
        if i == 1:  # Skip solution at infinity
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(P[i]), end=" ")
        else:
            print("{0:10.3f}".format(P[i]))

    print("T, K           ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(T[i]), end=" ")
        else:
            print("{0:10.3f}".format(T[i]))

    print("Density, kg/m^3", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(rho[i]), end=" ")
        else:
            print("{0:10.3f}".format(rho[i]))

    print("H, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.2f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.2f}".format(enthalpy[i]))

    print("U, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.2f}".format(energy[i]), end=" ")
        else:
            print("{0:10.2f}".format(energy[i]))

    print("G, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, kJ/kg-K     ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(entropy[i]), end=" ")
        else:
            print("{0:10.3f}".format(entropy[i]))

    print("M, (1/n)       ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(M_1n[i]), end=" ")
        else:
            print("{0:10.3f}".format(M_1n[i]))

    print("MW             ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(MW[i]), end=" ")
        else:
            print("{0:10.3f}".format(MW[i]))

    print("Cp_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_eq[i]))

    print("Cp_fr, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_fr[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_fr[i]))

    print("Cv_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(cv_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cv_eq[i]))

    print("Cv_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(cv_fr[i]), end=" ")
        else:
            print("{0:10.3f}".format(cv_fr[i]))

    print("Gamma_s        ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(gamma_s[i]), end=" ")
        else:
            print("{0:10.3f}".format(gamma_s[i]))

    print("Son. vel., m/s ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.2f}".format(v_sonic[i]), end=" ")
        else:
            print("{0:10.2f}".format(v_sonic[i]))

    print("Mach           ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(Mach[i]), end=" ")
        else:
            print("{0:10.3f}".format(Mach[i]))

    print()
    print("PERFORMANCE PARAMETERS")
    print()

    print("Ae/At          ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(ae_at[i]), end=" ")
        else:
            print("{0:10.3f}".format(ae_at[i]))

    print("C*, m/s        ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.2f}".format(c_star[i]), end=" ")
        else:
            print("{0:10.2f}".format(c_star[i]))

    print("Cf             ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(Cf[i]), end=" ")
        else:
            print("{0:10.3f}".format(Cf[i]))

    print("Isp, vac., m/s ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(Isp_vac[i]), end=" ")
        else:
            print("{0:10.3f}".format(Isp_vac[i]))

    print("Isp, m/s       ", end=" ")
    for i in range(num_pts):
        if i == 1:
            continue
        if i < num_pts-1:
            print("{0:10.3f}".format(Isp[i]), end=" ")
        else:
            print("{0:10.3f}".format(Isp[i]))

    print()

    print()
    print("MOLE FRACTIONS")
    print("")
    trace_species = []
    for prod in solution.mole_fractions:
        if np.any(solution.mole_fractions[prod] > 5e-6):
            print("{0:15s}".format(prod), end=" ")
            for j in range(len(solution.mole_fractions[prod])):
                if j == 1:
                    continue
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

This results in the following output to the terminal:

.. code-block:: console

    P, bar              53.317     44.531     28.221      5.332      0.533      0.053      0.189      0.075      0.044
    T, K              3383.845   3340.917   3179.445   2595.639   1787.832   1136.449   1469.622   1220.899   1089.829
    Density, kg/m^3      2.410      2.041      1.370      0.324      0.047      0.007      0.020      0.010      0.006
    H, kJ/kg          -1026.05   -1239.22   -2206.40   -5288.60   -8465.21  -10560.67   -9526.72  -10306.41  -10698.69
    U, kJ/kg          -3238.68   -3421.46   -4266.11   -6934.60   -9590.89  -11276.12  -10451.92  -11075.02  -11384.79
    G, kJ/kg          -64163.7   -63757.3   -61702.9   -53860.4   -41920.6   -31826.9   -37027.5   -33152.9   -31092.5
    S, kJ/kg-K          18.659     18.713     18.713     18.713     18.713     18.713     18.713     18.713     18.713
    M, (1/n)            12.716     12.729     12.835     13.111     13.205     13.207     13.207     13.207     13.207
    MW                  12.716     12.729     12.835     13.111     13.205     13.207     13.207     13.207     13.207
    Cp_eq, kJ/kg-K       8.325      8.296      7.582      5.027      3.459      2.978      3.225      3.043      2.943
    Cp_fr, kJ/kg-K       3.934      3.926      3.894      3.748      3.414      2.978      3.221      3.043      2.943
    Cv_eq, kJ/kg-K       7.130      7.111      6.517      4.282      2.827      2.349      2.595      2.413      2.313
    Cv_eq, kJ/kg-K       3.280      3.273      3.246      3.114      2.784      2.349      2.591      2.413      2.313
    Gamma_s              1.145      1.144      1.146      1.170      1.224      1.268      1.243      1.261      1.272
    Son. vel., m/s     1591.47    1580.21    1536.46    1387.56    1173.64     952.47    1072.25     984.43     934.26
    Mach                 0.000      0.413      1.000      2.104      3.287      4.585      3.845      4.376      4.708

    PERFORMANCE PARAMETERS

    Ae/At                0.000      1.580      1.000      2.226     11.523     64.689     25.000     50.000     75.000
    C*, m/s            2330.97    2330.97    2330.97    2330.97    2330.97    2330.97    2330.97    2330.97    2330.97
    Cf                   0.000      0.280      0.659      1.253      1.655      1.873      1.769      1.848      1.887
    Isp, vac., m/s       0.000   3995.119   2877.015   3483.517   4149.076   4530.669   4347.652   4486.620   4554.316
    Isp, m/s             0.000    652.938   1536.457   2919.774   3857.241   4366.833   4123.266   4308.214   4398.326


    MOLE FRACTIONS

    H                 0.033498   0.032976    0.02717  0.0089113 0.00024284 1.3423e-07 1.5313e-05 5.7955e-07 5.4001e-08
    H2                 0.29479    0.29453    0.29418    0.29676    0.30037    0.30052    0.30051    0.30052    0.30052
    H2O                0.63456    0.63671    0.65174    0.68892    0.69934    0.69948    0.69948    0.69948    0.69948
    H2O2            5.6145e-06 4.7882e-06 2.5763e-06 1.2856e-07 6.6819e-11 6.7889e-17 3.5584e-13 9.1967e-16 1.3584e-17
    HO2             1.4937e-05  1.303e-05 6.7146e-06 2.1669e-07  1.565e-11 1.5642e-19  1.586e-14 5.3494e-18 1.7518e-20
    O                0.0020678  0.0019707  0.0012568 9.2139e-05 2.2701e-08 8.3271e-16 4.0908e-11 2.3335e-14 1.0468e-16
    O2               0.0017218  0.0016517  0.0010892 8.7033e-05 2.3416e-08 1.0055e-15 4.4656e-11 2.7286e-14 1.2901e-16
    OH                0.033341   0.032151   0.024554  0.0052272 4.1616e-05 2.1685e-09 1.0839e-06 1.4722e-08 6.5905e-10

    TRACE SPECIES:
    O3              H2O(L)          H2O(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
