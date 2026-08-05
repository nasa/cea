Example 8 from RP-1311
======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example8.py` in the CEA repository.

Here we describe how to run example 8 from RP-1311 [1]_ using the Python API.
This is a rocket problem assuming an infinite-area combustor (IAC), using H\ :sub:`2`\ (L) and O\ :sub:`2`\ (L) as the reactants.

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
    subar = [1.58]  # Subsonic area ratio
    supar = [25.0, 50.0, 75.0]  # Supersonic area ratio

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

    solver.solve(solution, weights, pc, pi_p, subar=subar, supar=supar, hc=hc, iac=True)

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

    print("Density, kg/m^3", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(rho[i]), end=" ")
        else:
            print("{0:10.3f}".format(rho[i]))

    print("H, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.2f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.2f}".format(enthalpy[i]))

    print("U, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.2f}".format(energy[i]), end=" ")
        else:
            print("{0:10.2f}".format(energy[i]))

    print("G, kJ/kg       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, kJ/kg-K     ", end=" ")
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

    print("Cp_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_eq[i]))

    print("Cp_fr, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_fr[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_fr[i]))

    print("Cv_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cv_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cv_eq[i]))

    print("Cv_eq, kJ/kg-K ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cv_fr[i]), end=" ")
        else:
            print("{0:10.3f}".format(cv_fr[i]))

    print("Gamma_s        ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(gamma_s[i]), end=" ")
        else:
            print("{0:10.3f}".format(gamma_s[i]))

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

    print()
    print("PERFORMANCE PARAMETERS")
    print()

    print("Ae/At          ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(ae_at[i]), end=" ")
        else:
            print("{0:10.3f}".format(ae_at[i]))

    print("C*, m/s        ", end=" ")
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

    print("Isp, vac., m/s ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(Isp_vac[i]), end=" ")
        else:
            print("{0:10.3f}".format(Isp_vac[i]))

    print("Isp, m/s       ", end=" ")
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
        if np.any(solution.mole_fractions[prod] > 5e-6):
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

This results in the following output to the terminal:

.. code-block:: console

    P, bar              53.317     30.655      5.332      0.533      0.053     48.382      0.205      0.081      0.047
    T, K              3383.845   3185.673   2567.340   1759.893   1115.867   3348.687   1468.163   1219.613   1088.640
    Density, kg/m^3      2.410      1.486      0.328      0.048      0.008      2.214      0.022      0.011      0.007
    H, kJ/kg          -1026.05   -2208.65   -5428.75   -8561.51  -10621.81   -1239.65   -9531.43  -10310.32  -10702.19
    U, kJ/kg          -3238.68   -4271.01   -7055.36   -9669.56  -11324.30   -3425.26  -10455.71  -11078.12  -11387.53
    G, kJ/kg          -64163.7   -61648.7   -53331.6   -41398.6   -31442.3   -63721.3   -36925.2   -33066.5   -31014.6
    S, kJ/kg-K          18.659     18.659     18.659     18.659     18.659     18.659     18.659     18.659     18.659
    M, (1/n)            12.716     12.843     13.123     13.206     13.207     12.739     13.207     13.207     13.207
    MW                  12.716     12.843     13.123     13.206     13.207     12.739     13.207     13.207     13.207
    Cp_eq, kJ/kg-K       8.325      7.480      4.880      3.435      2.963      8.182      3.223      3.042      2.942
    Cp_fr, kJ/kg-K       3.934      3.895      3.739      3.399      2.963      3.927      3.220      3.042      2.942
    Cv_eq, kJ/kg-K       7.130      6.427      4.149      2.803      2.333      7.013      2.594      2.412      2.312
    Cv_eq, kJ/kg-K       3.280      3.248      3.106      2.769      2.333      3.275      2.590      2.412      2.312
    Gamma_s              1.145      1.147      1.172      1.225      1.270      1.145      1.243      1.261      1.272
    Son. vel., m/s     1591.47    1537.92    1380.99    1165.22     944.48    1581.88    1071.77     983.95     933.79
    Mach                 0.000      1.000      2.149      3.332      4.638      0.413      3.848      4.379      4.711

    PERFORMANCE PARAMETERS

    Ae/At                0.000      1.000      2.350     12.238     68.753      1.580     25.000     50.000     75.000
    C*, m/s            2332.34    2332.34    2332.34    2332.34    2332.34    2332.34    2332.34    2332.34    2332.34
    Cf                  0.0000     0.6594     1.2723     1.6645     1.8783     0.2802     1.7684     1.8476     1.8861
    Isp, vac., m/s       0.000   2878.925   3515.549   4167.551   4541.167   3997.593   4348.510   4487.303   4554.913
    Isp, m/s             0.000   1537.917   2967.386   3882.127   4380.811    653.592   4124.410   4309.122   4399.121

    MOLE FRACTIONS

    H                 0.033498   0.026523  0.0079332 0.00019077 8.6863e-08   0.032259 1.4438e-05 5.4374e-07  5.048e-08
    H2                 0.29479    0.29432    0.29711     0.3004    0.30052    0.29466    0.30051    0.30052    0.30052
    H2O                0.63456     0.6528     0.6903    0.69937    0.69948    0.63794    0.69948    0.69948    0.69948
    H2O2            5.6145e-06 2.6541e-06 1.0707e-07 4.5473e-11 3.3916e-17 4.9535e-06 3.4561e-13 8.8625e-16 1.3014e-17
    HO2             1.4937e-05 6.7089e-06 1.6718e-07 9.0802e-12 5.8316e-20  1.309e-05 1.4613e-14 4.8736e-18 1.5828e-20
    O                0.0020678  0.0012027 7.1141e-05 1.3226e-08 3.1274e-16  0.0018954 3.6204e-11 2.0422e-14 9.0861e-17
    O2               0.0017218  0.0010431 6.7407e-05 1.3699e-08 3.8095e-16  0.0015903 3.9534e-11  2.389e-14 1.1204e-16
    OH                0.033341   0.024095  0.0045231 3.0655e-05 1.2435e-09   0.031643 1.0169e-06 1.3724e-08  6.115e-10

    TRACE SPECIES:
    O3              H2O(L)          H2O(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
