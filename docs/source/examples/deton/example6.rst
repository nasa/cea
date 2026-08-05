Example 6 from RP-1311
======================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example6.py` in the CEA repository.

Here we describe how to run example 6 from RP-1311 [1]_ using the Python API.
This example is a detonation problem with H\ :sub:`2`\  and O\ :sub:`2`\  as reactants, and includes computing transport properties.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions.

Declare the reactant names.

.. code-block:: python

    reac_names = ["H2", "O2"]

Define the temperature and pressure of the initial reactant mixture, in SI units.

.. code-block:: python

    p1 = np.array([1.0, 20.0]) # Initial pressure (bar)
    T1 = np.array([298.15, 500.0])  # Initial temperature (K)

Define the amounts of each reactant; in this case, a chemical equivalence ratio `eq_ratios` is prescribed (:math:`r_{eq}` in the RP-1311 [2]_).
The arrays `fuel_weights` and `oxidant_weights` correspond to the `reac_names` list, and set the weight fraction of each that is part of the fuel and oxidant mixtures, respectively.
In this case, because we are using one fuel and one oxidizer, these are simply `1.0` to indicate which reactant is the fuel and which is the oxidizer.

.. code-block:: python

    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    eq_ratios = np.array([1.0])

Now instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

Now we instantiate the :class:`~cea.DetonationSolver` and :class:`~cea.DetonationSolution` objects.
Note that this is the point where we set the `transport` flag equal to `True`, so that transport properties will be computed when we solve the problem later.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.DetonationSolver(prod, reactants=reac, transport=True)
    solution = cea.DetonationSolution(solver)

Next, we convert the `chem_eq_ratios` to oxidant-to-fuel ratios.

.. code-block:: python

    of_ratios = len(eq_ratios)*[0.0]
    for i, eqrat in enumerate(eq_ratios):
        of_ratios[i] = reac.chem_eq_ratio_to_of_ratio(oxidant_weights,
                                                      fuel_weights,
                                                      eqrat)

We will now initialize an array to store each of the solution variables for printing the output later.

.. code-block:: python

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

Now we will loop over the arrays of temperature, pressure, and o/f ratio, and solve the detonation problem at each condition. Most of the lines of code here are storing the output for printing later.
Note that the :meth:`~cea.DetonationSolver.solve` method strictly requires an array of weight fractions, so we first compute that array using the reactant mixture object.

.. code-block:: python
    :emphasize-lines: 2, 5, 6

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

Finally, print everything out in a formatted manner consistent with the legacy CEA output format.

.. code-block:: python

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

This results in the following output to the terminal:

.. code-block:: console

    UNBURNED GAS

    o/f                  7.937      7.937      7.937      7.937
    P1, atm              0.987     19.738      0.987     19.738
    T1, K              298.150    298.150    500.000    500.000
    H1, cal/g           -0.000     -0.000    118.410    118.410
    M1 (1/n)            12.010     12.010     12.010     12.010
    Gamma1               1.402      1.402      1.386      1.386
    Son. Vel.1, m/s    537.868    537.868    692.592    692.592

    BURNED GAS

    P, atm              18.523    409.027     10.813    240.197
    T, K              3674.283   4283.382   3600.034   4209.491
    Density, g/cc    8.908e-04  1.775e-02  5.219e-04  1.042e-02
    H, cal/g           676.647    752.006    758.260    836.596
    U, cal/g           173.085    194.015    256.508    278.369
    G, cal/g          -14625.8   -15395.6   -14583.5   -15410.9
    S, cal/g-K           4.165      3.770      4.262      3.860
    M, (1/n)            14.500     15.255     14.258     14.985
    MW                  14.500     15.255     14.258     14.985
    Gamma_s             1.1288     1.1438     1.1266     1.1423
    Cp_eq, cal/g-K      3.8874     2.4544     4.3162     2.7219
    Cp_fr, cal/g-K      0.7784     0.7936     0.7765     0.7920
    Cv_eq, cal/g-K      3.1830     2.0235     3.5188     2.2326
    Cv_fr, cal/g-K      0.6413     0.6633     0.6372     0.6593
    Son. Vel., m/s   1542.1580  1634.1166  1537.8972  1633.4041

    Transport properties:

    Viscosity, mP       1.1412     1.2740     1.1244     1.2589

    with equilibrium reaction:

    Cp, cal/g-K         3.8874     2.4544     4.3162     2.7219
    Cv, cal/g-K         3.1830     2.0235     3.5188     2.2326
    Conductivity        9.1031     5.9556    10.0641     6.6597
    Prandtl number      0.4873     0.5250     0.4822     0.5145

    with frozen reaction:

    Cp, cal/g-K         0.7784     0.7936     0.7765     0.7920
    Cv, cal/g-K         0.6413     0.6633     0.6372     0.6593
    Conductivity        1.2903     1.4211     1.2823     1.4170
    Prandtl number      0.6884     0.7115     0.6809     0.7036

    DETONATION PARAMETERS

    P/P1               18.7688    20.7227    10.9569    12.1691
    T/T1               12.3236    14.3665     7.2001     8.4190
    M/M1                1.2073     1.2701     1.1872     1.2477
    rho/rho1            1.8387     1.8321     1.8066     1.8035
    Mach                5.2718     5.5661     4.0115     4.2533
    Velocity         2835.5321  2993.8052  2778.3142  2945.7966

    MOLE FRACTIONS

    H                 0.080079   0.047166   0.090963    0.05644
    H2                 0.16211    0.14397    0.16681    0.15214
    H2O                0.53189    0.60992    0.50727    0.57898
    H2O2             2.004e-05 0.00016339 1.3838e-05 0.00011628
    HO2             0.00018402  0.0006771 0.00014829 0.00056843
    O                 0.037436    0.02347   0.042144   0.027923
    O2                0.046851   0.037051   0.049053   0.039695
    OH                 0.14143    0.13758     0.1436    0.14414

    TRACE SPECIES:
    O3              H2O(L)          H2O(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
.. [2] Gordon, S., McBride, B.J., "Computer program for calculation of complex chemical equilibrium compositions and applications. Part 1: Analysis",
    NASA RP-1311, 1994. [NTRS](https://ntrs.nasa.gov/citations/19950013764)
