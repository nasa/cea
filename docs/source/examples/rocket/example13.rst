Example 13 from RP-1311
=======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example13.py` in the CEA repository.

Here we describe how to run example 13 from RP-1311 [1]_ using the Python API.

This is a rocket problem assuming an infinite-area combustor (IAC), using N\ :sub:`2`\ H\ :sub:`4`\ (L), Be(a), and H\ :sub:`2`\ O\ :sub:`2`\ (L) as reactants.
This problem demonstrates how to call the `insert` parameter for a rocket problem, adding BeO(L) to the initial guess for the equilibrium evaluation.
The resulting equilibrium mixtures also include non-negligible amounts of condensed species.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions and :data:`cea.R` inline when normalizing enthalpy.

Declare the reactants and set their amounts and initial temperatures.
The initial reactant temperatures, `reac_T`, will be used later to compute the chamber enthalpy.
The amounts of each are specified through the `fuel_weights`, `oxidant_weights`, and `pct_fuel` variables specified here.
Setting the `fuel_weights` array equal to `[0.8, 0.2, 0.0]` means that N\ :sub:`2`\ H\ :sub:`4`\ (L) constitutes 80\% of the fuel,
and Be(a) makes up the other 20\% of the fuel mixture by weight.
Similarly, setting the `oxidant_weights` array equal to `[0.0, 0.0, 1.0]` means that H\ :sub:`2`\ O\ :sub:`2`\ (L) constitutes 100\% of the oxidant.
These values will be used in conjunction later with `pct_fuel` to compute the overall weight fraction array of the reactant mixture.

.. code-block:: python

    reac_names = ["N2H4(L)", "Be(a)", "H2O2(L)"]
    fuel_weights = np.array([0.8, 0.2, 0.0])
    oxidant_weights = np.array([0.0, 0.0, 1.0])
    reac_T = np.array([298.15, 298.15, 298.15])  # Reactant temperatures (K)
    pct_fuel = 67.0

Next we set a variable for our trace parameter, defining the threshold below which species concentrations are assumed to be 0.

.. code-block:: python

    trace = 1e-10

Next, set some states for the rocket analysis. We will pass these values into the :class:`~cea.RocketSolver` later.

.. code-block:: python

    pc = cea.units.psi_to_bar(3000)  # Chamber pressure (bar)
    pi_p = [3.0, 10.0, 30.0, 300.0]   # Pressure ratio

Now we will instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.
We will also create the list of `insert` species that we will pass into the :class:`~cea.RocketSolver` next.

.. code-block:: python
    :emphasize-lines: 3

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    insert = ["BeO(L)"]

Now initialize the :class:`~cea.RocketSolver` and :class:`~cea.RocketSolution` class instances.
We pass in the optional arguments `insert` and `trace` to the :class:`~cea.RocketSolver` class initialization.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.RocketSolver(prod, reactants=reac, trace=trace, insert=insert)
    solution = cea.RocketSolution(solver)

Compute the array of weight fractions for the reactant mixture:

.. code-block:: python

    of_ratio = (100.0 - pct_fuel) / pct_fuel
    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

Compute the chamber enthalpy value based on the reactant weights and temperatures. We will pass this in later when we call :meth:`~cea.RocketSolver.solve`.
Note that this value is normalized by `R` here.

.. code-block:: python

    hc = reac.calc_property(cea.ENTHALPY, weights, reac_T)/cea.R

Now we can solve the :meth:`~cea.RocketSolver.solve` function:

.. code-block:: python

    solver.solve(solution, weights, pc, pi_p, iac=True, hc=hc)

Finally, query the solution variables and print them out:

.. code-block:: python

    num_pts = solution.num_pts
    T = solution.T
    P = cea.units.bar_to_atm(solution.P)
    rho = 1e-3*solution.density
    enthalpy = cea.units.joule_to_cal(solution.enthalpy)
    energy = cea.units.joule_to_cal(solution.energy)
    gibbs = cea.units.joule_to_cal(solution.gibbs_energy)
    entropy = cea.units.joule_to_cal(solution.entropy)
    n = solution.n
    M_1n = solution.M
    MW = solution.MW
    cp_eq = cea.units.joule_to_cal(solution.cp_eq)
    Mach = solution.Mach
    gamma_s = solution.gamma_s
    v_sonic = solution.sonic_velocity
    ae_at = solution.ae_at
    c_star = cea.units.m_per_s_to_ft_per_s(solution.c_star)
    Cf = solution.coefficient_of_thrust
    Isp = solution.Isp/9.80665
    Isp_vac = solution.Isp_vacuum/9.80665

    print("P, atm         ", end=" ")
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

    print("Density, g/cm^3", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3g}".format(rho[i]), end=" ")
        else:
            print("{0:10.3g}".format(rho[i]))

    print("H, cal/g       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.2f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.2f}".format(enthalpy[i]))

    print("U, cal/g       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.2f}".format(energy[i]), end=" ")
        else:
            print("{0:10.2f}".format(energy[i]))

    print("G, cal/g       ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, cal/g-K     ", end=" ")
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

    print("Cp, cal/g-K    ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_eq[i]))

    print("Gamma_s        ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(gamma_s[i]), end=" ")
        else:
            print("{0:10.4f}".format(gamma_s[i]))

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

    print("")
    print("PERFORMANCE PARAMETERS")
    print("")

    print("Ae/At          ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(ae_at[i]), end=" ")
        else:
            print("{0:10.4f}".format(ae_at[i]))

    print("C*, ft/s       ", end=" ")
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

    print("Ivac, lb-s/lb  ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(Isp_vac[i]), end=" ")
        else:
            print("{0:10.3f}".format(Isp_vac[i]))

    print("Isp, lb-s/lb   ", end=" ")
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
        if np.any(solution.mole_fractions[prod] > trace):
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

    P, atm             204.138    125.281     68.046     20.414      6.805      0.680
    T, K              3002.540   2851.000   2851.000   2449.698   2068.730   1398.106
    Density, g/cm^3     0.0138    0.00891    0.00483    0.00169   0.000669   9.91e-05
    H, cal/g           -234.01    -403.87    -611.80    -995.93   -1291.19   -1761.87
    U, cal/g           -592.87    -744.29    -952.71   -1287.95   -1537.45   -1928.24
    G, cal/g          -10072.6    -9745.9    -9953.9    -9023.0    -8069.9    -6343.1
    S, cal/g-K           3.277      3.277      3.277      3.277      3.277      3.277
    M, (1/n)            16.627     16.643     16.619     16.670     16.694     16.700
    MW                  13.431     13.420     13.405     13.378     13.377     13.378
    Cp, cal/g-K          0.974      0.000      0.000      0.813      0.747      0.665
    Gamma_s             1.1515     0.9979     0.9974     1.1791     1.1916     1.2180
    Son. vel., m/s     1314.87    1192.19    1192.73    1200.25    1108.03     920.76
    Mach                 0.000      1.000      1.491      2.104      2.684      3.883

    PERFORMANCE PARAMETERS

    Ae/At               0.0000     1.0000     1.2363     2.4856     5.3385    30.0004
    C*, ft/s           6386.75    6386.75    6386.75    6386.75    6386.75    6386.75
    Cf                  0.0000     0.6124     0.9133     1.2971     1.5279     1.8368
    Ivac, lb-s/lb        0.000    243.396    263.111    306.823    338.619    384.464
    Isp, lb-s/lb         0.000    121.571    181.306    257.481    303.294    364.613

    MOLE FRACTIONS

    Be              8.5544e-06 4.1033e-06 7.5534e-06 2.5505e-07 1.9172e-09 1.4393e-16
    Be(OH)2          0.0072678  0.0057344  0.0057215  0.0014065 0.00020281 4.5262e-07
    Be2O            1.1221e-06 4.0985e-07 7.5446e-07 5.5565e-09 5.3694e-12  5.606e-22
    Be2O2           1.8492e-07 9.0835e-08 1.6729e-07 3.1776e-09 9.8507e-12 3.6727e-20
    Be3O3           3.2166e-06 1.9697e-06 3.6277e-06  7.568e-08 2.5733e-10 1.1253e-18
    Be4O4           9.2887e-07 7.4489e-07 1.3719e-06  3.683e-08 1.6849e-10 1.5897e-18
    BeH             1.3243e-06 4.8345e-07 6.5487e-07 1.1009e-08 4.1984e-11 6.6186e-19
    BeH2            7.0102e-06 2.8949e-06 2.8856e-06 8.8386e-08 9.5162e-10 6.8426e-16
    BeN             4.3183e-08 1.3137e-08 1.7812e-08 1.6559e-10 2.9079e-13 4.1644e-22
    BeO              4.359e-07 1.8631e-07 3.4313e-07 8.1179e-09 3.7299e-11  6.566e-19
    BeOH            0.00022402 0.00012482 0.00016916 1.1082e-05  2.335e-07 7.1452e-13
    H                 0.007146  0.0055921  0.0075786  0.0028413 0.00062797 3.6751e-06
    H2                 0.51468    0.51518    0.51377    0.51519    0.51626    0.51664
    H2O               0.053371   0.054925   0.054802   0.059185   0.060446   0.060664
    H2O2            4.3791e-09 2.1323e-09 2.1285e-09  2.004e-10 7.8656e-12 3.4764e-16
    HNO             8.6849e-08 4.1263e-08 4.1208e-08 3.7031e-09 1.5209e-10 8.6972e-15
    HNO2            1.9242e-10 8.0962e-11 8.0894e-11 4.7618e-12 1.0221e-13 7.4123e-19
    HO2             2.0568e-09 9.3034e-10 1.2621e-09 7.7342e-11 1.3662e-12 2.9585e-18
    N               4.1872e-07 1.9115e-07  2.593e-07  1.691e-08 3.7932e-10 1.8929e-15
    N2                 0.22449    0.22435    0.22415    0.22374    0.22373    0.22376
    N2H2            7.2231e-09 2.7375e-09  1.481e-09 9.5737e-11 4.5049e-12 1.3131e-15
    N2O             7.9158e-09 3.8882e-09 3.8866e-09  3.857e-10 1.7835e-11 1.4552e-15
    N3H              1.411e-10 4.5333e-11 2.4548e-11 9.1327e-13 2.0287e-14 5.2666e-19
    NH              2.3329e-06 1.0891e-06 1.0871e-06 9.1901e-08 3.6387e-09 1.7282e-13
    NH2             1.8652e-05 9.8459e-06 7.2321e-06  1.105e-06 1.2065e-07 2.2644e-10
    NH2OH           2.6057e-09 1.0468e-09 5.6608e-10   4.32e-11 2.3068e-12 9.2313e-16
    NH3              0.0003067 0.00021088 0.00011398 4.9816e-05 2.7364e-05 1.3079e-05
    NO              1.6315e-05 1.0305e-05 1.3985e-05 2.5602e-06  2.032e-07 4.6237e-11
    O               2.0862e-06 1.1814e-06 2.1769e-06 2.3209e-07 7.2162e-09  5.494e-14
    O2              8.5536e-08 5.0299e-08 9.2727e-08  1.098e-08 3.6115e-10 3.0766e-15
    OH               0.0002671 0.00019071 0.00025858 7.0365e-05 9.3363e-06 1.0156e-08
    BeO(L)             0.19218    0.17345   0.031752          0          0          0
    BeO(a)                   0          0          0          0    0.19869    0.19892
    BeO(b)                   0   0.020211    0.16164     0.1975          0          0

    TRACE SPECIES:
    Be2             HNO3            N2H4            N2O3            N2O4            N2O5            N3              NH2NO2
    NO2             NO3             O3              Be(L)           Be(OH)2(b)      Be(a)           Be(b)           Be3N2(L)
    Be3N2(a)        Be3N2(b)        H2O(L)          H2O(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
