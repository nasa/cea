Example 12 from RP-1311
=======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example12.py` in the CEA repository.

Here we describe how to run example 12 from RP-1311 [1]_ using the Python API.

This is a rocket problem assuming an infinite-area combustor (IAC), with the composition frozen at the throat,
using CH\ :sub:`6`\ N\ :sub:`2`\ (L)  and N\ :sub:`2`\ O\ :sub:`4`\ (L)  as reactants.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions and :data:`cea.R` inline when normalizing enthalpy.

Declare the reactants and set their amounts.
The amounts of each are specified through the `fuel_weights`, `oxidant_weights`, and `of_ratio` variables specified here.
Setting the `fuel_weights` array equal to `[1.0, 0.0]` means that CH\ :sub:`6`\ N\ :sub:`2`\ (L)  constitutes 100\% of the fuel,
and similarly, setting the `oxidant_weights` array equal to `[0.0, 1.0]` means that N\ :sub:`2`\ O\ :sub:`4`\ (L)  constitutes 100\% of the oxidant.
These values will be used in conjunction later with `of_ratio` to compute the overall weight fraction array of the reactant mixture.

.. code-block:: python

    reac_names = ["CH6N2(L)", "N2O4(L)"]
    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    of_ratio = 2.5

Next, set some states for the rocket analysis. We will pass these values into the :class:`~cea.RocketSolver` later.
Note that here we define a variable `n_frz` that we will use later to declare which rocket station should begin the frozen composition.

.. code-block:: python
    :emphasize-lines: 4

    pc = cea.units.psi_to_bar(1000)  # Chamber pressure (bar)
    pi_p = 68.0457   # Pressure ratio
    supar = [5.0, 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0]  # Supersonic area ratio
    n_frz = 2

Now we will instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(["CO", "CO2", "H", "HNO", "HNO2", "HO2",
                        "H2", "H2O", "H2O2", "N", "NO", "NO2",
                        "N2", "N2O", "O", "OH", "O2", "HCO", "NH",
                        "CH4", "NH2", "NH3", "H2O(L)", "C(gr)"])

Now initialize the :class:`~cea.RocketSolver` and :class:`~cea.RocketSolution` class instances.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)

Compute the array of weight fractions for the reactant mixture:

.. code-block:: python

    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)

Compute the chamber enthalpy value based on the reactant weights and temperatures. We will pass this in later when we call :meth:`~cea.RocketSolver.solve`.
Note that this value is normalized by `R` here.

.. code-block:: python

    hc = reac.calc_property(cea.ENTHALPY, weights, 298.15)/cea.R

Now we can solve the :meth:`~cea.RocketSolver.solve` function:

.. code-block:: python

    solver.solve(solution, weights, pc, pi_p, supar=supar, iac=True, hc=hc, n_frz=n_frz)

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

    print("Cp, kJ/kg-K    ", end=" ")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.3f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp_eq[i]))

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

    print("")
    print("PERFORMANCE PARAMETERS")
    print("")

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
            print("{0:10.3f}".format(Cf[i]), end=" ")
        else:
            print("{0:10.3f}".format(Cf[i]))

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
    print("MASS FRACTIONS")
    print("")
    trace_species = []
    for prod in solution.mass_fractions:
        if np.any(solution.mass_fractions[prod] > 5e-6):
            print("{0:15s}".format(prod), end=" ")
            for j in range(len(solution.mass_fractions[prod])):
                if j < len(solution.mass_fractions[prod])-1:
                    print("{0:10.5g}".format(solution.mass_fractions[prod][j]), end=" ")
                else:
                    print("{0:10.5g}".format(solution.mass_fractions[prod][j]))
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

    P, bar              68.948     39.814      1.013      2.068      0.790      0.228      0.090      0.052      0.035      0.020      0.014
    T, K              3382.319   3203.673   1627.679   1868.526   1549.626   1202.295    983.296    870.192    796.165    700.162    637.715
    Density, kg/m^3      5.846      3.606      0.181      0.321      0.148      0.055      0.026      0.017      0.013      0.008      0.006
    H, kJ/kg            199.89    -426.81   -3389.58   -2960.26   -3525.96   -4113.12   -4463.73   -4637.96   -4749.30   -4890.46   -4980.28
    U, kJ/kg           -979.48   -1530.94   -3950.56   -3604.24   -4060.04   -4527.49   -4802.62   -4937.86   -5023.70   -5131.77   -5200.07
    G, kJ/kg          -36833.8   -35504.5   -21211.4   -23419.1   -20493.1   -17277.3   -15230.0   -14165.9   -13466.7   -12556.7   -11962.8
    S, kJ/kg-K          10.949     10.949     10.949     10.949     10.949     10.949     10.949     10.949     10.949     10.949     10.949
    M, (1/n)            23.845     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125
    MW                  23.845     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125     24.125
    Cp, kJ/kg-K          5.130      4.982      1.757      1.807      1.738      1.638      1.562      1.519      1.490      1.451      1.426
    Gamma_s              1.138      1.135      1.244      1.236      1.247      1.266      1.283      1.294      1.301      1.312      1.319
    Son. vel., m/s     1158.38    1119.55     835.41     892.06     816.21     724.39     659.43     622.86     597.49     562.57     538.37
    Mach                 0.000      1.000      3.207      2.818      3.344      4.054      4.631      4.994      5.266      5.672      5.979

    PERFORMANCE PARAMETERS

    Ae/At                0.000      1.000      8.342      5.000     10.000     25.000     50.000     75.000    100.000    150.000    200.000
    C*, m/s            1707.92    1707.92    1707.92    1707.92    1707.92    1707.92    1707.92    1707.92    1707.92    1707.92    1707.92
    Cf                   0.000      0.656      1.569      1.472      1.598      1.720      1.788      1.821      1.842      1.868      1.885
    Isp, vac., m/s       0.000   2105.782   2888.724   2770.175   2925.428   3078.093   3165.016   3206.992   3233.386   3266.348   3287.032
    Isp, m/s             0.000   1119.549   2679.355   2514.019   2729.780   2937.008   3054.053   3110.577   3146.170   3190.720   3218.748

    MASS FRACTIONS

    CO                0.076919   0.067456   0.067456   0.067456   0.067456   0.067456   0.067456   0.067456   0.067456   0.067456   0.067456
    CO2                0.15207    0.16694    0.16694    0.16694    0.16694    0.16694    0.16694    0.16694    0.16694    0.16694    0.16694
    H               0.00043944 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645 0.00033645
    HNO             1.8511e-05 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06 9.9871e-06
    HNO2             8.301e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06 4.7208e-06
    HO2             0.00014312 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05 9.0662e-05
    H2               0.0031271  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704  0.0026704
    H2O                0.28474    0.29296    0.29296    0.29296    0.29296    0.29296    0.29296    0.29296    0.29296    0.29296    0.29296
    H2O2            2.2385e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05 1.3449e-05
    N               5.0012e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06 2.5164e-06
    NO                0.021824   0.017486   0.017486   0.017486   0.017486   0.017486   0.017486   0.017486   0.017486   0.017486   0.017486
    NO2             4.8048e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05 3.1177e-05
    N2                 0.38097    0.38301    0.38301    0.38301    0.38301    0.38301    0.38301    0.38301    0.38301    0.38301    0.38301
    N2O             9.9328e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06  6.023e-06
    O                0.0050373  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866  0.0037866
    OH                0.034904    0.02888    0.02888    0.02888    0.02888    0.02888    0.02888    0.02888    0.02888    0.02888    0.02888
    O2                0.039711   0.036317   0.036317   0.036317   0.036317   0.036317   0.036317   0.036317   0.036317   0.036317   0.036317

    TRACE SPECIES:
    HCO             NH              CH4             NH2             NH3             H2O(L)          C(gr)


.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
