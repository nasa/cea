Example 11 from RP-1311
=======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example11.py` in the CEA repository.

Here we describe how to run example 11 from RP-1311 [1]_ using the Python API.
This is a rocket problem assuming an infinite-area combustor (IAC), using Li(cr) and F\ :sub:`2`\ (L)  as reactants.
The problem includes the use of ionized species, and computes transport properties of the resulting mixture.
The result of this problem also includes condensed species.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions and :data:`cea.R` inline when normalizing enthalpy.

Declare the reactants and set their amounts and initial temperatures.
The initial reactant temperatures, `T_reactant`, will be used later to compute the chamber enthalpy.
The amounts of each are specified through the `fuel_moles` and `oxidant_moles` variables specified here.

.. code-block:: python

    reac_names = ["Li(cr)", "F2(L)"]
    T_reactant = np.array([298.15, 85.02])  # Reactant temperatures (K)
    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 0.5556])

Next, define the chamber pressure, pressure ratio, and subsonic and supersonic area ratios:

.. code-block:: python

    pc = cea.units.psi_to_bar(1000.0)  # Chamber pressure (bar)
    pi_p = [68.0457]   # Pressure ratio
    subar = [10.0]
    supar = [10.0, 20.0, 100.0]  # Supersonic area ratio

Now we will instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.
Note that both of these initializations set `ions=True` to ensure that ionized species are included.

.. code-block:: python
    :emphasize-lines: 1, 2

    reac = cea.Mixture(reac_names, ions=True)
    prod = cea.Mixture(reac_names, products_from_reactants=True, ions=True)

Now initialize the :class:`~cea.RocketSolver` and :class:`~cea.RocketSolution` class instances.
We again set `ions=True` during the :class:`~cea.RocketSolver` class initialization,
and at this point, we also set `transport=True` which sets the flag to compute transport properties later
when we call :meth:`~cea.RocketSolver.solve`.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.RocketSolver(prod, reactants=reac, transport=True, ions=True)
    solution = cea.RocketSolution(solver)

Compute the array of weight fractions for the reactant mixture:

.. code-block:: python

    fuel_weights = reac.moles_to_weights(fuel_moles)
    oxidant_weights = reac.moles_to_weights(oxidant_moles)
    weights = fuel_weights + oxidant_weights

Compute the chamber enthalpy value based on the reactant weights and temperatures. We will pass this in later when we call :meth:`~cea.RocketSolver.solve`.
Note that this value is normalized by `R` here.

.. code-block:: python

    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R

Now we can solve the :meth:`~cea.RocketSolver.solve` function:

.. code-block:: python

    solver.solve(solution, weights, pc, pi_p, subar=subar, supar=supar, iac=True, hc=hc)

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
    visc = solution.viscosity
    cond_fr = solution.conductivity_fr
    cond_eq = solution.conductivity_eq
    prandtl_fr = solution.Pr_fr
    prandtl_eq = solution.Pr_eq
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
    print("TRANSPORT PROPERTIES")
    print("")

    print("Viscosity, mP   ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(visc[i]), end=" ")
        else:
            print("{0:10.4f}".format(visc[i]))

    print("")
    print("with equilibrium reaction:")
    print("")

    print("Cp, cal/g-K     ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_eq[i]))

    print("Cv, cal/g-K     ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cv_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cv_eq[i]))

    print("Conductivity    ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cond_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cond_eq[i]))

    print("Prandtl number  ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(prandtl_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(prandtl_eq[i]))

    print("")
    print("with frozen reaction:")
    print("")

    print("Cp, cal/g-K     ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cp_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_fr[i]))

    print("Cv, cal/g-K     ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cv_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(cv_fr[i]))

    print("Conductivity    ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(cond_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(cond_fr[i]))

    print("Prandtl number  ", end="")
    for i in range(num_pts):
        if i < num_pts-1:
            print("{0:10.4f}".format(prandtl_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(prandtl_fr[i]))

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

    P, bar              68.948     39.299      1.013     68.804      0.942      0.368      0.044
    T, K              5680.808   5332.680   3508.666   5679.461   3470.831   2926.665   1955.978
    Density, kg/m^3      3.218      1.994      0.087      3.212      0.082      0.038      0.007
    H, kJ/kg           -259.28   -1414.52   -7031.63    -263.73   -7115.90   -8109.52   -9765.47
    U, kJ/kg          -2402.02   -3385.05   -8196.38   -2405.80   -8266.77   -9072.36  -10394.18
    G, kJ/kg          -64725.7   -61930.3   -46848.3   -64714.9   -46503.3   -41321.6   -31962.1
    S, kJ/kg-K          11.348     11.348     11.348     11.348     11.348     11.348     11.348
    M, (1/n)            22.043     22.501     25.047     22.045     25.075     25.273     25.867
    MW                  22.043     22.501     25.047     22.045     25.075     25.273     25.867
    Cp, kJ/kg-K          6.742      6.570      2.645      6.741      2.529      1.631      3.122
    Gamma_s              1.178      1.173      1.195      1.178      1.199      1.267      1.196
    Son. vel., m/s     1588.67    1520.02    1179.71    1588.41    1174.84    1104.36     866.99
    Mach                 0.000      1.000      3.120      0.059      3.152      3.588      5.029

    TRANSPORT PROPERTIES

    Viscosity, mP       1.4414     1.3880     1.0810     1.4412     1.0734     0.9580     0.7308

    with equilibrium reaction:

    Cp, cal/g-K         6.7418     6.5700     2.6448     6.7414     2.5293     1.6315     3.1219
    Cv, cal/g-K         5.2923     5.2263     2.1940     5.2922     2.0924     1.2854     2.5518
    Conductivity       14.4426    13.6926     4.3206    14.4401     4.0891     2.2075     2.6678
    Prandtl number      0.6729     0.6660     0.6617     0.6728     0.6640     0.7080     0.8552

    with frozen reaction:

    Cp, cal/g-K         1.6854     1.6503     1.5156     1.6852     1.5137     1.4904     1.4584
    Cv, cal/g-K         1.3082     1.2807     1.1836     1.3081     1.1822     1.1614     1.1370
    Conductivity        3.1552     3.0168     2.2918     3.1547     2.2744     2.0066     1.4888
    Prandtl number      0.7700     0.7593     0.7148     0.7699     0.7144     0.7115     0.7159

    PERFORMANCE PARAMETERS

    Ae/At                0.000      1.000      9.468     10.000     10.000     20.000    100.000
    C*, m/s            2274.40    2274.40    2274.40    2274.40    2274.40    2274.40    2274.40
    Cf                   0.000      0.668      1.618      0.041      1.628      1.742      1.917
    Isp, vac., m/s       0.000   2816.400   3996.792  22791.160   4013.922   4205.378   4504.509
    Isp, m/s             0.000   1520.029   3680.312     94.377   3703.139   3962.382   4360.319

    MOLE FRACTIONS

    F                  0.20853    0.19341     0.1076    0.20848    0.10664    0.10062    0.10253
    F-               0.0045774  0.0035905 0.00028192  0.0045735 0.00025498 3.6579e-05 3.7622e-08
    F2              1.5866e-05 1.0131e-05 6.6385e-07  1.584e-05   6.47e-07 6.8415e-07   2.69e-06
    Li                 0.11766     0.1016  0.0082205     0.1176  0.0071628 0.00043266 2.0646e-08
    Li+              0.0073265  0.0058435 0.00037576  0.0073207 0.00033424 3.8991e-05 3.7631e-08
    Li-             2.3779e-05 1.2161e-05 6.1699e-09 2.3723e-05 4.4374e-09 7.1984e-12 1.4196e-21
    Li2             0.00020698 0.00010255 6.6038e-08 0.00020646 4.8707e-08 1.4744e-10 3.9217e-19
    Li2+            0.00016743 8.3248e-05 6.1651e-08 0.00016701  4.669e-08  2.969e-10 2.1021e-17
    Li2F2            0.0014414  0.0012625 0.00083528  0.0014407 0.00085229   0.001579   0.023952
    Li3F3           3.7935e-06 2.6207e-06 6.4793e-07 3.7882e-06 6.6756e-07   1.99e-06 0.00034821
    LiF                0.65715    0.69175    0.88259    0.65727    0.88467    0.89729    0.87316
    e-               0.0028928  0.0023241 9.3897e-05  0.0028905 7.9304e-05  2.412e-06 9.0735e-12

    TRACE SPECIES:
    F+              Li3+            Li(L)           Li(cr)          LiF(L)          LiF(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
