Example 1 from RP-1311
======================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example1.py` in the CEA repository.

Here we describe how to run example 1 from RP-1311 [1]_ using the Python API. This example is a TP equilibrium problem, with H\ :sub:`2`\  and Air as reactants.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions as needed.

.. code-block:: python

    # Unit conversions are handled with cea.units helpers.

Declare the product and reactant species names. Also note that `prod_names` is optional in general, but in this case, we explicitly define the list of species that we want to be in the final mixture (note for experienced CEA-users: this is akin to the `only` parameter in the legacy interface).

.. code-block:: python

    reac_names = ["H2", "Air"]
    prod_names = ["Ar",   "C",   "CO",  "CO2", "H",
                  "H2",   "H2O", "HNO", "HO2", "HNO2",
                  "HNO3", "N",   "NH",  "NO",  "N2",
                  "N2O3", "O",   "O2",  "OH",  "O3"]

Define the thermodynamic states at which we want to solve the equilibrium problem, in SI units.

.. code-block:: python

    pressures = cea.units.atm_to_bar(np.array([1.0, 0.1, 0.01]))
    temperatures = np.array([3000.0, 2000.0])

Define the amounts of each reactant; in this case, a chemical equivalence ratio `chem_eq_ratios` is prescribed (:math:`r_{eq}` in the RP-1311 [2]_).
The arrays `fuel_moles` and `oxidant_moles` correspond to the `reac_names` list, and set the mole fraction of each that is part of the fuel and oxidant mixtures, respectively. In this case, because we are using one fuel and one oxidizer, these are simply `1.0` to indicate which reactant is the fuel and which is the oxidizer.

.. code-block:: python

    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 1.0])
    chem_eq_ratios = np.array([1.0, 1.5])

Now having defined all of the relevant inputs to the problem, we can begin creating the required CEA objects, starting with the :class:`~cea.Mixture`.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(prod_names)

Next we instantiate the :class:`~cea.EqSolver` and :class:`~cea.EqSolution` objects.

.. code-block:: python

    solver = cea.EqSolver(prod, reactants=reac)
    solution = cea.EqSolution(solver)

Next, we convert the `chem_eq_ratios` to oxidant-to-fuel ratios. This also requires first converting the `fuel_moles` and `oxidant_moles` to `fuel_weights` and `oxidant_weights`, respectively.

.. code-block:: python

    fuel_weights = reac.moles_to_weights(fuel_moles)
    oxidant_weights = reac.moles_to_weights(oxidant_moles)
    of_ratios = len(chem_eq_ratios)*[0.0]
    for i, eqrat in enumerate(chem_eq_ratios):
        of_ratios[i] = reac.chem_eq_ratio_to_of_ratio(oxidant_weights,
                                                    fuel_weights,
                                                    eqrat)

We will now initialize an array to store each of the solution variables for printing the output later.

.. code-block:: python

    n = len(chem_eq_ratios)*len(pressures)*len(temperatures)
    of_ratio_out = np.zeros(n)
    T_out = np.zeros(n)
    P_out = np.zeros(n)
    rho = np.zeros(n)
    volume = np.zeros(n)
    enthalpy = np.zeros(n)
    energy = np.zeros(n)
    gibbs = np.zeros(n)
    entropy = np.zeros(n)
    molecular_weight_M = np.zeros(n)
    molecular_weight_MW = np.zeros(n)
    gamma_s = np.zeros(n)
    cp_eq = np.zeros(n)
    cp_fr = np.zeros(n)
    cv_eq = np.zeros(n)
    cv_fr = np.zeros(n)
    mole_fractions = {}
    i = 0

Finally, we can loop through the defined pressures, temperatures, and oxidant-to-fuel ratios to solve the equilibrium problem at each state. We will also retrieve the solution variables and store them in the arrays we just initialized, and convert some units before storing.
The key points to note here are:

1. The :meth:`~cea.EqSolver.solve` requires a list of reactant weights, which we compute using the `of_ratio_to_weights` method of the :class:`cea.Mixture` class.
2. The syntax of the :meth:`~cea.EqSolver.solve` method is `solver.solve(solution, cea.TP, temperature, pressure, weights)`, where `weights` is the list of reactant weights computed from the oxidant-to-fuel ratio.

.. code-block:: python

    for of_ratio in of_ratios:
        for p in pressures:
            for t in temperatures:
                weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
                solver.solve(solution, cea.TP, t, p, weights)

                # Store the output
                of_ratio_out[i] = of_ratio
                T_out[i] = t
                P_out[i] = cea.units.bar_to_atm(p)
                if solution.converged:
                    rho[i] = solution.density*1.e-3
                    volume[i] = solution.volume*1.e3
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

                if i == 0:
                    for prod in solution.mole_fractions:
                        mole_fractions[prod] = np.array([solution.mole_fractions[prod]])
                else:
                    for prod in mole_fractions:
                        mole_fractions[prod] = np.append(mole_fractions[prod], solution.mole_fractions[prod])

                i += 1

Finally, print everything out in a formatted manner consistent with the legacy CEA output format.

.. code-block:: python

    print("o/f             ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(of_ratio_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(of_ratio_out[i]))

    print("P, atm          ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(P_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(P_out[i]))

    print("T, K            ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(T_out[i]))

    print("Density, g/cc   ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3e}".format(rho[i]), end=" ")
        else:
            print("{0:10.3e}".format(rho[i]))

    print("Volume, cc/g    ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3e}".format(volume[i]), end=" ")
        else:
            print("{0:10.3e}".format(volume[i]))

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
    max_cols = 10
    nrows = (len(trace_species) + max_cols - 1) // max_cols
    for i in range(nrows):
        print(" ".join("{0:10s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))

This results in the following output to the terminal:

.. code-block:: console

    o/f                 34.296     34.296     34.296     34.296     34.296     34.296     22.853     22.853     22.853     22.853     22.853     22.853
    P, atm               1.000      1.000      0.100      0.100      0.010      0.010      1.000      1.000      0.100      0.100      0.010      0.010
    T, K              3000.000   2000.000   3000.000   2000.000   3000.000   2000.000   3000.000   2000.000   3000.000   2000.000   3000.000   2000.000
    Density, g/cc    9.178e-05  1.499e-04  8.085e-06  1.496e-05  6.614e-07  1.487e-06  8.122e-05  1.297e-04  7.116e-06  1.296e-05  5.672e-07  1.293e-06
    Volume, cc/g     1.090e+04  6.671e+03  1.237e+05  6.687e+04  1.512e+06  6.723e+05  1.231e+04  7.707e+03  1.405e+05  7.714e+04  1.763e+06  7.735e+05
    H, cal/g           663.714   -203.569   1369.546   -191.833   2647.099   -164.346    718.654   -120.698   1550.122   -116.204   3208.235   -101.779
    U, cal/g           399.856   -365.131   1070.016   -353.764   2280.923   -327.156    420.475   -307.350   1209.797   -303.012   2781.243   -289.088
    G, cal/g           -7974.5    -5290.4    -8616.7    -5662.8    -9381.0    -6036.5    -8818.8    -5830.6    -9545.4    -6260.6   -10424.5    -6691.2
    S, cal/g-K           2.879      2.543      3.329      2.735      4.009      2.936      3.179      2.855      3.698      3.072      4.544      3.295
    M, (1/n)            22.594     24.600     19.903     24.544     16.281     24.412     19.994     21.293     17.518     21.275     13.962     21.219
    MW                  22.594     24.600     19.903     24.544     16.281     24.412     19.994     21.293     17.518     21.275     13.962     21.219
    Gamma_s             1.1312     1.2258     1.1206     1.2026     1.1318     1.1668     1.1337     1.2529     1.1196     1.2389     1.1296     1.2051
    Cp_eq, cal/g-K      1.6816     0.4552     3.4398     0.5215     3.7168     0.6858     1.8253     0.4671     4.1929     0.5000     4.9230     0.6077
    Cp_fr, cal/g-K      0.4250     0.4008     0.4282     0.4008     0.4368     0.4010     0.4812     0.4520     0.4853     0.4521     0.4968     0.4522
    Cv_eq, cal/g-K      1.4371     0.3711     2.8458     0.4330     3.0554     0.5856     1.5585     0.3727     3.4492     0.4034     4.0091     0.5032
    Cv_fr, cal/g-K      0.3370     0.3200     0.3284     0.3199     0.3148     0.3196     0.3818     0.3587     0.3719     0.3587     0.3544     0.3586

    MOLE FRACTIONS

    Ar               0.0070982  0.0077283  0.0062528  0.0077107  0.0051148  0.0076691  0.0061933  0.0065959  0.0054263  0.0065904  0.0043249  0.0065728
    CO              0.00017073 1.0365e-05 0.00018417 2.0965e-05 0.00016776  4.081e-05 0.00017519 0.00015634 0.00016627  0.0001561  0.0001428 0.00015539
    CO2             7.1077e-05 0.00025289 2.8816e-05 0.00024169 6.4686e-06 0.00022042 3.5772e-05 6.8336e-05 1.8562e-05 6.8384e-05 4.5151e-06 6.8501e-05
    H                 0.040752 8.9523e-05    0.14289 0.00040886    0.31894  0.0018585    0.06025 0.00062121    0.18236  0.0019603    0.39312  0.0061596
    H2                0.067277  0.0030607   0.082718  0.0063841   0.041209   0.013191    0.14706    0.14737    0.13472    0.14676   0.062606     0.1449
    H2O                 0.2073    0.34207   0.095791    0.33714   0.011761    0.32637    0.22224    0.29508    0.11131     0.2945   0.014651    0.29261
    HO2             1.0257e-05 1.0068e-07 5.0804e-06 1.0266e-07  6.887e-07 1.0243e-07 3.6481e-06 2.2423e-10 3.3005e-06 7.1074e-10 5.7072e-07 2.2616e-09
    N               1.0577e-05 7.1874e-10 3.1348e-05   2.27e-09 8.9735e-05 7.1582e-09 9.9041e-06 6.6412e-10 2.9248e-05 2.0992e-09 8.2567e-05 6.6293e-09
    NO                0.012303 0.00048345   0.013705  0.0007215  0.0096685  0.0010659  0.0056508  8.003e-06  0.0091234 2.5353e-05  0.0072944 8.0572e-05
    N2                 0.58568    0.64414    0.51448    0.64255    0.42158    0.63891    0.51356    0.54995    0.44786    0.54948    0.35692    0.54799
    O                 0.015397 2.1253e-05   0.057868 0.00010042    0.14261  0.0004705  0.0075519 3.8075e-07   0.041288  3.816e-06    0.11694 3.8402e-05
    O2                0.018761  0.0010234   0.026501  0.0022849   0.016095  0.0050156  0.0045133 3.2847e-07   0.013491 3.2993e-06   0.010821 3.3412e-05
    OH                0.045174  0.0011285   0.059534  0.0024353   0.032747  0.0051864   0.032758 0.00014029   0.054207 0.00044369   0.033096   0.001403

    TRACE SPECIES:
    C          HNO        HNO2       HNO3       NH         N2O3       O3

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
.. [2] Gordon, S., McBride, B.J., "Computer program for calculation of complex chemical equilibrium compositions and applications. Part 1: Analysis",
    NASA RP-1311, 1994. [NTRS](https://ntrs.nasa.gov/citations/19950013764)
