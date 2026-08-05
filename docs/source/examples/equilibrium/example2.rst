Example 2 from RP-1311
=======================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example2.py` in the CEA repository.

Here we describe how to run example 2 from RP-1311 [1]_ using the Python API. This is a TV equilibrium problem, with H\ :sub:`2`\  and Air as reactants, and computes transport properties of the resulting equilibrium mixture.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions.

Next we define a flag to turn on transport properties; this could also be done in-line later in the code.

.. code-block:: python

    transport = True

Declare the product and reactant species names. Also note that `prod_names` is optional in general, but in this case, we explicitly define the list of species that we want to be in the final mixture (note for experienced CEA-users: this is akin to the `only` parameter in the legacy interface).

.. code-block:: python

    reac_names = ["H2", "Air"]
    prod_names = ["Ar",   "C",   "CO",  "CO2", "H",
                  "H2",   "H2O", "HNO", "HO2", "HNO2",
                  "HNO3", "N",   "NH",  "NO",  "N2",
                  "N2O3", "O",   "O2",  "OH",  "O3"]

Define the thermodynamic states at which we want to solve the equilibrium problem, in SI units.

.. code-block:: python

    densities = 1.0e3*np.array([9.1864e-5, 8.0877e-6, 6.6054e-7])  # kg/m^3
    temperatures = np.array([3000.0])

Define the amounts of each reactant; in this case, a weight equivalence ratio `phis` is prescribed (:math:`{\phi}` in the RP-1311 [2]_).
The arrays `fuel_moles` and `oxidant_moles` correspond to the `reac_names` list, and set the mole fraction of each that is part of the fuel and oxidant mixtures, respectively. In this case, because we are using one fuel and one oxidizer, these are simply `1.0` to indicate which reactant is the fuel and which is the oxidizer.

.. code-block:: python

    phis = np.array([1.0])
    fuel_moles = np.array([1.0, 0.0])
    oxidant_moles = np.array([0.0, 1.0])

Now having defined all of the relevant inputs to the problem, we can begin creating the required CEA objects, starting with the :class:`~cea.Mixture`.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(prod_names)

Next we instantiate the :class:`~cea.EqSolver` and :class:`~cea.EqSolution` objects. Note that the `transport` flag is passed at this point during the :class:`~cea.EqSolver` instantiation.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.EqSolver(prod, reactants=reac, transport=transport)
    solution = cea.EqSolution(solver)

Next, we convert the `phis` to oxidant-to-fuel ratios. This also requires first converting the `fuel_moles` and `oxidant_moles` to `fuel_weights` and `oxidant_weights`, respectively.

.. code-block:: python

    fuel_weights = reac.moles_to_weights(fuel_moles)
    oxidant_weights = reac.moles_to_weights(oxidant_moles)
    of_ratios = len(phis)*[0.0]
    for i, phi in enumerate(phis):
        of_ratios[i] = reac.weight_eq_ratio_to_of_ratio(oxidant_weights,
                                                        fuel_weights,
                                                        phi)

We will now initialize an array to store each of the solution variables for printing the output later.

.. code-block:: python

    n = len(phis)*len(densities)*len(temperatures)
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
    visc = np.zeros(n)
    cond_fr = np.zeros(n)
    cond_eq = np.zeros(n)
    prandtl_fr = np.zeros(n)
    prandtl_eq = np.zeros(n)
    mole_fractions = {}
    trace_species = []
    i = 0

Finally, we can loop through the defined pressures, temperatures, and oxidant-to-fuel ratios to solve the equilibrium problem at each state. We will also retrieve the solution variables and store them in the arrays we just initialized, and convert some units before storing.
The key points to note here are:

1. The :meth:`~cea.EqSolver.solve` requires a list of reactant weights, which we compute using the `of_ratio_to_weights` method of the :class:`cea.Mixture` class.
2. The syntax of the :meth:`~cea.EqSolver.solve` method is `solver.solve(solution, cea.TP, temperature, pressure, weights)`, where `weights` is the list of reactant weights computed from the oxidant-to-fuel ratio.

.. code-block:: python

    for of_ratio in of_ratios:
        for density in densities:
            for t in temperatures:
                weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
                solver.solve(solution, cea.TV, t, 1.0/density, weights)

                # Store the output
                of_ratio_out[i] = of_ratio
                T_out[i] = t
                if solution.converged:
                    rho[i] = solution.density*1.e-3
                    P_out[i] = cea.units.bar_to_atm(solution.P)
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

    o/f                 34.296     34.296     34.296
    P, atm               1.001      0.100      0.010
    T, K              3000.000   3000.000   3000.000
    Density, g/cc    9.186e-05  8.088e-06  6.605e-07
    Volume, cc/g     1.089e+04  1.236e+05  1.514e+06
    H, cal/g           663.554   1369.409   2647.694
    U, cal/g           399.704   1069.887   2281.487
    G, cal/g           -7974.3    -8616.6    -9381.4
    S, cal/g-K           2.879      3.329      4.010
    M, (1/n)            22.595     19.904     16.279
    MW                  22.595     19.904     16.279
    Gamma_s             1.1312     1.1206     1.1318

    Transport properties:

    Viscosity, mP       0.9358     0.9401     0.9482

    with equilibrium reaction:

    Cp, cal/g-K         1.6812     3.4395     3.7156
    Cv, cal/g-K         1.4368     2.8456     3.0544
    Conductivity        4.4411     9.6448     8.8614
    Prandtl number      0.3593     0.3377     0.4009

    with frozen reaction:

    Cp, cal/g-K         0.4250     0.4282     0.4368
    Cv, cal/g-K         0.3370     0.3284     0.3148
    Conductivity        0.6290     0.7265     0.8641
    Prandtl number      0.6322     0.5541     0.4793

    MOLE FRACTIONS

    Ar               0.0070984   0.006253  0.0051143
    CO              0.00017071 0.00018417 0.00016775
    CO2             7.1088e-05 2.8821e-05 6.4627e-06
    H                 0.040731    0.14287    0.31902
    H2                0.067266   0.082719   0.041181
    H2O                0.20733   0.095807   0.011743
    HO2             1.0258e-05 5.0811e-06 6.8767e-07
    N               1.0572e-05 3.1343e-05 8.9783e-05
    NO                0.012302   0.013705  0.0096653
    N2                 0.58569     0.5145    0.42155
    O                 0.015389    0.05786    0.14265
    O2                0.018757   0.026501   0.016086
    OH                0.045166   0.059534   0.032727

    TRACE SPECIES:
    C          HNO        HNO2       HNO3       NH         N2O3       O3

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
.. [2] Gordon, S., McBride, B.J., "Computer program for calculation of complex chemical equilibrium compositions and applications. Part 1: Analysis",
    NASA RP-1311, 1994. [NTRS](https://ntrs.nasa.gov/citations/19950013764)
