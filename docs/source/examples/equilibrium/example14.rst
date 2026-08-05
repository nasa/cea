Example 14 from RP-1311
========================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example14.py` in the CEA repository.

Here we describe how to run example 14 from RP-1311 [1]_ using the Python API. This example is a TP equilibrium problem with H\ :sub:`2`\ (L)  and O\ :sub:`2`\ (L)  as reactants.
This problem results in significant amounts of condensed species in the resulting equilibrium mixture, and is used to illustrate the effects of condensed species on volume and molecular weight.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions as needed.

.. code-block:: python

    # Unit conversions are handled with cea.units helpers.

Declare the reactant species names:

.. code-block:: python

    reac_names = ["H2(L)", "O2(L)"]

Define the thermodynamic states at which we want to solve the equilibrium problem, in SI units.

.. code-block:: python

    p = cea.units.atm_to_bar(0.05)
    temperatures = np.array([1000.0, 500.0, 350.0, 305.0, 304.3, 304.2, 304.0, 300.0])

Define the amounts of each reactant.

.. code-block:: python

    fuel_moles = np.array([100.0, 0.0])
    oxidant_moles = np.array([0.0, 60.0])

Create the :class:`~cea.Mixture` objects. The product mixture is created using the reactants species list, and setting `products_from_reactants=True` to obtain the full set of possible product species.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

Next we instantiate the :class:`~cea.EqSolver` and :class:`~cea.EqSolution` objects.

.. code-block:: python

    solver = cea.EqSolver(prod, reactants=reac)
    solution = cea.EqSolution(solver)

Now compute the weight fractions of the reactant mixture:

.. code-block:: python

    fuel_weights = reac.moles_to_weights(fuel_moles)
    oxidant_weights = reac.moles_to_weights(oxidant_moles)
    weights = fuel_weights + oxidant_weights

We will now initialize an array to store each of the solution variables for printing the output later.

.. code-block:: python

    n = len(temperatures)
    T_out = np.zeros(n)
    P_out = np.zeros(n)
    rho = np.zeros(n)
    enthalpy = np.zeros(n)
    energy = np.zeros(n)
    gibbs = np.zeros(n)
    entropy = np.zeros(n)
    molecular_weight_M = np.zeros(n)
    molecular_weight_MW = np.zeros(n)
    gamma_s = np.zeros(n)
    cp_eq = np.zeros(n)
    mole_fractions = {}
    i = 0

Finally, loop over the prescribed temperature values and store the output:

.. code-block:: python

    for t in temperatures:
        solver.solve(solution, cea.TP, t, p, weights)

        # Store the output
        T_out[i] = t
        P_out[i] = p
        if solution.converged:
            rho[i] = solution.density
            enthalpy[i] = solution.enthalpy
            energy[i] = solution.energy
            gibbs[i] = solution.gibbs_energy
            entropy[i] = solution.entropy
            molecular_weight_M[i] = solution.M
            molecular_weight_MW[i] = solution.MW
            gamma_s[i] = solution.gamma_s
            cp_eq[i] = solution.cp_eq

            if i == 0:
                for prod in solution.mole_fractions:
                    mole_fractions[prod] = np.array([solution.mole_fractions[prod]])
            else:
                for prod in mole_fractions:
                    mole_fractions[prod] = np.append(mole_fractions[prod], solution.mole_fractions[prod])

            i += 1

And print the output to the terminal:

.. code-block:: python

    print("P, bar          ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.5f}".format(P_out[i]), end=" ")
        else:
            print("{0:10.5f}".format(P_out[i]))

    print("T, K            ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(T_out[i]))

    print("Density, kg/m^3 ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3e}".format(rho[i]), end=" ")
        else:
            print("{0:10.3e}".format(rho[i]))

    print("H, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.2f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.2f}".format(enthalpy[i]))

    print("U, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.2f}".format(energy[i]), end=" ")
        else:
            print("{0:10.2f}".format(energy[i]))

    print("G, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, kJ/kg-K      ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(entropy[i]), end=" ")
        else:
            print("{0:10.4f}".format(entropy[i]))

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

    print("Cp, kJ/kg-K     ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_eq[i]))

    print("Gamma_s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(gamma_s[i]), end=" ")
        else:
            print("{0:10.4f}".format(gamma_s[i]))

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

    P, bar             0.05066    0.05066    0.05066    0.05066    0.05066    0.05066    0.05066    0.05066
    T, K              1000.000    500.000    350.000    305.000    304.300    304.200    304.000    300.000
    Density, kg/m^3  1.175e-02  2.350e-02  3.358e-02  3.853e-02  4.499e-02  4.715e-02  5.146e-02  1.304e-01
    H, kJ/kg         -10066.00  -11043.65  -11309.10  -11386.94  -11709.17  -11798.29  -11953.26  -12988.71
    U, kJ/kg         -10497.11  -11259.20  -11459.98  -11518.42  -11821.79  -11905.74  -12051.70  -13027.57
    G, kJ/kg          -23601.6   -17139.8   -15355.8   -14840.7   -14833.0   -14832.0   -14830.0   -14801.2
    S, kJ/kg-K         13.5356    12.1924    11.5619    11.3239    10.2655     9.9726     9.4630     6.0416
    M, (1/n)            19.287     19.287     19.287     19.287     22.466     23.541     25.676     64.179
    MW                  19.287     19.287     19.287     19.287     19.287     19.287     19.287     19.287
    Cp, kJ/kg-K         2.1108     1.8069     1.7370     1.7233   936.2315   848.8474   707.1961    95.8701
    Gamma_s             1.2567     1.3133     1.3301     1.3336     1.1080     1.1059     1.1018     1.0344

    MOLE FRACTIONS

    H2O                0.90909    0.90909    0.90909    0.90909    0.76756    0.72836    0.66025     0.2096
    O2                0.090909   0.090909   0.090909   0.090909   0.090909   0.090909   0.090909   0.090909
    H2O(L)                   0          0          0          0    0.14153    0.18073    0.24884    0.69949

    TRACE SPECIES:
    H          H2         H2O2       HO2        O          O3         OH         H2O(cr)

.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
