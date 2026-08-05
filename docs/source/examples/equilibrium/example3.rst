Example 3 from RP-1311
======================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example3.py` in the CEA repository.

Here we describe how to run example 3 from RP-1311 [1]_ using the Python API.
This is an HP equilibrium, with Air as the oxidant, and a mixture of C\ :sub:`7`\ H\ :sub:`8`\ (L), C\ :sub:`8`\ H\ :sub:`18`\ (L),n-octane as the fuel. This example also sets the optional `trace` parameter to include species concentrations as small as 1\ :sup:`-15`\  .

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

We also set a value for the `trace` parameter which we will use later. This sets the threshold below which species concentrations are assumed to be exactly 0.

.. code-block:: python

    trace = 1e-15

Declare the reactant species names, and a list of species to *omit* from the product mixture, `omit_names`.

.. code-block:: python

    reac_names = ["Air", "C7H8(L)", "C8H18(L),n-octa"]
    omit_names = ["CCN", "CNC", "C2N2", "C2O",
                  "C3H4,allene", "C3H4,propyne", "C3H4,cyclo-", "C3",
                  "C3H5,allyl", "C3H6,propylene", "C3H6,cyclo-", "C3H3,propargyl",
                  "C3H6O", "C3H7,n-propyl", "C3H7,i-propyl", "Jet-A(g)",
                  "C3O2", "C4", "C4H2", "C3H8O,2propanol",
                  "C4H4,1,3-cyclo-", "C4H6,butadiene", "C4H6,2-butyne", "C3H8O,1propanol",
                  "C4H8,tr2-butene", "C4H8,isobutene", "C4H8,cyclo-", "C4H6,cyclo-",
                  "(CH3COOH)2", "C4H9,n-butyl", "C4H9,i-butyl", "C4H8,1-butene",
                  "C4H9,s-butyl", "C4H9,t-butyl", "C4H10,isobutane", "C4H8,cis2-buten",
                  "C4H10,n-butane", "C4N2", "C5", "C3H8",
                  "C5H6,1,3cyclo-", "C5H8,cyclo-", "C5H10,1-pentene", "C10H21,n-decyl",
                  "C5H10,cyclo-", "C5H11,pentyl", "C5H11,t-pentyl", "C12H10,biphenyl",
                  "C5H12,n-pentane", "C5H12,i-pentane", "CH3C(CH3)2CH3", "C12H9,o-bipheny",
                  "C6H6", "C6H5OH,phenol", "C6H10,cyclo-", "C6H2",
                  "C6H12,1-hexene", "C6H12,cyclo-", "C6H13,n-hexyl", "C6H5,phenyl",
                  "C7H7,benzyl", "C7H8", "C7H8O,cresol-mx", "C6H5O,phenoxy",
                  "C7H14,1-heptene", "C7H15,n-heptyl", "C7H16,n-heptane", "C10H8,azulene",
                  "C8H8,styrene", "C8H10,ethylbenz", "C8H16,1-octene", "C10H8,napthlene",
                  "C8H17,n-octyl", "C8H18,isooctane", "C8H18,n-octane", "C9H19,n-nonyl",
                  "Jet-A(L)", "C6H6(L)", "H2O(s)", "H2O(L)"]

Define the thermodynamic states at which we want to solve the equilibrium problem, in SI units. Because this is an HP problem, we will compute the fixed enthalpy value from the mixture later, so only pressure is defined here.

.. code-block:: python

    pressures = np.array([100.0, 10.0, 1.0])

Define the amounts of each reactant; in this case, an oxidant-to-fuel ratio is prescribed (`of_ratios`).
The arrays `fuel_weights` and `oxidant_weights` correspond to the `reac_names` list, and set the weight fraction of each that is part of the fuel and oxidant mixtures, respectively.
These arrays correspond to an oxidant mixture of 100\% Air, and a fuel mixture which is 40\% C\ :sub:`7`\ H\ :sub:`8`\ (L) by weight, and 60\% C\ :sub:`8`\ H\ :sub:`18`\ (L),n-octane by weight.
Additionally, we need to define the temperature of each reactant (`T_reac`), in order to compute the reference enthalpy of the reactant mixture later.

.. code-block:: python

    fuel_weights = np.array([0.0, 0.4, 0.6])
    oxidant_weights = np.array([1.0, 0.0, 0.0])
    T_reac = np.array([700.0, 298.15, 298.15])
    of_ratios = np.array([17.0])

Now having defined all of the relevant inputs to the problem, we can begin creating the required CEA objects, starting with the :class:`~cea.Mixture`.
Note that when we create the products mixture, we pass the `reac_names` while setting the flag `products_from_reactants=True`, along with the list of products to ignore, `omit=omit_names`.
This will return the subset of species that could possibly result from the provided reactants, while ignoring any in `omit_names`.

.. code-block:: python
    :emphasize-lines: 2

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True, omit=omit_names)

Next we instantiate the :class:`~cea.EqSolver` and :class:`~cea.EqSolution` objects. Note that the `trace` flag is passed at this point during the :class:`~cea.EqSolver` instantiation.

.. code-block:: python
    :emphasize-lines: 1

    solver = cea.EqSolver(prod, reactants=reac, trace=trace)
    solution = cea.EqSolution(solver)

We will now initialize an array to store each of the solution variables for printing the output later.

.. code-block:: python

    n = len(of_ratios)*len(pressures)
    of_ratio_out = np.zeros(n)
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
    cp_fr = np.zeros(n)
    cv_eq = np.zeros(n)
    cv_fr = np.zeros(n)
    mole_fractions = {}
    trace_species = []
    i = 0

Finally, we can loop through the defined pressures and oxidant-to-fuel ratios to solve the equilibrium problem at each state. We will also retrieve the solution variables and store them in the arrays we just initialized, and convert some units before storing.
The key points to note here are:

1. The :meth:`~cea.EqSolver.solve` method requires a list of reactant weights, which we compute using the `of_ratio_to_weights` method of the :class:`cea.Mixture` class.
2. The :meth:`~cea.EqSolver.solve` method requires the fixed value of both thermodynamic states. Because this is an HP problem, we must first compute the mixture enthalpy, which depends on the `weights` array that we just computed.
3. When we pass the fixed enthalpy value to :meth:`~cea.EqSolver.solve`, we first normalize it by :data:`cea.R`.

.. code-block:: python
    :emphasize-lines: 2, 3, 5

    for of_ratio in of_ratios:
        weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
        h0 = reac.calc_property(cea.ENTHALPY, weights, T_reac)
        for p in pressures:
            solver.solve(solution, cea.HP, h0/cea.R, p, weights)

            # Store the output
            of_ratio_out[i] = of_ratio
            P_out[i] = p
            if solution.converged:
                T_out[i] = solution.T
                rho[i] = solution.density
                enthalpy[i] = solution.enthalpy
                energy[i] = solution.energy
                gibbs[i] = solution.gibbs_energy
                entropy[i] = solution.entropy
                molecular_weight_M[i] = solution.M
                molecular_weight_MW[i] = solution.MW
                gamma_s[i] = solution.gamma_s
                cp_eq[i] = solution.cp_eq
                cp_fr[i] = solution.cp_fr
                cv_eq[i] = solution.cv_eq
                cv_fr[i] = solution.cv_fr

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

    print("P, bar          ", end="")
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

    print("Density,kg/m^3  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3e}".format(rho[i]), end=" ")
        else:
            print("{0:10.3e}".format(rho[i]))

    print("H, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.3f}".format(enthalpy[i]))

    print("U, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(energy[i]), end=" ")
        else:
            print("{0:10.3f}".format(energy[i]))

    print("G, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, kJ/kg-K      ", end="")
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

    print("Cp_eq, kJ/kg-K  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_eq[i]))

    print("Cp_fr, kJ/kg-K  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cp_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_fr[i]))

    print("Cv_eq, kJ/kg-K  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cv_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cv_eq[i]))

    print("Cv_fr, kJ/kg-K  ", end="")
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
        if np.any(mole_fractions[prod] > trace):
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

    o/f                 17.000     17.000     17.000
    P, bar             100.000     10.000      1.000
    T, K              2418.660   2390.593   2338.840
    Density,kg/m^3   1.443e+01  1.457e+00  1.483e-01
    H, kJ/kg           317.838    317.838    317.838
    U, kJ/kg          -375.112   -368.523   -356.299
    G, kJ/kg          -19438.0   -20787.3   -21879.3
    S, kJ/kg-K           8.168      8.828      9.491
    M, (1/n)            29.021     28.959     28.846
    MW                  29.021     28.959     28.846
    Gamma_s             1.2257     1.2061     1.1800
    Cp_eq, kJ/kg-K      1.6094     1.8167     2.2066
    Cp_fr, kJ/kg-K      1.4176     1.4155     1.4117
    Cv_eq, kJ/kg-K      1.3122     1.5038     1.8640
    Cv_fr, kJ/kg-K      1.1311     1.1284     1.1234

    MOLE FRACTIONS

    Ar               0.0088617   0.008843  0.0088084
    CN              5.8514e-14  1.066e-13 1.1783e-13
    CO               0.0016773   0.004312   0.009186
    CO2                0.11537    0.11249    0.10716
    COOH             5.157e-08 2.3271e-08 8.6129e-09
    H               2.7554e-05 0.00012385 0.00045515
    H2              0.00025064 0.00066103  0.0014836
    H2O                0.10277    0.10136   0.099015
    H2O2            9.8795e-07 2.9479e-07 8.3673e-08
    HCHO,formaldehy 1.7116e-11 1.1607e-11 5.5535e-12
    HCN             1.0552e-11 1.1777e-11 8.6347e-12
    HCO             7.3805e-10 8.9173e-10 7.6001e-10
    HCOOH           6.4598e-09 1.6408e-09 3.4278e-10
    HNC             1.0156e-12 1.0934e-12 7.4827e-13
    HNCO            1.2316e-09 5.1426e-10 1.6422e-10
    HNO             4.2203e-07  2.079e-07 9.0989e-08
    HNO2            1.8498e-06  3.265e-07 5.7312e-08
    HNO3            1.1297e-09 6.6231e-11 4.0363e-12
    HO2             7.9912e-06 4.2569e-06 2.1612e-06
    N               1.1495e-08 2.7422e-08 5.0672e-08
    N2                 0.73547    0.73403    0.73136
    N2H2            5.1721e-13 1.1978e-13 2.0997e-14
    N2O             3.6449e-06 1.1134e-06 3.2962e-07
    N2O3            2.4376e-11 7.7057e-13 2.4337e-14
    N2O4            3.0966e-15 3.2818e-17 3.6635e-19
    N3              1.2336e-12 2.9997e-13 5.7498e-14
    N3H             3.8529e-13 5.2335e-14 5.5797e-15
    NCO             8.3116e-11 5.8222e-11 2.9535e-11
    NH              2.9342e-09  3.864e-09 3.8817e-09
    NH2             1.7088e-09 1.2781e-09 7.3713e-10
    NH2NO2          2.4362e-15 6.7282e-17  1.649e-18
    NH2OH           1.0212e-11  1.448e-12 1.6805e-13
    NH3             3.8906e-09 1.7182e-09 6.1263e-10
    NO               0.0067762  0.0065527  0.0061447
    NO2             2.3467e-05 7.5641e-06 2.4794e-06
    NO3             1.9246e-10 1.9492e-11 1.9954e-12
    O               0.00015505 0.00043112  0.0010664
    O2                0.026252   0.027365   0.029599
    O3              1.2096e-08 3.7201e-09 1.1154e-09
    OH               0.0023483   0.003816  0.0057214

    TRACE SPECIES:
    (HCOOH)2        C               C10H8,naphthale C2              C2H             C2H2,acetylene  C2H2,vinylidene C2H3,vinyl
    C2H4            C2H4O,ethylen-o C2H5            C2H5OH          C2H6            C3H3,1-propynl  C3H3,2-propynl  C3H6O,acetone
    C3H6O,propanal  C3H6O,propylox  C4H2,butadiyne  C4H6,1butyne    C4H6,2butyne    C6H14,n-hexane  C7H16,2-methylh CH
    CH2             CH2CO,ketene    CH2OH           CH3             CH3CHO,ethanal  CH3CN           CH3CO,acetyl    CH3COOH
    CH3N2CH3        CH3O            CH3O2CH3        CH3OCH3         CH3OH           CH3OOH          CH4             CNCOCN
    CNN             HCCN            HCCO            HO(CO)2OH       N2H4            N2O5            NCN             O(CH)2O
    OCCN            OHCH2COOH       C(gr)           H2O(cr)


.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
