Example 4 from RP-1311
======================
.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example4.py` in the CEA repository.

Here we describe how to run example 4 from RP-1311 [1]_ using the Python API.
This is a UV equilibrium problem, which is otherwise identical to :doc:`Example 3 <example3>`.
For that reason, the comments in this example will be limited to the minor differences.

.. code-block:: python

    import numpy as np
    import cea

    # Options
    trace = 1e-15

    # Species
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

Here, the fixed value of :math:`u/R` and :math:`{\rho}` (`density`) are pulled from the output of :doc:`Example 3 <example3>`.

.. code-block:: python
    :emphasize-lines: 2, 3

    # Thermo states
    u_R = -45.1343
    density = 14.428  # kg/m^3

    # Mixture states
    fuel_weights = np.array([0.0, 0.4, 0.6])
    oxidant_weights = np.array([1.0, 0.0, 0.0])
    T_reac = np.array([700.0, 298.15, 298.15])
    of_ratios = np.array([17.0])

    # Mixtures
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True, omit=omit_names)

    # Solver
    solver = cea.EqSolver(prod, reactants=reac, trace=trace)
    solution = cea.EqSolution(solver)

    # Initialize variable arrays
    n = len(of_ratios)
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

Note that the prescribed value of :math:`{\rho}` (`density`) is converted to specific volume before calling :meth:`~cea.EqSolver.solve`.

.. code-block:: python

    # Equilibrium solve
    for of_ratio in of_ratios:
        weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of_ratio)
        solver.solve(solution, cea.UV, u_R, 1.0/density, weights)

        # Store the output
        of_ratio_out[i] = of_ratio
        if solution.converged:
            P_out[i] = solution.P
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

    print("P, bar         ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(P_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(P_out[i]))

    print("T, K           ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(T_out[i]))

    print("Density, kg/m^3", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3e}".format(rho[i]), end=" ")
        else:
            print("{0:10.3e}".format(rho[i]))

    print("H, kJ/kg       ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(enthalpy[i]), end=" ")
        else:
            print("{0:10.3f}".format(enthalpy[i]))

    print("U, kJ/kg       ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(energy[i]), end=" ")
        else:
            print("{0:10.3f}".format(energy[i]))

    print("G, kJ/kg       ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.1f}".format(gibbs[i]), end=" ")
        else:
            print("{0:10.1f}".format(gibbs[i]))

    print("S, kJ/kg-K     ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(entropy[i]), end=" ")
        else:
            print("{0:10.3f}".format(entropy[i]))

    print("M, (1/n)       ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(molecular_weight_M[i]), end=" ")
        else:
            print("{0:10.3f}".format(molecular_weight_M[i]))

    print("MW             ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(molecular_weight_MW[i]), end=" ")
        else:
            print("{0:10.3f}".format(molecular_weight_MW[i]))

    print("Gamma_s        ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(gamma_s[i]), end=" ")
        else:
            print("{0:10.4f}".format(gamma_s[i]))

    print("Cp_eq, kJ/kg-K ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cp_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_eq[i]))

    print("Cp_fr, kJ/kg-K ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cp_fr[i]), end=" ")
        else:
            print("{0:10.4f}".format(cp_fr[i]))

    print("Cv_eq, kJ/kg-K ", end=" ")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(cv_eq[i]), end=" ")
        else:
            print("{0:10.4f}".format(cv_eq[i]))

    print("Cv_fr, kJ/kg-K ", end=" ")
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

Resulting in the output to the terminal:

.. code-block:: console

    P, bar              99.974
    T, K              2418.538
    Density, kg/m^3  1.443e+01
    H, kJ/kg           317.645
    U, kJ/kg          -375.270
    G, kJ/kg          -19437.2
    S, kJ/kg-K           8.168
    M, (1/n)            29.021
    MW                  29.021
    Gamma_s             1.2257
    Cp_eq, kJ/kg-K      1.6093
    Cp_fr, kJ/kg-K      1.4176
    Cv_eq, kJ/kg-K      1.3121
    Cv_fr, kJ/kg-K      1.1311

    MOLE FRACTIONS

    Ar               0.0088617
    CN                5.84e-14
    CO               0.0016763
    CO2                0.11537
    COOH            5.1528e-08
    H               2.7535e-05
    H2              0.00025052
    H2O                0.10277
    H2O2            9.8756e-07
    HCHO,formaldehy 1.7093e-11
    HCN             1.0535e-11
    HCO             7.3708e-10
    HCOOH           6.4546e-09
    HNC             1.0139e-12
    HNCO            1.2303e-09
    HNO             4.2176e-07
    HNO2            1.8492e-06
    HNO3            1.1294e-09
    HO2             7.9882e-06
    N               1.1482e-08
    N2                 0.73547
    N2H2            5.1653e-13
    N2O             3.6436e-06
    N2O3            2.4361e-11
    N2O4            3.0948e-15
    N3              1.2321e-12
    N3H              3.848e-13
    NCO             8.3006e-11
    NH              2.9308e-09
    NH2             1.7069e-09
    NH2NO2          2.4338e-15
    NH2OH           1.0202e-11
    NH3             3.8872e-09
    NO               0.0067747
    NO2             2.3462e-05
    NO3             1.9237e-10
    O               0.00015498
    O2                0.026252
    O3               1.209e-08
    OH               0.0023475

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
