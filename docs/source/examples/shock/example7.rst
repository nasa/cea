Example 7 from RP-1311
======================

.. note:: The Python script for this example is available at `source/bind/python/cea/samples/rp1311/example7.py` in the CEA repository.

Here we describe how to run example 7 from RP-1311 [1]_ using the Python API.
This is a shock problem example using a mixture of H\ :sub:`2`\ , O\ :sub:`2`\ , and Ar.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Use :mod:`cea.units` for unit conversions.

Declare the reactant names, and the moles of each:

.. code-block:: python

    reac_names = ["H2", "O2", "Ar"]
    moles = np.array([0.05, 0.05, 0.9])

Define the unshocked temperature and pressure in SI units, as well as the inital shock velocity in m/s.

.. code-block:: python

    p0 = cea.units.mmhg_to_bar(10.0)  # Unshocked pressure (bar)
    T0 = 300.0   # Unshocked temperature (K)
    u1 = np.array([1100, 1200, 1250, 1300, 1350, 1400])  # Initial velocity (m/s)

Now instantiate the reactant and product :class:`~cea.Mixture` objects.
To create the product :class:`~cea.Mixture`, we pass the list of reactant names along with the flag `products_from_reactants=True`, which will return the full set of possible product species.

.. code-block:: python

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

Initialize the :class:`~cea.ShockSolver` and :class:`~cea.ShockSolution` objects.
We want to compute properties across the incident and reflected shock waves, so we set `reflected=True` when we initilaize the :class:`~cea.ShockSolution` object.
This needs to be done at this point to properly size the variable arrays, but note that this option will need to be set again later when we call :meth:`~cea.ShockSolver.solve`.

.. code-block:: python
    :emphasize-lines: 2

    solver = cea.ShockSolver(prod, reactants=reac)
    solution = cea.ShockSolution(solver, reflected=True)

Next we convert the reactant moles to weights:

.. code-block:: python

    weights = reac.moles_to_weights(moles)

Initialize arrays for storing the solution variables:

.. code-block:: python

    n = len(u1)
    mach1 = np.zeros(n)
    u1_out = np.zeros(n)
    P1_out = np.zeros(n)
    T1_out = np.zeros(n)
    rho1 = np.zeros(n)
    H1 = np.zeros(n)
    U1 = np.zeros(n)
    G1 = np.zeros(n)
    S1 = np.zeros(n)
    M1 = np.zeros(n)
    cp1 = np.zeros(n)
    gamma1 = np.zeros(n)
    v_sonic1 = np.zeros(n)

    u2 = np.zeros(n)
    P2 = np.zeros(n)
    T2 = np.zeros(n)
    rho2 = np.zeros(n)
    H2 = np.zeros(n)
    U2 = np.zeros(n)
    G2 = np.zeros(n)
    S2 = np.zeros(n)
    M2 = np.zeros(n)
    cp2 = np.zeros(n)
    gamma2 = np.zeros(n)
    v_sonic2 = np.zeros(n)
    P2_P1 = np.zeros(n)
    T2_T1 = np.zeros(n)
    M2_M1 = np.zeros(n)
    rho2_rho1 = np.zeros(n)
    v2 = np.zeros(n)

    u5 = np.zeros(n)
    P5 = np.zeros(n)
    T5 = np.zeros(n)
    rho5 = np.zeros(n)
    H5 = np.zeros(n)
    U5 = np.zeros(n)
    G5 = np.zeros(n)
    S5 = np.zeros(n)
    M5 = np.zeros(n)
    cp5 = np.zeros(n)
    gamma5 = np.zeros(n)
    v_sonic5 = np.zeros(n)
    P5_P2 = np.zeros(n)
    T5_T2 = np.zeros(n)
    M5_M2 = np.zeros(n)
    rho5_rho2 = np.zeros(n)
    u5_p_v2 = np.zeros(n)
    mole_fractions_incd = {}
    mole_fractions_refl = {}

Loop over the initial shock velocities and solve the shock problem, noting again that we set `reflected=True` to compute the equilibrium reflected shock solution.
This example still gates downstream use on ``solution.converged``. If a solve
returns ``solution.last_error == cea.LAST_VALID_SOLUTION``, the retained state
is still not a converged shock point.

.. code-block:: python
    :emphasize-lines: 4

    i = 0
    for u in u1:
        # Solve the shock problem
        solver.solve(solution, weights, T0, p0, u1=u, reflected=True)

        # Store the output
        u1_out[i] = u
        P1_out[i] = p0
        T1_out[i] = T0
        if solution.converged:
            mach1[i] = solution.Mach[0]
            rho1[i] = solution.density[0]
            H1[i] = solution.enthalpy[0]
            U1[i] = solution.energy[0]
            G1[i] = solution.gibbs_energy[0]
            S1[i] = solution.entropy[0]
            M1[i] = solution.M[0]
            cp1[i] = solution.cp_eq[0]
            gamma1[i] = solution.gamma_s[0]
            v_sonic1[i] = solution.sonic_velocity[0]

            u2[i] = solution.velocity[1]
            P2[i] = solution.P[1]
            T2[i] = solution.T[1]
            rho2[i] = solution.density[1]
            H2[i] = solution.enthalpy[1]
            U2[i] = solution.energy[1]
            G2[i] = solution.gibbs_energy[1]
            S2[i] = solution.entropy[1]
            M2[i] = solution.M[1]
            cp2[i] = solution.cp_eq[1]
            gamma2[i] = solution.gamma_s[1]
            v_sonic2[i] = solution.sonic_velocity[1]
            P2_P1[i] = solution.P21
            T2_T1[i] = solution.T21
            M2_M1[i] = solution.M21
            rho2_rho1[i] = 1.0/solution.rho12
            v2[i] = solution.v2

            u5[i] = solution.velocity[2]
            P5[i] = solution.P[2]
            T5[i] = solution.T[2]
            rho5[i] = solution.density[2]
            H5[i] = solution.enthalpy[2]
            U5[i] = solution.energy[2]
            G5[i] = solution.gibbs_energy[2]
            S5[i] = solution.entropy[2]
            M5[i] = solution.M[2]
            cp5[i] = solution.cp_eq[2]
            gamma5[i] = solution.gamma_s[2]
            v_sonic5[i] = solution.sonic_velocity[2]
            P5_P2[i] = solution.P52
            T5_T2[i] = solution.T52
            M5_M2[i] = solution.M52
            rho5_rho2[i] = solution.rho52
            u5_p_v2[i] = solution.u5_p_v2

            if i == 0:
                for prod in solution.mole_fractions:
                        mole_fractions_incd[prod] = np.array([solution.mole_fractions[prod][1]])
                        mole_fractions_refl[prod] = np.array([solution.mole_fractions[prod][2]])
            else:
                for prod in mole_fractions_incd:
                    mole_fractions_incd[prod] = np.append(mole_fractions_incd[prod], solution.mole_fractions[prod][1])
                for prod in mole_fractions_refl:
                    mole_fractions_refl[prod] = np.append(mole_fractions_refl[prod], solution.mole_fractions[prod][2])

        i += 1

Print out the results:

.. code-block:: python

    print()
    print("INITIAL GAS (1)")
    print()

    print("Mach Number1    ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(mach1[i]), end=" ")
        else:
            print("{0:10.3f}".format(mach1[i]))

    print("u1, m/s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(u1_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(u1_out[i]))

    print("P, bar          ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(P1_out[i]), end=" ")
        else:
            print("{0:10.4f}".format(P1_out[i]))

    print("T, K            ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T1_out[i]), end=" ")
        else:
            print("{0:10.3f}".format(T1_out[i]))

    print("Density, kg/m^3 ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4e}".format(rho1[i]), end=" ")
        else:
            print("{0:10.4e}".format(rho1[i]))

    print("H, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(H1[i]), end=" ")
        else:
            print("{0:10.3f}".format(H1[i]))

    print("U, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(U1[i]), end=" ")
        else:
            print("{0:10.3f}".format(U1[i]))

    print("G, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(G1[i]), end=" ")
        else:
            print("{0:10.3f}".format(G1[i]))

    print("S, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(S1[i]), end=" ")
        else:
            print("{0:10.3f}".format(S1[i]))

    print("M, (1/n)        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(M1[i]), end=" ")
        else:
            print("{0:10.3f}".format(M1[i]))

    print("Cp, kJ/kg-K     ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(cp1[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp1[i]))

    print("Gamma_s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(gamma1[i]), end=" ")
        else:
            print("{0:10.3f}".format(gamma1[i]))

    print("Son. Vel., m/s  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(v_sonic1[i]), end=" ")
        else:
            print("{0:10.3f}".format(v_sonic1[i]))

    print()
    print("SHOCKED GAS (2) - Incident, Equilibrium")
    print()
    print("u2, m/s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(u2[i]), end=" ")
        else:
            print("{0:10.3f}".format(u2[i]))

    print("P, bar          ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(P2[i]), end=" ")
        else:
            print("{0:10.4f}".format(P2[i]))

    print("T, K            ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T2[i]), end=" ")
        else:
            print("{0:10.3f}".format(T2[i]))

    print("Density, kg/m^3 ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4e}".format(rho2[i]), end=" ")
        else:
            print("{0:10.4e}".format(rho2[i]))

    print("H, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(H2[i]), end=" ")
        else:
            print("{0:10.3f}".format(H2[i]))

    print("U, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(U2[i]), end=" ")
        else:
            print("{0:10.3f}".format(U2[i]))

    print("G, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(G2[i]), end=" ")
        else:
            print("{0:10.3f}".format(G2[i]))

    print("S, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(S2[i]), end=" ")
        else:
            print("{0:10.3f}".format(S2[i]))

    print("M, (1/n)        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(M2[i]), end=" ")
        else:
            print("{0:10.3f}".format(M2[i]))

    print("Cp, kJ/kg-K     ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(cp2[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp2[i]))

    print("Gamma_s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(gamma2[i]), end=" ")
        else:
            print("{0:10.3f}".format(gamma2[i]))

    print("Son. Vel., m/s  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(v_sonic2[i]), end=" ")
        else:
            print("{0:10.3f}".format(v_sonic2[i]))

    print()
    print("P2/P1           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(P2_P1[i]), end=" ")
        else:
            print("{0:10.3f}".format(P2_P1[i]))

    print("T2/T1           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T2_T1[i]), end=" ")
        else:
            print("{0:10.3f}".format(T2_T1[i]))

    print("M2/M1           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(M2_M1[i]), end=" ")
        else:
            print("{0:10.3f}".format(M2_M1[i]))

    print("rho2/rho1       ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(rho2_rho1[i]), end=" ")
        else:
            print("{0:10.3f}".format(rho2_rho1[i]))

    print("v2, m/s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(v2[i]), end=" ")
        else:
            print("{0:10.3f}".format(v2[i]))

    print()
    print("MOLE FRACTIONS")
    print()
    trace_species = []
    for prod in mole_fractions_incd:
        if np.any(mole_fractions_incd[prod] > 5e-6):
            print("{0:15s}".format(prod), end=" ")
            for j in range(n):
                if j < n-1:
                    print("{0:10.5g}".format(mole_fractions_incd[prod][j]), end=" ")
                else:
                    print("{0:10.5g}".format(mole_fractions_incd[prod][j]))
        else:
            trace_species.append(prod)

    if len(trace_species) > 0:
        print()
        print("TRACE SPECIES:")
        max_cols = 8
        nrows = (len(trace_species) + max_cols - 1) // max_cols
        for i in range(nrows):
            print(" ".join("{0:15s}".format(trace_species[j]) for j in range(i * max_cols, min((i + 1) * max_cols, len(trace_species)))))

    print()
    print("SHOCKED GAS (5) - Reflected, Equilibrium")
    print()
    print("u5, m/s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(u5[i]), end=" ")
        else:
            print("{0:10.3f}".format(u5[i]))

    print("P, bar          ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4f}".format(P5[i]), end=" ")
        else:
            print("{0:10.4f}".format(P5[i]))

    print("T, K            ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T5[i]), end=" ")
        else:
            print("{0:10.3f}".format(T5[i]))

    print("Density, kg/m^3 ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.4e}".format(rho5[i]), end=" ")
        else:
            print("{0:10.4e}".format(rho5[i]))

    print("H, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(H5[i]), end=" ")
        else:
            print("{0:10.3f}".format(H5[i]))

    print("U, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(U5[i]), end=" ")
        else:
            print("{0:10.3f}".format(U5[i]))

    print("G, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(G5[i]), end=" ")
        else:
            print("{0:10.3f}".format(G5[i]))

    print("S, kJ/kg        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(S5[i]), end=" ")
        else:
            print("{0:10.3f}".format(S5[i]))

    print("M, (1/n)        ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(M5[i]), end=" ")
        else:
            print("{0:10.3f}".format(M5[i]))

    print("Cp, kJ/kg-K     ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(cp5[i]), end=" ")
        else:
            print("{0:10.3f}".format(cp5[i]))

    print("Gamma_s         ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(gamma5[i]), end=" ")
        else:
            print("{0:10.3f}".format(gamma5[i]))

    print("Son. Vel., m/s  ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(v_sonic5[i]), end=" ")
        else:
            print("{0:10.3f}".format(v_sonic5[i]))

    print()
    print("P5/P2           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(P5_P2[i]), end=" ")
        else:
            print("{0:10.3f}".format(P5_P2[i]))

    print("T5/T2           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(T5_T2[i]), end=" ")
        else:
            print("{0:10.3f}".format(T5_T2[i]))

    print("M5/M2           ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(M5_M2[i]), end=" ")
        else:
            print("{0:10.3f}".format(M5_M2[i]))

    print("rho5/rho2       ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(rho5_rho2[i]), end=" ")
        else:
            print("{0:10.3f}".format(rho5_rho2[i]))

    print("u5+v2, m/s      ", end="")
    for i in range(n):
        if i < n-1:
            print("{0:10.3f}".format(u5_p_v2[i]), end=" ")
        else:
            print("{0:10.3f}".format(u5_p_v2[i]))

    print()
    print("MOLE FRACTIONS")
    print()
    trace_species = []
    for prod in mole_fractions_refl:
        if np.any(mole_fractions_refl[prod] > 5e-6):
            print("{0:15s}".format(prod), end=" ")
            for j in range(n):
                if j < n-1:
                    print("{0:10.5g}".format(mole_fractions_refl[prod][j]), end=" ")
                else:
                    print("{0:10.5g}".format(mole_fractions_refl[prod][j]))
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

    INITIAL GAS (1)

    Mach Number1         3.353      3.658      3.810      3.962      4.115      4.267
    u1, m/s           1100.000   1200.000   1250.000   1300.000   1350.000   1400.000
    P, bar              0.0133     0.0133     0.0133     0.0133     0.0133     0.0133
    T, K               300.000    300.000    300.000    300.000    300.000    300.000
    Density, kg/m^3 2.0126e-02 2.0126e-02 2.0126e-02 2.0126e-02 2.0126e-02 2.0126e-02
    H, kJ/kg             1.062      1.062      1.062      1.062      1.062      1.062
    U, kJ/kg           -65.182    -65.182    -65.182    -65.182    -65.182    -65.182
    G, kJ/kg         -1556.265  -1556.265  -1556.265  -1556.265  -1556.265  -1556.265
    S, kJ/kg             5.191      5.191      5.191      5.191      5.191      5.191
    M, (1/n)            37.654     37.654     37.654     37.654     37.654     37.654
    Cp, kJ/kg-K          0.574      0.574      0.574      0.574      0.574      0.574
    Gamma_s              1.625      1.625      1.625      1.625      1.625      1.625
    Son. Vel., m/s     328.087    328.087    328.087    328.087    328.087    328.087

    SHOCKED GAS (2) - Incident, Equilibrium

    u2, m/s            666.380    575.193    558.801    546.820    536.956    527.905
    P, bar              0.1093     0.1642     0.1872     0.2104     0.2343     0.2591
    T, K              1528.516   1816.425   1930.785   2040.819   2147.535   2249.341
    Density, kg/m^3 3.3222e-02 4.1988e-02 4.5020e-02 4.7847e-02 5.0600e-02 5.3374e-02
    H, kJ/kg           384.032    555.644    626.184    696.556    768.446    841.973
    U, kJ/kg            54.945    164.500    210.327    256.838    305.453    356.577
    G, kJ/kg         -8121.784  -9579.928 -10165.744 -10731.423 -11281.059 -11805.155
    S, kJ/kg             5.565      5.580      5.589      5.600      5.611      5.623
    M, (1/n)            38.618     38.612     38.603     38.589     38.566     38.530
    Cp, kJ/kg-K          0.588      0.610      0.629      0.659      0.703      0.765
    Gamma_s              1.578      1.550      1.529      1.499      1.463      1.422
    Son. Vel., m/s     720.739    778.703    797.279    811.953    822.963    830.836

    P2/P1                8.200     12.318     14.043     15.781     17.572     19.432
    T2/T1                5.095      6.055      6.436      6.803      7.158      7.498
    M2/M1                1.026      1.025      1.025      1.025      1.024      1.023
    rho2/rho1            1.651      2.086      2.237      2.377      2.514      2.652
    v2, m/s            433.620    624.807    691.199    753.180    813.044    872.095

    MOLE FRACTIONS

    Ar                 0.92305    0.92289     0.9227    0.92236    0.92179    0.92093
    H               1.2167e-07  7.121e-06 2.5615e-05 7.6407e-05 0.00019702 0.00044351
    H2              2.6593e-06 4.9279e-05 0.00012308 0.00026816 0.00052464  0.0009274
    H2O               0.051235    0.05093   0.050603   0.050071   0.049255   0.048091
    O               2.8709e-06 5.5898e-05 0.00014199 0.00031455 0.00062675  0.0011323
    O2                0.025619    0.02549   0.025366   0.025185   0.024937   0.024628
    OH               8.642e-05 0.00057729  0.0010447  0.0017282  0.0026641  0.0038441

    TRACE SPECIES:
    H2O2            HO2             O3              H2O(L)          H2O(cr)

    SHOCKED GAS (5) - Reflected, Equilibrium

    u5, m/s            252.707    303.716    315.779    325.157    333.059    339.920
    P, bar              0.2590     0.4832     0.5832     0.6880     0.8009     0.9248
    T, K              2108.375   2576.766   2712.530   2829.429   2936.955   3038.970
    Density, kg/m^3 5.7006e-02 8.6377e-02 9.8543e-02 1.1083e-01 1.2352e-01 1.3693e-01
    H, kJ/kg           740.707   1120.094   1267.269   1411.248   1557.667   1707.473
    U, kJ/kg           286.301    560.643    675.433    790.470    909.312   1032.110
    G, kJ/kg        -11015.930 -13316.978 -13970.349 -14527.742 -15038.011 -15520.148
    S, kJ/kg             5.576      5.603      5.617      5.633      5.651      5.669
    M, (1/n)            38.578     38.296     38.107     37.896     37.664     37.413
    Cp, kJ/kg-K          0.680      1.078      1.277      1.462      1.628      1.768
    Gamma_s              1.481      1.310      1.274      1.253      1.239      1.232
    Son. Vel., m/s     820.448    855.925    868.334    881.790    896.418    912.129

    P5/P2                2.369      2.942      3.115      3.270      3.418      3.570
    T5/T2                1.379      1.419      1.405      1.386      1.368      1.351
    M5/M2                0.999      0.992      0.987      0.982      0.977      0.971
    rho5/rho2            1.716      2.057      2.189      2.316      2.441      2.566
    u5+v2, m/s         180.914    321.091    375.421    428.022    479.986    532.175

    MOLE FRACTIONS

    Ar                 0.92209    0.91534    0.91084     0.9058    0.90023    0.89425
    H               0.00012694  0.0029053  0.0055041  0.0088064   0.012802   0.017411
    H2              0.00038641  0.0033392  0.0049636  0.0064716   0.007792  0.0088491
    H2O               0.049672   0.041577   0.037066     0.0325   0.027947   0.023548
    O               0.00045802  0.0045891  0.0074745   0.010813   0.014638   0.018923
    O2                0.025057   0.023285    0.02251   0.021717   0.020846   0.019865
    OH               0.0022103  0.0089647   0.011639   0.013892   0.015742   0.017153

    TRACE SPECIES:
    H2O2            HO2             O3              H2O(L)          H2O(cr)


.. [1] McBride, B.J., Gordon, S., "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description: Users Manual and Program Description - 2",
    NASA RP-1311, 1996. [NTRS](https://ntrs.nasa.gov/citations/19960044559)
