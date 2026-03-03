Example 5 from RP-1311
======================
.. note:: The python script for this example is available in the `source/bind/python/cea/samples` directory of the CEA repository.

Here we describe how to run example 5 from RP-1311 [1]_ using the Python API.
This is an HP equilibrium problem for a solid propellant blend with five reactants,
including one custom reactant (``CHOS-Binder``) that is not present in ``thermo.lib``.

First import the required libraries:

.. code-block:: python

    import numpy as np
    import cea

Define the pressure schedule and reactant composition by weight fraction:

.. code-block:: python

    pressures = np.array([34.473652, 17.236826, 8.618413, 3.447365, 0.344737])
    weights = np.array([0.7206, 0.1858, 0.09, 0.002, 0.0016], dtype=np.float64)
    T_reac = np.array([298.15, 298.15, 298.15, 298.15, 298.15], dtype=np.float64)

Create a :class:`~cea.Reactant` for ``CHOS-Binder`` using SI values.
In the Python API, custom reactant ``enthalpy`` is in J/kg and ``temperature`` is in K.
Use :mod:`cea.units` helpers for pre-conversion as needed.

.. code-block:: python

    chos_binder_mw_kg_per_mol = 14.6652984484e-3
    chos_binder_h_si = cea.units.cal_to_joule(-2999.082) / chos_binder_mw_kg_per_mol

    reactants = [
        "NH4CLO4(I)",
        cea.Reactant(
            name="CHOS-Binder",
            formula={"C": 1.0, "H": 1.86955, "O": 0.031256, "S": 0.008415},
            molecular_weight=14.6652984484,
            enthalpy=chos_binder_h_si,
            temperature=298.15,
        ),
        "AL(cr)",
        "MgO(cr)",
        "H2O(L)",
    ]

Then build reactant/product mixtures, solve the HP problem, and print the full table output:

.. literalinclude:: ../../../../source/bind/python/cea/samples/rp1311/example5.py
   :language: python

.. rubric:: References

.. [1] B. J. McBride and S. Gordon, *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications*, NASA RP-1311, 1996.
