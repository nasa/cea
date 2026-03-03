Python Interface
****************

The Python API provides a convenient way to interact with the program through Python, while maintaining the performance of the underlying Fortran code.
No unit conversions are performed by default; inputs and outputs are in the documented CEA units, and users are responsible for converting as needed.
Use the conversion factors in :mod:`cea.units` when working across unit systems.

Mixture
-------
The :class:`~cea.Mixture` class is used to define a mixture of product or reactant species. It allows the user to specify the composition of the mixture and provides methods to compute thermodynamic curve fit properties.
The instances of this class are then passed as inputs to the available solver classes (e.g., :class:`~cea.EqSolver`, :class:`~cea.RocketSolver`, :class:`~cea.ShockSolver`, or :class:`~cea.DetonationSolver`).
Custom reactants can be provided through :class:`~cea.Reactant` objects (including mixed lists of strings and Reactant objects).
For :class:`~cea.Reactant`, ``enthalpy`` and ``temperature`` are SI-only in the Python API (J/kg and K, respectively); use :mod:`cea.units` helpers for pre-conversion.

.. autoclass:: cea.Reactant
   :members:

.. autoclass:: cea.Mixture
   :members:

EqSolver
--------

.. autoclass:: cea.EqSolver
   :members:

EqSolution
----------

.. autoclass:: cea.EqSolution
   :members:

RocketSolver
------------

.. autoclass:: cea.RocketSolver
   :members:

RocketSolution
--------------

.. autoclass:: cea.RocketSolution
   :members:

ShockSolver
-----------

.. autoclass:: cea.ShockSolver
   :members:

ShockSolution
-------------

.. autoclass:: cea.ShockSolution
   :members:

DetonationSolver
----------------

.. autoclass:: cea.DetonationSolver
   :members:

DetonationSolution
------------------

.. autoclass:: cea.DetonationSolution
   :members:
