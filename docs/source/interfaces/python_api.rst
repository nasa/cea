Python Interface
****************

The Python API provides a convenient way to interact with the program through Python, while maintaining the performance of the underlying Fortran code.
No unit conversions are performed by default; inputs and outputs are in the documented CEA units, and users are responsible for converting as needed.
Use the conversion factors in :mod:`cea.units` when working across unit systems.

MATLAB-Oriented Wrapper Module
------------------------------

The :mod:`cea.matlab` module provides a MATLAB-friendly compatibility layer on
top of the Python package. It is intended for use through MATLAB's Python
bridge, not as a separate compiled MATLAB extension.

Import :mod:`cea` for enums/constants such as ``cea.TP`` and for
:mod:`cea.units`, then import :mod:`cea.matlab` for the solve entry points:

* :func:`cea.matlab.eq_solve` for equilibrium problems
* :func:`cea.matlab.rocket_solve` for rocket performance problems
* :func:`cea.matlab.shock_solve` for incident/reflected shock problems
* :func:`cea.matlab.detonation_solve` for detonation problems

These wrappers return flat namespace-style results composed of Python scalars,
NumPy arrays, and dictionaries so MATLAB can consume the outputs without
working directly with the compiled solution classes. Example scripts live in
``source/bind/matlab/samples/``.

.. autofunction:: cea.matlab.eq_solve

.. autofunction:: cea.matlab.rocket_solve

.. autofunction:: cea.matlab.shock_solve

.. autofunction:: cea.matlab.detonation_solve

Mixture
-------
The :class:`~cea.Mixture` class is used to define a mixture of product or reactant species. It allows the user to specify the composition of the mixture and provides methods to compute thermodynamic curve fit properties.
The instances of this class are then passed as inputs to the available solver classes (e.g., :class:`~cea.EqSolver`, :class:`~cea.RocketSolver`, :class:`~cea.ShockSolver`, or :class:`~cea.DetonationSolver`).
Custom reactants can be provided through :class:`~cea.Reactant` objects (including mixed lists of strings and Reactant objects).
For :class:`~cea.Reactant`, ``temperature`` is in K, and ``enthalpy`` must be paired with explicit
``enthalpy_units`` (for example ``"J/kg"`` or ``"kJ/mol"``). Weight-based enthalpy units require
``molecular_weight`` for conversion, while molar enthalpy units do not.

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

Shock solves can emit a ``RuntimeWarning`` with ``LAST_VALID_SOLUTION`` when the
incident-equilibrium path retains a last valid state without converging the
shock iteration. In that case ``ShockSolution.last_error`` records the soft
failure, while ``ShockSolution.converged`` remains the authoritative strict
convergence signal.

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
