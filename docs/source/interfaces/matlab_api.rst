MATLAB
******

The :mod:`cea.matlab` module provides a MATLAB-friendly compatibility layer on
top of the Python package. It is intended for use through MATLAB's Python
bridge, not as a separate compiled MATLAB extension.

Install the Python package first, point MATLAB at that Python interpreter with
``pyenv(...)``, then import both :mod:`cea` and :mod:`cea.matlab`. Import
:mod:`cea` for enums/constants such as ``cea.TP`` and for :mod:`cea.units`,
then import :mod:`cea.matlab` for the solve entry points.

These wrappers return flat namespace-style results composed of Python scalars,
NumPy arrays, and dictionaries so MATLAB can consume the outputs without
working directly with the compiled solution classes. Example scripts live in
``source/bind/matlab/samples/``.

Available wrappers:

* :func:`cea.matlab.eq_solve` for equilibrium problems
* :func:`cea.matlab.rocket_solve` for rocket performance problems
* :func:`cea.matlab.shock_solve` for incident/reflected shock problems
* :func:`cea.matlab.detonation_solve` for detonation problems

.. autofunction:: cea.matlab.eq_solve

.. autofunction:: cea.matlab.rocket_solve

.. autofunction:: cea.matlab.shock_solve

.. autofunction:: cea.matlab.detonation_solve
