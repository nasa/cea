CEA MATLAB Binding
==================
MATLAB support is provided through the Python package namespace.
The supported entry points for MATLAB callers live in ``cea.matlab``.

Recommended approach:
- Import the root package for constants and units:
  ``py.importlib.import_module('cea')``
- Import the MATLAB-oriented wrapper module for solve entry points:
  ``py.importlib.import_module('cea.matlab')``
- Use ``cea.matlab.eq_solve(...)`` for MATLAB-facing equilibrium solves; it
  returns a flat Python namespace of scalars, arrays, and dictionaries rather
  than a raw ``EqSolution`` object.

Legacy note:
- ``source/bind/matlab/ceam.pyx`` is an old experimental artifact and is not the
  supported interface path.

If you need a native MATLAB interface, please open an issue or contribute fixes.
