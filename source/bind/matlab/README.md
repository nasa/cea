CEA MATLAB Binding
==================
MATLAB support is provided through the Python package namespace rather than a
separate compiled MATLAB extension. Point MATLAB at a Python interpreter with
the installed ``cea`` package, then import the supported wrapper entry points
from ``cea.matlab``.

Recommended approach:
- Import the root package for constants and units:
  ``py.importlib.import_module('cea')``
- Import the MATLAB-oriented wrapper module for solve entry points:
  ``py.importlib.import_module('cea.matlab')``
- Use the MATLAB-facing wrapper entry points:
  ``cea.matlab.eq_solve(...)``
  ``cea.matlab.rocket_solve(...)``
  ``cea.matlab.shock_solve(...)``
  ``cea.matlab.detonation_solve(...)``
- Each wrapper returns a flat Python namespace of scalars, arrays, and
  dictionaries rather than a raw ``EqSolution`` / ``RocketSolution`` /
  ``ShockSolution`` / ``DetonationSolution`` object.
- The wrapper module is pure Python. The ``CEA_ENABLE_BIND_MATLAB`` CMake
  option is a compatibility knob for MATLAB-via-Python workflows; it does not
  build a separate supported ``ceam`` extension.
- Example scripts live in ``source/bind/matlab/samples/``:
  ``equilibrium_example.m``, ``rocket_example.m``, ``shock_example.m``,
  ``detonation_example.m``.

Legacy note:
- ``source/bind/matlab/ceam.pyx`` is an old experimental artifact and is not the
  supported interface path.

If you need a native MATLAB interface, please open an issue or contribute fixes.
