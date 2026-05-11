C API
======

Notes:

- ``cea_mixture_get_species_name`` and ``cea_mixture_get_species_names`` return heap-allocated strings and are deprecated.
  Use ``cea_mixture_get_species_name_buf`` / ``cea_mixture_get_species_names_buf`` instead. If you use the deprecated APIs,
  free results with ``cea_string_free`` or ``cea_string_array_free``.
- Input string arrays passed into the API (species, reactants, omit, insert) are copied by the C/Fortran layer during
  the call. Callers retain ownership and may free or reuse their buffers after the call returns.
- ``cea_rocket_solver_get_eqsolver`` is no longer exposed; use ``cea_rocket_solver_get_size`` and other RocketSolver APIs.
- Rocket, shock, and detonation solution property getters now require a ``len`` argument and will return
  ``CEA_INVALID_SIZE`` if ``len`` is smaller than the internal number of points.
- ``cea_shock_solver_solve`` may return ``CEA_LAST_VALID_SOLUTION`` when the incident-equilibrium solver retains
  a last valid state for inspection without converging the shock iteration. In that case the returned properties
  remain readable, but ``cea_shock_solution_get_converged`` still reports false.

.. doxygenfile:: cea.h
   :project: cea

The C-API used the Fortran bind-c interface to call Fortran routines from C. This provides a way to use Fortran libraries in C programs, as well as an interface for other languages that can call C functions, including Python.
