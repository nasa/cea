Quick Start
***********

This guide walks through the bare essentials: installing CEA, running the
command-line solver on the supplied sample problems, and calling the Python API.

Prerequisites
-------------

* macOS, Linux, or Windows with a Fortran 2008 compiler (``gfortran`` ≥ 10 or
  Intel ``ifort`` 2021+).
* `CMake <https://cmake.org>`_ ≥ 3.19 and a build tool (Ninja or Make).
* Python ≥ 3.11 if you plan to use the Python binding.

Installation
------------

1. Install from PyPI (recommended for most users)::

       python -m pip install cea

   This provides the Python API immediately (``import cea``) without a local
   source build.
   Prebuilt release assets are also available from the GitHub Releases page,
   including platform executables and shared libraries (``.so``, ``.dll``,
   ``.lib``, ``.dll.a``).

2. If you need the CLI executable (``cea``) or native build artifacts, build
   from source.

   Clone the repository::

       git clone https://github.com/nasa/cea
       cd cea

3. Configure and build (the ``dev`` preset enables the command-line executable,
   the libraries, and the Python binding on Linux/macOS)::

       cmake --preset dev
       cmake --build build-dev

   On Windows, prefer an explicit generator instead of the GNU-based preset::

       cmake -S . -B build-win -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Debug -DCEA_ENABLE_BIND_PYTHON=ON
       cmake --build build-win --config Debug

4. Install into a staging directory (defaults to ``build-dev/install`` unless
   you set ``CMAKE_INSTALL_PREFIX``)::

       cmake --install build-dev

   For the Windows build directory above, install from ``build-win`` instead::

       cmake --install build-win --config Debug

   After installation, add ``<install-prefix>/bin`` to your ``PATH`` so the
   ``cea`` executable is discoverable. On Windows (PowerShell)::

       setx PATH "$env:PATH;<install-prefix>\\bin"

   On macOS/Linux::

       export PATH="<install-prefix>/bin:$PATH"

Running the Sample Problems
---------------------------

CEA ships with the NASA RP-1311 example suite in ``samples/``.  Once the
executable is in your ``PATH`` you can run every problem in one shot::

    cea samples/rp1311_examples.inp

The solver writes results to ``samples/rp1311_examples.out``; open that file in
your editor to inspect the species tables and property profiles.  To run a
single problem, point the executable at a specific ``.inp`` file::

    cea samples/example1.inp

Use the ``-h`` flag to view additional CLI options for controlling verbosity and
output naming.

Quick Python Example
--------------------

The editable install emitted by ``cmake --build`` exposes the Python API.  From
the repository root::

    pip install -e .
    python - <<'PY'
    import numpy as np
    from cea import EqSolver, TP, TEMPERATURE, PRESSURE

    solver = EqSolver()
    solver.set_problem_type(TP)
    solver.set_state_property(TEMPERATURE, 3500.0)  # Kelvin
    solver.set_state_property(PRESSURE, 1.0e6)      # Pa
    solver.set_reactant_fraction("H2", 2.0)
    solver.set_reactant_fraction("O2", 1.0)
    solution = solver.solve()
    print("Equilibrium temperature:", solution.temperature, "K")
    print("Major species:", dict(zip(solution.species_names,
                                     np.round(solution.mass_fractions, 3))))
    PY

CEA does not perform unit conversions by default; inputs and outputs are in the
documented CEA units, and users must convert as needed. Use the conversion
factors in :mod:`cea.units` when working across unit systems.

The ``EqSolver`` and its siblings ``RocketSolver``, ``ShockSolver``, and
``DetonationSolver`` expose the same properties as the Fortran core.  See
:doc:`interfaces/python_api` for the full API reference.

Reporting Issues
----------------

If you encounter a bug or surprising result, please open a GitHub issue so we
can track it.  A solid report includes:

* A short summary plus expected vs. actual behavior.
* Steps to reproduce, including the command line or script you ran.
* Your platform, compiler, and CEA version or Git commit.
* Any relevant input/output artifacts.  If you find a discrepancy from the
  legacy code, include the ``problem.inp`` file so we can reproduce it.


Next Steps
----------

* :doc:`installation` – deeper coverage of build options, database generation,
  and platform-specific notes.
* :doc:`developer_guide` – workflows for contributors and advanced users.
* :doc:`example_index` – detailed documentation of every RP-1311 example.
