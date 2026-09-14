Quick Start
***********

This guide walks through the bare essentials: installing CEA, running the
command-line solver on the supplied sample problems, and calling the Python
and MATLAB APIs.

Prerequisites
-------------

* macOS, Linux, or Windows with a Fortran 2008 compiler (``gfortran`` ≥ 10 or
  Intel ``ifort`` 2021+).
* `CMake <https://cmake.org>`_ ≥ 3.19 and a build tool (Ninja or Make).
* Python ≥ 3.11 if you plan to use the Python binding.
* MATLAB, if you plan to use the MATLAB binding — no prior Python
  experience needed. The Quick MATLAB Example below walks through
  everything, including installing Python itself if you don't already have
  a copy.

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

After installing the Python binding, solve a stoichiometric H\ :sub:`2`/O\ :sub:`2`
constant-enthalpy, constant-pressure (HP) combustion problem as follows::

    import numpy as np
    import cea

    reactants = cea.Mixture(["H2", "O2"])
    products = cea.Mixture(["H2", "O2"], products_from_reactants=True)
    solver = cea.EqSolver(products, reactants=reactants)
    solution = cea.EqSolution(solver)

    moles = np.array([2.0, 1.0])
    weights = reactants.moles_to_weights(moles)
    initial_enthalpy = reactants.calc_property(cea.ENTHALPY, weights, 298.15)
    pressure = cea.units.atm_to_bar(1.0)

    solver.solve(solution, cea.HP, initial_enthalpy / cea.R, pressure, weights)
    print(f"Adiabatic flame temperature: {solution.T:.1f} K")

CEA does not perform unit conversions by default; inputs and outputs are in the
documented CEA units, and users must convert as needed. Use the conversion
factors in :mod:`cea.units` when working across unit systems.

The ``EqSolver`` and its siblings ``RocketSolver``, ``ShockSolver``, and
``DetonationSolver`` expose the same properties as the Fortran core.  See
:doc:`interfaces/python_api` for the full API reference.

Quick MATLAB Example
--------------------

CEA doesn't ship a native MATLAB toolbox. Instead, MATLAB calls CEA through
a small bridge to Python, using MATLAB's built-in ``pyenv`` feature. You
don't need to know any Python to use it — follow the steps below once, then
the two commands under "Every MATLAB Session" are all you'll retype.

If you already have a working Python installation with ``cea`` installed,
skip to "Every MATLAB Session" below.

One-Time Setup
~~~~~~~~~~~~~~

1. Install Python, if you don't already have it. Download the Windows
   installer for Python 3.12 from
   `python.org <https://www.python.org/downloads/>`_ and run it. On the
   first installer screen, check **"Add python.exe to PATH"** before
   clicking "Install Now" — this lets you type ``python`` in a Command
   Prompt window. (Python 3.12 is used here because it works with every
   current MATLAB release; if MATLAB later refuses to load it, see the
   troubleshooting note at the end of this step.)

2. Open a Command Prompt window — a plain text window for typing commands,
   separate from MATLAB. Click the Start menu (or press the Windows key),
   type ``cmd``, and press Enter, or click "Command Prompt" in the search
   results. If you had a Command Prompt window open before you installed
   Python, close it and open a new one — it won't see the update otherwise.
   Then install ``cea``::

       python -m pip install cea

   This downloads a ready-to-use package — no compiler, no conda, nothing
   else to build.

   *Troubleshooting:* if MATLAB later reports that this Python version
   isn't supported, your MATLAB release may need an older or newer Python
   than 3.12. Check MathWorks' `Python compatibility table
   <https://www.mathworks.com/support/requirements/python-compatibility.html>`_
   for your release, install that version from python.org instead (same
   steps as above), and run ``python -m pip install cea`` again using that
   version.

3. In the same Command Prompt window, find the full path to the Python you
   just installed — you'll paste it into MATLAB below::

       where python

   This prints one or more paths ending in ``python.exe``; copy the one
   under the Python version you just installed (e.g.
   ``C:\Users\<you>\AppData\Local\Programs\Python\Python312\python.exe``).

Every MATLAB Session
~~~~~~~~~~~~~~~~~~~~

Paste these lines into the MATLAB Command Window, using the path from step 3
above, before doing anything else with ``cea``::

    pyenv('Version', 'C:\path\to\python.exe');

    cea = py.importlib.import_module('cea');
    ceam = py.importlib.import_module('cea.matlab');

This only needs to run once per MATLAB session — running ``pyenv`` a second
time after these lines have already run will error, so if you need to
change the Python path, restart MATLAB first.

*Tip:* save these three lines as a MATLAB script, e.g. ``setup_cea.m``, so
each session you just type ``setup_cea`` instead of retyping them.

Solving a Problem
~~~~~~~~~~~~~~~~~

With the session set up, solve a stoichiometric H\ :sub:`2`/O\ :sub:`2`
constant-enthalpy, constant-pressure (HP) combustion problem — the
adiabatic flame temperature of hydrogen burning in oxygen::

    reactants = py.list({'H2', 'O2'});
    pressure = cea.units.atm_to_bar(1.0);

    solution = ceam.eq_solve(cea.HP, reactants, ...
        fuel_amounts=py.numpy.array([2.0, 0.0]), ...
        oxid_amounts=py.numpy.array([0.0, 1.0]), ...
        moles=true, ...
        T_reac=298.15, ...
        P=pressure);

    fprintf('Adiabatic flame temperature: %.1f K\n', solution.T);

This should print ``Adiabatic flame temperature: 3074.5 K``.
``fuel_amounts`` and ``oxid_amounts`` each list one amount per entry in
``reactants``: ``[2.0, 0.0]`` is 2 mol of ``H2`` and 0 mol of ``O2`` on the
fuel side, ``[0.0, 1.0]`` is 0 mol ``H2`` and 1 mol ``O2`` on the oxidizer
side — together, 2 mol H2 to 1 mol O2.

``solution`` holds the result as plain numbers and arrays you can read
directly with dot notation, the same as any other MATLAB struct — no
further conversion needed. ``solution.T`` above is the temperature in K;
see :doc:`interfaces/matlab_api` for the full list of result fields and the
other three solver functions (rocket, shock, and detonation problems).

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
