Developer Guide
===============

The modern CEA code base is intentionally modular so that contributors can work
on the core Fortran solver, the command line executable, the language bindings,
or the documentation in isolation.  This page summarizes the workflows that we
use internally when developing the project.

.. contents::
   :local:
   :depth: 2

Audience and Expectations
-------------------------

This guide targets contributors who plan to build CEA from source or extend the
solver.  It assumes you are comfortable with Git, CMake, and either Fortran or
Python development.  If you only need to *use* CEA, start with the :doc:`quick
start <quick_start>` and :doc:`installation <installation>` guides instead.

Project Layout
--------------

.. list-table::
   :header-rows: 1

   * - Path
     - Description
   * - ``source/``
     - Fortran 2008 implementation of the equilibrium, rocket, shock, and
       detonation solvers plus the command-line driver.
   * - ``source/bind/``
     - Language bindings (Python, C, and planned MATLAB) that wrap the Fortran
       core via Cython/C-APIs.
   * - ``extern/fbasics/``
     - Helper functions vendored in-tree (no top-level git submodules).
   * - ``data/`` & ``samples/``
     - Reference thermodynamic/transport databases and the RP‑1311 example
       problems that we use for regression testing.
   * - ``test/``
     - Python-based integration regression harness plus supporting utilities.
   * - ``docs/``
     - Sphinx documentation sources (``docs/source``); generated snapshots live
       in ``docs/_build/html`` and ``docs/_build/doctrees``.

Prerequisites
-------------

* A Fortran 2008 compiler (``gfortran`` ≥ 10 or recent Intel ``ifort``).
* `CMake <https://cmake.org>`_ ≥ 3.19 and either Ninja or Make.
* Python ≥ 3.11 with ``pip`` and ``virtualenv``/``conda`` (used for the Python
  binding and integration tests).
* Sphinx + Breathe for documentation builds (install separately; not bundled
  with the default Python or conda environment).
* Optional but strongly recommended: `Doxygen <https://www.doxygen.nl/>`_ to
  regenerate the XML that feeds the Sphinx API docs, and `PFUnit
  <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`_ support libraries for
  unit testing (vendored in ``extern/GFE``).

We provide a ``environment.yml`` file for quickly provisioning a ``conda`` /
``mamba`` environment with the complete toolchain::

    mamba env create -f environment.yml
    conda activate cea-dev

Environment Setup
-----------------

1. Clone the repository::

      git clone https://github.com/nasa/cea
      cd cea

2. Create or activate the Python environment (see ``environment.yml`` above) and
   install the binding in editable mode so ``import cea`` resolves to your
   working tree::

      pip install -e .

3. (Optional) export helpful environment variables:

   * ``FC``, ``CC``, ``CXX`` – override compilers if you installed multiple.
   * ``PFUNIT_DIR`` – path to an existing PFUnit installation if you prefer not
     to build the vendored copy.
   * ``CEA_DATA_DIR`` – alternate location for custom ``thermo.lib`` /
     ``trans.lib`` databases when testing bespoke chemistry sets.

Building the Core Solver
------------------------

We use `CMake presets <CMakePresets.json>`_ to capture the most common
configurations.  The ``dev`` preset enables all bindings and produces a debug
build in ``build-dev``::

    cmake --preset dev
    cmake --build build-dev

This compiles the ``cea`` executable (``build-dev/source/cea``) plus the shared
library ``libcea``.  A few tips:

* Set ``CEA_BUILD_TESTING=ON`` (enabled automatically for top-level builds) to
  compile the PFUnit suite and register CTest targets.
* To build with a different compiler family, use the ``dev-intel`` preset or
  create a custom preset that inherits from ``gnu``/``intel`` in
  ``CMakePresets.json``.
* To install into a staging directory, configure with
  ``-DCMAKE_INSTALL_PREFIX=/path/to/install`` and invoke ``cmake --install
  build-dev`` when the build completes.

Python Binding Workflow
-----------------------

The Python package is backed by `scikit-build-core
<https://scikit-build-core.readthedocs.io/>`_ and lives under
``source/bind/python``.  The editable install described above is usually enough
for day-to-day development.  Additional commands that often come in handy:

* Build a fresh wheel: ``pip wheel --no-deps -w dist .``.
* Smoke-test the binding after rebuilding the Fortran core::

      python - <<'PY'
      import cea
      solver = cea.EqSolver()
      print("Loaded CEA", cea.__version__, "with", solver.num_elements, "elements")
      PY

Because the binding links against the ``libcea`` artifacts produced by CMake,
rebuild the project whenever you touch Fortran sources before rerunning Python
tests.

Testing
-------

Unit and regression tests complement each other.  Please add or update tests
whenever you modify solver behavior.

* **PFUnit suite** – Enabled when ``CEA_BUILD_TESTING`` is ``ON``.  The target
  ``cea_core_test`` lives inside ``build-*/source``.  Run all tests via::

      ctest --output-on-failure

  (execute from the corresponding build directory).

  PFUnit is a maintainer-only dependency; to run these tests locally, install
  PFUnit (https://github.com/Goddard-Fortran-Ecosystem/pFUnit) and set ``PFUNIT_DIR``
  to its install prefix.

* **CLI regression test** – ``ctest`` also registers ``cea_main_test`` which
  runs every RP‑1311 example input in ``samples/`` and checks the exit status of
  helper commands such as ``cea -h``.

* **Integration harness** – ``test/main_interface/test_main.py`` re-runs the
  original RP‑1311 problems and compares the generated ``*.out`` files against
  reference outputs.  Ensure the ``cea`` executable is on ``PATH`` or adjust the
  ``run_dir`` variable near the top of the script before running::

      python test/main_interface/test_main.py

Documentation Workflow
----------------------

The published documentation combines handwritten Sphinx pages with API content
generated from the Fortran sources via Doxygen/Breathe.

1. Rebuild the Doxygen XML whenever you touch public Fortran modules::

       doxygen Doxyfile

   The XML will be emitted under ``docs/doxygen/xml``.

2. Build the HTML docs::

       cd docs
       make html

   The rendered site is in ``docs/_build/html``.  This repository also carries a
   snapshot of ``docs/_build/html`` and ``docs/_build/doctrees`` for release publishing and
   quick browsing; those files are generated and can lag behind.  For normal
   documentation updates, commit only ``docs/source`` and regenerate the HTML
   snapshots when preparing a release or when explicitly requested.

.. note::
   Always regenerate the Doxygen XML before running ``make html`` so the C and
   Fortran API pages stay in sync; stale XML is the most common cause of build
   warnings such as duplicate declarations.

Coding Guidelines
-----------------

* **Fortran**

  - Target Fortran 2008 (the compiler flags enforce ``implicit none`` and long
    source lines).
  - Prefer module procedures within the existing ``cea_*`` files; introduce new
    modules when functionality does not fit an existing concern.
  - Keep solver-facing interfaces deterministic: no I/O or global state outside
    of the logging facilities defined in ``cea.f90``.
  - Add PFUnit tests in ``source/*_test.pf`` whenever you touch a solver module.

* **Python**

  - Follow `PEP 8 <https://peps.python.org/pep-0008/>`_ style with type hints.
  - Minimize Python-level loops; all heavy lifting should be delegated to the
    Fortran library via the Cython wrappers in ``CEA.pyx``.
  - Keep the public API documented.  Sphinx pulls docstrings directly from the
    binding, so write docstrings as you go and include usage examples when it
    helps users.

Reporting Issues
----------------

If you encounter a bug or unexpected results, please open a GitHub issue so the
team can track it.  A good report usually includes:

* A short summary and a clear statement of expected vs. actual behavior.
* Minimal, reproducible steps (command line or script).
* Platform details plus compiler and CEA version or Git commit.
* Input and output artifacts that reproduce the issue.  If you find a
  discrepancy from the legacy code, include the ``problem.inp`` file so we can
  reproduce it.

Contribution Workflow
---------------------

1. Fork the repository (export control rules still apply) and open a feature
   branch named ``feat/<short-description>`` or ``fix/<ticket-id>``.
2. Make focused commits with descriptive messages.  Keep unrelated refactors and
   whitespace churn out of functional changes to simplify review.
3. Run the relevant tests and document verification steps in your pull request.
   At minimum this means ``ctest`` plus the integration harness when you change
   solver numerics, and ``make html`` if you touched documentation.
4. Update ``docs/source`` or ``samples/`` whenever the user-facing behavior
   changes.  Small inline notes are preferred over relying solely on release
   notes.
5. Submit a pull request.  A reviewer will work with you to
   iterate until the change is ready to merge.

Release engineering (tagging, packaging, uploading artifacts) is handled by the
core maintainers, but a pull request is more likely to merge quickly when it
contains a clear description, tests, and documentation updates.
