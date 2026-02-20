Install
*******
CEA can be installed on Windows, MacOS, and Linux systems.

Recommended: Install from PyPI
------------------------------

For most users, the recommended installation path is the published Python
package:

::

    python -m pip install cea

This installs the Python API (``import cea``) without requiring a local source
build.

If you need the command-line executable (``cea``), C/Fortran artifacts, or
custom compiler/toolchain control, follow the source build instructions below.

Prebuilt release assets are also available on the GitHub Releases page,
including platform executables and shared libraries (``.so``, ``.dll``,
``.lib``, ``.dll.a``).

Prerequisites
-------------

Before building, make sure the following tools are available on your system:

* A Fortran 2008 compiler (``gfortran`` ≥ 10 or Intel ``ifort`` 2021+).  CEA is
  routinely tested with GNU and Intel toolchains.
* `CMake <https://cmake.org>`_ ≥ 3.19 plus a build backend (Ninja or Make).
* Python ≥ 3.11 with ``pip`` for the optional binding, docs, and integration
  tests.
* Git, if you plan to clone the repository and build from source.

Clone the repository::

    git clone https://github.com/nasa/cea
    cd cea

Build and Install from Source
-----------------------------

The CEA software package is compiled and installed using CMake. The basic
installation process is as follows:

::

    cd <cea_source_dir>
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<cea_install_dir> -DCEA_BUILD_TESTING=OFF ..
    cmake --build .
    cmake --install .

This will build and install the `cea` executable, `libcea` library, default
thermodynamic and transport property databases, documentation, and sample
problems to the user-specified `cea_install_dir`.

Upon installation, all that is required to use the `cea` applications is to add
the CEA install directory to the users `PATH` environment variables, e.g.:

::

    export PATH="<cea_install_dir>/bin:$PATH"

Once properly configured, you should be able to run the provided sample problems
from any working directory as follows:

::

    cea <cea_source_dir>/samples/rp1311_examples.inp

.. tip::

   If you regularly switch between compilers or build types, consider using the
   predefined `CMake presets <../CMakePresets.json>`_ (``cmake --preset dev``)
   which configure a debug build with all bindings enabled under ``build-dev``.
   The ``dev`` preset targets GNU toolchains and is not intended for Windows; on
   Windows, use the explicit generator flow shown in the Windows Notes section.

Minimal Builds
--------------

If you want a Fortran-only build or a Fortran+C build without Python/Cython/NumPy
dependencies, use the presets below.

Fortran-only (no C/Python bindings)::

    cmake --preset core
    cmake --build --preset core
    cmake --install build-core

Fortran + C (no Python bindings)::

    cmake --preset core-c
    cmake --build --preset core-c
    cmake --install build-core-c

If you are not using presets, set ``-DCEA_ENABLE_BIND_PYTHON=OFF`` and also
disable the MATLAB wrapper (it forces Python on). For Fortran-only,
add ``-DCEA_ENABLE_BIND_C=OFF``.


Selecting Compilers and Generators
----------------------------------

CEA honors standard CMake environment variables and arguments. Useful options:

* Linux/macOS (CMake): set ``CC``, ``CXX``, and ``FC`` in your shell before
  configuring. Example::

    CC=gcc CXX=g++ FC=gfortran cmake -DCMAKE_BUILD_TYPE=Release ..

* Linux/macOS (pip/scikit-build): use ``CMAKE_ARGS`` to pass CMake settings::

    CMAKE_ARGS="-DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_GENERATOR=Ninja" pip install .

* Windows: run from a developer prompt and select a generator explicitly::

    cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release ..

  For Intel oneAPI, set ``FC=ifort`` (or ``ifx``) before configuring.

GNU Toolchain Defaults (Optional)
---------------------------------

If you want to use the historical GNU defaults (gcc/g++/gfortran + flags) you
can opt into the provided toolchain file::

    CMAKE_ARGS="-GNinja -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/gnu.cmake" pip install .

For wheel builds, the convenience script uses the same defaults::

    bash scripts/build_wheel.sh

Running Tests
-------------

Set ``CEA_BUILD_TESTING=ON`` when configuring to compile the PFUnit test suite
and register ctest entries::

    cmake -DCMAKE_BUILD_TYPE=Release -DCEA_BUILD_TESTING=ON ..
    cmake --build .
    ctest --output-on-failure

The top-level ``cea_main_test`` target runs every RP-1311 sample via the CLI.
The ``cea_core_test`` PFUnit executable exercises the Fortran modules directly.

Python Binding
--------------

The new python binding provides direct access to compiled CEA routines.
The basic installation process is as follows:

::

    cd <cea_source_dir>
    pip install .

For iterative development against a local checkout, use editable mode instead::

    pip install -e .

A binary wheel distributions can also be generated with the following:

::

    cd <cea_source_dir>
    pip wheel --no-deps -w dist .

This will build a stand-alone compiled binary wheel distribution in ./dist
directory. This distribution can then be installed on compatible local hosts
with:

::

    pip install path/to/wheel/<wheel-file-name>

The python binding to CEA has been succesfully compiled and executed on MacOS,
Linux, and Windows systems.

Manual Database Generation
--------------------------
CEA requires thermodynamic and transport property databases. When using the
provided CMake build system, these databases are generated automatically during
the build and installed along side the `cea` executable.  However, in many
applications it is necessary to perform calculations with modified versions of
the provided databases.

To generate modified databases, run the `cea` program in compilation mode with
the modified `thermo.inp` and `trans.inp` files. See the `data` directory for
baseline versions of these files. This will produce `thermo.lib` and
`trans.lib` in the current directory.

::

    ./cea --compile-thermo path/to/thermo.inp
    ./cea --compile-trans path/to/trans.inp

To use the customized databases, copy them into the working directory where you
will be executing the `cea` program (ususally the same directory as the `.inp`
problem definition file). Database files in the working directory will take
precedence over the installed database files.

Windows Notes
-------------

* When building with Microsoft Visual Studio or Intel oneAPI on Windows, run the
  commands from a developer prompt so the compiler environment variables are set
  correctly.
* Use PowerShell syntax for environment variables::

      setx PATH "$env:PATH;<cea_install_dir>\\bin"

* If Ninja is unavailable, add ``-G "Visual Studio 17 2022"`` (or a similar
  generator) to the ``cmake`` configure command.
* Windows Subsystem for Linux (WSL) provides a native Linux environment.  Install
  a recent Ubuntu distribution from the Microsoft Store, then inside the WSL
  shell install build tools with::

      sudo apt update
      sudo apt install gfortran ninja-build cmake python3 python3-pip git

  After that follow the standard Linux instructions in this guide.  Files you
  compile inside WSL live in the Linux filesystem (``\\wsl$``) and can be run
  directly from the WSL prompt.
