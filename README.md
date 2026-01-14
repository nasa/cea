[![Build and docs](https://github.com/nasa/cea/actions/workflows/docs.yml/badge.svg)](https://github.com/nasa/cea/actions/workflows/docs.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![Release](https://img.shields.io/github/release/nasa/cea.svg)

# CEA (Chemical Equilibrium with Applications)
<img src="docs/source/images/logo.png" alt="CEA logo" width="200">

A modernized version of NASA's Chemical Equilibrium with Applications.

Online documentation and examples are located at <https://nasa.github.io/cea/>

## Overview
The NASA software package CEA (Chemical Equilibrium with Applications) enables
the rapid solution of chemical equilibrium problems for complex mixtures.
The core solver computes equilibrium product concentrations given a set
of reactants and thermodynamic states. These product concentrations are then
used to compute the thermodynamic and transport properties of the equilibrium
mixture. Applications include estimation of theoretical rocket performance,
Chapman-Jouguet detonation characteristics, and shock-tube parameters for
incident and reflected shocks. Associated with the program are independent
databases with transport and thermodynamic properties of individual species. Over
2000 species are contained in the thermodynamic database.

This software repository is a complete re-implementation of the original CEA
software, with initial development supported by the NASA Engineering & Safety
Center (NESC). The software represents the latest evolution of a series of
computer programs that developed at the NASA Glenn (formerly Lewis) Research
Center since the 1950s. The primary goals of the re-implementation were to modernize
the CEA code base to adopt modern software engineering practices and improve
CEA's ability to interface with other software packages and analysis environments
via well-defined programming APIs in multiple languages.

## Build and Install

The CEA software package is compiled and installed using CMake v3.19+. The core
software has no external dependencies, and is known to build successfully on a
wide range of platforms and using the Intel and GNU Fortran compilers. The basic
installation process is as follows:

    cd <cea_source_dir>
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<cea_install_dir> -DCEA_BUILD_TESTING=OFF ..
    cmake --build .
    cmake --install .

This will build and install the `cea` executable, `libcea` library, default
thermodynamic and transport property databases, documentation, and sample
problems to the user-specified `cea_install_dir`.

Upon installation, all that is required to use the `cea` applications is to add
the CEA install directory to the user's `PATH` environment variable, e.g.:

    export PATH="<cea_install_dir>/bin:$PATH"

Once properly configured, you should be able to run the provided sample problems
from any working directory as follows:

    cea <cea_source_dir>/samples/rp1311_examples.inp

### Build Prerequisites

To build the Python bindings from source, Ninja is required (scikit-build-core
uses the Ninja generator). Ensure `ninja` is available on your `PATH` before
running `pip install .` or `pip install -e .`.

### Minimal Builds

If you want a Fortran-only build or a Fortran+C build without Python/Cython/NumPy
dependencies, use the presets below.

Fortran-only (no C/Python bindings):

    cmake --preset core
    cmake --build build-core
    cmake --install build-core

Fortran + C (no Python bindings):

    cmake --preset core-c
    cmake --build build-core-c
    cmake --install build-core-c

If you are not using presets, set `-DCEA_ENABLE_BIND_PYTHON=OFF` and also disable
the MATLAB wrapper (it forces Python on). For Fortran-only, also set
`-DCEA_ENABLE_BIND_C=OFF`.

### Python Binding

The new Python binding provides direct access to compiled CEA routines.
The basic installation process is as follows:

    cd <cea_source_dir>
    pip install .

A binary wheel distribution can also be generated with the following:

    cd <cea_source_dir>
    pip wheel --no-deps -w dist .

This will build a standalone binary wheel distribution in the `./dist`
directory. This distribution can then be installed on compatible local hosts
with:

    pip install path/to/wheel/<wheel-file-name>

The Python binding to CEA has been successfully compiled and executed on macOS,
Linux, and Windows systems.

## Examples

Legacy CLI (classic `.inp` deck - run this from the `build/source` directory):

    ./cea ../samples/example1

Python example (runs the H2/O2 case after installing the Python bindings):

    python source/bind/python/cea/samples/h2_02.py


## Manual Database Generation
CEA requires thermodynamic and transport property databases. When using the
provided CMake build system, these databases are generated automatically during
the build and installed alongside the `cea` executable.  However, in many
applications it is necessary to perform calculations with modified versions of
the provided databases.

To generate modified databases, run the `cea` program in compilation mode with
the modified `thermo.inp` and `trans.inp` files. See the `data` directory for
baseline versions of these files. This will produce `thermo.lib` and
`trans.lib` in the current directory.

    ./cea --compile-thermo path/to/thermo.inp
    ./cea --compile-trans path/to/trans.inp

To use the customized databases, copy them into the working directory where you
will be executing the `cea` program (usually the same directory as the `.inp`
problem definition file). Database files in the working directory will take
precedence over the installed database files.


## References
 1. McBride, B.J., Zehe, M. J., Gordon, S., "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species",
    NASA TP-2002-211556, 2002. [NTRS](https://ntrs.nasa.gov/citations/20020036214)
 2. McBride, B.J., Gordon, S., and Reno, M.A., "Thermodynamic Data for Fifty Reference Elements",
    NASA TP-3287/REV1, 2001. [NTRS](https://ntrs.nasa.gov/citations/20010021116)
 3. Gordon, S., McBride, B.J., "Thermodynamic Data to 20 000 K for Monatomic Gases",
    NASA TP-1999-208523, 1999. [NTRS](https://ntrs.nasa.gov/citations/19990063361)
 4. Svehla, R.A., "Transport Coefficients for the NASA Lewis Chemical Equilibrium Program",
    NASA TM-4647, 1995. [NTRS](https://ntrs.nasa.gov/citations/19950021761)
 5. McBride, B.J., and Gordon, S., "Computer Program for Calculating and Fitting Thermodynamic Functions",
    NASA RP-1271, 1992. [NTRS](https://ntrs.nasa.gov/citations/19930003779)
