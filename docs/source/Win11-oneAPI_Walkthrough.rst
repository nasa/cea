Installation Walkthrough: Win w/Intel oneAPI
********************************************
This walkthrough provides an example of installation on Windows using
the Intel oneAPI compilers. It was developed using oneAPI 2025.3.

The example includes installation of the core CEA application, as well
as the C and Python bindings.

Prerequisites
-------------

This example assumes that the following tools are installed on your system:

* Visual Studio 2019 or later.
* Intel oneAPI Base Toolkit. As an alternate, Intel C++ Essentials should 
  work but is untested.
* Intel oneAPI HPC Toolkit. As an alternate, Intel Fortran Essentials should
  work but is untested.
* Python ≥ 3.11 with ``pip`` for the optional binding, docs, and integration
  tests.
* Git, if you plan to clone the repository and build from source.
* Doxygen, if you plan to generate a local copy of the documentation.

Preparations
------------
Start an Intel oneAPI command prompt so that the compiler variables and paths
are correctly setup. Create a directory to work in and cd into it::

    mkdir \ceabuild
    cd \ceabuild
	
Prepare the Python environment. For example, if using Anaconda::

	conda create --name ceaEnv pip
	conda activate ceaEnv
	
If using vanilla Python::

	py -m venv ceaEnv
	ceaEnv\Scripts\activate

Install Python packages::

	pip install cython numpy setuptools pandas pytest

Clone the repository::

    git clone https://github.com/nasa/cea.git


Build 
-----

The CEA software package is compiled and installed using CMake. The basic
installation process is as follows::

    cd cea
    mkdir build-intel
    cd build-intel
    cmake -DCMAKE_INSTALL_PREFIX=<cea_install_dir> -DCEA_BUILD_TESTING=OFF -DCEA_ENABLE_BIND_PYTHON=ON --preset intel-ifx ..
    cmake --build .
    
Test the Core CEA Executable
----------------------------
The cea/test directory tree contains example problems. It includes Python
code that runs the examples, compares results to reference output, and 
reports whether or not the outputs match.

::

	cd source
	python ..\..\test\main_interface\test_main.py
	
The expected result is 13 out of 14 tests pass, with a small error on 
example 11::

	                    Reference    | Test         | Rel. Error 
	--------------------------------------------------------------
	F-                :   2.5000e-04 |   2.6000e-04 |       4.000%

Different results could indicate that something is wrong with the build.

Create Documentation (Optional)
-------------------------------
Prerequisites:

* The Doxygen documentation generator.  
* The Sphinx Python package, installable with pip.
* The Breathe Python package, installable with pip.

CEA html documentation can be created by installing the Python binding in
edit mode, running Doxygen to create docs for the API, and then using Sphinx
to create the html files. Change directories to the **cea root** (the directory that 
contains Doxyfile) and execute the following commands::

	pip install -e .
	doxygen
	cd docs
	sphinx-build -M html source .
	
Documentation will be generated in the docs/html directory.

Install CEA
-----------
To install CEA, run the following from the build-intel directory::

	cmake --install .
	
This will install the ``cea`` executable, ``libcea`` library, C library, 
and the default thermodynamic and transport property databases to the 
``cea_install_dir`` that was specified when configuring CMake. If you
elected to generate documentation it can be manually copied to the
installation directory.

After installing the executable, install and test the Python binding::

	cd ..
	pip install .
	pytest source\bind\python\tests

All tests should run (no skips!) and pass.

After the Python interface is tested and confirmed to be working, you can make a
wheel for further distribution. Runtime libraries are statically linked so the 
wheel will be self-contained (but it still needs numpy)::

	pip wheel --no-deps -w dist .
	
and the wheel will end up in the cea/dist directory. This directory can be moved
into the CEA installation directory so that the wheel is convienently available
for installation into new Python virtual environments using::

	pip install <path-to-wheel>\<wheel-name>

Setup User Environment
----------------------

After installation, all that is required to use the `cea` executable is to add
the CEA install directory to the user's `PATH` environment variable, e.g.:

::

    setx PATH=<cea_install_dir>\bin;%PATH%

Once properly configured, you should be able to run the provided sample problems
from any working directory as follows:

::

    cea <cea_source_dir>\samples\rp1311_examples.inp


Minimal Builds
--------------

If you want a Fortran-only build or a Fortran+C build without Python/Cython/NumPy
dependencies, use the following options when configuring CMake:

Fortran-only (no C/Python bindings)::

    -DCEA_ENABLE_BIND_PYTHON=OFF -DCEA_ENABLE_BIND_C=OFF -DCEA_ENABLE_BIND_MATLAB=OFF

Fortran + C (no Python bindings)::

    -DCEA_ENABLE_BIND_PYTHON=OFF -DCEA_ENABLE_BIND_C=ON -DCEA_ENABLE_BIND_MATLAB=OFF


Intel oneAPI on Linux
---------------------

This procedure can also be used to build CEA on Linux using Intel oneAPI. To build on Linux,
add an additional option when configuring CMake:
	
	-G "Unix Makefiles"


