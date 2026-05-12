Interfaces
==========

CEA offers several interfaces to interact with the software, allowing users to choose the most suitable method for their workflow. Currently supported interfaces include:

Python API
----------
A modern interface that allows users to interact with CEA using Python, making
it easier to integrate into Python-based workflows and scripts. This interface
also includes the ``cea.matlab`` compatibility wrappers used from MATLAB via
its Python bridge.

.. toctree::
   :maxdepth: 1

   interfaces/python_api

Fortran API
-----------
Fortran subroutine interface.

.. toctree::
   :maxdepth: 1

   interfaces/fortran_api

C API
-----
C-language binding of the Fortran API.

.. toctree::
   :maxdepth: 1

   interfaces/c_api

Legacy Interface
----------------
A command-line interface designed for backward compatibility with legacy versions of CEA (i.e. CEA2). It uses a text-based file format for input and output, making it suitable for users familiar with older versions of the software.

.. toctree::
   :maxdepth: 1

   interfaces/legacy_interface
