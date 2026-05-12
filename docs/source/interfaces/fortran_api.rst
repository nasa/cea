Fortran
*******

The Fortran API provides a subroutine interface to the functionality of the
library. This serves as the basis for all other language bindings, such as C
and Python. The sections below are generated directly from the Fortran source
files so the complete set of procedures remains searchable.

Core Library
------------

.. doxygenfile:: cea.f90
   :project: cea
   :no-link:

Inputs and Units
----------------

.. doxygenfile:: input.f90
   :project: cea
   :no-link:

.. doxygenfile:: units.f90
   :project: cea
   :no-link:

.. doxygenfile:: fits.f90
   :project: cea
   :no-link:

Mixture and Property Evaluation
-------------------------------

.. doxygenfile:: atomic_data.f90
   :project: cea
   :no-link:

.. doxygenfile:: mixture.f90
   :project: cea
   :no-link:

.. doxygenfile:: thermo.f90
   :project: cea
   :no-link:

.. doxygenfile:: transport.f90
   :project: cea
   :no-link:

Solvers
-------

.. doxygenfile:: equilibrium.f90
   :project: cea
   :no-link:

.. doxygenfile:: rocket.f90
   :project: cea
   :no-link:

.. doxygenfile:: shock.f90
   :project: cea
   :no-link:

.. doxygenfile:: detonation.f90
   :project: cea
   :no-link:
