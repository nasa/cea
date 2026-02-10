
CEA (Chemical Equilibrium with Applications)
============================================

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

Getting Started
===============
.. toctree::
   :maxdepth: 1

   installation
   quick_start
   Win11 Install Example <Win11-oneAPI_Walkthrough>

Interfaces
===========
.. toctree::
   :maxdepth: 2

   interface_index

Examples
========
.. toctree::
   :maxdepth: 3

   example_index

References
==========
.. toctree::
   :maxdepth: 1

   faq
   developer_guide
   theory
   ref_list

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
