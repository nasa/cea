FAQs
****

This page collects the most common questions from new users and contributors.

.. _faq-build-from-source:

How do I build CEA from source?
===============================

Configure one of the provided CMake presets, build, and then install.  From the
repository root run::

    cmake --preset dev
    cmake --build build-dev

This creates the ``cea`` executable plus the ``libcea`` library inside
``build-dev``.  See the :doc:`developer guide <developer_guide>` for more
details on compiler requirements and how to customize presets or installation
prefixes.

.. _faq-custom-databases:

Where should I install or customize thermodynamic databases?
============================================================

When you install CEA with ``cmake --install``, the default ``thermo.lib`` and
``trans.lib`` are copied into ``<install-prefix>/data``.  To use modified
databases, place your custom files next to the `.inp` case you plan to run;
files in the working directory take precedence over the installed versions.  The
:doc:`manual database generation workflow <installation>` describes how to
regenerate those files from ``thermo.inp`` and ``trans.inp``.

.. _faq-docs-build:

How do I regenerate the documentation after editing ``docs/source``?
====================================================================

Install the Python requirements for the binding and the docs: ``pip install -e
.`` plus ``python -m pip install sphinx breathe``.  Then run ``make html`` from
the ``docs`` directory.  The rendered site is placed in ``docs/_build/html``.  If you
change public Fortran interfaces, rerun ``doxygen Doxyfile`` first so Breathe
can pick up the latest XML API descriptions.

.. _faq-support:

Where do I report bugs or ask for help?
=======================================

Please open a GitHub issue and include a concise summary, expected vs. actual
behavior, and the exact steps to reproduce.  Helpful details include platform,
compiler, CEA version or Git commit, and any logs or outputs.  If you find a
discrepancy from the legacy code, include the ``problem.inp`` file so we can
reproduce it quickly.  If the issue contains export-controlled content,
coordinate with the maintainers over the approved NASA channels instead of
posting publicly.  Pull requests are always welcome when you already have a fix
in hand.
