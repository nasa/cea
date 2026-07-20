Excel
*****

The Excel interface provides worksheet user-defined functions (UDFs) for CEA
from a macro-enabled workbook. The worksheet functions are implemented in VBA
and call a native ``cea_excel.dll`` that wraps the C and Fortran CEA APIs.

This interface is beta/experimental. It is included as a complete optional
interface for this release path, but future releases may refine workbook
conventions or UDF signatures as users exercise the interface.

For most spreadsheet-oriented workflows, prefer the Python interface with
Pandas: run CEA from Python, collect results in a ``pandas.DataFrame``, and
export with ``DataFrame.to_csv()`` or ``DataFrame.to_excel()``. Use the
Excel/VBA interface when you specifically need live worksheet formulas backed
by CEA inside 64-bit Windows Excel.

Platform Support
----------------

The Excel interface should only be attempted with 64-bit Windows and 64-bit
Excel/VBA. It is not an XLL add-in and it is not an Office.js add-in. macOS
Excel/VBA and 32-bit Office are not supported.

The build option defaults to ``OFF`` so ordinary CEA builds are unchanged.

Build the Native Library
------------------------

Run the build from a Windows developer prompt with C, C++, and Fortran
compilers available. Configure with the Excel binding enabled, then build the
``cea_excel`` target. For example, with Visual Studio:

.. code-block:: console

   cmake -S . -B build-excel -G "Visual Studio 17 2022" -DCEA_ENABLE_BIND_EXCEL=ON
   cmake --build build-excel --target cea_excel --config Release

If your shell already selects the desired CMake generator and compilers, the
generator argument may be omitted:

.. code-block:: console

   cmake -S . -B build-excel -DCEA_ENABLE_BIND_EXCEL=ON
   cmake --build build-excel --target cea_excel --config Release

The build writes ``cea_excel.dll`` below ``source/bind/excel/lib``. With a
multi-configuration generator, the DLL may be written below
``source/bind/excel/lib/Release``.

Workbook Setup
--------------

Use ``source/bind/excel/cea_template.xlsm`` as the starting workbook, or import
the VBA modules into another macro-enabled workbook:

* ``source/bind/excel/vba/modCEADeclare.bas`` loads and declares the native DLL.
* ``source/bind/excel/vba/modCEAUDF.bas`` exposes the worksheet UDFs.
* ``source/bind/excel/vba/modCEATest.bas`` provides smoke-test macros.

Enable macros according to your organization's Excel security policy, then save
the workbook before using the UDFs. The loader searches for
``cea_excel.dll`` relative to the workbook folder in these locations:

* next to the workbook
* ``lib\`` below the workbook folder
* ``lib\Release\`` below the workbook folder
* ``Release\`` below the workbook folder

Keep required native DLL dependencies in the same folder as ``cea_excel.dll``.
Run the ``TestCEAVersion`` or ``TestCEAAdd`` macros first if formulas cannot
load the native library.

Reactant Table
--------------

All solve and helper UDFs take a reactant table range. The first row must
contain these headers:

.. list-table::
   :header-rows: 1

   * - name
     - role
     - base_amount
     - basis
   * - ``H2(L)``
     - ``fuel``
     - ``1``
     - ``weight``
   * - ``O2(L)``
     - ``oxidizer``
     - ``1``
     - ``weight``

``role`` accepts ``fuel``, ``oxidizer``, and ``inert``. ``basis`` accepts
``weight`` or ``mole``. The species names must match CEA species names in the
thermodynamic database.

For ratio-based solves, the table defines the fixed fuel, oxidizer, and inert
vectors. The worksheet formula supplies the per-row amount mode, and the native
wrapper computes the actual solve weights.

Solve UDFs
----------

Each solve UDF spills one horizontal result row by default. The first three
columns are:

* ``status``: native CEA status code for the call.
* ``converged``: ``TRUE`` when the CEA solution reports strict convergence.
* ``message``: ``OK`` or a diagnostic message from the wrapper/native solver.

Pass ``include_headers=TRUE`` to spill a two-row header/value block. This is
useful while building a sheet or when selecting custom properties/species.

Equilibrium solves:

.. code-block:: text

   CEA_TP_SOLVE(reactants, T, P, [of_ratio], [phi], [r_eq],
                [pct_fuel], [weights], [T_reac], [properties],
                [species], [transport], [ions], [include_headers])
   CEA_HP_SOLVE(reactants, H, P, ...)
   CEA_SP_SOLVE(reactants, S, P, ...)
   CEA_TV_SOLVE(reactants, T, V, ...)
   CEA_UV_SOLVE(reactants, U, V, ...)
   CEA_SV_SOLVE(reactants, S, V, ...)

Rocket solves:

.. code-block:: text

   CEA_ROCKET_IAC_SOLVE(reactants, pc, [pi_p], [subar], [supar],
                        [of_ratio], [phi], [r_eq], [pct_fuel],
                        [weights], [T_reac], [n_frz], [hc], [tc],
                        [tc_est], [omit], [insert], [properties],
                        [species], [transport], [ions], [include_headers])

   CEA_ROCKET_FAC_SOLVE(reactants, pc, [pi_p], [subar], [supar],
                        [of_ratio], [phi], [r_eq], [pct_fuel],
                        [weights], [T_reac], [n_frz], [hc], [tc],
                        [mdot], [ac_at], [tc_est], [omit], [insert],
                        [properties], [species], [transport], [ions],
                        [include_headers])

Shock and detonation solves:

.. code-block:: text

   CEA_SHOCK_SOLVE(reactants, T0, p0, [u1], [Mach1], [of_ratio],
                   [phi], [r_eq], [pct_fuel], [weights], [reflected],
                   [incident_frozen], [reflected_frozen], [omit], [insert],
                   [properties], [species], [transport], [ions],
                   [include_headers])

   CEA_DETONATION_SOLVE(reactants, T1, p1, [of_ratio], [phi], [r_eq],
                        [pct_fuel], [weights], [frozen], [omit], [insert],
                        [properties], [species], [transport], [ions],
                        [include_headers])

Amount Modes
------------

Provide exactly one amount mode to each solve:

* ``of_ratio``: oxidizer/fuel weight ratio.
* ``phi``: fuel-air equivalence ratio.
* ``r_eq``: CEA equivalence ratio.
* ``pct_fuel``: percent fuel.
* ``weights``: explicit per-reactant weights whose length matches the reactant
  table.

For example, with a reactant table in ``A1:D3``:

.. code-block:: text

   =CEA_TP_SOLVE($A$1:$D$3, 3000, 10, of_ratio:=6, include_headers:=TRUE)

Optional ``T_reac`` values may be a scalar or a vector matching the reactant
table, depending on the solve type.

Properties and Species
----------------------

If ``properties`` is omitted, each solver returns a compact default property
set. Pass a worksheet range of property names to control the result columns.
Transport properties are accepted only when ``transport=TRUE``.

Supported property names are:

* Equilibrium: ``T``, ``P``, ``volume``, ``density``, ``M``, ``MW``,
  ``enthalpy``, ``energy``, ``entropy``, ``gibbs_energy``, ``gamma_s``,
  ``cp_fr``, ``cv_fr``, ``cp_eq``, ``cv_eq``, ``viscosity``,
  ``conductivity_fr``, ``conductivity_eq``, ``Pr_fr``, ``Pr_eq``.
* Rocket: all equilibrium names plus ``Mach``, ``sonic_velocity``, ``ae_at``,
  ``c_star``, ``coefficient_of_thrust``, ``Isp``, and ``Isp_vacuum``.
* Shock: ``T``, ``P``, ``velocity``, ``Mach``, ``sonic_velocity``, ``rho12``,
  ``rho52``, ``P21``, ``P52``, ``T21``, ``T52``, ``M21``, ``M52``, ``v2``,
  ``u5_p_v2``, ``volume``, ``density``, ``M``, ``MW``, ``enthalpy``,
  ``energy``, ``entropy``, ``gibbs_energy``, ``gamma_s``, ``cp_fr``,
  ``cv_fr``, ``cp_eq``, ``cv_eq``, ``viscosity``, ``conductivity_fr``,
  ``conductivity_eq``, ``Pr_fr``, ``Pr_eq``.
* Detonation: ``P1``, ``T1``, ``H1``, ``M1``, ``gamma1``,
  ``sonic_velocity1``, ``P``, ``T``, ``density``, ``enthalpy``, ``energy``,
  ``gibbs_energy``, ``entropy``, ``Mach``, ``velocity``,
  ``sonic_velocity``, ``gamma_s``, ``P_P1``, ``T_T1``, ``M_M1``,
  ``rho_rho1``, ``cp_fr``, ``cv_fr``, ``cp_eq``, ``cv_eq``, ``M``, ``MW``,
  ``viscosity``, ``conductivity_fr``, ``conductivity_eq``, ``Pr_fr``,
  ``Pr_eq``.

Species are omitted by default to keep sweep rows narrow. Pass a ``species``
range to append selected mass and mole fractions. Equilibrium and detonation
species columns are labeled like ``mass_H2O`` and ``mole_H2O``. Rocket and
shock species columns include station labels, such as ``mass_H2O_throat`` or
``mole_H2O_incident``.

Helper UDFs
-----------

The helper UDFs use the same reactant table convention:

* ``CEA_WEIGHTS_FROM_OF`` computes per-reactant weights from an O/F ratio.
* ``CEA_OF_FROM_EQUIVALENCE`` and ``CEA_OF_FROM_PHI`` convert ratio inputs to
  O/F.
* ``CEA_EQUIVALENCE_FROM_OF`` and ``CEA_PHI_FROM_OF`` convert O/F to ratio
  outputs.
* ``CEA_WEIGHTS_FROM_MOLES`` and ``CEA_MOLES_FROM_WEIGHTS`` convert reactant
  amounts.
* ``CEA_PER_MOLE_FROM_PER_WEIGHT`` and ``CEA_PER_WEIGHT_FROM_PER_MOLE`` convert
  per-reactant extensive properties.
* ``CEA_CALC_THERMO`` evaluates mixture thermodynamic properties for a set of
  weights and temperatures.
* ``CEA_TO_SI`` and the ``CEA_PRESSURE_SI``, ``CEA_TEMPERATURE_SI``,
  ``CEA_ENERGY_SI``, ``CEA_DENSITY_SI``, and ``CEA_VOLUME_SI`` aliases convert
  supported unit strings to SI values.
* ``CEA_VERSION`` returns the native wrapper version.
* ``CEA_LAST_ERROR`` returns the last native wrapper error message.

Diagnostics and Tests
---------------------

If a workbook reports a DLL load failure, check that the workbook is saved and
that ``cea_excel.dll`` and its dependencies are in one of the loader search
folders listed above. Then run ``TestCEAVersion`` or ``TestCEAAdd`` from
``modCEATest.bas``.

For native validation, configure a Windows build with
``CEA_ENABLE_BIND_EXCEL=ON`` and run:

.. code-block:: console

   ctest -R "cea_excel_(version|api)_test" --output-on-failure
