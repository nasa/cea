CEA Excel/VBA Wrapper
=====================

This directory contains the 64-bit Windows Excel/VBA interface for CEA. The
public interface is worksheet UDFs backed by a native `cea_excel.dll`; VBA
parses worksheet ranges and the native wrapper owns CEA setup, solve calls, and
property extraction.

Status and recommended workflow
-------------------------------
This interface is beta/experimental. It is included as a complete optional
interface, but future releases may refine workbook conventions or UDF
signatures as users exercise the interface.

For most spreadsheet-oriented workflows, prefer the Python interface with
Pandas. Run CEA from Python, collect results in a `pandas.DataFrame`, and export
with `DataFrame.to_csv()` or `DataFrame.to_excel()` when the goal is to analyze
or share results in spreadsheet form. Use this Excel/VBA interface when you
specifically need live worksheet formulas backed by CEA inside 64-bit Windows
Excel.

Platform support
----------------
- Supported: 64-bit Windows with 64-bit Excel/VBA only.
- Not supported: macOS Excel/VBA, 32-bit Office, XLL, or Office.js.
- Do not attempt to use the workbook interface outside 64-bit Windows Excel.

Installation overview
---------------------
1. Build the native `cea_excel.dll` on 64-bit Windows.
2. Place the DLL and its native dependencies beside the macro-enabled workbook,
   or in one of the supported loader subfolders listed below.
3. Use `cea_template.xlsm`, or import the VBA modules into another `.xlsm`
   workbook.
4. Save the workbook and enable macros according to your organization's Excel
   security policy.
5. Run `TestCEAVersion` or `TestCEAAdd` before building formulas.

Build commands
--------------
Run the build from a Windows developer prompt with C, C++, and Fortran
compilers available. For example, with Visual Studio:

```bash
cmake -S . -B build-excel -G "Visual Studio 17 2022" -DCEA_ENABLE_BIND_EXCEL=ON
cmake --build build-excel --target cea_excel --config Release
```

If your shell already selects the desired CMake generator and compilers, the
generator argument may be omitted:

```bash
cmake -S . -B build-excel -DCEA_ENABLE_BIND_EXCEL=ON
cmake --build build-excel --target cea_excel --config Release
```

The DLL is written below `source/bind/excel/lib/`; multi-config generators may
place it in `source/bind/excel/lib/Release/`.

Workbook setup
--------------
1. Build `cea_excel.dll`.
2. Copy `cea_excel.dll` next to the macro-enabled workbook, or into a `lib`,
   `lib\Release`, or `Release` folder below the workbook folder.
3. Keep native DLL dependencies in the same folder as `cea_excel.dll`.
4. Use `cea_template.xlsm`, or import these VBA modules into another workbook:
   - `vba/modCEADeclare.bas`
   - `vba/modCEAUDF.bas`
   - `vba/modCEATest.bas` for diagnostics
5. Save the workbook before calling the UDFs. The loader resolves DLL paths
   relative to the saved workbook location.

Reactant table
--------------
Solve and helper UDFs expect a table whose first row contains these headers:

| name | role | base_amount | basis |
| ---- | ---- | ----------- | ----- |
| H2(L) | fuel | 1 | weight |
| O2(L) | oxidizer | 1 | weight |

`role` may be `fuel`, `oxidizer`, or `inert`. `basis` may be `weight` or
`mole`. The table is intended to be fixed per sheet, while per-row formulas
vary O/F, phi, equivalence ratio, explicit weights, pressure, temperature, or
rocket schedules.

Worksheet solve UDFs
--------------------
Equilibrium:
- `CEA_TP_SOLVE`
- `CEA_HP_SOLVE`
- `CEA_SP_SOLVE`
- `CEA_TV_SOLVE`
- `CEA_UV_SOLVE`
- `CEA_SV_SOLVE`

Rocket, shock, and detonation:
- `CEA_ROCKET_IAC_SOLVE`
- `CEA_ROCKET_FAC_SOLVE`
- `CEA_SHOCK_SOLVE`
- `CEA_DETONATION_SOLVE`

Each solve UDF spills one horizontal row by default. The first three columns
are `status`, `converged`, and `message`; remaining columns are selected
properties and optional selected species fractions. Pass `include_headers=TRUE`
to spill a two-row header/value block for setup and debugging.

Amount handling
---------------
Provide exactly one amount mode per solve:
- `of_ratio`
- `phi`
- `r_eq`
- `pct_fuel`
- `weights`

For ratio modes, the native wrapper converts the fixed reactant table into fuel
and oxidizer weight vectors and computes the per-row solve weights internally.
For explicit weights, pass a range whose length matches the reactant table.
CEA's existing pressure convention is bar; the unit helper maps Pa/kPa/MPa,
atm, psi, and mmHg into that convention.

Property and species selection
------------------------------
The optional `properties` range limits property output. If omitted, each solver
uses a compact default set. Transport properties are returned only when
`transport=TRUE`.

The optional `species` range controls species-fraction output. Species are not
emitted by default to keep sweep rows narrow. Selected species add stable
columns such as `mass_H2O`, `mole_H2O`, or station-specific rocket/shock
variants.

Helper UDFs
-----------
- `CEA_WEIGHTS_FROM_OF`
- `CEA_OF_FROM_EQUIVALENCE`
- `CEA_OF_FROM_PHI`
- `CEA_WEIGHTS_FROM_MOLES`
- `CEA_MOLES_FROM_WEIGHTS`
- `CEA_EQUIVALENCE_FROM_OF`
- `CEA_PHI_FROM_OF`
- `CEA_PER_MOLE_FROM_PER_WEIGHT`
- `CEA_PER_WEIGHT_FROM_PER_MOLE`
- `CEA_CALC_THERMO`
- `CEA_TO_SI`, plus `CEA_PRESSURE_SI`, `CEA_TEMPERATURE_SI`,
  `CEA_ENERGY_SI`, `CEA_DENSITY_SI`, and `CEA_VOLUME_SI`
- `CEA_VERSION`
- `CEA_LAST_ERROR`

Diagnostics
-----------
`modCEATest.bas` still provides the smoke-test macros `TestCEAAdd`,
`TestCEAVersion`, and `WriteCEAVersionToReadme`. Use these first when a workbook
cannot load the native DLL.

Native tests
------------
From a Windows build directory configured with `CEA_ENABLE_BIND_EXCEL=ON`, run:

```bash
ctest -R "cea_excel_(version|api)_test" --output-on-failure
```
