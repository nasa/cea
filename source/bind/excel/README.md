CEA Excel/VBA Wrapper (First Pass)
==================================
This directory contains a minimal native wrapper that lets Excel/VBA verify the
native bridge by calling exported C functions from 64-bit Windows Excel/VBA.

What this first pass does
-------------------------
- Builds a small shared library named `cea_excel`.
- Exports `int cea_excel_version(char *version, int version_len)`.
- Exports `int cea_excel_test_add(int a, int b)` for trivial smoke testing.
- Returns the existing CEA library version as `major.minor.patch`.
- Provides starter VBA modules for 64-bit Windows smoke testing.

What this first pass does not do
--------------------------------
- It does not build an XLL.
- It does not build an Office.js add-in.
- It does not expose the full CEA C API to VBA.
- It does not install or register the library for Excel automatically.
- It does not embed VBA into `cea_template.xlsm`.

Build option
------------
- CMake option: `CEA_ENABLE_BIND_EXCEL`
- Default: `OFF`
- Current repository presets leave it disabled.

Build commands
--------------
Explicit configure/build:

```bash
cmake -S . -B build-excel -DCEA_ENABLE_BIND_EXCEL=ON
cmake --build build-excel --target cea_excel --config Release
```

Native library output
---------------------
- Single-config generators place the wrapper in `source/bind/excel/lib/`.
- Multi-config generators may place it in a configuration subdirectory such as
  `source/bind/excel/lib/Release/`.

Expected filenames
------------------
- Windows: `cea_excel.dll`

Platform support
----------------
- Supported: 64-bit Windows Excel/VBA
- Not supported: macOS Excel/VBA
- Not supported: 32-bit Office

VBA smoke test files
--------------------
- `vba/modCEADeclare.bas`
- `vba/modCEATest.bas`

Manual workbook setup
---------------------
This first pass does not modify `cea_template.xlsm`. To create a workbook for
manual testing:

1. Create or open a macro-enabled workbook in Excel.
2. Add a sheet named `README`.
3. Import `vba/modCEADeclare.bas`.
4. Import `vba/modCEATest.bas`.
5. Add a button on the `README` sheet and assign
   `WriteCEAVersionToReadme`.
6. Save the workbook as `cea_template.xlsm` if desired.

Testing from Excel on Windows
-----------------------------
1. Build `cea_excel.dll`.
2. Copy `cea_excel.dll` next to the macro-enabled workbook. For local build
   tree testing, `lib/cea_excel.dll`, `lib/Release/cea_excel.dll`, and
   `Release/cea_excel.dll` below the workbook folder are also checked.
3. Keep any native DLL dependencies in the same folder as `cea_excel.dll`.
4. In Excel, run `TestCEAAdd` first, then `TestCEAVersion` or
   `WriteCEAVersionToReadme`.

Known limitations
-----------------
- The native library must be discoverable by Excel/VBA.
- The VBA loader resolves the native library relative to the saved workbook
  location; it does not require editing absolute paths into the VBA modules.
- The wrapper is intentionally kept minimal; this first pass only exposes a
  trivial add function and the version string bridge.
- On Windows, this first pass targets 64-bit Excel.
- 32-bit Office is not supported or planned.
- macOS direct dylib loading from VBA is not supported in this wrapper.
