# Changelog

All notable user-visible changes to this project are documented here.

## [Unreleased]

### Changed

### Fixed
- Corrected downstream shock Mach numbers to use finalized gas speeds and
  reported sound speeds, including frozen and retained last-valid states.
  Documented shock reference frames and the legacy reflected-frozen sound-speed
  basis; `M21`/`M52` remain molecular-weight ratios, with no thermochemical or
  legacy output changes.
- Fixed premature convergence of minor elements such as carbon in database Air
  by checking element-relative balances and gas updates. This can add Newton
  iterations and slightly change equilibrium results.
- Hard-mode equilibrium sensitivities now use the species threshold saved at
  convergence, rather than the final reporting threshold. Combined with the
  convergence correction, this resolves the H2/Air CO/CO2 pressure discrepancy.
- Preserved double precision for both state inputs to Python `EqSolver.solve`,
  removing single-precision rounding from finite-difference perturbations.

### Added
- Database-Air derivative regression cases and a reproducible paper/optimization
  audit; see [the investigation](docs/validation/database_air_derivatives.md),
  including the remaining limitations at hard species-activation boundaries.

## [3.3.3] - 2026-08-24

### Changed
- Reworked equilibrium total-derivative post-processing as a vector-valued
  response Jacobian and reused one Jacobian factorization across all direct
  sensitivity right-hand sides. Smooth-mode derivatives now analytically
  differentiate the convergence threshold instead of using finite-difference
  stabilization for small reported species.

### Fixed

### Added

## [3.3.2] - 2026-08-05

### Changed
- Updated the RP-1311 Python documentation examples and quick start to use the
  current API and match the maintained sample scripts.

### Fixed
- Fixed volume and density property queries for liquid and other condensed
  reactant mixtures, including consistent SI conversions for scalar and
  per-species temperature inputs.

### Added
- Added Fortran, C, and Python APIs for querying a database reactant's inclusive
  valid temperature range in Kelvin before constructing a mixture.

## [3.3.1] - 2026-07-29

### Changed

### Fixed
- Frozen rocket calculations now switch from equilibrium to frozen composition
  after `n_frz` when the boundary falls within a multi-point exit schedule,
  restoring CEA2 behavior for both IAC and FAC calculations.
- Equilibrium solves now support gas temperatures through the CEA2
  high-temperature fit range instead of reporting non-convergence above 6600 K.

### Added

## [3.3.0] - 2026-07-20

### Changed
- Excel bindings are now an optional beta/experimental interface that defaults
  to disabled in standard builds and is documented as supported only for
  64-bit Windows Excel/VBA.

### Fixed

### Added
- Added a beta/experimental Excel/VBA interface with worksheet UDFs backed by a
  native `cea_excel.dll`, including equilibrium, rocket, shock, detonation,
  unit-conversion, mixture-ratio, and thermodynamic helper functions.
- Added the Excel native wrapper API, VBA modules, macro-enabled template
  workbook, CMake build target, and native smoke/API tests.
- Added Excel interface documentation, installation guidance, workbook setup
  instructions, and recommendations to use the Python interface with Pandas
  when spreadsheet-formatted exports are sufficient.

## [3.2.1] - 2026-06-23

### Changed
- Legacy CLI `outp long` now enables extra-precision numeric output formatting
  across equilibrium, rocket, shock, and detonation reports.
- Python `cea.eq_solve` remains a deprecated compatibility shim, but now
  returns the legacy object-based `EqSolution` result again via the compiled
  binding instead of forwarding to the MATLAB wrapper.
- Python initialization is quiet at the default log level; database load paths
  are printed only when the log level is `LOG_INFO` or `LOG_DEBUG`.
- MATLAB-facing docs now point to the Python package namespace instead of an
  active separate MATLAB extension.

### Fixed
- C and Python binding calls now recover fatal Fortran `abort` paths as
  `CEA_FORTRAN_ABORT` / Python `RuntimeError` with the recovered message instead
  of terminating the caller process.
- Legacy input parsing now treats bare `h` and `u` reactant energy units as
  `J/mole`, matching the intended molar default.
- Legacy SI transport output now labels conductivity as
  `MILLIWATTS/(CM)(K)` instead of the non-SI millicalorie label.
- Rocket throat states are now marked complete so throat properties are printed
  even when no exit-cone condition is requested.
- Corrected the malformed Paraffin reactant record in `data/thermo.inp`.

### Added
- Added the public `CEA_FORTRAN_ABORT` error code and C helpers for retrieving
  the last recovered fatal Fortran error message.
- Added a pure-Python `cea.matlab` module with a Matlab-oriented `eq_solve`
  entry point that returns a flat Python namespace of scalars, arrays, and
  dictionaries.
- Added Matlab-oriented `cea.matlab.rocket_solve`, `shock_solve`, and
  `detonation_solve` wrapper entry points, plus MATLAB sample scripts and
  docs for those problem types.

## [3.2.0] - 2026-05-12

### Changed
- Python `cea.Reactant` now supports explicit `enthalpy_units` for custom reactant enthalpy input (`J/kg`, `kJ/kg`, `cal/kg`, `kcal/kg`, `J/mol`, `kJ/mol`, `cal/mol`, `kcal/mol`; `/mol` and `/mole` spellings accepted).
- Python `cea.Reactant` now requires explicit `enthalpy_units` whenever `enthalpy` is provided; omitting the units now raises `ValueError` instead of falling back to the legacy implicit `J/kg` default.
- Weight-based enthalpy units (`*/kg`) now require `Reactant.molecular_weight` for conversion, while molar units (`*/mol`) do not require molecular weight for enthalpy interpretation.
- Python docs and examples were updated to pass explicit `enthalpy_units` for custom reactants (including RP-1311 Example 5) instead of relying on implicit defaults.
- Smooth-truncation derivative assembly in the equilibrium solver now uses an effective Jacobian/chain-rule formulation across TP/HP/SP/TV/UV/SV solves so reported sensitivities track the smoothed species-activity model more consistently near truncation thresholds.

### Fixed
- Smooth-truncation total derivatives for reported gas species and derived thermodynamic outputs were aligned with the final reported-species mapping, improving behavior for trace-species sensitivities near truncation cutoffs.
- Corrected user-facing viscosity unit labels in the equilibrium core comments and Python bindings to report millipoise consistently.

### Added
- Added Python regression coverage for `Reactant` enthalpy-unit handling, including required explicit units, explicit molar/weight units, molecular-weight requirements for weight units, and accepted unit spellings.
- Added pfunit regression coverage for smooth-truncation derivative chain-rule behavior, targeted finite-difference agreement across TP/HP/SP/TV/UV/SV constraints, and reported-entropy consistency checks.

## [3.1.2] - 2026-03-23

### Changed
- Legacy CLI `.out` files now echo each problem’s original input block before the corresponding equilibrium, rocket, shock, or detonation output section.
- Frozen rocket failure handling was aligned with the established partial-output behavior: successful upstream stations are retained, output truncates at the last valid station, and the existing warning text is emitted after the partial report.
- Frozen rocket stop checks now use the active thermo-fit interval for condensed species plus the existing lower-temperature frozen-fit guard, and frozen retained points now populate transport properties before output.

### Fixed
- Non-converging frozen rocket cases no longer produce empty output, overrun into invalid exit columns, or crash while writing partial results.
- Difficult frozen rocket chamber solves can now retain a usable reduced-component state after repeated singular matrices, restoring partial output for reproduced failure cases.

### Added

## [3.1.1] - 2026-03-18

### Changed
- Python `Mixture` input validation now accepts `str` and `cea.Reactant` entries (including mixed lists) and no longer accepts raw `bytes` species names (`#53`).
- Added SI-focused custom-reactant handling at the Python API layer: `Reactant.temperature` is specified in K and `Reactant.enthalpy` in J/kg (converted internally for core input) (`#53`).
- Legacy input parsing now supports repeated `outp` dataset keywords (including multiline forms) by merging successive `outp` entries during dataset assembly (`#52`).
- FAC rocket chamber-closure iteration logic in `RocketSolver_solve_fac` was updated toward CEA2 parity: Option-1 pressure correction direction now follows legacy semantics, the Option-1 convergence check is normalized to assigned injector pressure, the fixed 4-pass outer loop was replaced with tolerance-driven iteration plus a bounded safety guard, and FAC combustor-end reseeding now refreshes from the current infinity state each chamber iteration (`#54`).
- Legacy SHCK branch sequencing and output selection were reworked toward CEA2 parity: explicit incident requests now report only the requested incident branch (with reflected state attached when requested), while reflected-only requests continue exploring the mixed equilibrium/frozen permutations.
- Shock reports now include incident/reflected transport tables (`Visc`, `Cp`, `Conductivity`, `Prandtl Number`) when `outp tran` is enabled and transport data are available.
- Python test/main-interface helper workflows no longer require `pandas`; CSV result output is now written with the standard library and the installation guidance/environment were updated accordingly.

### Fixed
- Legacy CLI equilibrium/rocket/shock workflows now propagate `include_ions` into generated product mixtures so ionized products are retained when requested (`#52`).
- Reflected-shock velocity outputs were corrected for reflected equilibrium/frozen workflows: the legacy CEA2 `u5`/`u5+v2` relations effectively used `rho5/rho2` where the reflected-shock mass balance requires `rho5/rho2 - 1`, so the old equations underpredicted reflected velocities and did not satisfy the reflected conservation equations. CEA now reports conservation-consistent `u5`, `u5+v2`, sonic speed, reflected ratios, and reflected-frozen header thermodynamic values.
- Shock failure handling now stops at the last valid state instead of emitting partially invalid reflected results; failed downstream shock states are explicitly zeroed for consistent output.
- Shock transport-property calculations now preserve and restore stable transport seeds across singular-recovery paths, fixing incorrect or missing transport values in difficult reflected-shock cases.

### Added
- Added C and Python support for custom reactant data (including species not present in `thermo.lib`) in parity with the main interface workflow used by RP-1311 Example 5 (`#53`).
- Added new C-API constructors for generating product mixtures from input-reactant payloads: `cea_mixture_create_products_from_input_reactants` and `cea_mixture_create_products_from_input_reactants_w_ions` (`#53`).
- Added a shared bindc parser path for `cea_reactant_input -> ReactantInput` conversion to reduce duplicated C-binding logic (`#53`).
- Added Python `cea.Reactant` and mixed-input `Mixture(...)` support in the Cython binding (`#53`).
- Added shock regression coverage for transport output population across equilibrium/frozen branches and for singular-recovery reflected-shock cases.
- Added a Python `pytest` reflected-shock conservation check that verifies incident/reflected mass, momentum, and energy balances and confirms the corrected `u5+v2` relation.

## [3.1.0] - 2026-03-02

### Changed
- Added equilibrium total-derivative capabilities across core solver paths (TP/HP/TV/UV/SP/SV), including chain-rule state/property derivatives and related thermo/property plumbing updates (`#50`).
- Exposed derivative workflows in the C and Python interfaces, including derivative object lifecycle, analytic-vs-finite-difference result accessors, and central finite-difference verification options (`#50`).
- Introduced optional smooth species truncation controls (with configurable width) through solver options in Fortran/C/Python interfaces to improve derivative behavior near composition cutoffs (`#50`).
- C binding headers and the RP-1311 C sample were updated to resolve compiler warnings (`#49`).

### Fixed

### Added
- Added extensive regression coverage for derivative correctness and stability, including new pfunit derivative suites and Python tests for smooth-truncation behavior and sample execution (`#50`).
- Added derivative-focused Python examples and documentation (including SP and multi-state derivative sample scripts and run helpers) (`#50`).

## [3.0.4] - 2026-02-27

### Changed
- Command-line input parsing now accepts explicit `.inp` filenames (`#44`).
- Python RP-1311 sample scripts were moved into `source/bind/python/cea/samples/rp1311/` for clearer organization (`#47`).
- Expanded `reac` input compatibility with CEA2-style forms, including case-insensitive keywords, exploded formulas with implicit coefficients, `den` density aliases, molecular-weight aliases, and stricter mixed mole/weight basis validation (`#48`).
- User-specified reactant enthalpy input is now applied as a runtime override for reactant thermo initialization (including database species) (`#48`).

### Fixed
- Fixed a crashing output case and restored missing output values (`#45`).
- Reusing `EqSolution` across solve calls now resets transient iteration state and recovers from prior non-converged attempts using the last stable warm-start seed, preventing reuse-related non-convergence regressions (`#47`).

### Added
- Added missing Python test dependencies to improve out-of-the-box test runs (`#41`).
- Added Fortran and Python regression tests covering `EqSolution` reuse and detonation/equilibrium convergence behavior (`#47`).
- Added `reac` parser regression tests for custom species inputs, molecular-weight aliases, density aliases/default units, case-insensitive tokens, and implicit formula coefficients (`#48`).

## [3.0.3] - 2026-02-20

### Changed
- Hardened PyPI publishing CI with OIDC preflight checks and manual `workflow_dispatch` target selection (`testpypi`/`pypi`); installation docs now lead with `python -m pip install cea` and mention GitHub Releases binary assets. (`#39`, `#40`).

### Fixed

### Added

## [3.0.2] – 2026-02-20

### Added
- Initial Microsoft Visual C++ (MSVC) support in the C bindings/build pipeline (`#29`).

### Changed
- `RocketSolver_solve` now treats `pi_p` as optional for better API compatibility and easier caller usage (`#36`).

### Fixed
- Corrected detonation-wave upstream enthalpy calculations and upstream molecular weight handling in the legacy interface (`#33`).
- Resolved failing transport-property test cases to restore expected regression behavior (`#37`).

## [3.0.1] – 2026-02-11

### Added
- Python bindings and examples/tests for the CEA interface (`Feat/icx python binding`, #23).
- `rp1311_examples.inp` and updated C sample output including `M (1/n)` (`#24`).
- Windows 11 Intel oneAPI walkthrough documentation.
- GitHub issue and pull request templates for contribution workflow.
- A local `findloc` implementation in `extern/fbasics` for legacy compiler support.

### Changed
- Core equilibrium/rocket/shock solver iteration logic and convergence rules were refined, including area-ratio loops, throat updates, and frozen-flow gating/criteria.
- Condensed-species handling was tightened across initialization and iteration paths (ordering, active-element indexing, phase checks, and first-iteration thermo evaluation).
- Build and CI configuration was expanded and hardened (compiler flags/presets, broader build/test coverage, and workflow updates).
- Precompiled thermo/transport database blobs were removed from the repository/package test assets in favor of build-time generation.
- Documentation, contribution guidance, and sample scripts were updated for current workflows.

### Fixed
- Singular-matrix recovery paths in equilibrium condensed-species solves (including element-row singularities and related index handling).
- Multiple convergence and initialization edge cases in shock/frozen/throat calculations (including denominator guards and corrected initial values).
- Several compiler-compatibility issues (including Intel ifx workarounds and initialization safety fixes in mixture/transport paths).
- Intermittent CTest instability from an extraneous C format specifier.
- Sample/README correctness issues (including reversed example arguments and text typos).

## [3.0.0] – 2025-12-31
- Initial public release.
