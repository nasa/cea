# Changelog

All notable user-visible changes to this project are documented here.

## [Unreleased]

### Changed
- Python `Mixture` input validation now accepts `str` and `cea.Reactant` entries (including mixed lists) and no longer accepts raw `bytes` species names (`#53`).
- Added SI-focused custom-reactant handling at the Python API layer: `Reactant.temperature` is specified in K and `Reactant.enthalpy` in J/kg (converted internally for core input) (`#53`).
- Legacy input parsing now supports repeated `outp` dataset keywords (including multiline forms) by merging successive `outp` entries during dataset assembly (`#52`).

### Fixed
- Legacy CLI equilibrium/rocket/shock workflows now propagate `include_ions` into generated product mixtures so ionized products are retained when requested (`#52`).

### Added
- Added C and Python support for custom reactant data (including species not present in `thermo.lib`) in parity with the main interface workflow used by RP-1311 Example 5 (`#53`).
- Added new C-API constructors for generating product mixtures from input-reactant payloads:
  - `cea_mixture_create_products_from_input_reactants` (`#53`).
  - `cea_mixture_create_products_from_input_reactants_w_ions` (`#53`).
- Added a shared bindc parser path for `cea_reactant_input -> ReactantInput` conversion to reduce duplicated C-binding logic (`#53`).
- Added Python `cea.Reactant` and mixed-input `Mixture(...)` support in the Cython binding (`#53`).

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
