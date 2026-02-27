# Changelog

All notable user-visible changes to this project are documented here.

## [Unreleased]

### Changed

### Fixed

### Added

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
