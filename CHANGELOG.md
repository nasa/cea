# Changelog

All notable user-visible changes to this project are documented here.

## [Unreleased]

### Changed

### Fixed

### Added

## [3.0.1] – Unreleased

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
