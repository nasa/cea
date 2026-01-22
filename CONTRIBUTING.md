Contributing to CEA (Chemical Equilibrium with Applications)
============================================================

Thank you for your interest in contributing to CEA (Chemical Equilibrium with Applications). CEA is a scientific computing library with a long legacy and a broad user base across industry, academia, and government. Contributions are welcome, but correctness, stability, and clarity are paramount. This document describes how to contribute effectively and responsibly.

## Guiding Principles

CEA is governed by the following priorities:

- Numerical correctness > performance > style
- Backward compatibility matters
- Explicit is better than clever
- Small, well-scoped changes are preferred
- Scientific transparency is required

Changes that improve readability, documentation, or interfaces are often welcome. Changes that alter numerical behavior require careful justification and review.

## Getting Started

### Setting Up Your Development Environment

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/your-username/cea.git
   cd cea
   ```

2. **Choose a development preset**
   CEA provides CMake presets for development. For contributors, use one of:
   - `dev` – GNU/gfortran with all bindings (recommended for most contributors)
   - `dev-intel` – Intel ifort with all bindings
   - `dev-ub-hunt` – GNU with undefined behavior detection (for advanced debugging)

3. **Configure and build**
   ```bash
   cmake --preset dev
   cmake --build build-dev
   ```

4. **Run tests to verify your setup**
   ```bash
   cd build-dev
   ctest
   ```

For Python development, see the [Python Binding Rebuilds](#python-binding-rebuilds-editable-installs) section below.

See [README.md](README.md) for detailed build instructions and requirements.

## Codebase Overview

CEA consists of three primary layers:

- `source/` – Core Fortran solvers and scientific logic
- `source/bind/c/` – C ABI bindings
- `source/bind/python/` – Python interface layer

In addition, CEA relies on validated thermodynamic and transport datasets:

- `thermo.inp`
- `trans.inp`

⚠️ These data files must not be modified without explicit coordination with the maintainers.

## Types of Contributions

We welcome contributions in the following areas:

- Documentation improvements
- Bug fixes (with regression tests when applicable)
- Interface improvements (C and Python bindings)
- Build and packaging improvements
- New tests or validation cases

Large refactors or algorithmic changes should be discussed before implementation.

## Contribution Workflow

1. **Open an Issue (Recommended)**
   Before starting work, please open an issue describing:
   - What you want to change
   - Why the change is needed
   - Whether numerical behavior is expected to change
   This helps avoid duplicated effort and ensures alignment.
   Use the issue templates for bug reports and feature requests:
   - `.github/ISSUE_TEMPLATE/bug_report.md`
   - `.github/ISSUE_TEMPLATE/feature_request.md`

2. **Make Small, Focused Changes**
   Each pull request should ideally address one concern:
   - One bug
   - One feature
   - One documentation improvement
   Avoid bundling unrelated changes.

3. **Respect Layer Boundaries**
   | Layer | Guidelines |
   | --- | --- |
   | Fortran core (`source`) | Do not change algorithms unless explicitly approved |
   | C bindings (`source/bind/c`) | Thin adapters only; explicit ownership rules |
   | Python bindings (`source/bind/python`) | Pythonic APIs, explicit validation, no solver changes |

   Do not mix solver changes with interface or documentation changes in the same PR.

## Numerical Changes & Scientific Integrity

Changes that may affect numerical behavior require extra care.

Required for numerical changes:

- Clear explanation of why the change is needed
- Identification of impacted equations, tolerances, or invariants
- Comparison against previous results (when possible)
- Regression tests or validation evidence

If your change modifies:

- Loop ordering
- Floating-point reductions
- Tolerances or convergence logic
- Thermodynamic assumptions

Call this out explicitly in your PR description.

## Testing Expectations

- New features should include tests when feasible
- Bug fixes should include regression tests when possible

At minimum, contributors should verify that:

- Existing tests pass
- No unintended behavior changes occur

### Running Tests

CEA uses pFUnit (part of the Goddard Fortran Ecosystem) for unit testing. Tests are enabled via the `CEA_BUILD_TESTING` CMake option.

**Prerequisites:**
1. Clone GFE into the `extern/` directory:
   ```bash
   cd extern
   git clone https://github.com/Goddard-Fortran-Ecosystem/GFE.git
   cd GFE
   git submodule update --init
   ```

2. Configure CEA with testing enabled (the `dev*` presets enable this by default):
   ```bash
   cmake --preset dev
   cmake --build build-dev
   ```

**Running all tests:**
```bash
cd build-dev
ctest
```

**Running specific tests:**
```bash
# Run core CEA unit tests with verbose output
ctest -R cea_core_test -V

# Run all main interface tests
ctest -R cea_main_test

# Run C binding tests
ctest -R cea_bindc
```

**Available test suites:**
- `cea_core_test` – Core Fortran solver unit tests
- `cea_main_test*` – Legacy CLI interface tests
- `cea_bindc_*` – C binding validation tests

For Python tests, see the [Python Binding Rebuilds](#python-binding-rebuilds-editable-installs) section.

## Python Binding Rebuilds (Editable Installs)

The Python package uses scikit-build-core with Ninja. Editable installs do not rebuild on import,
so you must rebuild after changing `*.pyx`, `*.pxd`, or native sources.

Recommended workflow:

- Ensure Ninja is installed and available on `PATH`.
- Ensure `numpy` is installed in your active Python environment.
- Rebuild the editable install (uses `--no-build-isolation`): `make py-rebuild`
- Run Python tests: `pytest source/bind/python/tests`

Python database lookup follows the same precedence as the CLI: current working
directory, `CEA_DATA_DIR` (if set), packaged `cea/data`, then the repo `data/`
directory when running from a source checkout. Set `CEA_DATA_DIR` to point to
custom databases if needed.

If you are unsure how to test a change, ask in the issue or PR.

## Documentation Contributions

Documentation improvements are strongly encouraged.

Acceptable documentation formats include:

- Markdown (.md)
- reStructuredText (.rst)
- Plain text (.txt)

Documentation should:

- Be technically accurate
- Avoid marketing language
- Clearly describe assumptions and limitations
- Use consistent terminology with existing CEA documentation

## Style and Formatting

CEA favors:

- Clear, readable code
- Consistent naming
- Minimal diffs

There is no strict formatting tool enforced, but contributors should follow existing conventions in the surrounding code.

## Commit and Pull Request Guidelines

### Branch Naming

Branch names must use one of the following prefixes:

- `feat/` – New features
- `fix/` – Bug fixes
- `docs/` – Documentation changes
- `chore/` – Build, tooling, dependencies
- `refactor/` – Code restructuring (no behavior change)
- `test/` – Test additions or improvements
- `exp/` – Experimental work
- `wip/` – Work in progress

Prefer short, descriptive topic names (e.g., `fix/python-import-error`, `docs/build-instructions`).

### Commits

- Use clear, descriptive commit messages
- Avoid "drive-by" formatting or whitespace changes

### Pull Requests

Your PR description should include:

- Summary of changes
- Motivation
- Impacted components
- Whether numerical behavior changes (yes/no)

The PR template provides a concise checklist:
- `.github/pull_request_template.md`

## What Not to Change

Please do not:

- Modify `thermo.inp` or `trans.inp`
- Rewrite solver algorithms without prior discussion
- Introduce undocumented behavior changes
- Add dependencies without justification

## Review Process

All contributions are reviewed with an emphasis on:

- Scientific correctness
- Maintainability
- Backward compatibility
- Clarity for future contributors

Review feedback is part of the process and is not a rejection of the contribution.

## License and Attribution

By contributing to CEA, you agree that your contributions will be licensed under the project’s license and may be redistributed as part of the project.

## Questions?

If you are unsure about:

- Whether a change is appropriate
- How to structure a contribution
- How to validate a result

Please open an issue or ask before proceeding. Early communication is encouraged.

Thank you for helping improve CEA.
