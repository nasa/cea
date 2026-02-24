# AGENTS

This repository is a scientific computing library (CEA). Numerical correctness and backward compatibility are the top priorities.

## Priorities

1. Numerical correctness
2. Backward compatibility for legacy users
3. Clarity and maintainability
4. Performance

## Hard Constraints

- Do NOT modify `thermo.inp` or `trans.inp` unless explicitly instructed.
- Avoid algorithmic changes unless explicitly requested.
- Keep changes small and focused; avoid drive-by formatting and whitespace churn.
- In Fortran files (`.f90`, `.pf`), do not write lines longer than 132 characters; use continuation with `&` when needed.
- Preserve existing scientific behavior and bitwise-sensitive execution where feasible.

## Where To Look First

Use these files as source of truth before making changes:

- `CONTRIBUTING.md` for contribution expectations and testing policy
- `docs/source/developer_guide.rst` for architecture and workflows
- `README.md` for build/install and database behavior
- `CMakePresets.json` for supported build configurations
- `test/main_interface/test_main.py` for legacy CLI regression harness
- `source/bind/python/tests/` for Python binding test coverage

## Architecture Boundaries

- Fortran core in `source/`: scientific logic and solver behavior
- C bindings in `source/bind/c/`: thin ABI adapters only
- Python bindings in `source/bind/python/`: Python API ergonomics only, not solver behavior changes

Rules:

- Do not move scientific logic into C/Python binding layers.
- Do not mix solver refactors and binding refactors in one change unless explicitly requested.
- Maintain backward-compatible public behavior unless the user requests otherwise.

## Task Routing Matrix

### Fortran Core Changes (`source/`)

- Read: `CONTRIBUTING.md`, `docs/source/developer_guide.rst`
- Validate: build + relevant `ctest` targets, plus regression checks when behavior risk exists

### C Binding Changes (`source/bind/c/`)

- Read: `docs/source/interfaces/c_api.rst`, `CONTRIBUTING.md`
- Validate: C binding tests and affected integration tests

### Python Binding Changes (`source/bind/python/`)

- Read: `docs/source/interfaces/python_api.rst`, `CONTRIBUTING.md`
- Validate: `pytest source/bind/python/tests`
- Note: rebuild editable installs after `*.pyx`/`*.pxd` or native changes

### Docs-Only Changes

- Read: `docs/source/developer_guide.rst` and nearby interface/example docs
- Validate: build docs when needed and ensure command examples match current repo workflows

## Numerical Change Protocol

Apply this protocol if a change touches any of the following:

- Floating-point operation ordering or reductions
- Loop ordering that affects accumulation/results
- Convergence logic or tolerances
- Thermodynamic assumptions, species handling, or database access behavior

Required actions:

1. Explicitly label the change as numerical-risk in your summary/PR notes.
2. Describe the expected numerical impact.
3. Run regression/validation comparisons against prior behavior where feasible.
4. Add or update tests that lock in intended behavior.

## Definition Of Done

Before considering a task complete:

1. Change scope is minimal and focused on one concern.
2. Layer boundaries are respected.
3. Required build/tests for touched layers pass.
4. Numerical-risk changes include explicit impact notes and validation evidence.
5. No protected data file edits (`thermo.inp`, `trans.inp`, and `data/` variants) unless explicitly requested.
6. If behavior or workflow changed, relevant docs/tests were updated.

## Preferred Commands

Use established repo workflows:

- Configure/build: `cmake --preset dev` then `cmake --build build-dev`
- Run tests: `ctest` (from build dir)
- Python tests: `pytest source/bind/python/tests`
- Python binding rebuild (editable): `make py-rebuild`

## Guidance Maintenance

Keep this file concise and operational. If a recurring rule is discovered during implementation, add it to the correct source-of-truth doc (`CONTRIBUTING.md` or developer docs) and keep `AGENTS.md` as a high-signal entrypoint.
