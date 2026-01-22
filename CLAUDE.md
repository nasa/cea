# CLAUDE.md - AI Agent Instructions for CEA

This repository contains CEA (Chemical Equilibrium with Applications), a NASA scientific computing library for chemical equilibrium calculations. The codebase consists of validated Fortran solvers with C and Python bindings.

## Core Principles

1. **Numerical correctness is paramount** - Preserve bitwise results and scientific behavior
2. **Backward compatibility matters** - Legacy user base depends on stable behavior
3. **Clarity over cleverness** - Prefer readable code over micro-optimizations
4. **Small, focused changes** - One concern per change; avoid bundling unrelated work

## Critical Rules

### Data Files (DO NOT MODIFY)
- **NEVER** modify `thermo.inp` or `trans.inp` unless explicitly instructed by the user
- **NEVER** modify `data/thermo.inp` or `data/trans.inp`
- These files contain validated thermodynamic and transport property databases

### Algorithmic Changes
- **AVOID** changing solver algorithms unless explicitly requested
- **AVOID** modifying convergence logic, tolerances, or loop ordering
- **CALL OUT** any changes that might affect numerical behavior
- **ADD** validation tests when numerical behavior might change

### Code Quality
- **AVOID** drive-by formatting or whitespace changes
- **AVOID** restructuring code for style alone
- **KEEP** diffs minimal and focused
- **PRESERVE** existing naming conventions and patterns

## Architecture & Layer Boundaries

```
source/                 - Fortran core (scientific solvers)
source/bind/c/          - C ABI bindings (thin adapters)
source/bind/python/     - Python interface (Cython/NumPy)
data/                   - Thermodynamic and transport databases
```

**Respect layer boundaries:**
- Do not mix Fortran solver changes with binding changes
- Do not change scientific logic from binding layers
- Keep C bindings thin (no business logic)
- Python layer may add Pythonic conveniences but not alter solver behavior

## Build System

- **Build tool**: CMake 3.19+
- **Languages**: Fortran (core), C (bindings), Python (bindings)
- **Compilers**: Intel and GNU Fortran are primary targets
- **Python build**: scikit-build-core with Ninja generator
- **Presets available**: `core` (Fortran only), `core-c` (Fortran+C), default (all bindings)

### Python Binding Development
When modifying Python bindings:
1. Changes to `*.pyx` or `*.pxd` files require rebuild
2. Editable installs do not auto-rebuild
3. Rebuild command: `make py-rebuild` (requires Ninja on PATH)
4. Test command: `pytest source/bind/python/tests`
5. Ensure `numpy` is installed in the active Python environment

## Testing

- **Run existing tests** before and after changes
- **Add regression tests** for bug fixes when possible
- **Add validation tests** for new features
- **Verify** no unintended behavior changes occur

Test locations:
- Fortran: Built with `-DCEA_BUILD_TESTING=ON`
- Python: `pytest source/bind/python/tests`
- Legacy CLI validation: `test/main_interface/test_main.py`
- Samples: `samples/` directory contains reference cases

## Common Scenarios

### Bug Fixes
- Identify root cause before fixing
- Add regression test if feasible
- Keep fix minimal and focused
- Document what was broken and how it's fixed

### Documentation
- Improvements always welcome
- Use Markdown or reStructuredText
- Be technically accurate
- Avoid marketing language
- Document assumptions and limitations

### Interface Improvements (C/Python bindings)
- Allowed and encouraged
- Do not change solver behavior
- Add clear docstrings/comments
- Follow language idioms (Pythonic for Python, etc.)

### Refactoring
- Discuss with user first for large refactors
- Prove numerical equivalence
- Provide before/after comparison
- Small, incremental refactors preferred

## What to Call Out

**Always explicitly mention if your change involves:**
- Floating-point arithmetic reordering
- Loop ordering changes
- Tolerance or convergence criterion modifications
- Changes to thermodynamic assumptions
- Database compilation or access patterns
- Changes that might affect bitwise reproducibility

## Dependencies

**Core** (Fortran): No external dependencies
**C bindings**: No additional dependencies
**Python bindings**: Requires:
- Cython
- NumPy
- scikit-build-core
- Ninja (build time)

## References

When in doubt, refer to:
- [CONTRIBUTING.md](CONTRIBUTING.md) - Detailed contribution guidelines
- [README.md](README.md) - Build and usage instructions
- `docs/` - Full documentation

## Summary

Priority order: **Numerical correctness > Backward compatibility > Performance > Style**

When uncertain about a change, ask the user first.
