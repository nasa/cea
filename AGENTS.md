# AGENTS

This repository is a scientific computing library (CEA). Follow these rules:

- Do NOT modify `data/thermo.inp` or `data/trans.inp` unless explicitly instructed.
- Numerical correctness is paramount; preserve bitwise results and scientific behavior.
- Avoid algorithmic changes unless explicitly requested.
- Prefer clarity over cleverness or micro-optimizations.
- Keep changes small and focused; avoid drive-by formatting or whitespace churn.
- In Fortran source files, keep lines within the default free-form GNU limit (target `<= 132`
  characters) so CI builds that do not pass `-ffree-line-length-none` do not fail.
- Maintain backward compatibility for the legacy user base.
- If numerical behavior might change, call it out and add validation or tests when possible.
- Respect layer boundaries: Fortran core in `source/`, C bindings in `source/bind/c/`,
  Python bindings in `source/bind/python/`.
- If the public-facing Python API in `source/bind/python/CEA.pyx` changes, update the
  corresponding type-hinting files in `source/bind/python/cea/` (for example `.pyi` stubs).
- Run minimal validation for touched areas:
  - Core/CMake/Fortran/C changes: run relevant `ctest` targets from the build directory.
  - Python binding changes: run `make py-rebuild` then `pytest source/bind/python/tests`.
  - Documentation-only changes: check referenced commands/paths against `README.md` and `CONTRIBUTING.md`.
- Prefer established workflows:
  - Configure/build with `cmake --preset dev` and `cmake --build build-dev` when appropriate.
  - Use focused tests first, then broader test runs when numerical behavior may be affected.
