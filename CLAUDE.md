# AGENTS

This repository is a scientific computing library (CEA). Follow these rules:

- Do NOT modify `thermo.inp` or `trans.inp` unless explicitly instructed.
- Numerical correctness is paramount; preserve bitwise results and scientific behavior.
- Avoid algorithmic changes unless explicitly requested.
- Prefer clarity over cleverness or micro-optimizations.
- Keep changes small and focused; avoid drive-by formatting or whitespace churn.
- Maintain backward compatibility for the legacy user base.
- If numerical behavior might change, call it out and add validation or tests when possible.
- Respect layer boundaries: Fortran core in `source/`, C bindings in `source/bind/c/`,
  Python bindings in `source/bind/python/`.
