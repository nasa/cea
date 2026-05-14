CEA Language Bindings
=====================
This directory contains the language bindings that wrap the CEA core.

Subdirectories:
- `c/` — C ABI shim and headers (see `c/README.md`)
- `python/` — Python binding (see `python/README.md`)
- `matlab/` — MATLAB usage docs and legacy experimental artifacts (see `matlab/README.md`)
- `excel/` — experimental Excel interface (see `excel/README.md`)

Build configuration:
- Use `cmake --preset dev` to enable the current development bindings.
- Or set `CEA_ENABLE_BIND_C`, `CEA_ENABLE_BIND_PYTHON`,
  `CEA_ENABLE_BIND_MATLAB`, and `CEA_ENABLE_BIND_EXCEL` individually.
- `CEA_ENABLE_BIND_EXCEL` defaults to `OFF` so standard builds remain unchanged.
- The Excel wrapper is currently supported on 64-bit Windows only.
- The MATLAB option is a compatibility knob that forces Python bindings on for
  MATLAB-via-Python workflows; it does not build a separate supported MATLAB
  extension module.
