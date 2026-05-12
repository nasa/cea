CEA Python Binding
==================
This directory contains the Cython-based Python interface to the CEA core
library. The Python package lives under `source/bind/python/cea`.

Build/install:
- From the repo root: `pip install .` (or `pip install -e .` for editable).
- CMake presets (e.g., `cmake --preset dev`) also enable the binding.

Usage:
- See `docs/source/interfaces/python_api.rst`.
- MATLAB callers should use MATLAB's Python bridge with `import cea` and
  `import cea.matlab`; the MATLAB-oriented wrappers are documented in
  `source/bind/matlab/README.md`.
- Sample scripts live under `source/bind/python/cea/samples`.
