# DistGeom nanobind migration

## Module converted
`Code/DistGeom` — Python bindings for basic distance geometry operations
(triangle smoothing and bounds matrix embedding).

## Files created/modified

- `Code/DistGeom/nbWrap/DistGeom.cpp` — nanobind version of the bindings
- `Code/DistGeom/nbWrap/CMakeLists.txt` — updated to use `rdkit_nanobind_extension`
- `Code/DistGeom/CMakeLists.txt` — added `if(RDK_BUILD_NANOBIND_WRAPPERS)` block

## Conversion notes

- `#include <RDBoost/python.h>`, `<RDBoost/Wrap.h>`, and `<RDBoost/pyint_api.h>`
  replaced with `<nanobind/nanobind.h>` and relevant nanobind STL headers.
- `namespace python = boost::python` replaced with `namespace nb = nanobind`
  and `using namespace nb::literals` for `"arg"_a` syntax.
- `BOOST_PYTHON_MODULE` replaced with `NB_MODULE`.
- `python::scope().attr("__doc__")` replaced with `m.doc()`.
- `throw_value_error(...)` replaced with `throw nb::value_error(...)`.
- `python::list` parameter in `embedBoundsMatrix` replaced with `nb::list`.
- The weights loop previously used raw CPython sequence/integer APIs
  (`PySequence_GetItem`, `PyInt_AsLong`, `PyFloat_AsDouble`, manual `Py_DecRef`);
  replaced with `nb::cast<nb::sequence>` and `nb::cast<int/double>` for
  cleaner, reference-counted access.
- Docstrings converted from C-string concatenation with line continuations
  to `R"DOC(...)DOC"` raw string literals.
- Copyright updated to include 2026 and "other RDKit contributors".
- `rough_test.py` excluded from nbWrap (test file references the Boost
  extension by name; the `add_pytest` in nbWrap/CMakeLists.txt points back
  to `../Wrap/rough_test.py` so the existing test still runs against
  the nanobind build).
