//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  Python wrapper for the BulkSimilarity library. The result is returned as
//  a 2D NumPy array of double-precision Tanimoto coefficients.

#define PY_ARRAY_UNIQUE_SYMBOL rdbulksimilarity_array_API

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>
#include "numpy/arrayobject.h"

#include <DataStructs/BitVects.h>
#include <DataStructs/BulkSimilarity/BulkSimilarity.h>
#include <DataStructs/ExplicitBitVect.h>

#include <cstdint>
#include <string>
#include <vector>

namespace python = boost::python;

namespace {

//! RAII handle for a PyArrayObject. Decrefs on scope exit unless release()
//! is called. Used to make sure half-built result arrays don't leak when
//! the C++ kernel throws.
class PyArrayGuard {
 public:
  explicit PyArrayGuard(PyArrayObject *arr) : d_arr(arr) {}
  PyArrayGuard(const PyArrayGuard &) = delete;
  PyArrayGuard &operator=(const PyArrayGuard &) = delete;
  ~PyArrayGuard() {
    if (d_arr != nullptr) {
      Py_DECREF(reinterpret_cast<PyObject *>(d_arr));
    }
  }
  PyArrayObject *get() const { return d_arr; }
  PyArrayObject *release() {
    PyArrayObject *out = d_arr;
    d_arr = nullptr;
    return out;
  }

 private:
  PyArrayObject *d_arr{nullptr};
};

std::vector<const ExplicitBitVect *> extractFps(python::object seq,
                                                const char *which) {
  const Py_ssize_t length = python::len(seq);
  std::vector<const ExplicitBitVect *> out;
  out.reserve(static_cast<std::size_t>(length));
  for (Py_ssize_t i = 0; i < length; ++i) {
    python::extract<const ExplicitBitVect *> ex(seq[i]);
    if (!ex.check()) {
      const std::string msg = std::string("expected ExplicitBitVect in `") +
                              which + "` at index " + std::to_string(i);
      throw_value_error(msg.c_str());
    }
    out.push_back(ex());
  }
  return out;
}

PyObject *tanimotoMatrixPy(python::object probes, python::object targets) {
  const auto probeFps = extractFps(probes, "probes");
  const auto targetFps = extractFps(targets, "targets");

  npy_intp dims[2];
  dims[0] = static_cast<npy_intp>(probeFps.size());
  dims[1] = static_cast<npy_intp>(targetFps.size());
  PyArrayGuard guard(
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNew(2, dims, NPY_DOUBLE)));
  if (guard.get() == nullptr) {
    throw_value_error("failed to allocate Tanimoto matrix");
  }
  if (probeFps.empty() || targetFps.empty()) {
    return reinterpret_cast<PyObject *>(guard.release());
  }

  std::size_t probeBits = 0;
  std::size_t targetBits = 0;
  auto packedProbes =
      RDKit::BulkSimilarity::packFingerprints(probeFps, probeBits);
  auto packedTargets =
      RDKit::BulkSimilarity::packFingerprints(targetFps, targetBits);
  if (probeBits != targetBits) {
    throw_value_error(
        "probe and target fingerprints must have the same bit length");
  }
  const std::size_t words = RDKit::BulkSimilarity::wordsForBits(probeBits);
  double *out = static_cast<double *>(PyArray_DATA(guard.get()));

  RDKit::BulkSimilarity::tanimotoMatrix(packedProbes.data(), probeFps.size(),
                                        packedTargets.data(), targetFps.size(),
                                        words, out);
  return reinterpret_cast<PyObject *>(guard.release());
}

}  // namespace

BOOST_PYTHON_MODULE(cBulkSimilarity) {
  rdkit_import_array();

  python::scope().attr("__doc__") =
      "Bulk Tanimoto similarity over many ExplicitBitVect fingerprints.\n"
      "\n"
      "BulkTanimotoMatrix(probes, targets) returns an (M, N) NumPy array of\n"
      "Tanimoto coefficients between every probe and every target. All\n"
      "fingerprints in a single call must share the same bit length, and\n"
      "that length must be a positive multiple of 64.\n"
      "\n"
      "The kernel parallelises across probe rows over the hardware thread\n"
      "count and uses an AVX-512 popcount path at runtime on CPUs that\n"
      "support VPOPCNTDQ (Intel Ice Lake and later, AMD Zen 4 and later).\n";

  python::def("BulkTanimotoMatrix", tanimotoMatrixPy,
              (python::arg("probes"), python::arg("targets")),
              "Compute an M x N Tanimoto similarity matrix between two\n"
              "iterables of ExplicitBitVect fingerprints. Returns a NumPy\n"
              "float64 ndarray of shape (len(probes), len(targets)).");

  python::def("BulkSimilarityActiveKernel",
              &RDKit::BulkSimilarity::activeKernel,
              "Name of the popcount kernel selected at startup: either\n"
              "'scalar' or 'avx512vpopcntdq'.");
}
