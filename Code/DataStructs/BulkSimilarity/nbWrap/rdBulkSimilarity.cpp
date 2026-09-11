//
//  Copyright (C) 2026 Clay Moore and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  nanobind wrapper for the BulkSimilarity library. This mirrors the
//  Boost.Python extension in ../Wrap/ so that both binding backends expose
//  the same rdkit.DataStructs entry points; the two are never built at the
//  same time (RDK_BUILD_BOOST_PYTHON_WRAPPERS and RDK_BUILD_NANOBIND_WRAPPERS
//  are mutually exclusive) and they share the test file in ../Wrap/.

#include <DataStructs/BitVects.h>
#include <DataStructs/BulkSimilarity/BulkSimilarity.h>
#include <DataStructs/ExplicitBitVect.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

//! The extracted pointers reference C++ objects owned by Python. Because
//! the input is an nb::iterable it may be a generator, whose items are not
//! kept alive by any container, so we hold an owning nb::object for each
//! element for as long as the pointers are in use.
struct FingerprintList {
  std::vector<nb::object> keepAlive;
  std::vector<const ExplicitBitVect *> ptrs;
};

FingerprintList extractFps(const nb::iterable &seq, const char *which) {
  FingerprintList out;
  std::size_t idx = 0;
  for (auto item : seq) {
    nb::object owned = nb::borrow<nb::object>(item);
    try {
      out.ptrs.push_back(&nb::cast<const ExplicitBitVect &>(owned));
    } catch (const nb::cast_error &) {
      const std::string msg = std::string("expected ExplicitBitVect in `") +
                              which + "` at index " + std::to_string(idx);
      throw nb::value_error(msg.c_str());
    }
    out.keepAlive.push_back(std::move(owned));
    ++idx;
  }
  return out;
}

nb::ndarray<nb::numpy, double, nb::ndim<2>> tanimotoMatrixPy(
    const nb::iterable &probes, const nb::iterable &targets) {
  const FingerprintList probeFps = extractFps(probes, "probes");
  const FingerprintList targetFps = extractFps(targets, "targets");

  const std::size_t numProbes = probeFps.ptrs.size();
  const std::size_t numTargets = targetFps.ptrs.size();

  // Held by unique_ptr until the very end so that an exception from
  // packFingerprints (e.g. mismatched bit lengths) cannot leak the buffer.
  std::unique_ptr<double[]> buf(new double[numProbes * numTargets]());

  if (numProbes != 0 && numTargets != 0) {
    std::size_t probeBits = 0;
    std::size_t targetBits = 0;
    const auto packedProbes =
        RDKit::BulkSimilarity::packFingerprints(probeFps.ptrs, probeBits);
    const auto packedTargets =
        RDKit::BulkSimilarity::packFingerprints(targetFps.ptrs, targetBits);
    if (probeBits != targetBits) {
      throw nb::value_error(
          "probe and target fingerprints must have the same bit length");
    }
    const std::size_t words = RDKit::BulkSimilarity::wordsForBits(probeBits);
    RDKit::BulkSimilarity::tanimotoMatrix(packedProbes.data(), numProbes,
                                          packedTargets.data(), numTargets,
                                          words, buf.get());
  }

  // Ownership transfers to NumPy via the capsule.
  // Pattern from
  // https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
  double *raw = buf.release();
  nb::capsule owner(raw, [](void *p) noexcept {
    delete[] reinterpret_cast<double *>(p);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(
      /* data  = */ raw,
      /* shape = */ {numProbes, numTargets},
      /* owner = */ owner);
}

}  // namespace

NB_MODULE(cBulkSimilarity, m) {
  m.doc() =
      "Bulk Tanimoto similarity over many ExplicitBitVect fingerprints.\n"
      "\n"
      "BulkTanimotoMatrix(probes, targets) returns an (M, N) NumPy array of\n"
      "Tanimoto coefficients between every probe and every target. All\n"
      "fingerprints in a single call must share the same bit length, and\n"
      "that length must be a positive multiple of 64.\n"
      "\n"
      "The kernel parallelises across probe rows over the hardware thread\n"
      "count and uses an AVX-512 popcount path at runtime on CPUs that\n"
      "support VPOPCNTDQ (Intel Ice Lake and later, AMD Zen 4 and later).";

  m.def("BulkTanimotoMatrix", &tanimotoMatrixPy, "probes"_a, "targets"_a,
        R"DOC(Compute an M x N Tanimoto similarity matrix between two
iterables of ExplicitBitVect fingerprints. Returns a NumPy float64
ndarray of shape (len(probes), len(targets)).)DOC");

  m.def("BulkSimilarityActiveKernel", &RDKit::BulkSimilarity::activeKernel,
        R"DOC(Name of the popcount kernel selected at startup: either
'scalar' or 'avx512vpopcntdq'.)DOC");
}
