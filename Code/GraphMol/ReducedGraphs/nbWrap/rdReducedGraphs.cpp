//
//  Copyright (C) 2013-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/ReducedGraphs/ReducedGraphs.h>
#include <Numerics/Vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::ndarray<nb::numpy, double, nb::ndim<1>> dvToNumpyArray(
    RDNumeric::DoubleVector *dv) {
  size_t n = dv->size();
  double *data = new double[n];
  memcpy(data, dv->getData(), n * sizeof(double));
  delete dv;
  nb::capsule owner(data, [](void *p) noexcept {
    delete[] reinterpret_cast<double *>(p);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
}

}  // namespace

NB_MODULE(rdReducedGraphs, m) {
  m.doc() =
      "Module containing functions to generate and work with reduced graphs";

  m.def(
      "GenerateMolExtendedReducedGraph",
      [](const RDKit::ROMol &mol, nb::object atomTypes) -> RDKit::ROMol * {
        if (!atomTypes.is_none()) {
          throw std::invalid_argument(
              "specification of atom types not yet supported");
        }
        return RDKit::ReducedGraphs::generateMolExtendedReducedGraph(mol);
      },
      "mol"_a, "atomTypes"_a = nb::none(), nb::rv_policy::take_ownership,
      "Returns the reduced graph for a molecule");

  m.def(
      "GenerateErGFingerprintForReducedGraph",
      [](const RDKit::ROMol &mol, nb::object atomTypes, double fuzzIncrement,
         int minPath, int maxPath) {
        if (!atomTypes.is_none()) {
          throw std::invalid_argument(
              "specification of atom types not yet supported");
        }
        auto *dv =
            RDKit::ReducedGraphs::generateErGFingerprintForReducedGraph(
                mol, nullptr, fuzzIncrement, minPath, maxPath);
        return dvToNumpyArray(dv);
      },
      "mol"_a, "atomTypes"_a = nb::none(), "fuzzIncrement"_a = 0.3,
      "minPath"_a = 1, "maxPath"_a = 15,
      "Returns the ErG fingerprint vector for a reduced graph");

  m.def(
      "GetErGFingerprint",
      [](const RDKit::ROMol &mol, nb::object atomTypes, double fuzzIncrement,
         int minPath, int maxPath) {
        if (!atomTypes.is_none()) {
          throw std::invalid_argument(
              "specification of atom types not yet supported");
        }
        auto *dv = RDKit::ReducedGraphs::getErGFingerprint(
            mol, nullptr, fuzzIncrement, minPath, maxPath);
        return dvToNumpyArray(dv);
      },
      "mol"_a, "atomTypes"_a = nb::none(), "fuzzIncrement"_a = 0.3,
      "minPath"_a = 1, "maxPath"_a = 15,
      "Returns the ErG fingerprint vector for a molecule");
}
