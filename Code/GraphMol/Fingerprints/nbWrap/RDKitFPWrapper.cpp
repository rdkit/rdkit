//
//  Copyright (C) 2018-2026 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace RDKitFPWrapper {

template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    unsigned int minPath, unsigned int maxPath, bool useHs, bool branchedPaths,
    bool useBondOrder, bool countSimulation, nb::object py_countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
    nb::object py_atomInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;

  if (!py_atomInvGen.is_none()) {
    atomInvariantsGenerator =
        nb::cast<AtomInvariantsGenerator *>(py_atomInvGen)->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  if (!py_countBounds.is_none()) {
    countBounds.clear();
    for (auto item : py_countBounds) {
      countBounds.push_back(nb::cast<std::uint32_t>(item));
    }
  }

  return RDKitFP::getRDKitFPGenerator<OutputType>(
      minPath, maxPath, useHs, branchedPaths, useBondOrder,
      atomInvariantsGenerator, countSimulation, countBounds, fpSize,
      numBitsPerFeature, true);
}

void exportRDKit(nb::module_ &m) {
  nb::class_<RDKitFP::RDKitFPArguments, FingerprintArguments>(
      m, "RDKitFingerprintOptions")
      .def_rw("minPath", &RDKitFP::RDKitFPArguments::d_minPath,
              "minimum path length (in bonds) to be included")
      .def_rw("maxPath", &RDKitFP::RDKitFPArguments::d_maxPath,
              "maximum path length (in bonds) to be included")
      .def_rw("useHs", &RDKitFP::RDKitFPArguments::df_useHs,
              "use explicit Hs in the paths (if molecule has explicit Hs)")
      .def_rw("branchedPaths", &RDKitFP::RDKitFPArguments::df_branchedPaths,
              "generate branched subgraphs, not just linear ones")
      .def_rw("useBondOrder", &RDKitFP::RDKitFPArguments::df_useBondOrder,
              "include bond orders in the path hashes");

  m.def(
      "GetRDKitFPGenerator", &getRDKitFPGenerator<std::uint64_t>,
      "minPath"_a = 1, "maxPath"_a = 7, "useHs"_a = true,
      "branchedPaths"_a = true, "useBondOrder"_a = true,
      "countSimulation"_a = false, "countBounds"_a = nb::none(),
      "fpSize"_a = 2048, "numBitsPerFeature"_a = 2,
      "atomInvariantsGenerator"_a = nb::none(),
      R"DOC(Get an RDKit fingerprint generator

ARGUMENTS:
    - minPath: the minimum path length (in bonds) to be included
    - maxPath: the maximum path length (in bonds) to be included
    - useHs: toggles inclusion of Hs in paths (if the molecule has explicit Hs)
    - branchedPaths: toggles generation of branched subgraphs, not just linear paths
    - useBondOrder: toggles inclusion of bond orders in the path hashes
    - countSimulation: if set, use count simulation while generating the fingerprint
    - countBounds: boundaries for count simulation, corresponding bit
      will be set if the count is higher than the number provided for that spot
    - fpSize: size of the generated fingerprint, does not affect the sparse versions
    - numBitsPerFeature: the number of bits set per path/subgraph found
    - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

This generator supports the following AdditionalOutput types:
    - atomToBits: which bits each atom is involved in
    - atomCounts: how many bits each atom sets
    - bitPaths: map from bitId to vectors of bond indices for the individual subgraphs

RETURNS: FingerprintGenerator
)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetRDKitAtomInvGen",
      []() -> AtomInvariantsGenerator * {
        return new RDKitFP::RDKitFPAtomInvGenerator();
      },
      R"DOC(Get an RDKit atom invariants generator

RETURNS: AtomInvariantsGenerator
)DOC",
      nb::rv_policy::take_ownership);
}
}  // namespace RDKitFPWrapper

}  // namespace RDKit
