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
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace AtomPairWrapper {

template <typename OutputType>
FingerprintGenerator<OutputType> *getAtomPairGenerator(
    unsigned int minDistance, unsigned int maxDistance, bool includeChirality,
    bool use2D, bool countSimulation, nb::object py_countBounds,
    std::uint32_t fpSize, nb::object py_atomInvGen) {
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

  return AtomPair::getAtomPairGenerator<OutputType>(
      minDistance, maxDistance, includeChirality, use2D,
      atomInvariantsGenerator, countSimulation, fpSize, countBounds, true);
}

void exportAtompair(nb::module_ &m) {
  nb::class_<AtomPair::AtomPairArguments, FingerprintArguments>(
      m, "AtomPairFingerprintOptions")
      .def_rw("use2D", &AtomPair::AtomPairArguments::df_use2D,
              "use 2D distances")
      .def_rw("minDistance", &AtomPair::AtomPairArguments::d_minDistance,
              "minimum distance to be included")
      .def_rw("maxDistance", &AtomPair::AtomPairArguments::d_maxDistance,
              "maximum distance to be included");

  m.def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint64_t>,
      "minDistance"_a = 1,
      "maxDistance"_a = (unsigned int)(AtomPair::maxPathLen - 1),
      "includeChirality"_a = false, "use2D"_a = true,
      "countSimulation"_a = true, "countBounds"_a = nb::none(),
      "fpSize"_a = 2048, "atomInvariantsGenerator"_a = nb::none(),
      R"DOC(Get an atom pair fingerprint generator

ARGUMENTS:
    - minDistance: minimum distance between atoms to be considered in a
      pair, default is 1 bond
    - maxDistance: maximum distance between atoms to be considered in a
      pair, default is maxPathLen-1 bonds
    - includeChirality: if set, chirality will be used in the atom
      invariants, this is ignored if atomInvariantsGenerator is provided
    - use2D: if set, the 2D (topological) distance matrix will be used
    - countSimulation: if set, use count simulation while generating the fingerprint
    - countBounds: boundaries for count simulation, corresponding bit
      will be set if the count is higher than the number provided for that spot
    - fpSize: size of the generated fingerprint, does not affect the sparse versions
    - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

This generator supports the following AdditionalOutput types:
    - atomToBits: which bits each atom is involved in
    - atomCounts: how many bits each atom sets
    - bitInfoMap: map from bitId to (atomId, radius) pairs

RETURNS: FingerprintGenerator
)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetAtomPairAtomInvGen",
      [](const bool includeChirality) -> AtomInvariantsGenerator * {
        return new AtomPairAtomInvGenerator(includeChirality);
      },
      "includeChirality"_a = false,
      R"DOC(Get an atom pair atom-invariant generator

ARGUMENTS:
    - includeChirality: if set, chirality will be taken into account for invariants
RETURNS: AtomInvariantsGenerator
)DOC",
      nb::rv_policy::take_ownership);
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
