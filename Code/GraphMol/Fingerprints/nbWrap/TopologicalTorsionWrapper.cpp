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
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace TopologicalTorsionWrapper {

template <typename OutputType>
FingerprintGenerator<OutputType> *getTopologicalTorsionFPGenerator(
    const bool includeChirality, const uint32_t torsionAtomCount,
    const bool countSimulation, nb::object py_countBounds,
    const std::uint32_t fpSize, nb::object py_atomInvGen) {
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

  return TopologicalTorsion::getTopologicalTorsionGenerator<OutputType>(
      includeChirality, torsionAtomCount, atomInvariantsGenerator,
      countSimulation, fpSize, countBounds, false);
}

void exportTopologicalTorsion(nb::module_ &m) {
  // Topological torsion fingerprint does not support 32 bit output yet
  nb::class_<TopologicalTorsion::TopologicalTorsionArguments,
             FingerprintArguments>(m, "TopologicalTorsionFingerprintOptions")
      .def_rw(
          "torsionAtomCount",
          &TopologicalTorsion::TopologicalTorsionArguments::d_torsionAtomCount,
          "number of atoms to be included in the paths")
      .def_rw(
          "onlyShortestPaths",
          &TopologicalTorsion::TopologicalTorsionArguments::df_onlyShortestPaths,
          R"DOC(whether or not to only include paths which are the shortest path between the start and end atoms)DOC");

  m.def(
      "GetTopologicalTorsionGenerator",
      &getTopologicalTorsionFPGenerator<std::uint64_t>,
      "includeChirality"_a = false, "torsionAtomCount"_a = 4,
      "countSimulation"_a = true, "countBounds"_a = nb::none(),
      "fpSize"_a = 2048, "atomInvariantsGenerator"_a = nb::none(),
      R"DOC(Get an atom pair fingerprint generator

ARGUMENTS:
    - includeChirality: includeChirality argument for both the default
      atom invariants generator and the fingerprint arguments
    - torsionAtomCount: the number of atoms to include in the "torsions"
    - countSimulation: if set, use count simulation while generating the fingerprint
    - countBounds: boundaries for count simulation, corresponding bit
      will be set if the count is higher than the number provided for that spot
    - fpSize: size of the generated fingerprint, does not affect the sparse versions
    - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

This generator supports the following AdditionalOutput types:
    - atomToBits: which bits each atom is involved in
    - atomCounts: how many bits each atom sets
    - bitPaths: map from bitId to vectors of atom indices

RETURNS: FingerprintGenerator
)DOC",
      nb::rv_policy::take_ownership);
}
}  // namespace TopologicalTorsionWrapper

}  // namespace RDKit
