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
#include <GraphMol/Fingerprints/MorganGenerator.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace MorganWrapper {

template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    unsigned int radius, bool countSimulation, bool includeChirality,
    bool useBondTypes, bool onlyNonzeroInvariants,
    bool,  // includeRingMembership
    nb::object py_countBounds, std::uint32_t fpSize, nb::object py_atomInvGen,
    nb::object py_bondInvGen, bool includeRedundantEnvironments) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  BondInvariantsGenerator *bondInvariantsGenerator = nullptr;

  if (!py_atomInvGen.is_none()) {
    atomInvariantsGenerator =
        nb::cast<AtomInvariantsGenerator *>(py_atomInvGen)->clone();
  }

  if (!py_bondInvGen.is_none()) {
    bondInvariantsGenerator =
        nb::cast<BondInvariantsGenerator *>(py_bondInvGen)->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  if (!py_countBounds.is_none()) {
    countBounds.clear();
    for (auto item : py_countBounds) {
      countBounds.push_back(nb::cast<std::uint32_t>(item));
    }
  }

  return MorganFingerprint::getMorganGenerator<OutputType>(
      radius, countSimulation, includeChirality, useBondTypes,
      onlyNonzeroInvariants, includeRedundantEnvironments,
      atomInvariantsGenerator, bondInvariantsGenerator, fpSize, countBounds,
      true, true);
}

void exportMorgan(nb::module_ &m) {
  nb::class_<MorganFingerprint::MorganArguments, FingerprintArguments>(
      m, "MorganFingerprintOptions")
      .def_rw(
          "onlyNonzeroInvariants",
          &MorganFingerprint::MorganArguments::df_onlyNonzeroInvariants,
          "use include atoms which have nonzero invariants")
      .def_rw("radius", &MorganFingerprint::MorganArguments::d_radius,
              "the radius of the fingerprints to generate")
      .def_rw(
          "includeRedundantEnvironments",
          &MorganFingerprint::MorganArguments::df_includeRedundantEnvironments,
          "include redundant environments in the fingerprint");

  m.def(
      "GetMorganGenerator", getMorganGenerator<std::uint64_t>,
      "radius"_a = 3, "countSimulation"_a = false,
      "includeChirality"_a = false, "useBondTypes"_a = true,
      "onlyNonzeroInvariants"_a = false, "includeRingMembership"_a = true,
      "countBounds"_a = nb::none(), "fpSize"_a = 2048,
      "atomInvariantsGenerator"_a = nb::none(),
      "bondInvariantsGenerator"_a = nb::none(),
      "includeRedundantEnvironments"_a = false,
      R"DOC(Get a morgan fingerprint generator

ARGUMENTS:
    - radius: the number of iterations to grow the fingerprint
    - countSimulation: if set, use count simulation while generating the fingerprint
    - includeChirality: if set, chirality information will be added to
      the generated fingerprint
    - useBondTypes: if set, bond types will be included as a part of
      the default bond invariants
    - countBounds: boundaries for count simulation, corresponding bit
      will be set if the count is higher than the number provided for that spot
    - fpSize: size of the generated fingerprint, does not affect the sparse versions
    - atomInvariantsGenerator: atom invariants to be used during fingerprint generation

This generator supports the following AdditionalOutput types:
    - atomToBits: which bits each atom is the center of
    - atomCounts: how many bits each atom sets
    - bitInfoMap: map from bitId to (atomId1, radius) pairs

RETURNS: FingerprintGenerator
)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetMorganAtomInvGen",
      [](const bool includeRingMembership) -> AtomInvariantsGenerator * {
        return new MorganFingerprint::MorganAtomInvGenerator(
            includeRingMembership);
      },
      "includeRingMembership"_a,
      R"DOC(Get a morgan atom invariants generator

ARGUMENTS:
    - includeRingMembership: if set, whether or not the atom is in a ring
      will be used in the invariant list

RETURNS: AtomInvariantsGenerator
)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetMorganFeatureAtomInvGen",
      [](nb::object py_patterns) -> AtomInvariantsGenerator * {
        if (py_patterns.is_none()) {
          return new MorganFingerprint::MorganFeatureAtomInvGenerator(nullptr);
        }
        std::vector<const ROMol *> patterns;
        for (auto item : py_patterns) {
          patterns.push_back(nb::cast<const ROMol *>(item));
        }
        return new MorganFingerprint::MorganFeatureAtomInvGenerator(&patterns);
      },
      "patterns"_a = nb::none(),
      R"DOC(Get a morgan feature atom invariants generator

ARGUMENTS:
    - patterns: if provided should contain the queries used to assign
      atom-types. if not provided, feature definitions adapted from
      reference: Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will
      be used for Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.

RETURNS: AtomInvariantsGenerator
)DOC",
      nb::rv_policy::take_ownership);

  m.def(
      "GetMorganBondInvGen",
      [](const bool useBondTypes,
         const bool useChirality) -> BondInvariantsGenerator * {
        return new MorganFingerprint::MorganBondInvGenerator(useBondTypes,
                                                             useChirality);
      },
      "useBondTypes"_a = true,
      "useChirality"_a = false,
      R"DOC(Get a morgan bond invariants generator

ARGUMENTS:
    - useBondTypes: if set, bond types will be included as a part of the bond invariants
    - useChirality: if set, chirality information will be included as a
      part of the bond invariants

RETURNS: BondInvariantsGenerator
)DOC",
      nb::rv_policy::take_ownership);
}
}  // namespace MorganWrapper

}  // namespace RDKit
