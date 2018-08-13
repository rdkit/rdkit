//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

using namespace RDKit;
namespace python = boost::python;

namespace RDKit {
namespace MorganWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants, const bool includeRingMembership,
    python::object &py_countBounds, const std::uint32_t fpSize,
    python::object &py_atomInvGen, python::object &py_bondInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  BondInvariantsGenerator *bondInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check() && atomInvGen()) {
    atomInvariantsGenerator = atomInvGen()->clone();
  }

  python::extract<BondInvariantsGenerator *> bondInvGen(py_bondInvGen);
  if (bondInvGen.check() && bondInvGen()) {
    bondInvariantsGenerator = bondInvGen()->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  python::extract<std::vector<std::uint32_t>> countBoundsE(py_countBounds);
  if (countBoundsE.check() && !countBoundsE().empty()) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return MorganFingerprint::getMorganGenerator<OutputType>(
      radius, countSimulation, includeChirality, useBondTypes,
      onlyNonzeroInvariants, atomInvariantsGenerator, bondInvariantsGenerator,
      fpSize, countBounds, true, true);
}

AtomInvariantsGenerator *getMorganAtomInvGen(const bool includeRingMembership) {
  return new MorganFingerprint::MorganAtomInvGenerator(includeRingMembership);
}

AtomInvariantsGenerator *getMorganFeatureAtomInvGen(
    python::object &py_patterns) {
  std::vector<const ROMol *> patterns;
  python::extract<std::vector<const ROMol *>> patternsE(py_patterns);
  if (patternsE.check()) {
    patterns = patternsE();
    return new MorganFingerprint::MorganFeatureAtomInvGenerator(&patterns);
  } else {
    return new MorganFingerprint::MorganFeatureAtomInvGenerator(nullptr);
  }
}

BondInvariantsGenerator *getMorganBondInvGen(const bool useBondTypes,
                                             const bool useChirality) {
  return new MorganFingerprint::MorganBondInvGenerator(useBondTypes,
                                                       useChirality);
}

void exportMorgan() {
  /*python::def(
      "GetMorganGenerator32", getMorganGenerator<std::uint32_t>,
      (python::arg("radius") = 3, python::arg("useCountSimulation") = true,
       python::arg("includeChirality") = false,
       python::arg("useBondTypes") = true,
       python::arg("onlyNonzeroInvariants") = false,
       python::arg("includeRingMembership") = true,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object(),
       python::arg("bondInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());*/

  python::def(
      "GetMorganGenerator", getMorganGenerator<std::uint64_t>,
      (python::arg("radius") = 3, python::arg("useCountSimulation") = true,
       python::arg("includeChirality") = false,
       python::arg("useBondTypes") = true,
       python::arg("onlyNonzeroInvariants") = false,
       python::arg("includeRingMembership") = true,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object(),
       python::arg("bondInvariantsGenerator") = python::object()),
      "Get a morgan fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - radius:  the number of iterations to grow the fingerprint\n"
      "    - useCountSimulation: if set, use count simulation while generating "
      "the fingerprint\n"
      "    - includeChirality: if set, chirality information will be added to "
      "the generated fingerprint\n"
      "    - useBondTypes: if set, bond types will be included as a part of "
      "the default bond invariants\n"
      "    - useCountSimulation:  if set, use count simulation while  "
      "generating the fingerprint\n"
      "    - countBounds: boundaries for count simulation, corresponding bit "
      "will be  set if the count is higher than the number provided for that "
      "spot\n"
      "    - fpSize: size of the generated fingerprint, does not affect the "
      "sparse versions\n"
      "    - atomInvariantsGenerator: atom invariants to be used during "
      "fingerprint generation\n\n"
      "  RETURNS: FingerprintGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  python::def("GetMorganAtomInvGen", &getMorganAtomInvGen,
              (python::arg("includeRingMembership") = false),
              "Get a morgan atom invariants generator\n\n"
              "  ARGUMENTS:\n"
              "    - includeRingMembership: if set, whether or not the atom is "
              "in a ring will be used in the invariant list\n\n"
              "  RETURNS: AtomInvariantsGenerator\n\n",
              python::return_value_policy<python::manage_new_object>());

  python::def(
      "GetMorganFeatureAtomInvGen", &getMorganFeatureAtomInvGen,
      (python::arg("patterns") = python::object()),
      "Get a morgan feature atom invariants generator\n\n"
      "  ARGUMENTS:\n"
      "    - patterns: if provided should contain the queries used to assign "
      "atom-types. if not provided, feature definitions adapted from "
      "reference: Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will "
      "be used for Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.\n\n"
      "  RETURNS: AtomInvariantsGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  python::def(
      "GetMorganBondInvGen", &getMorganBondInvGen,
      (python::arg("useBondTypes") = true, python::arg("useChirality") = false),
      "Get a morgan bond invariants generator\n\n"
      "  ARGUMENTS:\n"
      "    - useBondTypes: if set, bond types will be included as a part of "
      "the bond invariants\n"
      "    - useChirality: if set, chirality information will be included as a "
      "part of the bond invariants\n\n"
      "  RETURNS: BondInvariantsGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace MorganWrapper

}  // namespace RDKit
