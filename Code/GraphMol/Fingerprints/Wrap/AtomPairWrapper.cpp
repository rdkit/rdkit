//
//  Copyright (C) 2018-2021 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <RDBoost/Wrap.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace python = boost::python;

namespace RDKit {
namespace AtomPairWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getAtomPairGenerator(
    unsigned int minDistance, unsigned int maxDistance, bool includeChirality,
    bool use2D, bool countSimulation, python::object &py_countBounds,
    std::uint32_t fpSize, python::object &py_atomInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check() && atomInvGen()) {
    atomInvariantsGenerator = atomInvGen()->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  if (py_countBounds) {
    auto tmp = pythonObjectToVect<std::uint32_t>(py_countBounds);
    countBounds = *tmp;
  }

  return AtomPair::getAtomPairGenerator<OutputType>(
      minDistance, maxDistance, includeChirality, use2D,
      atomInvariantsGenerator, countSimulation, fpSize, countBounds, true);
}

AtomInvariantsGenerator *getAtomPairAtomInvGen(const bool includeChirality) {
  return new AtomPairAtomInvGenerator(includeChirality);
}

void exportAtompair() {
  python::def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint64_t>,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("countSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      "Get an atom pair fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - minDistance: minimum distance between atoms to be considered in a "
      "pair, default is 1 bond\n"
      "    - maxDistance: maximum distance between atoms to be considered in a "
      "pair, default is maxPathLen-1 bonds\n"
      "    - includeChirality: if set, chirality will be used in the atom  "
      "invariants, this is ignored if atomInvariantsGenerator is provided\n"
      "    - use2D: if set, the 2D (topological) distance matrix  will be "
      "used\n"
      "    - countSimulation:  if set, use count simulation while  "
      "generating the fingerprint\n"
      "    - countBounds: boundaries for count simulation, corresponding bit "
      "will be  set if the count is higher than the number provided for that "
      "spot\n"
      "    - fpSize: size of the generated fingerprint, does not affect the "
      "sparse versions\n"
      "    - atomInvariantsGenerator: atom invariants to be used during "
      "fingerprint generation\n\n"
      "This generator supports the following AdditionalOutput types:\n"
      "    - atomToBits: which bits each atom is involved in\n"
      "    - atomCounts: how many bits each atom sets\n"
      "    - bitInfoMap: map from bitId to (atomId, radius) pairs\n\n"
      "  RETURNS: FingerprintGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  python::def("GetAtomPairAtomInvGen", &getAtomPairAtomInvGen,
              (python::arg("includeChirality") = false),
              "Get an atom pair atom-invariant generator\n\n"
              "  ARGUMENTS:\n"
              "    - includeChirality: if set, chirality will be taken into "
              "account for invariants\n"
              "  RETURNS: AtomInvariantsGenerator\n\n",
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
