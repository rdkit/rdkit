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
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <RDBoost/Wrap.h>

using namespace RDKit;
namespace python = boost::python;

namespace RDKit {
namespace RDKitFPWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    unsigned int minPath, unsigned int maxPath, bool useHs, bool branchedPaths,
    bool useBondOrder, bool countSimulation, python::object &py_countBounds,
    std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
    python::object &py_atomInvGen) {
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

  return RDKitFP::getRDKitFPGenerator<OutputType>(
      minPath, maxPath, useHs, branchedPaths, useBondOrder,
      atomInvariantsGenerator, countSimulation, countBounds, fpSize,
      numBitsPerFeature, true);
}

AtomInvariantsGenerator *getRDKitAtomInvGen() {
  return new RDKitFP::RDKitFPAtomInvGenerator();
}

void exportRDKit() {
  python::def(
      "GetRDKitFPGenerator", &getRDKitFPGenerator<std::uint64_t>,
      (python::arg("minPath") = 1, python::arg("maxPath") = 7,
       python::arg("useHs") = true, python::arg("branchedPaths") = true,
       python::arg("useBondOrder") = true,
       python::arg("countSimulation") = false,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048, python::arg("numBitsPerFeature") = 2,
       python::arg("atomInvariantsGenerator") = python::object()),
      "Get an RDKit fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - minPath: the minimum path length (in bonds) to be included\n"
      "    - maxPath: the maximum path length (in bonds) to be included\n"
      "    - useHs: toggles inclusion of Hs in paths (if the molecule has "
      "explicit Hs)\n"
      "    - branchedPaths: toggles generation of branched subgraphs, not just "
      "linear paths\n"
      "    - useBondOrder: toggles inclusion of bond orders in the path "
      "hashes\n"
      "    - countSimulation:  if set, use count simulation while  "
      "generating the fingerprint\n"
      "    - countBounds: boundaries for count simulation, corresponding bit "
      "will be  set if the count is higher than the number provided for that "
      "spot\n"
      "    - fpSize: size of the generated fingerprint, does not affect the "
      "sparse versions\n"
      "    - numBitsPerFeature: the number of bits set per path/subgraph "
      "found\n"
      "    - atomInvariantsGenerator: atom invariants to be used during "
      "fingerprint generation\n\n"
      "This generator supports the following AdditionalOutput types:\n"
      "    - atomToBits: which bits each atom is involved in\n"
      "    - atomCounts: how many bits each atom sets\n"
      "    - bitPaths: map from bitId to vectors of bond indices for the "
      "individual subgraphs\n\n"
      "  RETURNS: FingerprintGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  python::def("GetRDKitAtomInvGen", &getRDKitAtomInvGen,
              "Get an RDKit atom invariants generator\n\n"
              "  RETURNS: AtomInvariantsGenerator\n\n",
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace RDKitFPWrapper

}  // namespace RDKit
