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
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>

using namespace RDKit;
namespace python = boost::python;

namespace RDKit {
namespace RDKitFPWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    const unsigned int minPath, const unsigned int maxPath, const bool useHs,
    const bool branchedPaths, const bool useBondOrder,
    const bool countSimulation, python::object &py_countBounds,
    const std::uint32_t fpSize, python::object &py_atomInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check() && atomInvGen()) {
    atomInvariantsGenerator = atomInvGen()->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  python::extract<std::vector<std::uint32_t>> countBoundsE(py_countBounds);
  if (countBoundsE.check() && !countBoundsE().empty()) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return RDKitFP::getRDKitFPGenerator<OutputType>(
      minPath, maxPath, useHs, branchedPaths, useBondOrder,
      atomInvariantsGenerator, countSimulation, countBoundsC, fpSize, true);
}

AtomInvariantsGenerator *getRDKitAtomInvGen() {
  return new RDKitFP::RDKitFPAtomInvGenerator();
}

void exportRDKit() {
  /*python::def("GetRDKitFPGenerator32", &getRDKitFPGenerator<std::uint32_t>,
              (python::arg("minPath") = 1, python::arg("maxPath") = 7,
               python::arg("useHs") = true, python::arg("branchedPaths") = true,
               python::arg("useBondOrder") = true,
               python::arg("countSimulation") = true,
               python::arg("countBounds") = python::object(),
               python::arg("fpSize") = 2048,
               python::arg("atomInvariantsGenerator") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());*/

  python::def(
      "GetRDKitFPGenerator", &getRDKitFPGenerator<std::uint64_t>,
      (python::arg("minPath") = 1, python::arg("maxPath") = 7,
       python::arg("useHs") = true, python::arg("branchedPaths") = true,
       python::arg("useBondOrder") = true,
       python::arg("countSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048,
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

  python::def("GetRDKitAtomInvGen", &getRDKitAtomInvGen,
              "Get an RDKit atom invariants generator\n\n"
              "  RETURNS: AtomInvariantsGenerator\n\n",
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace RDKitFPWrapper

}  // namespace RDKit
