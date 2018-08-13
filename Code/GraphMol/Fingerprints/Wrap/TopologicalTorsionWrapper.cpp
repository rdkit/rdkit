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
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>

using namespace RDKit;
namespace python = boost::python;

namespace RDKit {
namespace TopologicalTorsionWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getTopologicalTorsionFPGenerator(
    const bool includeChirality, const uint32_t torsionAtomCount,
    const bool countSimulation, python::object &py_countBounds,
    const std::uint32_t fpSize, python::object &py_atomInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  BondInvariantsGenerator *bondInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check() && atomInvGen()) {
    atomInvariantsGenerator = atomInvGen();
    atomInvariantsGenerator = atomInvariantsGenerator->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  python::extract<std::vector<std::uint32_t>> countBoundsE(py_countBounds);
  if (countBoundsE.check() && !countBoundsE().empty()) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return TopologicalTorsion::getTopologicalTorsionGenerator<OutputType>(
      includeChirality, torsionAtomCount, atomInvariantsGenerator,
      countSimulation, countBoundsC, fpSize, false);
}

void exportTopologicalTorsion() {
  // Topological torsion fingerprint does not support 32 bit output yet

  python::def(
      "GetTopologicalTorsionGenerator",
      &getTopologicalTorsionFPGenerator<std::uint64_t>,
      (python::arg("includeChirality") = false,
       python::arg("torsionAtomCount") = 4,
       python::arg("countSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("fpSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      "Get an atom pair fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - includeChirality: includeChirality argument for both the default "
      "atom invariants generator and the fingerprint arguments\n"
      "    - torsionAtomCount: the number of atoms to include in the "
      "\"torsions\"\n"
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

  return;
}
}  // namespace TopologicalTorsionWrapper

}  // namespace RDKit
