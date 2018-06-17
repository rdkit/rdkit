
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
using namespace RDKit::FingerprintWrapper;
namespace python = boost::python;

namespace RDKit {
namespace AtomPairWrapper {
FingerprintGeneratorWrapper *getAtomPairGeneratorWrapped(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation,
    AtomInvGeneratorWrapper &atomInvGeneratorWrapper,
    BondInvGeneratorWrapper &bondInvGeneratorWrapper) {
  AtomInvariantsGenerator *atomInvariantsGenerator =
      atomInvGeneratorWrapper.dp_atomInvariantsGenerator;
  BondInvariantsGenerator *bondInvariantsGenerator =
      bondInvGeneratorWrapper.dp_bondInvariantsGenerator;
  AtomEnvironmentGenerator *atomPairEnvGenerator =
      new AtomPair::AtomPairEnvGenerator();
  FingerprintArguments *atomPairArguments = new AtomPair::AtomPairArguments(
      useCountSimulation, includeChirality, use2D, minDistance, maxDistance);

  bool ownsAtomInvGenerator = false;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new AtomPairAtomInvGenerator(includeChirality);
    ownsAtomInvGenerator = true;
  }

  FingerprintGenerator *fingerprintGenerator = new FingerprintGenerator(
      atomPairEnvGenerator, atomPairArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, false);

  FingerprintGeneratorWrapper *wrapped = new FingerprintGeneratorWrapper();
  wrapped->dp_fingerprintGenerator = fingerprintGenerator;
  return wrapped;
}

AtomInvGeneratorWrapper *getAtomPairAtomInvGen(const bool includeChirality) {
  AtomInvGeneratorWrapper *atomInvGeneratorWrapper =
      new AtomInvGeneratorWrapper();

  atomInvGeneratorWrapper->dp_atomInvariantsGenerator =
      new AtomPairAtomInvGenerator(includeChirality);

  return atomInvGeneratorWrapper;
}

void exportAtompair() {
  std::string docString = "";
  python::def(
      "getAtomPairGeneratorWrapped", &getAtomPairGeneratorWrapped,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("atomInvGeneratorWrapper") = AtomInvGeneratorWrapper(),
       python::arg("BondInvGeneratorWrapper") = BondInvGeneratorWrapper()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString = "";
  python::def("getAtomPairAtomInvGen", &getAtomPairAtomInvGen,
              (python::arg("includeChirality") = false), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
