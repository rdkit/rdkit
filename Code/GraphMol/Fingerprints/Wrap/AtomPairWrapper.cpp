
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace python = boost::python;

namespace RDKit {
namespace AtomPairWrapper {
FingerprintGenerator<std::uint32_t> *getAtomPairGeneratorWrapped(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation, python::object &py_atomInvGen,
    python::object &py_bondInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  BondInvariantsGenerator *bondInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check()) {
    atomInvariantsGenerator = atomInvGen();
  }

  python::extract<BondInvariantsGenerator *> bondInvGen(py_bondInvGen);
  if (bondInvGen.check()) {
    bondInvariantsGenerator = bondInvGen();
  }

  return AtomPair::getAtomPairGenerator(
      minDistance, maxDistance, includeChirality, use2D, useCountSimulation,
      atomInvariantsGenerator, bondInvariantsGenerator);
}

AtomInvariantsGenerator *getAtomPairAtomInvGen(const bool includeChirality) {
  return new AtomPairAtomInvGenerator(includeChirality);
}

void exportAtompair() {
  std::string docString = "";
  python::def(
      "GetAtomPairGenerator", &getAtomPairGeneratorWrapped,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("atomInvariantsGenerator") = python::object(),
       python::arg("bondInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<
          python::manage_new_object,
          python::return_internal_reference<
              6, python::return_internal_reference<7>>>());

  docString = "";
  python::def("GetAtomPairAtomInvGen", &getAtomPairAtomInvGen,
              (python::arg("includeChirality") = false), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
