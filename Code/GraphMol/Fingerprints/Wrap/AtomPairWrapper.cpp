
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace python = boost::python;

namespace RDKit {
namespace AtomPairWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getAtomPairGenerator(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation, python::object &py_countBounds,
    const std::uint32_t foldedSize, python::object &py_atomInvGen) {
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

  return AtomPair::getAtomPairGenerator<OutputType>(
      minDistance, maxDistance, includeChirality, use2D,
      atomInvariantsGenerator, useCountSimulation, foldedSize, countBoundsC,
      true);
}

AtomInvariantsGenerator *getAtomPairAtomInvGen(const bool includeChirality) {
  return new AtomPairAtomInvGenerator(includeChirality);
}

void exportAtompair() {
  std::string docString = "";
  /*python::def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint32_t>,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("foldedSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());*/

  docString = "";
  python::def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint64_t>,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("foldedSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString = "";
  python::def("GetAtomPairAtomInvGen", &getAtomPairAtomInvGen,
              (python::arg("includeChirality") = false), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
