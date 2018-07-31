
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
    const std::uint32_t foldedSize, python::object &py_atomInvGen,
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

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  python::extract<std::vector<std::uint32_t>> countBoundsE(py_countBounds);
  if (countBoundsE.check()) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return TopologicalTorsion::getTopologicalTorsionGenerator<OutputType>(
      includeChirality, torsionAtomCount, atomInvariantsGenerator,
      bondInvariantsGenerator, countSimulation, countBoundsC, foldedSize);
}

void exportTopologicalTorsion() {
  std::string docString = "";

  // Topological torsion fingerprint does not support 32 bit output yet
  /*python::def("GetTopologicalTorsionGenerator32",
              &getTopologicalTorsionFPGenerator<std::uint32_t>,
              (python::arg("includeChirality") = false,
               python::arg("torsionAtomCount") = 4,
               python::arg("countSimulation") = true,
               python::arg("countBounds") = python::object(),
               python::arg("foldedSize") = 2048,
               python::arg("atomInvariantsGenerator") = python::object(),
               python::arg("bondInvariantsGenerator") = python::object()),
              docString.c_str(),
              python::return_value_policy<
                  python::manage_new_object,
                  python::return_internal_reference<
                      6, python::return_internal_reference<7>>>());*/

  docString = "";

  python::def("GetTopologicalTorsionGenerator64",
              &getTopologicalTorsionFPGenerator<std::uint64_t>,
              (python::arg("includeChirality") = false,
               python::arg("torsionAtomCount") = 4,
               python::arg("countSimulation") = true,
               python::arg("countBounds") = python::object(),
               python::arg("foldedSize") = 2048,
               python::arg("atomInvariantsGenerator") = python::object(),
               python::arg("bondInvariantsGenerator") = python::object()),
              docString.c_str(),
              python::return_value_policy<
                  python::manage_new_object,
                  python::return_internal_reference<
                      6, python::return_internal_reference<7>>>());

  return;
}
}  // namespace TopologicalTorsionWrapper

}  // namespace RDKit
