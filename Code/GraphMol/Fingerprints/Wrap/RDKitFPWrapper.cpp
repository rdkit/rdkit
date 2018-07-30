
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

  return RDKitFP::getRDKitFPGenerator<OutputType>(
      minPath, maxPath, useHs, branchedPaths, useBondOrder,
      atomInvariantsGenerator, bondInvariantsGenerator, countSimulation,
      countBoundsC, foldedSize);
}

AtomInvariantsGenerator *getRDKitAtomInvGen() {
  return new RDKitFP::RDKitFPAtomInvGenerator();
}

void exportRDKit() {
  std::string docString = "";

  python::def("GetRDKitFPGenerator32", &getRDKitFPGenerator<std::uint32_t>,
              (python::arg("minPath") = 1, python::arg("maxPath") = 7,
               python::arg("useHs") = true, python::arg("branchedPaths") = true,
               python::arg("useBondOrder") = true,
               python::arg("countSimulation") = true,
               python::arg("countBounds") = python::object(),
               python::arg("foldedSize") = 2048,
               python::arg("atomInvariantsGenerator") = python::object(),
               python::arg("bondInvariantsGenerator") = python::object()),
              docString.c_str(),
              python::return_value_policy<
                  python::manage_new_object,
                  python::return_internal_reference<
                      7, python::return_internal_reference<8>>>());

  python::def("GetRDKitFPGenerator64", &getRDKitFPGenerator<std::uint64_t>,
              (python::arg("minPath") = 1, python::arg("maxPath") = 7,
               python::arg("useHs") = true, python::arg("branchedPaths") = true,
               python::arg("useBondOrder") = true,
               python::arg("countSimulation") = true,
               python::arg("countBounds") = python::object(),
               python::arg("foldedSize") = 2048,
               python::arg("atomInvariantsGenerator") = python::object(),
               python::arg("bondInvariantsGenerator") = python::object()),
              docString.c_str(),
              python::return_value_policy<
                  python::manage_new_object,
                  python::return_internal_reference<
                      7, python::return_internal_reference<8>>>());

  docString = "";
  // todo add inv generator creator

  return;
}
}  // namespace RDKitFPWrapper

}  // namespace RDKit
