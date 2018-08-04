
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
    python::object &py_countBounds, const std::uint32_t foldedSize,
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
  if (countBoundsE.check() && !countBoundsE().empty() ) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return MorganFingerprint::getMorganGenerator<OutputType>(
      radius, countSimulation, includeChirality, useBondTypes,
      onlyNonzeroInvariants, atomInvariantsGenerator, bondInvariantsGenerator,
      foldedSize, countBounds, true, true);
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
  std::string docString = "";

  /*python::def(
      "GetMorganGenerator32", getMorganGenerator<std::uint32_t>,
      (python::arg("radius") = 3, python::arg("useCountSimulation") = true,
       python::arg("includeChirality") = false,
       python::arg("useBondTypes") = true,
       python::arg("onlyNonzeroInvariants") = false,
       python::arg("includeRingMembership") = true,
       python::arg("countBounds") = python::object(),
       python::arg("foldedSize") = 2048,
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
       python::arg("foldedSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object(),
       python::arg("bondInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString = "";
  python::def("GetMorganAtomInvGen", &getMorganAtomInvGen,
              (python::arg("includeRingMembership") = false), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "";
  python::def("GetMorganFeatureAtomInvGen", &getMorganFeatureAtomInvGen,
              (python::arg("patterns") = python::object()), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "";
  python::def(
      "GetMorganBondInvGen", &getMorganBondInvGen,
      (python::arg("useBondTypes") = true, python::arg("useChirality") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace MorganWrapper

}  // namespace RDKit
