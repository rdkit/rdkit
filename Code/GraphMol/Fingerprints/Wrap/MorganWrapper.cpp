
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

using namespace RDKit;
using namespace RDKit::FingerprintWrapper;
namespace python = boost::python;

namespace RDKit {
namespace MorganWrapper {
FingerprintGeneratorWrapper *getMorganGeneratorWrapped(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants, const bool includeRingMembership,
    python::object &py_atomInvGen, python::object &py_bondInvGen) {
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

  AtomEnvironmentGenerator *morganEnvGenerator =
      new MorganFingerprint::MorganEnvGenerator();
  FingerprintArguments *morganArguments =
      new MorganFingerprint::MorganArguments(radius, countSimulation,
                                             includeChirality, useBondTypes,
                                             onlyNonzeroInvariants);

  bool ownsAtomInvGenerator = false;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator =
        new MorganFingerprint::MorganAtomInvGenerator(includeRingMembership);
    ownsAtomInvGenerator = true;
  }

  FingerprintGenerator *fingerprintGenerator = new FingerprintGenerator(
      morganEnvGenerator, morganArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, false);

  FingerprintGeneratorWrapper *wrapped = new FingerprintGeneratorWrapper();
  wrapped->dp_fingerprintGenerator = fingerprintGenerator;
  return wrapped;
}

AtomInvariantsGenerator *getMorganAtomInvGen(const bool includeRingMembership) {
  return new MorganFingerprint::MorganAtomInvGenerator(includeRingMembership);
}

void exportMorgan() {
  std::string docString = "";

  python::def(
      "GetMorganGenerator", &getMorganGeneratorWrapped,
      (python::arg("radius") = 3, python::arg("useCountSimulation") = true,
       python::arg("includeChirality") = false,
       python::arg("useBondTypes") = true,
       python::arg("onlyNonzeroInvariants") = false,
       python::arg("includeRingMembership") = true,
       python::arg("atomInvariantsGenerator") = python::object(),
       python::arg("bondInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<
          python::manage_new_object,
          python::return_internal_reference<
              7, python::return_internal_reference<8>>>());

  docString = "";
  python::def("GetAMorganAtomInvGen", &getMorganAtomInvGen,
              (python::arg("includeRingMembership") = false), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace MorganWrapper

}  // namespace RDKit
