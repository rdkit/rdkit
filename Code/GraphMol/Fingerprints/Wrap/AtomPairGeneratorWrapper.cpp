#include <iostream>
#include <string>
#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <cstdint>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace python = boost::python;

class FingerprintGeneratorWrapper {
 public:
  FingerprintGenerator *dp_fingerprintGenerator;

  SparseIntVect<std::uint32_t> *getFingerprint(const ROMol &mol,
                                               python::object py_fromAtoms,
                                               python::object py_ignoreAtoms,
                                               const int confId) {
    std::vector<std::uint32_t> *fromAtoms = nullptr;
    if (py_fromAtoms) {
      unsigned int len =
          python::extract<unsigned int>(py_fromAtoms.attr("__len__")());
      if (len) {
        fromAtoms = new std::vector<std::uint32_t>();
        for (unsigned int i = 0; i < len; ++i) {
          fromAtoms->push_back(python::extract<std::uint32_t>(py_fromAtoms[i]));
        }
      }
    }

    std::vector<std::uint32_t> *ignoreAtoms = nullptr;
    if (py_ignoreAtoms) {
      unsigned int len =
          python::extract<unsigned int>(py_ignoreAtoms.attr("__len__")());
      if (len) {
        ignoreAtoms = new std::vector<std::uint32_t>();
        for (unsigned int i = 0; i < len; ++i) {
          ignoreAtoms->push_back(
              python::extract<std::uint32_t>(py_ignoreAtoms[i]));
        }
      }
    }

    return dp_fingerprintGenerator->getFingerprint(mol, fromAtoms, ignoreAtoms,
                                                   confId, nullptr);
  };

  ~FingerprintGeneratorWrapper() { delete dp_fingerprintGenerator; }
};

FingerprintGeneratorWrapper *getAtomPairGeneratorWrapped(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;
  BondInvariantsGenerator *bondInvariantsGenerator = nullptr;
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

std::string docString = "";

BOOST_PYTHON_MODULE(rdAtomPairGenerator) {
  python::class_<FingerprintGeneratorWrapper>(
      "FingerprintGeneratorWrapper",
      python::init<FingerprintGeneratorWrapper>())
      .def("getFingerprint", &FingerprintGeneratorWrapper::getFingerprint,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1),
           docString.c_str(),
           python::return_value_policy<python::manage_new_object>());

  python::def(
      "getAtomPairGeneratorWrapped", &getAtomPairGeneratorWrapped,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
}
