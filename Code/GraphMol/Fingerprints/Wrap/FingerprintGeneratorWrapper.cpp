#include <iostream>
#include <string>
#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.h>
#include <GraphMol/Fingerprints/Wrap/AtomPairWrapper.cpp>
#include <cstdint>

namespace python = boost::python;

namespace RDKit {
namespace FingerprintWrapper {

SparseIntVect<std::uint32_t> *FingerprintGeneratorWrapper::getFingerprint(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) {
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

FingerprintGeneratorWrapper::FingerprintGeneratorWrapper() {
  dp_fingerprintGenerator = nullptr;
}

FingerprintGeneratorWrapper::~FingerprintGeneratorWrapper() {
  delete dp_fingerprintGenerator;
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

  AtomPairWrapper::exportAtompair();
}

}  // namespace FingerprintWrapper
}  // namespace RDKit