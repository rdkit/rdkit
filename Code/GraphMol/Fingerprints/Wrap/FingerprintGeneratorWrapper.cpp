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

struct ArgumentHelper {
  std::vector<std::uint32_t> *dp_fromAtoms;
  std::vector<std::uint32_t> *dp_ignoreAtoms;
};

ArgumentHelper convertPyArguments(python::object py_fromAtoms,
                                  python::object py_ignoreAtoms) {
  ArgumentHelper argumentHelper;

  argumentHelper.dp_fromAtoms = nullptr;
  argumentHelper.dp_ignoreAtoms = nullptr;

  if (!py_fromAtoms.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_fromAtoms.attr("__len__")());
    if (len) {
      argumentHelper.dp_fromAtoms = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < len; ++i) {
        argumentHelper.dp_fromAtoms->push_back(
            python::extract<std::uint32_t>(py_fromAtoms[i]));
      }
    }
  }

  if (!py_ignoreAtoms.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_ignoreAtoms.attr("__len__")());
    if (len) {
      argumentHelper.dp_ignoreAtoms = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < len; ++i) {
        argumentHelper.dp_ignoreAtoms->push_back(
            python::extract<std::uint32_t>(py_ignoreAtoms[i]));
      }
    }
  }

  return argumentHelper;
}

SparseIntVect<std::uint32_t> *FingerprintGeneratorWrapper::getFingerprint(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  ArgumentHelper argumentHelper =
      convertPyArguments(py_fromAtoms, py_ignoreAtoms);

  SparseIntVect<std::uint32_t> *result =
      dp_fingerprintGenerator->getFingerprint(mol, argumentHelper.dp_fromAtoms,
                                              argumentHelper.dp_ignoreAtoms,
                                              confId, nullptr);

  delete argumentHelper.dp_fromAtoms;
  delete argumentHelper.dp_ignoreAtoms;

  return result;
}

SparseBitVect *FingerprintGeneratorWrapper::getFingerprintAsBitVect(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  ArgumentHelper argumentHelper =
      convertPyArguments(py_fromAtoms, py_ignoreAtoms);

  SparseBitVect *result = dp_fingerprintGenerator->getFingerprintAsBitVect(
      mol, argumentHelper.dp_fromAtoms, argumentHelper.dp_ignoreAtoms, confId,
      nullptr);

  delete argumentHelper.dp_fromAtoms;
  delete argumentHelper.dp_ignoreAtoms;

  return result;
}

SparseIntVect<std::uint32_t> *FingerprintGeneratorWrapper::getFoldedFingerprint(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  ArgumentHelper argumentHelper =
      convertPyArguments(py_fromAtoms, py_ignoreAtoms);

  SparseIntVect<std::uint32_t> *result =
      dp_fingerprintGenerator->getFoldedFingerprint(
          mol, argumentHelper.dp_fromAtoms, argumentHelper.dp_ignoreAtoms,
          confId, nullptr);

  delete argumentHelper.dp_fromAtoms;
  delete argumentHelper.dp_ignoreAtoms;

  return result;
}

ExplicitBitVect *FingerprintGeneratorWrapper::getFoldedFingerprintAsBitVect(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  ArgumentHelper argumentHelper =
      convertPyArguments(py_fromAtoms, py_ignoreAtoms);

  ExplicitBitVect *result =
      dp_fingerprintGenerator->getFoldedFingerprintAsBitVect(
          mol, argumentHelper.dp_fromAtoms, argumentHelper.dp_ignoreAtoms,
          confId, nullptr);

  delete argumentHelper.dp_fromAtoms;
  delete argumentHelper.dp_ignoreAtoms;

  return result;
}

FingerprintGeneratorWrapper::FingerprintGeneratorWrapper() {
  dp_fingerprintGenerator = nullptr;
}

FingerprintGeneratorWrapper::~FingerprintGeneratorWrapper() {
  delete dp_fingerprintGenerator;
}

BOOST_PYTHON_MODULE(rdAtomPairGenerator) {
  std::string docString = "";

  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  docString = "";
  python::class_<FingerprintGeneratorWrapper>(
      "FingerprintGenerator", python::init<FingerprintGeneratorWrapper>())
      .def("GetFingerprint", &FingerprintGeneratorWrapper::getFingerprint,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1),
           docString.c_str(),
           python::return_value_policy<python::manage_new_object>())
      .def("GetFingerprintAsBitVect",
           &FingerprintGeneratorWrapper::getFingerprintAsBitVect,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1),
           docString.c_str(),
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprint",
           &FingerprintGeneratorWrapper::getFoldedFingerprint,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1),
           docString.c_str(),
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprintAsBitVect",
           &FingerprintGeneratorWrapper::getFoldedFingerprintAsBitVect,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1),
           docString.c_str(),
           python::return_value_policy<python::manage_new_object>());

  AtomPairWrapper::exportAtompair();
}

}  // namespace FingerprintWrapper
}  // namespace RDKit