#include <iostream>
#include <string>
#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.h>
#include <GraphMol/Fingerprints/Wrap/AtomPairWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/MorganWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/RDKitFPWrapper.cpp>
#include <cstdint>

namespace python = boost::python;

namespace RDKit {
namespace FingerprintWrapper {

void convertPyArguments(python::object py_fromAtoms,
                        python::object py_ignoreAtoms,
                        std::vector<std::uint32_t> *&fromAtoms,
                        std::vector<std::uint32_t> *&ignoreAtoms) {
  if (!py_fromAtoms.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_fromAtoms.attr("__len__")());
    if (len) {
      fromAtoms = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < len; ++i) {
        fromAtoms->push_back(python::extract<std::uint32_t>(py_fromAtoms[i]));
      }
    }
  }

  if (!py_ignoreAtoms.is_none()) {
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

  return;
}

SparseIntVect<std::uint32_t> *FingerprintGeneratorWrapper::getFingerprint(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, fromAtoms, ignoreAtoms);

  SparseIntVect<std::uint32_t> *result =
      dp_fingerprintGenerator->getFingerprint(mol, fromAtoms, ignoreAtoms,
                                              confId, nullptr);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

SparseBitVect *FingerprintGeneratorWrapper::getFingerprintAsBitVect(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, fromAtoms, ignoreAtoms);

  SparseBitVect *result = dp_fingerprintGenerator->getFingerprintAsBitVect(
      mol, fromAtoms, ignoreAtoms, confId, nullptr);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

SparseIntVect<std::uint32_t> *FingerprintGeneratorWrapper::getFoldedFingerprint(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, fromAtoms, ignoreAtoms);

  SparseIntVect<std::uint32_t> *result =
      dp_fingerprintGenerator->getFoldedFingerprint(mol, fromAtoms, ignoreAtoms,
                                                    confId, nullptr);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

ExplicitBitVect *FingerprintGeneratorWrapper::getFoldedFingerprintAsBitVect(
    const ROMol &mol, python::object py_fromAtoms,
    python::object py_ignoreAtoms, const int confId) const {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, fromAtoms, ignoreAtoms);

  ExplicitBitVect *result =
      dp_fingerprintGenerator->getFoldedFingerprintAsBitVect(
          mol, fromAtoms, ignoreAtoms, confId, nullptr);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

FingerprintGeneratorWrapper::FingerprintGeneratorWrapper() {
  dp_fingerprintGenerator = nullptr;
}

FingerprintGeneratorWrapper::~FingerprintGeneratorWrapper() {
  delete dp_fingerprintGenerator;
}

BOOST_PYTHON_MODULE(rdFingerprintGenerator) {
  std::string docString = "";

  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  docString = "";
  // todo remove the wrapper class if possible
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
  MorganWrapper::exportMorgan();
  RDKitFPWrapper::exportRDKit();
}

}  // namespace FingerprintWrapper
}  // namespace RDKit