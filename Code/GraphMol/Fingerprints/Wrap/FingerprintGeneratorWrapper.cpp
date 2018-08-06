#include <iostream>
#include <string>
#include <boost/python.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/Wrap/AtomPairWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/MorganWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/RDKitFPWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/TopologicalTorsionWrapper.cpp>
#include <cstdint>

namespace python = boost::python;

namespace RDKit {
namespace FingerprintWrapper {

void convertPyArguments(python::object py_fromAtoms,
                        python::object py_ignoreAtoms,
                        python::object py_atomInvs, python::object py_bondInvs,
                        std::vector<std::uint32_t> *&fromAtoms,
                        std::vector<std::uint32_t> *&ignoreAtoms,
                        std::vector<std::uint32_t> *&customAtomInvariants,
                        std::vector<std::uint32_t> *&customBondInvariants) {
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

  if (!py_atomInvs.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_atomInvs.attr("__len__")());
    if (len) {
      customAtomInvariants = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < len; ++i) {
        customAtomInvariants->push_back(
            python::extract<std::uint32_t>(py_atomInvs[i]));
      }
    }
  }

  if (!py_bondInvs.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_bondInvs.attr("__len__")());
    if (len) {
      customBondInvariants = new std::vector<std::uint32_t>();
      for (unsigned int i = 0; i < len; ++i) {
        customBondInvariants->push_back(
            python::extract<std::uint32_t>(py_bondInvs[i]));
      }
    }
  }

  return;
}

template <typename OutputType>
SparseIntVect<OutputType> *getFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;

  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);

  SparseIntVect<OutputType> *result =
      fpGen->getFingerprint(mol, fromAtoms, ignoreAtoms, confId, nullptr,
                            customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseBitVect *getFingerprintAsBitVect(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);

  SparseBitVect *result = fpGen->getFingerprintAsBitVect(
      mol, fromAtoms, ignoreAtoms, confId, nullptr, customAtomInvariants,
      customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseIntVect<OutputType> *getFoldedFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);

  SparseIntVect<OutputType> *result =
      fpGen->getFoldedFingerprint(mol, fromAtoms, ignoreAtoms, confId, nullptr,
                                  customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
ExplicitBitVect *getFoldedFingerprintAsBitVect(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);

  ExplicitBitVect *result = fpGen->getFoldedFingerprintAsBitVect(
      mol, fromAtoms, ignoreAtoms, confId, nullptr, customAtomInvariants,
      customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

BOOST_PYTHON_MODULE(rdFingerprintGenerator) {

  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  python::class_<FingerprintGenerator<std::uint32_t>, boost::noncopyable>(
      "FingerprintGenerator32", python::no_init)
      .def("GetFingerprint", getFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFingerprintAsBitVect", getFingerprintAsBitVect<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a fingerprint as SparseBitVect\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprint", getFoldedFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a folded fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprintAsBitVect",
           getFoldedFingerprintAsBitVect<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a folded fingerprint as ExplicitBitVect\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>());

  python::class_<FingerprintGenerator<std::uint64_t>, boost::noncopyable>(
      "FingerprintGenerator64", python::no_init)
      .def("GetFingerprint", getFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFingerprintAsBitVect", getFingerprintAsBitVect<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a fingerprint as SparseBitVect\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprint", getFoldedFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a folded fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFoldedFingerprintAsBitVect",
           getFoldedFingerprintAsBitVect<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a folded fingerprint as ExplicitBitVect\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indicies of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indicies of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>());

  AtomPairWrapper::exportAtompair();
  MorganWrapper::exportMorgan();
  RDKitFPWrapper::exportRDKit();
  TopologicalTorsionWrapper::exportTopologicalTorsion();
}

}  // namespace FingerprintWrapper
}  // namespace RDKit