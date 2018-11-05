//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
SparseIntVect<OutputType> *getSparseCountFingerprint(
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

  SparseIntVect<OutputType> *result = fpGen->getSparseCountFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, nullptr, customAtomInvariants,
      customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseBitVect *getSparseFingerprint(
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

  SparseBitVect *result =
      fpGen->getSparseFingerprint(mol, fromAtoms, ignoreAtoms, confId, nullptr,
                                  customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseIntVect<std::uint32_t> *getCountFingerprint(
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

  SparseIntVect<std::uint32_t> *result =
      fpGen->getCountFingerprint(mol, fromAtoms, ignoreAtoms, confId, nullptr,
                                 customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
ExplicitBitVect *getFingerprint(const FingerprintGenerator<OutputType> *fpGen,
                                const ROMol &mol, python::object py_fromAtoms,
                                python::object py_ignoreAtoms, const int confId,
                                python::object py_atomInvs,
                                python::object py_bondInvs) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);

  ExplicitBitVect *result =
      fpGen->getFingerprint(mol, fromAtoms, ignoreAtoms, confId, nullptr,
                            customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
std::string getInfoString(const FingerprintGenerator<OutputType> *fpGen) {
  return std::string(fpGen->infoString());
}

const std::vector<const ROMol *> convertPyArgumentsForBulk(
    const python::list &py_molVect) {
  std::vector<const ROMol *> molVect;
  if (!py_molVect.is_none()) {
    unsigned int len =
        python::extract<unsigned int>(py_molVect.attr("__len__")());
    if (len) {
      for (unsigned int i = 0; i < len; i++) {
        molVect.push_back(python::extract<const ROMol *>(py_molVect[i]));
      }
    }
  }
  return molVect;
}

python::list getSparseCountFPBulkPy(python::list &py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getSparseCountFPBulk(molVect, fPType);
  python::list result;

  for (auto it = tempResult->begin(); it != tempResult->end(); ++it) {
    result.append((boost::shared_ptr<SparseIntVect<std::uint64_t>>)*it);
  }
  delete tempResult;
  return result;
}

python::list getSparseFPBulkPy(python::list &py_molVect, FPType fpType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getSparseFPBulk(molVect, fpType);
  python::list result;

  for (auto it = tempResult->begin(); it != tempResult->end(); ++it) {
    // todo every other bulk method casts results to boost::shared_ptr, except
    // this one. It should also be boost::shared_ptr
    result.append(*it);
  }
  delete tempResult;
  return result;
}

python::list getCountFPBulkPy(python::list &py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getCountFPBulk(molVect, fPType);
  python::list result;

  for (auto it = tempResult->begin(); it != tempResult->end(); ++it) {
    result.append((boost::shared_ptr<SparseIntVect<std::uint32_t>>)*it);
  }
  delete tempResult;
  return result;
}

python::list getFPBulkPy(python::list &py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getFPBulk(molVect, fPType);
  python::list result;

  for (auto it = tempResult->begin(); it != tempResult->end(); ++it) {
    result.append((boost::shared_ptr<ExplicitBitVect>)*it);
  }
  delete tempResult;
  return result;
}

BOOST_PYTHON_MODULE(rdFingerprintGenerator) {
  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  python::class_<FingerprintGenerator<std::uint32_t>, boost::noncopyable>(
      "FingerprintGenerator32", python::no_init)
      .def("GetSparseCountFingerprint",
           getSparseCountFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a sparse count fingerprint\n\n"
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
      .def("GetSparseFingerprint", getSparseFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a sparse fingerprint\n\n"
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
      .def("GetCountFingerprint", getCountFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a count fingerprint\n\n"
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
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetInfoString", getInfoString<std::uint32_t>,
           "Returns a string containing information about the fingerprint "
           "generator\n\n"
           "  RETURNS: an information string\n\n");

  python::class_<FingerprintGenerator<std::uint64_t>, boost::noncopyable>(
      "FingerprintGenerator64", python::no_init)
      .def("GetSparseCountFingerprint",
           getSparseCountFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a sparse count fingerprint\n\n"
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
      .def("GetSparseFingerprint", getSparseFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a sparse fingerprint\n\n"
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
      .def("GetCountFingerprint", getCountFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list()),
           "Generates a count fingerprint\n\n"
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
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetInfoString", getInfoString<std::uint64_t>,
           "Returns a string containing information about the fingerprint "
           "generator\n\n"
           "  RETURNS: an information string\n\n");

  python::enum_<FPType>("FPType")
      .value("RDKitFP", FPType::RDKitFP)
      .value("MorganFP", FPType::MorganFP)
      .value("AtomPairFP", FPType::AtomPairFP)
      .value("TopologicalTorsionFP", FPType::TopologicalTorsionFP)
      .export_values();
  ;

  python::def("GetSparseCountFPs", &getSparseCountFPBulkPy,
              (python::arg("molecules") = python::list(),
               python::arg("fpType") = FPType::MorganFP),
              "");

  python::def("GetSparseFPs", &getSparseFPBulkPy,
              (python::arg("molecules") = python::list(),
               python::arg("fpType") = FPType::MorganFP),
              "");

  python::def("GetCountFPs", &getCountFPBulkPy,
              (python::arg("molecules") = python::list(),
               python::arg("fpType") = FPType::MorganFP),
              "");

  python::def("GetFPs", &getFPBulkPy,
              (python::arg("molecules") = python::list(),
               python::arg("fpType") = FPType::MorganFP),
              "");

  AtomPairWrapper::exportAtompair();
  MorganWrapper::exportMorgan();
  RDKitFPWrapper::exportRDKit();
  TopologicalTorsionWrapper::exportTopologicalTorsion();
}

}  // namespace FingerprintWrapper
}  // namespace RDKit
