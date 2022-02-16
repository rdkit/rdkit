//
//  Copyright (C) 2018-2021 Boran Adas and other RDKit contributors
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
#include <RDBoost/boost_numpy.h>
#include <numpy/npy_common.h>
#include <RDBoost/import_array.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/Wrap/AtomPairWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/MorganWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/RDKitFPWrapper.cpp>
#include <GraphMol/Fingerprints/Wrap/TopologicalTorsionWrapper.cpp>
#include <cstdint>

namespace python = boost::python;
namespace np = boost::python::numpy;

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
      fromAtoms->reserve(len);
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
      ignoreAtoms->reserve(len);
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
      customAtomInvariants->reserve(len);
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
      customBondInvariants->reserve(len);
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
    const int confId, python::object py_atomInvs, python::object py_bondInvs,
    python::object py_additionalOutput) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;

  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = python::extract<AdditionalOutput *>(py_additionalOutput);
  }

  SparseIntVect<OutputType> *result = fpGen->getSparseCountFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseBitVect *getSparseFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs,
    python::object py_additionalOutput) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = python::extract<AdditionalOutput *>(py_additionalOutput);
  }

  SparseBitVect *result = fpGen->getSparseFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
SparseIntVect<std::uint32_t> *getCountFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs,
    python::object py_additionalOutput) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = python::extract<AdditionalOutput *>(py_additionalOutput);
  }

  SparseIntVect<std::uint32_t> *result = fpGen->getCountFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
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
                                python::object py_bondInvs,
                                python::object py_additionalOutput) {
  std::vector<std::uint32_t> *fromAtoms = nullptr;
  std::vector<std::uint32_t> *ignoreAtoms = nullptr;
  std::vector<std::uint32_t> *customAtomInvariants = nullptr;
  std::vector<std::uint32_t> *customBondInvariants = nullptr;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = python::extract<AdditionalOutput *>(py_additionalOutput);
  }

  ExplicitBitVect *result = fpGen->getFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants);

  delete fromAtoms;
  delete ignoreAtoms;

  return result;
}

template <typename OutputType>
python::object getNumPyFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs,
    python::object py_additionalOutput) {
  std::unique_ptr<ExplicitBitVect> ebv{
      getFingerprint(fpGen, mol, py_fromAtoms, py_ignoreAtoms, confId,
                     py_atomInvs, py_bondInvs, py_additionalOutput)};

  npy_intp size[1] = {static_cast<npy_intp>(ebv->size())};
  PyObject *arr = PyArray_ZEROS(1, size, NPY_UINT8, 0);
  PyObject *one = PyInt_FromLong(1);
  for (auto i = 0u; i < ebv->size(); ++i) {
    if ((*ebv)[i]) {
      PyArray_SETITEM(
          (PyArrayObject *)arr,
          static_cast<char *>(PyArray_GETPTR1((PyArrayObject *)arr, i)), one);
    }
  }
  Py_DECREF(one);
  python::handle<> res(arr);
  return python::object(res);
}

template <typename OutputType>
python::object getNumPyCountFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    python::object py_fromAtoms, python::object py_ignoreAtoms,
    const int confId, python::object py_atomInvs, python::object py_bondInvs,
    python::object py_additionalOutput) {
  std::unique_ptr<SparseIntVect<uint32_t>> fp{
      getCountFingerprint(fpGen, mol, py_fromAtoms, py_ignoreAtoms, confId,
                          py_atomInvs, py_bondInvs, py_additionalOutput)};

  npy_intp size[1] = {static_cast<npy_intp>(fp->size())};
  PyObject *arr = PyArray_ZEROS(1, size, NPY_UINT32, 0);
  for (auto i = 0u; i < fp->size(); ++i) {
    auto v = (*fp)[i];
    if (v) {
      PyObject *val = PyInt_FromLong(v);
      PyArray_SETITEM(
          (PyArrayObject *)arr,
          static_cast<char *>(PyArray_GETPTR1((PyArrayObject *)arr, i)), val);
      Py_DECREF(val);
    }
  }
  python::handle<> res(arr);
  return python::object(res);
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

  for (auto &it : *tempResult) {
    result.append(boost::shared_ptr<SparseIntVect<std::uint64_t>>(it));
  }
  delete tempResult;
  return result;
}

python::list getSparseFPBulkPy(python::list &py_molVect, FPType fpType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getSparseFPBulk(molVect, fpType);
  python::list result;

  for (auto &it : *tempResult) {
    // todo every other bulk method casts results to boost::shared_ptr, except
    // this one. It should also be boost::shared_ptr
    result.append(boost::shared_ptr<SparseBitVect>(it));
  }
  delete tempResult;
  return result;
}

python::list getCountFPBulkPy(python::list &py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getCountFPBulk(molVect, fPType);
  python::list result;

  for (auto &it : *tempResult) {
    result.append(boost::shared_ptr<SparseIntVect<std::uint32_t>>(it));
  }
  delete tempResult;
  return result;
}

python::list getFPBulkPy(python::list &py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getFPBulk(molVect, fPType);
  python::list result;

  for (auto &it : *tempResult) {
    result.append(boost::shared_ptr<ExplicitBitVect>(it));
  }
  delete tempResult;
  return result;
}

python::object getAtomCountsHelper(const AdditionalOutput &ao) {
  if (!ao.atomCounts) {
    return python::object();
  }
  python::list res;
  for (const auto v : *ao.atomCounts) {
    res.append(v);
  }
  return python::tuple(res);
}
python::object getAtomToBitsHelper(const AdditionalOutput &ao) {
  if (!ao.atomToBits) {
    return python::object();
  }
  python::list res;
  for (const auto &lst : *ao.atomToBits) {
    python::list local;
    for (const auto v : lst) {
      local.append(v);
    }
    res.append(python::tuple(local));
  }
  return python::tuple(res);
}
python::object getBitPathsHelper(const AdditionalOutput &ao) {
  if (!ao.bitPaths) {
    return python::object();
  }
  python::dict res;
  for (const auto &pr : *ao.bitPaths) {
    python::list local;
    for (const auto &lst : pr.second) {
      python::list inner;
      for (const auto v : lst) {
        inner.append(v);
      }
      local.append(python::tuple(inner));
    }
    res[pr.first] = python::tuple(local);
  }
  return std::move(res);
}
python::object getBitInfoMapHelper(const AdditionalOutput &ao) {
  if (!ao.bitInfoMap) {
    return python::object();
  }
  python::dict res;
  for (const auto &pr : *ao.bitInfoMap) {
    python::list local;
    for (const auto &v : pr.second) {
      python::tuple inner = python::make_tuple(v.first, v.second);
      local.append(inner);
    }
    res[pr.first] = python::tuple(local);
  }
  return std::move(res);
}

BOOST_PYTHON_MODULE(rdFingerprintGenerator) {
  rdkit_import_array();
  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  python::class_<AdditionalOutput, boost::noncopyable>("AdditionalOutput")
      .def("AllocateAtomToBits", &AdditionalOutput::allocateAtomToBits)
      .def("AllocateBitInfoMap", &AdditionalOutput::allocateBitInfoMap)
      .def("AllocateBitPaths", &AdditionalOutput::allocateBitPaths)
      .def("AllocateAtomCounts", &AdditionalOutput::allocateAtomCounts)
      .def("GetAtomToBits", &getAtomToBitsHelper)
      .def("GetBitInfoMap", &getBitInfoMapHelper)
      .def("GetBitPaths", &getBitPathsHelper)
      .def("GetAtomCounts", &getAtomCountsHelper);

  python::class_<FingerprintGenerator<std::uint32_t>, boost::noncopyable>(
      "FingerprintGenerator32", python::no_init)
      .def("GetSparseCountFingerprint",
           getSparseCountFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a sparse count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetSparseFingerprint", getSparseFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a sparse fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetCountFingerprint", getCountFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetCountFingerprintAsNumPy",
           getNumPyCountFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a numpy array containing the fingerprint\n\n")
      .def("GetFingerprint", getFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFingerprintAsNumPy", getNumPyFingerprint<std::uint32_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a numpy array containing the fingerprint\n\n")
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
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a sparse count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetSparseFingerprint", getSparseFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a sparse fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetCountFingerprint", getCountFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a SparseIntVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetCountFingerprintAsNumPy",
           getNumPyCountFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a count fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a numpy array containing the fingerprint\n\n")
      .def("GetFingerprint", getFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a ExplicitBitVect containing fingerprint\n\n",
           python::return_value_policy<python::manage_new_object>())
      .def("GetFingerprintAsNumPy", getNumPyFingerprint<std::uint64_t>,
           (python::arg("mol"), python::arg("fromAtoms") = python::list(),
            python::arg("ignoreAtoms") = python::list(),
            python::arg("confId") = -1,
            python::arg("customAtomInvariants") = python::list(),
            python::arg("customBondInvariants") = python::list(),
            python::arg("additionalOutput") = python::object()),
           "Generates a fingerprint\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - fromAtoms: indices of atoms to use while generating the "
           "fingerprint\n"
           "    - ignoreAtoms: indices of atoms to exclude while generating "
           "the fingerprint\n"
           "    - confId: 3D confirmation to use, only used by AtomPair "
           "fingerprint\n"
           "    - customAtomInvariants: custom atom invariants to be used, "
           "overrides invariants from the invariant generator\n"
           "    - customBondInvariants: custom bond invariants to be used, "
           "overrides invariants from the invariant generator\n\n"
           "    - additionalOutput: AdditionalOutput instance used to return "
           "extra information about the bits\n\n"
           "  RETURNS: a numpy array containing fingerprint\n\n")
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
