//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/Wrap.h>
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

#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

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

  FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
                                additionalOutput, customAtomInvariants,
                                customBondInvariants);
  auto result = fpGen->getSparseCountFingerprint(mol, args);

  delete fromAtoms;
  delete ignoreAtoms;

  return result.release();
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

  FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
                                additionalOutput, customAtomInvariants,
                                customBondInvariants);
  auto result = fpGen->getSparseFingerprint(mol, args);

  delete fromAtoms;
  delete ignoreAtoms;

  return result.release();
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

  FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
                                additionalOutput, customAtomInvariants,
                                customBondInvariants);
  auto result = fpGen->getCountFingerprint(mol, args);

  delete fromAtoms;
  delete ignoreAtoms;

  return result.release();
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

  FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
                                additionalOutput, customAtomInvariants,
                                customBondInvariants);
  auto result = fpGen->getFingerprint(mol, args);

  delete fromAtoms;
  delete ignoreAtoms;

  return result.release();
}

template <typename ReturnType, typename FuncType>
python::tuple mtgetFingerprints(FuncType func, python::object mols,
                                int numThreads) {
  unsigned int nmols = python::extract<unsigned int>(mols.attr("__len__")());
  std::vector<const ROMol *> tmols;
  for (auto i = 0u; i < nmols; ++i) {
    tmols.push_back(python::extract<const ROMol *>(mols[i])());
  }

  typename decltype(std::function{
      func})::result_type fps;  // sometimes you really have to love C++ syntax
  {
    NOGIL gil;
    fps = std::move(func(tmols, numThreads));
  }
  python::list result;
  for (auto &fp : fps) {
    result.append(boost::shared_ptr<ReturnType>(fp.release()));
  }
  return python::tuple(result);
}

template <typename OutputType>
python::tuple getFingerprints(const FingerprintGenerator<OutputType> *fpGen,
                              python::object mols, int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<ExplicitBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                              numThreads);
}
template <typename OutputType>
python::tuple getCountFingerprints(
    const FingerprintGenerator<OutputType> *fpGen, python::object mols,
    int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getCountFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseIntVect<std::uint32_t>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
}
template <typename OutputType>
python::tuple getSparseFingerprints(
    const FingerprintGenerator<OutputType> *fpGen, python::object mols,
    int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getSparseFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                            numThreads);
}
template <typename OutputType>
python::tuple getSparseCountFingerprints(
    const FingerprintGenerator<OutputType> *fpGen, python::object mols,
    int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getSparseCountFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseIntVect<OutputType>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
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

template <typename OutputType>
FingerprintArguments *getOptions(FingerprintGenerator<OutputType> *fpGen) {
  return fpGen->getOptions();
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
  return res;
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
  return res;
}

namespace {
template <typename T>
void wrapGenerator(const std::string &nm) {
  python::class_<FingerprintGenerator<T>, boost::noncopyable>(nm.c_str(),
                                                              python::no_init)
      .def("GetSparseCountFingerprint", getSparseCountFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetSparseFingerprint", getSparseFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetCountFingerprint", getCountFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetCountFingerprintAsNumPy", getNumPyCountFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetFingerprint", getFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetFingerprintAsNumPy", getNumPyFingerprint<T>,
           ((python::arg("self"), python::arg("mol")),
            python::arg("fromAtoms") = python::list(),
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
      .def("GetFingerprints", getFingerprints<T>,
           ((python::arg("self"), python::arg("mols")),
            python::arg("numThreads") = 1),
           "Generates fingerprints for a sequence of molecules\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - numThreads: number of threads to use\n\n"
           "  RETURNS: a tuple of ExplicitBitVects\n\n")
      .def("GetCountFingerprints", getCountFingerprints<T>,
           ((python::arg("self"), python::arg("mols")),
            python::arg("numThreads") = 1),
           "Generates count fingerprints for a sequence of molecules\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - numThreads: number of threads to use\n\n"
           "  RETURNS: a tuple of SparseIntVects\n\n")
      .def("GetSparseFingerprints", getSparseFingerprints<T>,
           ((python::arg("self"), python::arg("mols")),
            python::arg("numThreads") = 1),
           "Generates sparse fingerprints for a sequence of molecules\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - numThreads: number of threads to use\n\n"
           "  RETURNS: a tuple of SparseBitVects\n\n")
      .def("GetSparseCountFingerprints", getSparseCountFingerprints<T>,
           ((python::arg("self"), python::arg("mols")),
            python::arg("numThreads") = 1),
           "Generates sparse count fingerprints for a sequence of molecules\n\n"
           "  ARGUMENTS:\n"
           "    - mol: molecule to be fingerprinted\n"
           "    - numThreads: number of threads to use\n\n"
           "  RETURNS: a tuple of SparseIntVects\n\n")
      .def("GetInfoString", getInfoString<T>, python::args("self"),
           "Returns a string containing information about the fingerprint "
           "generator\n\n"
           "  RETURNS: an information string\n\n")
      .def("GetOptions", getOptions<T>,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>(),
           python::args("self"), "return the fingerprint options object");
}

void setCountBoundsHelper(FingerprintArguments &opts, python::object bounds) {
  pythonObjectToVect(bounds, opts.d_countBounds);
}
}  // namespace

BOOST_PYTHON_MODULE(rdFingerprintGenerator) {
  rdkit_import_array();
  python::class_<AtomInvariantsGenerator, boost::noncopyable>(
      "AtomInvariantsGenerator", python::no_init);

  python::class_<BondInvariantsGenerator, boost::noncopyable>(
      "BondInvariantsGenerator", python::no_init);

  python::class_<AdditionalOutput, boost::noncopyable>("AdditionalOutput")
      .def("AllocateAtomToBits", &AdditionalOutput::allocateAtomToBits,
           python::args("self"), "synonym for CollectAtomToBits()")
      .def("AllocateBitInfoMap", &AdditionalOutput::allocateBitInfoMap,
           python::args("self"), "synonym for CollectBitInfoMap()")
      .def("AllocateBitPaths", &AdditionalOutput::allocateBitPaths,
           python::args("self"), "synonym for CollectBitPaths()")
      .def("AllocateAtomCounts", &AdditionalOutput::allocateAtomCounts,
           python::args("self"), "synonym for CollectAtomCounts()")
      .def(
          "CollectAtomToBits", &AdditionalOutput::allocateAtomToBits,
          python::args("self"),
          "toggle collection of information mapping each atom to the bits it is involved in.")
      .def(
          "CollectBitInfoMap", &AdditionalOutput::allocateBitInfoMap,
          python::args("self"),
          "toggles collection of information mapping each atom to more detail about the atom environment (not available from all fingerprints)")
      .def(
          "CollectBitPaths", &AdditionalOutput::allocateBitPaths,
          python::args("self"),
          "toggles collection of information matching each atom to information about the paths it is involved in (not available from all fingerprints).")
      .def(
          "CollectAtomCounts", &AdditionalOutput::allocateAtomCounts,
          python::args("self"),
          "toggles collection of information about the number of bits each atom is involved in")
      .def("GetAtomToBits", &getAtomToBitsHelper, python::args("self"))
      .def("GetBitInfoMap", &getBitInfoMapHelper, python::args("self"))
      .def("GetBitPaths", &getBitPathsHelper, python::args("self"))
      .def("GetAtomCounts", &getAtomCountsHelper, python::args("self"));

  python::class_<FingerprintArguments, boost::noncopyable>("FingerprintOptions",
                                                           python::no_init)
      .def_readwrite("countSimulation",
                     &FingerprintArguments::df_countSimulation,
                     "use count simulation")
      .def_readwrite(
          "includeChirality", &FingerprintArguments::df_includeChirality,
          "include chirality in atom invariants (not for all fingerprints)")
      .def_readwrite("fpSize", &FingerprintArguments::d_fpSize,
                     "size of the fingerprints created")
      .def_readwrite("numBitsPerFeature",
                     &FingerprintArguments::d_numBitsPerFeature,
                     "number of bits to set for each feature")
      .def("SetCountBounds", &setCountBoundsHelper,
           python::args("self", "bounds"), "set the bins for the count bounds");

  wrapGenerator<std::uint32_t>("FingeprintGenerator32");
  wrapGenerator<std::uint64_t>("FingeprintGenerator64");

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
