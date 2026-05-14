//
//  Copyright (C) 2018-2026 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <string>
#include <numpy/npy_common.h>
#include <numpy/ndarrayobject.h>
#include <RDBoost/import_array.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/nbWrap/AtomPairWrapper.cpp>
#include <GraphMol/Fingerprints/nbWrap/MorganWrapper.cpp>
#include <GraphMol/Fingerprints/nbWrap/RDKitFPWrapper.cpp>
#include <GraphMol/Fingerprints/nbWrap/TopologicalTorsionWrapper.cpp>
#include <cstdint>

#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace FingerprintWrapper {

void convertPyArguments(
    nb::object py_fromAtoms, nb::object py_ignoreAtoms,
    nb::object py_atomInvs, nb::object py_bondInvs,
    std::unique_ptr<std::vector<std::uint32_t>> &fromAtoms,
    std::unique_ptr<std::vector<std::uint32_t>> &ignoreAtoms,
    std::unique_ptr<std::vector<std::uint32_t>> &customAtomInvariants,
    std::unique_ptr<std::vector<std::uint32_t>> &customBondInvariants) {
  auto convert = [](nb::object obj,
                    std::unique_ptr<std::vector<std::uint32_t>> &vec) {
    if (!obj.is_none()) {
      size_t len = nb::len(obj);
      if (len) {
        vec.reset(new std::vector<std::uint32_t>());
        vec->reserve(len);
        for (auto item : obj) {
          vec->push_back(nb::cast<std::uint32_t>(item));
        }
      }
    }
  };
  convert(py_fromAtoms, fromAtoms);
  convert(py_ignoreAtoms, ignoreAtoms);
  convert(py_atomInvs, customAtomInvariants);
  convert(py_bondInvs, customBondInvariants);
}

template <typename OutputType>
SparseIntVect<OutputType> *getSparseCountFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    nb::object py_fromAtoms, nb::object py_ignoreAtoms, const int confId,
    nb::object py_atomInvs, nb::object py_bondInvs,
    nb::object py_additionalOutput) {
  std::unique_ptr<std::vector<std::uint32_t>> fromAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> ignoreAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> customAtomInvariants;
  std::unique_ptr<std::vector<std::uint32_t>> customBondInvariants;

  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = nb::cast<AdditionalOutput *>(py_additionalOutput);
  }

  FingerprintFuncArguments args(fromAtoms.get(), ignoreAtoms.get(), confId,
                                additionalOutput, customAtomInvariants.get(),
                                customBondInvariants.get());
  auto result = fpGen->getSparseCountFingerprint(mol, args);

  return result.release();
}

template <typename OutputType>
SparseBitVect *getSparseFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    nb::object py_fromAtoms, nb::object py_ignoreAtoms, const int confId,
    nb::object py_atomInvs, nb::object py_bondInvs,
    nb::object py_additionalOutput) {
  std::unique_ptr<std::vector<std::uint32_t>> fromAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> ignoreAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> customAtomInvariants;
  std::unique_ptr<std::vector<std::uint32_t>> customBondInvariants;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = nb::cast<AdditionalOutput *>(py_additionalOutput);
  }

  FingerprintFuncArguments args(fromAtoms.get(), ignoreAtoms.get(), confId,
                                additionalOutput, customAtomInvariants.get(),
                                customBondInvariants.get());
  auto result = fpGen->getSparseFingerprint(mol, args);

  return result.release();
}

template <typename OutputType>
SparseIntVect<std::uint32_t> *getCountFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    nb::object py_fromAtoms, nb::object py_ignoreAtoms, const int confId,
    nb::object py_atomInvs, nb::object py_bondInvs,
    nb::object py_additionalOutput) {
  std::unique_ptr<std::vector<std::uint32_t>> fromAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> ignoreAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> customAtomInvariants;
  std::unique_ptr<std::vector<std::uint32_t>> customBondInvariants;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = nb::cast<AdditionalOutput *>(py_additionalOutput);
  }

  FingerprintFuncArguments args(fromAtoms.get(), ignoreAtoms.get(), confId,
                                additionalOutput, customAtomInvariants.get(),
                                customBondInvariants.get());
  auto result = fpGen->getCountFingerprint(mol, args);

  return result.release();
}

template <typename OutputType>
ExplicitBitVect *getFingerprint(const FingerprintGenerator<OutputType> *fpGen,
                                const ROMol &mol, nb::object py_fromAtoms,
                                nb::object py_ignoreAtoms, const int confId,
                                nb::object py_atomInvs, nb::object py_bondInvs,
                                nb::object py_additionalOutput) {
  std::unique_ptr<std::vector<std::uint32_t>> fromAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> ignoreAtoms;
  std::unique_ptr<std::vector<std::uint32_t>> customAtomInvariants;
  std::unique_ptr<std::vector<std::uint32_t>> customBondInvariants;
  convertPyArguments(py_fromAtoms, py_ignoreAtoms, py_atomInvs, py_bondInvs,
                     fromAtoms, ignoreAtoms, customAtomInvariants,
                     customBondInvariants);
  AdditionalOutput *additionalOutput = nullptr;
  if (!py_additionalOutput.is_none()) {
    additionalOutput = nb::cast<AdditionalOutput *>(py_additionalOutput);
  }

  FingerprintFuncArguments args(fromAtoms.get(), ignoreAtoms.get(), confId,
                                additionalOutput, customAtomInvariants.get(),
                                customBondInvariants.get());
  auto result = fpGen->getFingerprint(mol, args);

  return result.release();
}

template <typename ReturnType, typename FuncType>
nb::tuple mtgetFingerprints(FuncType func, nb::object mols, int numThreads) {
  std::vector<const ROMol *> tmols;
  for (auto item : mols) {
    tmols.push_back(nb::cast<const ROMol *>(item));
  }

  typename decltype(std::function{func})::result_type fps;
  {
    nb::gil_scoped_release release;
    fps = std::move(func(tmols, numThreads));
  }
  nb::list result;
  for (auto &fp : fps) {
    result.append(nb::cast(fp.release(), nb::rv_policy::take_ownership));
  }
  return nb::tuple(result);
}

template <typename OutputType>
nb::tuple getFingerprints(const FingerprintGenerator<OutputType> *fpGen,
                          nb::object mols, int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<ExplicitBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                              numThreads);
}

template <typename OutputType>
nb::tuple getCountFingerprints(const FingerprintGenerator<OutputType> *fpGen,
                               nb::object mols, int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getCountFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseIntVect<std::uint32_t>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
}

template <typename OutputType>
nb::tuple getSparseFingerprints(const FingerprintGenerator<OutputType> *fpGen,
                                nb::object mols, int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getSparseFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseBitVect, decltype(fpfunc)>(fpfunc, mols,
                                                            numThreads);
}

template <typename OutputType>
nb::tuple getSparseCountFingerprints(
    const FingerprintGenerator<OutputType> *fpGen, nb::object mols,
    int numThreads) {
  auto fpfunc = [&fpGen](const std::vector<const ROMol *> &tmols,
                         int numThreads) {
    return fpGen->getSparseCountFingerprints(tmols, numThreads);
  };
  return mtgetFingerprints<SparseIntVect<OutputType>, decltype(fpfunc)>(
      fpfunc, mols, numThreads);
}

template <typename OutputType>
nb::object getNumPyFingerprint(const FingerprintGenerator<OutputType> *fpGen,
                               const ROMol &mol, nb::object py_fromAtoms,
                               nb::object py_ignoreAtoms, const int confId,
                               nb::object py_atomInvs, nb::object py_bondInvs,
                               nb::object py_additionalOutput) {
  std::unique_ptr<ExplicitBitVect> ebv{
      getFingerprint(fpGen, mol, py_fromAtoms, py_ignoreAtoms, confId,
                     py_atomInvs, py_bondInvs, py_additionalOutput)};

  npy_intp size[1] = {static_cast<npy_intp>(ebv->size())};
  PyObject *arr = PyArray_ZEROS(1, size, NPY_UINT8, 0);
  PyObject *one = PyLong_FromLong(1);
  for (auto i = 0u; i < ebv->size(); ++i) {
    if ((*ebv)[i]) {
      PyArray_SETITEM(
          (PyArrayObject *)arr,
          static_cast<char *>(PyArray_GETPTR1((PyArrayObject *)arr, i)), one);
    }
  }
  Py_DECREF(one);
  return nb::steal<nb::object>(arr);
}

template <typename OutputType>
nb::object getNumPyCountFingerprint(
    const FingerprintGenerator<OutputType> *fpGen, const ROMol &mol,
    nb::object py_fromAtoms, nb::object py_ignoreAtoms, const int confId,
    nb::object py_atomInvs, nb::object py_bondInvs,
    nb::object py_additionalOutput) {
  std::unique_ptr<SparseIntVect<uint32_t>> fp{
      getCountFingerprint(fpGen, mol, py_fromAtoms, py_ignoreAtoms, confId,
                          py_atomInvs, py_bondInvs, py_additionalOutput)};

  npy_intp size[1] = {static_cast<npy_intp>(fp->size())};
  PyObject *arr = PyArray_ZEROS(1, size, NPY_UINT32, 0);
  for (auto i = 0u; i < fp->size(); ++i) {
    auto v = (*fp)[i];
    if (v) {
      PyObject *val = PyLong_FromLong(v);
      PyArray_SETITEM(
          (PyArrayObject *)arr,
          static_cast<char *>(PyArray_GETPTR1((PyArrayObject *)arr, i)), val);
      Py_DECREF(val);
    }
  }
  return nb::steal<nb::object>(arr);
}

const std::vector<const ROMol *> convertPyArgumentsForBulk(
    nb::object py_molVect) {
  std::vector<const ROMol *> molVect;
  if (!py_molVect.is_none()) {
    for (auto item : py_molVect) {
      molVect.push_back(nb::cast<const ROMol *>(item));
    }
  }
  return molVect;
}

nb::list getSparseCountFPBulkPy(nb::object py_molVect, FPType fPType) {
  const auto molVect = convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getSparseCountFPBulk(molVect, fPType);
  nb::list result;

  for (auto &it : *tempResult) {
    result.append(nb::cast(it, nb::rv_policy::take_ownership));
  }
  delete tempResult;
  return result;
}

nb::list getSparseFPBulkPy(nb::object py_molVect, FPType fpType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getSparseFPBulk(molVect, fpType);
  nb::list result;

  for (auto &it : *tempResult) {
    result.append(nb::cast(it, nb::rv_policy::take_ownership));
  }
  delete tempResult;
  return result;
}

nb::list getCountFPBulkPy(nb::object py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getCountFPBulk(molVect, fPType);
  nb::list result;

  for (auto &it : *tempResult) {
    result.append(nb::cast(it, nb::rv_policy::take_ownership));
  }
  delete tempResult;
  return result;
}

nb::list getFPBulkPy(nb::object py_molVect, FPType fPType) {
  const std::vector<const ROMol *> molVect =
      convertPyArgumentsForBulk(py_molVect);
  auto tempResult = getFPBulk(molVect, fPType);
  nb::list result;

  for (auto &it : *tempResult) {
    result.append(nb::cast(it, nb::rv_policy::take_ownership));
  }
  delete tempResult;
  return result;
}

nb::object getAtomCountsHelper(const AdditionalOutput &ao) {
  if (!ao.atomCounts) {
    return nb::none();
  }
  nb::list res;
  for (const auto v : *ao.atomCounts) {
    res.append(v);
  }
  return nb::tuple(res);
}

nb::object getAtomToBitsHelper(const AdditionalOutput &ao) {
  if (!ao.atomToBits) {
    return nb::none();
  }
  nb::list res;
  for (const auto &lst : *ao.atomToBits) {
    nb::list local;
    for (const auto v : lst) {
      local.append(v);
    }
    res.append(nb::tuple(local));
  }
  return nb::tuple(res);
}

nb::object getBitPathsHelper(const AdditionalOutput &ao) {
  if (!ao.bitPaths) {
    return nb::none();
  }
  nb::dict res;
  for (const auto &pr : *ao.bitPaths) {
    nb::list local;
    for (const auto &lst : pr.second) {
      nb::list inner;
      for (const auto v : lst) {
        inner.append(v);
      }
      local.append(nb::tuple(inner));
    }
    res[nb::int_(pr.first)] = nb::tuple(local);
  }
  return res;
}

nb::object getBitInfoMapHelper(const AdditionalOutput &ao) {
  if (!ao.bitInfoMap) {
    return nb::none();
  }
  nb::dict res;
  for (const auto &pr : *ao.bitInfoMap) {
    nb::list local;
    for (const auto &v : pr.second) {
      local.append(nb::make_tuple(v.first, v.second));
    }
    res[nb::int_(pr.first)] = nb::tuple(local);
  }
  return res;
}

nb::object getAtomsPerBitHelper(const AdditionalOutput &ao) {
  if (!ao.atomsPerBit) {
    return nb::none();
  }
  nb::dict res;
  for (const auto &pr : *ao.atomsPerBit) {
    nb::list local;
    for (const auto &lst : pr.second) {
      nb::list inner;
      for (const auto v : lst) {
        inner.append(v);
      }
      local.append(nb::tuple(inner));
    }
    res[nb::int_(pr.first)] = nb::tuple(local);
  }
  return res;
}

namespace {
template <typename T>
void wrapGenerator(nb::module_ &m, const std::string &nm) {
  nb::class_<FingerprintGenerator<T>>(m, nm.c_str())
      .def("GetSparseCountFingerprint", getSparseCountFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a sparse count fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a SparseIntVect containing fingerprint
)DOC",
           nb::rv_policy::take_ownership)
      .def("GetSparseFingerprint", getSparseFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a sparse fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a SparseBitVect containing fingerprint
)DOC",
           nb::rv_policy::take_ownership)
      .def("GetCountFingerprint", getCountFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a count fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a SparseIntVect containing fingerprint
)DOC",
           nb::rv_policy::take_ownership)
      .def("GetCountFingerprintAsNumPy", getNumPyCountFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a count fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a numpy array containing the fingerprint
)DOC")
      .def("GetFingerprint", getFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a ExplicitBitVect containing fingerprint
)DOC",
           nb::rv_policy::take_ownership)
      .def("GetFingerprintAsNumPy", getNumPyFingerprint<T>,
           "mol"_a, "fromAtoms"_a = nb::none(), "ignoreAtoms"_a = nb::none(),
           "confId"_a = -1, "customAtomInvariants"_a = nb::none(),
           "customBondInvariants"_a = nb::none(),
           "additionalOutput"_a = nb::none(),
           R"DOC(Generates a fingerprint

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - fromAtoms: indices of atoms to use while generating the fingerprint
    - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
    - confId: 3D confirmation to use, only used by AtomPair fingerprint
    - customAtomInvariants: custom atom invariants to be used,
      overrides invariants from the invariant generator
    - customBondInvariants: custom bond invariants to be used,
      overrides invariants from the invariant generator
    - additionalOutput: AdditionalOutput instance used to return extra information about the bits

RETURNS: a numpy array containing the fingerprint
)DOC")
      .def("GetFingerprints", getFingerprints<T>, "mols"_a,
           "numThreads"_a = 1,
           R"DOC(Generates fingerprints for a sequence of molecules

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - numThreads: number of threads to use

RETURNS: a tuple of ExplicitBitVects
)DOC")
      .def("GetCountFingerprints", getCountFingerprints<T>, "mols"_a,
           "numThreads"_a = 1,
           R"DOC(Generates count fingerprints for a sequence of molecules

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - numThreads: number of threads to use

RETURNS: a tuple of SparseIntVects
)DOC")
      .def("GetSparseFingerprints", getSparseFingerprints<T>, "mols"_a,
           "numThreads"_a = 1,
           R"DOC(Generates sparse fingerprints for a sequence of molecules

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - numThreads: number of threads to use

RETURNS: a tuple of SparseBitVects
)DOC")
      .def("GetSparseCountFingerprints", getSparseCountFingerprints<T>,
           "mols"_a, "numThreads"_a = 1,
           R"DOC(Generates sparse count fingerprints for a sequence of molecules

ARGUMENTS:
    - mol: molecule to be fingerprinted
    - numThreads: number of threads to use

RETURNS: a tuple of SparseIntVects
)DOC")
      .def("GetInfoString",
           [](const FingerprintGenerator<T> *fpGen) {
             return std::string(fpGen->infoString());
           },
           R"DOC(Returns a string containing information about the fingerprint generator

RETURNS: an information string
)DOC")
      .def("GetOptions",
           [](FingerprintGenerator<T> *fpGen) { return fpGen->getOptions(); },
           nb::rv_policy::reference_internal, "return the fingerprint options object")
      .def("ToJSON", &generatorToJSON<T>,
           "Serialize a FingerprintGenerator to JSON");
}

}  // namespace

NB_MODULE(rdFingerprintGenerator, m) {
  rdkit_import_array();

  nb::class_<AtomInvariantsGenerator>(m, "AtomInvariantsGenerator");
  nb::class_<BondInvariantsGenerator>(m, "BondInvariantsGenerator");

  nb::class_<AdditionalOutput>(m, "AdditionalOutput")
      .def(nb::init<>())
      .def("AllocateAtomToBits", &AdditionalOutput::allocateAtomToBits,
           "synonym for CollectAtomToBits()")
      .def("AllocateBitInfoMap", &AdditionalOutput::allocateBitInfoMap,
           "synonym for CollectBitInfoMap()")
      .def("AllocateBitPaths", &AdditionalOutput::allocateBitPaths,
           "synonym for CollectBitPaths()")
      .def("AllocateAtomCounts", &AdditionalOutput::allocateAtomCounts,
           "synonym for CollectAtomCounts()")
      .def("AllocateAtomsPerBit", &AdditionalOutput::allocateAtomsPerBit,
           "synonym for CollectAtomsPerBit()")
      .def("CollectAtomToBits", &AdditionalOutput::allocateAtomToBits,
           R"DOC(toggle collection of information mapping each atom to the bits it is involved in.)DOC")
      .def("CollectBitInfoMap", &AdditionalOutput::allocateBitInfoMap,
           R"DOC(toggles collection of information mapping each atom to more detail about the atom environment (not available from all fingerprints))DOC")
      .def("CollectBitPaths", &AdditionalOutput::allocateBitPaths,
           R"DOC(toggles collection of information matching each atom to information about the paths it is involved in (not available from all fingerprints).)DOC")
      .def("CollectAtomCounts", &AdditionalOutput::allocateAtomCounts,
           R"DOC(toggles collection of information about the number of bits each atom is involved in)DOC")
      .def("CollectAtomsPerBit", &AdditionalOutput::allocateAtomsPerBit,
           R"DOC(toggles collection of information about all atoms involved in setting each bit)DOC")
      .def("GetAtomToBits", &getAtomToBitsHelper)
      .def("GetBitInfoMap", &getBitInfoMapHelper)
      .def("GetBitPaths", &getBitPathsHelper)
      .def("GetAtomCounts", &getAtomCountsHelper)
      .def("GetAtomsPerBit", &getAtomsPerBitHelper);

  nb::class_<FingerprintArguments>(m, "FingerprintOptions")
      .def_rw("countSimulation", &FingerprintArguments::df_countSimulation,
              "use count simulation")
      .def_rw(
          "includeChirality", &FingerprintArguments::df_includeChirality,
          "include chirality in atom invariants (not for all fingerprints)")
      .def_rw("fpSize", &FingerprintArguments::d_fpSize,
              "size of the fingerprints created")
      .def_rw("numBitsPerFeature", &FingerprintArguments::d_numBitsPerFeature,
              "number of bits to set for each feature")
      .def("SetCountBounds",
           [](FingerprintArguments &opts, nb::object bounds) {
             opts.d_countBounds.clear();
             for (auto item : bounds) {
               opts.d_countBounds.push_back(nb::cast<std::uint32_t>(item));
             }
           },
           "bounds"_a,
           "set the bins for the count bounds");

  wrapGenerator<std::uint32_t>(m, "FingerprintGenerator32");
  wrapGenerator<std::uint64_t>(m, "FingerprintGenerator64");

  nb::enum_<FPType>(m, "FPType")
      .value("RDKitFP", FPType::RDKitFP)
      .value("MorganFP", FPType::MorganFP)
      .value("AtomPairFP", FPType::AtomPairFP)
      .value("TopologicalTorsionFP", FPType::TopologicalTorsionFP)
      .export_values();

  m.def("GetSparseCountFPs", &getSparseCountFPBulkPy,
        "molecules"_a = nb::none(), "fpType"_a = FPType::MorganFP, "");

  m.def("GetSparseFPs", &getSparseFPBulkPy, "molecules"_a = nb::none(),
        "fpType"_a = FPType::MorganFP, "");

  m.def("GetCountFPs", &getCountFPBulkPy, "molecules"_a = nb::none(),
        "fpType"_a = FPType::MorganFP, "");

  m.def("GetFPs", &getFPBulkPy, "molecules"_a = nb::none(),
        "fpType"_a = FPType::MorganFP, "");

  m.def("FingerprintGeneratorFromJSON",
        [](const std::string &jsonStr) {
          return generatorFromJSON(jsonStr).release();
        },
        "jsonString"_a,
        "Deserialize a FingerprintGenerator from a JSON string",
        nb::rv_policy::take_ownership);

  AtomPairWrapper::exportAtompair(m);
  MorganWrapper::exportMorgan(m);
  RDKitFPWrapper::exportRDKit(m);
  TopologicalTorsionWrapper::exportTopologicalTorsion(m);
}

}  // namespace FingerprintWrapper
}  // namespace RDKit
