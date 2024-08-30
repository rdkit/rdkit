//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLPROCESSING_H
#define RD_MOLPROCESSING_H

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <DataStructs/BitVects.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/GeneralFileReader.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <mutex>
#endif

namespace RDKit {
namespace MolProcessing {

auto defaultSupplierOptions = GeneralMolSupplier::SupplierOptions();

namespace details {
inline std::unique_ptr<FileParsers::MolSupplier> getSupplier(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options) {
#ifdef RDK_BUILD_THREADSAFE_SSS
  static std::once_flag flag;
  std::call_once(flag, []() { defaultSupplierOptions.numWriterThreads = 4; });
#endif
  return GeneralMolSupplier::getSupplier(fileName, options);
}

}  // namespace details

#ifdef RDK_BUILD_THREADSAFE_SSS
namespace {
std::mutex &fp_mutex_get() {
  // create on demand
  static std::mutex _mutex;
  return _mutex;
}

void fp_mutex_create() {
  std::mutex &mutex = fp_mutex_get();
  std::lock_guard<std::mutex> test_lock(mutex);
}

std::mutex &get_fp_mutex() {
  static std::once_flag flag;
  std::call_once(flag, fp_mutex_create);
  return fp_mutex_get();
}
}  // namespace
#endif
//! \brief Get fingerprints for all of the molecules in a file
/*!
   \param fileName the name of the file to read
   \param options options controlling how the file is read, if not provided
           four threads will be used whegn reading the file
   \param generator the fingerprint generator to use, if not provided,
           Morgan fingerprints with radius of 3 will be used.

   \return an ExplicitBitVect,bitset pair, the first containing the
           fingerprints and the second a bitset indicating which molecules were
           successfully read
*/
template <typename OutputType = std::uint32_t>
std::vector<std::unique_ptr<ExplicitBitVect>> getFingerprintsForMolsInFile(
    const std::string &fileName,
    const GeneralMolSupplier::SupplierOptions &options = defaultSupplierOptions,
    FingerprintGenerator<OutputType> *generator = nullptr) {
  auto suppl = details::getSupplier(fileName, options);

  std::unique_ptr<FingerprintGenerator<OutputType>> morgan;
  if (generator == nullptr) {
    morgan.reset(MorganFingerprint::getMorganGenerator<OutputType>(3));
    generator = morgan.get();
  }
  std::map<unsigned int, std::unique_ptr<ExplicitBitVect>> fingerprints;

  // if we are using a multi-threaded supplier then we can register a write
  // callback to do our processing multi-threaded too
#ifdef RDK_BUILD_THREADSAFE_SSS
  auto tsuppl =
      dynamic_cast<v2::FileParsers::MultithreadedMolSupplier *>(suppl.get());
  if (tsuppl) {
    auto fpfunc = [&](RWMol &mol, const std::string &, unsigned int recordId) {
      auto fp = generator->getFingerprint(mol);
      {
        std::lock_guard<std::mutex> lock(get_fp_mutex());
        fingerprints[recordId].reset(fp);
      }
    };
    tsuppl->setWriteCallback(fpfunc);
    while (!tsuppl->atEnd()) {
      auto mol = tsuppl->next();
    }
    auto maxv = 0u;
    for (const auto &pr : fingerprints) {
      maxv = std::max(maxv, pr.first);
    }
    std::vector<std::unique_ptr<ExplicitBitVect>> fp_res(maxv);
    for (auto &fp : fingerprints) {
      fp_res[fp.first - 1] = std::move(fp.second);
    }
    return fp_res;
  } else {
#else
  {
#endif
    // otherwise we just loop through the molecules
    std::vector<std::unique_ptr<ExplicitBitVect>> fp_res;
    while (!suppl->atEnd()) {
      auto mol = suppl->next();
      if (mol) {
        auto fp = generator->getFingerprint(*mol);
        fp_res.emplace_back(fp);
      } else {
        fp_res.emplace_back(nullptr);
      }
    }
    return fp_res;
  }
}

}  // namespace MolProcessing
}  // namespace RDKit
#endif
